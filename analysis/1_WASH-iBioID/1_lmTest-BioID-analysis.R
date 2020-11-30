#!/usr/bin/env Rscript

# title: WASH iBioID Proteomics Analysis
# description: Preprocessing and statistical analysis of WASH1 (Washc1) iBioID.
# author: Tyler W A Bradshaw

## ---- inputs

FDR_alpha = 0.1 # FDR significance threshold for protein enrichment
enrichment_threshold = log2(2.0) # enrichment threshold

## Input data in root/data/BioID.zip/
zipfile = "BioID.zip"
datafile = "BioID_raw_protein.csv" 

## Output in root/tables:
# * WASH_BioID_Results.xlsx

## Output in root/data:
# * wash_interactome.rda # The WASH iBioID proteome.

## ---- renv

root <- "~/projects/SwipProteomics"
renv::load(root)

# library(SwipProteomics)
devtools::load_all(root)

# JC's bioid annotations 
data(bioid_anno)


## ---- Functions 

lquote <- function(string, single = TRUE) {
  # wrap a string in single or double quotes
  single_quote <- "'"
  double_quote <- "\""
  if (single) {
    # Single quote.
    return(paste0(single_quote, string, single_quote))
  } else {
    # Double quote.
    return(paste0(double_quote, string, double_quote))
  }
}


check_group_CV <- function(tidy_prot) {
  # Check CV of WASH, Control, and SPQC groups.
  # Calculate the Coefficient of Covarition (CV).
	cv <- function(x) { sd(x)/mean(x) }
	# Annotate with QC and biological replicates groups.
	tmp_dt <- tidy_prot %>% as.data.table()
	tmp_dt$Group <- sapply(strsplit(tmp_dt$Sample," "),"[",2)
	# Calculate protein-wise CV.
	tmp_res <- tmp_dt %>% group_by(Accession, Group) %>% 
		dplyr::summarize(n = length(Intensity), 
				 CV=cv(Intensity),.groups="drop")
	# NOTE: NA can arise if protein was not quantified in all replicates.
	cv_summary <- tmp_res %>% dplyr::filter(!is.na(CV)) %>% 
		group_by(Group) %>% 
		summarize('Mean(CV)' = mean(CV,na.rm=TRUE),
			  'SD(CV)' = sd(CV, na.rm=TRUE),
			  'n' = length(CV),.groups="drop")
	return(cv_summary)
}


reformatDEP <- function(prot_df, gene_map) {
  # reformat the data for DEP
  # NOTE: this function only works with expected input
  # returns a summarized experiment object formated for DEP
  # NOTE: to access data within a SummarizedExperiment object, use:
  # library(SummarizedExperiment)
  # rowDat(se)
  # colData(se)
  # assay(se)
  # cast data into wide format, the data.frame should contain 
  # unique protein IDs in 'ID' column
  dep_prot <- prot_df %>% 
  	dcast(Accession ~ Sample, value.var = "Intensity") %>%
  		as.matrix(rownames="Accession") %>% 
  		as.data.table(keep.rownames="ID")
  # IDs should not be duplicate or missing
  stopifnot(!any(duplicated(dep_prot$ID)))
  stopifnot(!any(is.na(dep_prot$ID)))
  # protein data.frame should contain unique protein names in 'name' column
  # we will annotate proteins with gene 'symbol's
  dep_prot$name <- gene_map$symbol[match(dep_prot$ID,gene_map$uniprot)]
  stopifnot(!any(is.na(dep_prot$name)))
  # Multiple Uniprot Accession IDs may be mapped to the same gene Symbol.
  # If any duplicated, make names unique base::make.unique.
  dep_prot$name <- make.unique(dep_prot$name,sep="-")
  stopifnot(!any(duplicated(dep_prot$name)))
  ## prepare exp_design for DEP
  # build experiment design data.frame from SWIP TMT sample data
  # NOTE: exp_design should contain colums for 'label', 'condition' 
  # and 'replicate' information
  # NOTE: 'condition' is important, it must be the contrast of interest
  # NOTE: DEP currently does not support more complicated experimental designs.
  # create exp design
  samples <- unique(prot_df$Sample)
  condition <- sapply(strsplit(samples," "),"[",2)
  replicate <- sapply(strsplit(samples," "),"[",3)
  exp_design <- data.frame(label = samples, condition, replicate) %>%
  	dplyr::filter(condition != "QC")
  stopifnot(all(c("label","condition","replicate") %in% colnames(exp_design)))
  ## build SummarizedExperiment (se) object
  # use DEP helper function to build a SE object
  # specify the column indices containing the numeric data (idy)
  idy <- which(colnames(dep_prot) %notin% c("ID","name"))
  dep_se <- DEP::make_se(dep_prot,columns=idy,exp_design)
  return(dep_se)
}


tidyProt <- function(raw_data,id.vars,species=NULL,
		     samples=NULL,summary=FALSE){
	# description: a function to tidy proteomics data.
	suppressPackageStartupMessages({
		library(data.table)
		library(tibble)
		library(dplyr)
	})
	dt <- as.data.table(raw_data) %>%
		melt(id.vars = id.vars,
		     variable.name="Sample",
		     value.name="Intensity",
		     variable.factor=FALSE) # Don't coerce to factor.
	# Remove proteins that do not coorespond to species of interest.
	idx <- grepl(paste0("OS=",species),dt$Description)
	if (!all(idx)) {
		n_out <- length(unique(dt$Accession[!idx]))
		msg <- paste(n_out,"proteins are not from", lquote(species),
			     "and will be removed.")
		warning(msg,call.=FALSE)
		dt <- dt %>% filter(grepl(paste0("OS=",species),Description))
	}
	# Insure zeros are NA.
	dt$Intensity[dt$Intensity==0] <- NA
	# Return tidy df.
	return(as.data.table(dt))
}


normSL <- function(tp, groupBy=NULL){
	suppressPackageStartupMessages({
		library(dplyr)
		library(data.table)
	})
	# Create an expression that will be evaluated to dynamically 
	# group data based on user provided grouping factors. Default's
	# to Sample -- every replicate.
	tp <- ungroup(as.data.frame(tp))
	cmd <- paste0("group_by(tp, ",paste(groupBy,collapse=", "),")")
	# Calculate column sums for grouped samples.
	# FIXME: if .groups is set to 'drop' to avoid Warning msg, then error?
	data_list <- eval(parse(text=cmd)) %>% 
		summarize(Total=sum(Intensity,na.rm=TRUE),.groups="drop") %>%
		group_split()
	# Calculate normalization factors.
	data_SL <- lapply(data_list,function(x) {
				     x$"Mean(Total)" <- mean(x$Total,na.rm=TRUE)
				     x$NormFactor <- x$"Mean(Total)"/x$Total
				     x$Norm <- x$Total*x$NormFactor
			      return(x) }) %>% bind_rows()
	# Collect normalization factors as named vector.
	SL_factors <- data_SL$NormFactor
	names(SL_factors) <- as.character(data_SL[["Sample"]])
	# Normalize measurements by sample factors.
	tp$Intensity <- tp$Intensity * SL_factors[as.character(tp[["Sample"]])]
	return(as.data.table(tp))
}


normSP <- function(tp, pool){
	# perform normalization to pooled QC samples
	# Store a copy of the data.
	tp <- ungroup(tp)
	tp_copy <- tp
	# Group pooled together.
	tp$Group <- as.numeric(grepl(paste(pool,collapse="|"),tp$Sample))
	tp_list <- tp %>% group_by(Accession,Group) %>% 
		dplyr::summarize(Mean_Intensity=mean(Intensity,na.rm=TRUE),
			  n = length(Intensity), .groups="drop") %>%
	as.data.table() %>% 
	arrange(Accession,Group) %>% 
	group_by(Accession) %>% group_split()
	# Loop to calculate normalization factors.
        new_list <- list()
	for (i in 1:length(tp_list)){
		x <- tp_list[[i]]
		x$NormFactor <- c(x$Mean_Intensity[2]/x$Mean_Intensity[1],1)
		x$Norm_Mean_Intensity <- x$Mean_Intensity * x$NormFactor
		new_list[[i]] <- x
	}
	tp_list <- new_list
	# Collect in a df.
	df <- do.call(rbind,tp_list) %>% 
		dplyr::select(Accession,Group,NormFactor)
	# Merge with input data.
	tp_norm <- left_join(tp,df,by=c("Accession","Group"))
	# Perform normalization step.
	tp_norm$Intensity <- tp_norm$Intensity * tp_norm$NormFactor
	tp_norm <- tp_norm %>% dplyr::select(colnames(tp_copy))
	tp_norm <- as.data.table(tp_norm)
	# Return the normalized data.
	return(tp_norm)
}


imputeKNNprot <- function(tidy_prot,ignore="QC",k=10,rowmax=0.5,colmax=0.8,
			  quiet=TRUE){
	# Determine how many missing values each protein has,
	# and then determine the permisible number of missing values 
	# for any given protein.
	# Imports.
	suppressPackageStartupMessages({
		library(impute)
		library(dplyr)
		library(tibble)
		library(data.table)
	})
	# Store a copy of the input data.
	tp <- tp_in <- tidy_prot
	# How many missing values are there?
	tp$Intensity[tp$Intensity==0] <- NA
	N_missing <- sum(is.na(tp$Intensity))
	if (N_missing ==0) { 
		message(paste("There are no missing values.",
			   "Returning untransformed data."))
		return(tp_in)
	}
	# Separate data to be imputed.
	if (!is.null(ignore)) {
		# Do not impute:
		tp_ignore <- tp %>% filter(grepl(ignore,Sample)) %>% 
			dplyr::select(Accession,Sample,Intensity) %>% 
			as.data.table
		# Data to be imputed:
		tp_impute <- tp %>% filter(!grepl(ignore,Sample)) %>% 
			as.data.table
	} else {
		tp_ignore <- NULL
	}
	# Cast the data into a matrix.
	dm <- tp_impute %>%
		dcast(Accession ~ Sample, value.var="Intensity") %>% 
		as.matrix(rownames=TRUE)
	# Don't impute rows (proteins) with too many missing values.
	n_missing <- apply(dm,1,function(x) sum(is.na(x)))
	limit <- ncol(dm) * rowmax
	rows_to_ignore <- n_missing > limit
	n_ignore <- sum(rows_to_ignore)
	if (n_ignore > 0) {
		msg <- paste(n_ignore,"proteins have more than",limit,
			     "missing values, these values cannot be imputed,",
			     "and will be ignored.")
		warning(msg,call.=FALSE)
	}
	# Total number of missing values to be imputed.
	n_imputed <- sum(is.na(dm[!rows_to_ignore,]))
	if (!quiet){
		message(paste("There are",n_imputed, "missing",
		      "values that will be replaced by imputing."))
	}
	# Perform KNN imputing.
		# 
	if (quiet) {
		# Suppress output from impute.knn.
		silence({
			data_knn <- impute.knn(log2(dm[!rows_to_ignore,]),
					       k=k,colmax=colmax,rowmax=rowmax)
		})
	} else {
		data_knn <- impute.knn(log2(dm[!rows_to_ignore,]),
				       k=k,colmax=colmax,rowmax=rowmax)
	}
	# Collect the imputed data.
	dm_knn <- dm
	dm_knn[!rows_to_ignore,] <- 2^data_knn$data
	# Melt into tidy df.
	dt_knn <- as.data.table(dm_knn,keep.rownames="Accession")
	tp_imputed <- melt(dt_knn,id.vars="Accession",
			  variable.name="Sample",value.name="Intensity")
	# Combine with any samples that were ignored.
	tp_imputed <- rbind(tp_ignore,tp_imputed)
	# Combine with input meta data.
	tp_in$Intensity <- NULL
	tp_in$Sample <- as.character(tp_in$Sample)
	tp_imputed$Sample <- as.character(tp_imputed$Sample)
	tp_out <- left_join(tp_in,tp_imputed,by=c("Sample","Accession")) %>%
		as.data.table
	return(tp_out)
} #EOF


## Prepare the workspace -------------------------------------------------------

# imports
suppressPackageStartupMessages({
	library(dplyr) # For manipulating the data
	library(getPPIs) # For mapping gene identifiers
	library(geneLists) # For a list of mito proteins
	library(data.table) # For working with data.tables
})


# project directories:
datadir <- file.path(root,"data") # key datasets
rdatdir <- file.path(root,"rdata") # temp data files
tabsdir <- file.path(root,"tables") # final xlsx tables
downdir <- file.path(root,"downloads") # misc/temp files

# Create dirs if they dont exist.
if (!dir.exists(datadir)){ dir.create(datadir) }
if (!dir.exists(rdatdir)){ dir.create(rdatdir) }
if (!dir.exists(tabsdir)){ dir.create(tabsdir) }
if (!dir.exists(downdir)){ dir.create(downdir) }


## Load the raw proteomics data -----------------------------------------------

# extract the raw data from zipped file
myfile <- file.path(datadir, zipfile)
unzip(myfile) # unzip 

# read into R with data.table::fread
here <- getwd()
myfile <- file.path(here, tools::file_path_sans_ext(zipfile), datafile)
raw_prot <- fread(myfile)

# clean-up
myfile <- file.path(downdir,tools::file_path_sans_ext(zipfile))
unlink(myfile,recursive=TRUE)
unlink("./BioID", recursive=TRUE)

# tidy-up the data
message("\nLoading raw Swip BioID protein data.")

tidy_prot <- tidyProt(raw_prot,species="Mus musculus",
		      id.vars=c("Accession","Description","Peptides"))

# insure that keratins have been removed
idx <- grepl("Keratin|keratin",tidy_prot$Description)
keratins <- tidy_prot %>% dplyr::filter(idx) %>% dplyr::select(Accession) %>% 
	unlist() %>% unique()

warning(paste(length(keratins),
	      "Keratin proteins remain, and  will be removed."))
tidy_prot <- tidy_prot %>% dplyr::filter(Accession %notin% keratins)

# Load mitochondrial protein list from twesleyb/geneLists
data(mitocarta2)

# map entrez to uniprot
mito_entrez <- unlist(mitocarta2,use.names=FALSE)
mito_prot <- getIDs(mito_entrez,from="entrez",to="uniprot",species="mouse")

nMito <- sum(mito_prot %in% tidy_prot$Accession)
warning(paste(nMito,"mitochondrial proteins will be removed as contaminants."))

tidy_prot <- tidy_prot %>% dplyr::filter(Accession %notin% mito_prot)
nProt <- length(unique(tidy_prot$Accession))

message(paste0("\nTotal number of proteins quantified: ",
	      formatC(nProt,big.mark=","),"."))

# cast the tidy raw protein data into a matrix -- we will save this as part of
# an excel workbook with the results
tidy_dm <- tidy_prot %>% as.data.table() %>%
	dcast(Accession + Description + Peptides ~ Sample, 
	      value.var = "Intensity")

# Status.
message("\nSummary of group CVs:")
knitr::kable(check_group_CV(tidy_prot))


## create gene map ------------------------------------------------------------

# map uniprot to entrez
uniprot <- unique(tidy_prot$Accession)
entrez <- mgi_batch_query(uniprot)

# fix missing
entrez[is.na(entrez)]
missing_entrez <- c("P10853" = 319180)
entrez[names(missing_entrez)] <- missing_entrez

# map entrez to symbols
symbol <- getIDs(entrez,from="Entrez",to="Symbol",species="mouse")

# there should be no missing entrez
stopifnot(!any(is.na(entrez)))

# there should be no missing symbols
stopifnot(!any(is.na(symbol)))

gene_map <- data.table(uniprot, entrez, symbol)


## Sample loading normalization -----------------------------------------------

message("\nPerforming sample loading normalization.")

tidy_prot <- tidy_prot %>% mutate(Condition = sapply(strsplit(Sample,"\\ "),"[",2))
tidy_prot <- tidy_prot %>% mutate(Run = sapply(strsplit(Sample,"\\ "),"[",1))
tidy_prot <- tidy_prot %>% mutate(BioReplicate = sapply(strsplit(Sample,"\\ "),"[",3))
tidy_prot <- tidy_prot %>% mutate(Subject = interaction(Condition,BioReplicate))

SL_prot <- normSL(tidy_prot,groupBy="Sample")

# Check, column sums should now be equal:
message("Total intensity sums are equal after sample loading normalization:")
df <- SL_prot %>% group_by(Sample) %>% 
	summarize("Total Intensity"=sum(Intensity,na.rm=TRUE),.groups="drop")
knitr::kable(df)


## Sample pool normalization --------------------------------------------------
# Perform normalization to QC samples

message("\nPerforming sample pool normalization used pooled QC samples.")

SPN_prot <- normSP(SL_prot,pool="QC")


## Protein level filtering ----------------------------------------------------

message("\nFiltering proteins...")

# Remove one hit wonders -- proteins identified by a single peptide.
one_hit_wonders <- unique(SPN_prot$Accession[SL_prot$Peptides == 1])
n_ohw <- length(one_hit_wonders)
filt_prot <- SPN_prot %>% dplyr::filter(Peptides>1)

message(paste("... Number of one-hit-wonders:",n_ohw))

# Remove proteins with unreliable QC (more than 1 missing values).
df <- filt_prot %>% dplyr::filter(grepl("QC",Sample)) %>% group_by(Accession) %>%
	summarize(N=length(Intensity),
		  n_missing=sum(is.na(Intensity)),.groups="drop")
out <- unique(df$Accession[df$n_missing>0])
filt_prot <- filt_prot %>% dplyr::filter(Accession %notin% out)
n_out <- length(out)

message(paste("... Number of proteins with missing QC data:",n_out))

# Remove proteins with more than 50% missingness as these cannot be imputed.
df <- filt_prot %>% dplyr::filter(!grepl("QC",Sample)) %>% group_by(Accession) %>%
	summarize(N=length(Intensity),
		  n_missing=sum(is.na(Intensity)),.groups="drop")
out <- unique(df$Accession[df$n_missing>0.5*df$N])
filt_prot <- filt_prot %>% dplyr::filter(Accession %notin% out)
n_out <- length(out)

message(paste("... Number of proteins with too many missing values:",n_out))

# Insure that we are still working with a data.table.
filt_prot <- as.data.table(filt_prot)

# Status:
prots <- unique(filt_prot$Accession)
n_prot <- length(prots)
message("\nFinal number of quantifiable proteins: ", 
	formatC(n_prot,big.mark=","))


## Impute missing values ------------------------------------------------------
# Proteins with missing values are less abundant than those without. 
# This is evidence that missing values are MNAR and can be imputed with
# the KNN algorithm.

message("\nImputing missing protein values using the KNN algorithm (k=10).")

imp_prot <- imputeKNNprot(filt_prot,k=10,rowmax=0.5,colmax=0.8,quiet=FALSE)

# Status:
knitr::kable(check_group_CV(imp_prot))


## ---- statistical testing

wash_bioid <- imp_prot

swip <- gene_map$uniprot[match("Washc4",gene_map$symbol)]

fx <- log2(Intensity) ~ 0 + Condition
fm <- lm(fx, wash_bioid %>% subset(Accession == swip))

LT <- getContrast(fm,"WASH","Control")
lmTestContrast(fm, LT) %>% mutate(Pvalue=formatC(Pvalue)) %>% knitr::kable()

## loop to perform test for all proteins

fm_list <- list()
s2_list <- list()
df_list <- list()

# fit the models
proteins <- unique(wash_bioid$Accession)
for (prot in proteins) {
	fm <- lm(fx, wash_bioid %>% subset(Accession == prot))
	fm_list[[prot]] <- fm
	s2_list[[prot]] <- sigma(fm)^2
	df_list[[prot]] <- fm$df.residual
}

# t-statistic moderation
eb_var <- limma::squeezeVar(unlist(s2_list), unlist(df_list))
df_prior <- eb_var$df.prior
s2_prior <- eb_var$s2.prior
if (is.null(s2_prior)) { s2_prior <- 0 }

# with moderation
fx <- log2(Intensity) ~ 0 + Condition
fm <- lm(fx, wash_bioid %>% subset(Accession == swip))
lmTestContrast(fm,LT,s2_prior,df_prior) %>% mutate(Pvalue=formatC(Pvalue)) %>% knitr::kable()

# loop to perform moderated comparisons
result_list <- list()

for (prot in proteins) {
	fm <- fm_list[[prot]]
	result_list[[prot]] <- lmTestContrast(fm,LT,s2_prior,df_prior)
}

# collect results
results_df <- bind_rows(result_list,.id="Protein")

# p.adjust
results_df <- results_df %>% mutate(FDR = p.adjust(Pvalue,method="bonferroni"))

# sig and up
results_df <- results_df %>% 
	mutate(up = log2FC > enrichment_threshold) %>% 
	mutate(sig = FDR < FDR_alpha) %>%
	mutate(candidate = up & sig)


idx <- match(results_df$Protein,gene_map$uniprot)
results_df <- results_df %>% 
	tibble::add_column(Symbol = gene_map$symbol[idx], .after="Protein") %>%
	tibble::add_column(Entrez = gene_map$entrez[idx], .after="Symbol")


# final sort
results_df <- results_df %>% 
	arrange(desc(sig), desc(up), Pvalue, desc(log2FC)) %>%
	dplyr::select(-sig,-up)


## ---- save results
myfile <- file.path(root,"tables","S1_WASH-BioID_Results.xlsx")
write_excel(list("WASH-BioID"=results_df), myfile)


## check 
#sig_prots <- unique(results_df$Protein[results_df$candidate])
#anno_prots <- unique(bioid_anno$Protein)
#sum(anno_prots %in% sig_prots)/length(anno_prots)
#mapID(anno_prots[anno_prots %notin% sig_prots],"uniprot","symbol")
