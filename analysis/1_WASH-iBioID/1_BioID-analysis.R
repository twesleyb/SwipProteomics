#!/usr/bin/env Rscript

# title: WASH iBioID Proteomics Analysis
# description: Preprocessing and statistical analysis of WASH1 (Washc1) iBioID.
# author: Tyler W A Bradshaw

## User parameters to change:
FDR_alpha = 0.05 # FDR significance threshold for protein enrichment.
enrichment_threshold = log2(3.0) # enrichment threshold.

## Input data in root/data/BioID.zip/
zipfile = "BioID.zip"
datafile = "BioID_raw_protein.csv" 

## Output in root/tables:
# * WASH_BioID_Results.xlsx

## Output in root/data:
# * wash_interactome.rda # The WASH iBioID proteome.


## Functions ------------------------------------------------------------------

# Get the repository's root directory.
getrd <- function(here=getwd(), dpat= ".git") {
	in_root <- function(h=here, dir=dpat) { 
		check <- any(grepl(dir,list.dirs(h,recursive=FALSE))) 
		return(check)
	}
	# Loop to find root.
	while (!in_root(here)) { 
		here <- dirname(here) 
	}
	root <- here
	return(root)
}


# Check CV of WASH, Control, and SPQC groups.
check_group_CV <- function(tidy_prot) {
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


# reformat the data for DEP
reformatDEP <- function(tidy_prot,gene_map) {
  # NOTE: this function only works with expected input
  # returns a summarized experiment object formated for DEP
  # NOTE: you can access the data contained in a sumamrized experiment object with 
  # library(SummarizedExperiment)
  # rowDat(se)
  # colData(se)
  # assay(se)
  # cast data into wide format, the data.frame should contain 
  # unique protein IDs in 'ID' column
  dep_prot <- tidy_prot %>% 
  	dplyr::filter(!grepl("QC",Sample)) %>% # drop QC
  	dcast(Accession ~ Sample, value.var = "Intensity") %>%
  		as.matrix(rownames="Accession") %>% 
  		as.data.table(keep.rownames="ID")
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
  samples <- unique(tidy_prot$Sample)
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


## Prepare the workspace -------------------------------------------------------

# Load renv
root <- getrd()
renv::load(root)

# imports
suppressPackageStartupMessages({
	library(dplyr) # For manipulating the data
	#library(edgeR) # For statitical comparisons
	suppressWarnings({
	  library(getPPIs) # For mapping gene identifiers
	})
	library(geneLists) # For a list of mito proteins
	library(data.table) # For working with data.tables
})

# Load any additional project specific functions and data.
suppressWarnings({
  suppressMessages({ devtools::load_all() })
})

# Project directories.
datadir <- file.path(root,"data") # key datasets
rdatdir <- file.path(root,"rdata") # temp data files
tabsdir <- file.path(root,"tables") # final xlsx tables
downdir <- file.path(root,"downloads") # misc/temp files

# Create dirs if they dont exist.
if (!dir.exists(tabsdir)){ dir.create(tabsdir) }
if (!dir.exists(rdatdir)){ dir.create(rdatdir) }
if (!dir.exists(downdir)){ dir.create(downdir) }


## Load the raw proteomics data -----------------------------------------------

# Extract the raw data from zipped file
myfile <- file.path(datadir,zipfile)
unzip(myfile) # unzip 

# Read into R with data.table::fread
myfile <- file.path(getwd(),tools::file_path_sans_ext(zipfile),datafile)
raw_prot <- fread(myfile)

# Clean-up
myfile <- file.path(downdir,tools::file_path_sans_ext(zipfile))
unlink(myfile,recursive=TRUE)
unlink("./BioID", recursive=TRUE)

# Tidy the data
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
mito_list <- geneLists("mitocarta2")
if (length(mito_list)==0) { 
	stop("Problem loading mitochrondrial contaiminants.")
}
data(list=mito_list)
mito_entrez <- unlist(mitocarta2,use.names=FALSE)
mito_prot <- getIDs(mito_entrez,from="entrez",to="uniprot",species="mouse")

# Status
nMito <- sum(mito_prot %in% tidy_prot$Accession)
warning(paste(nMito,"mitochondrial proteins will be removed as contaminants."))

# Remove mitochondrial proteins as contaminants.
tidy_prot <- tidy_prot %>% dplyr::filter(Accession %notin% mito_prot)

# Summary:
nProt <- length(unique(tidy_prot$Accession))
message(paste0("\nTotal number of proteins quantified: ",
	      formatC(nProt,big.mark=","),"."))

# Cast tidy prot into a matrix.
tidy_dm <- tidy_prot %>% as.data.table() %>%
	dcast(Accession + Description + Peptides ~ Sample, 
	      value.var = "Intensity")

# Status.
message("\nSummary of group CVs:")
knitr::kable(check_group_CV(tidy_prot))


## create gene map ------------------------------------------------------------

# map uniprot to entrez
uniprot <- unique(tidy_prot$Accession)
entrez <- mgi_batch_query(uniprot,quiet=TRUE)

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


## Statistical testing with DEP::limma ----------------------------------------

# coerce data into dep se object
dep_se <- reformatDEP(imp_prot,gene_map)

# normalize the data
message("\nPerforming VSN normalization.")

# suppress meanSdPlot message
suppressMessages({ vsn_se <- DEP::normalize_vsn(dep_se) })

# analysis of protein-level differential abundance based on linear models and
# empirical Bayes statistics with DEP (wrapped around limma)
dep_results <- DEP::test_diff(vsn_se, type = "manual", test = "WASH_vs_Control")

# NOTE: from limma::topTable() adjust.method = BH 'Benjamini Hochberg'
dep_results <- DEP::add_rejections(dep_results, alpha = FDR_alpha)

# generate a results table
results_df <- DEP::get_results(dep_results)

## Clean-up results
results_df$significant <- NULL
idy <- apply(results_df,2,function(x) all(is.na(x)))
colnames(results_df) <- gsub("WASH_vs_Control_","",colnames(results_df))
results_df <- results_df[,!idy]
colnames(results_df)[colnames(results_df) == "p.val"] <- "PValue"
colnames(results_df)[colnames(results_df) == "p.adj"] <- "FDR"
colnames(results_df)[colnames(results_df) == "ratio"] <- "logFC"
colnames(results_df)[colnames(results_df) == "name"] <- "Symbol"
colnames(results_df)[colnames(results_df) == "ID"] <- "Accession"
results_df <- results_df %>% arrange(PValue,logFC)
results_df$significant <- NULL

# convert logFC column to percent control
pc <- 100*(2^results_df$logFC)
results_df <- tibble::add_column(results_df,`Percent Control (%)`=pc,
				  .after="logFC")

# Summary:
sig <- results_df$FDR < FDR_alpha
up <- results_df$logFC > enrichment_threshold
nsig <- sum(sig & up)
message(paste0("\nNumber of significantly enriched proteins ",
	      "(log2FC > ",round(enrichment_threshold,2),
	      "; FDR < ",FDR_alpha,"): "), nsig,".")

# add gene identifiers to data table
idx <- match(results_df$Accession,gene_map$uniprot)
results_df <- tibble::add_column(results_df,
				  "Entrez"=gene_map$entrez[idx],
			          .after="Accession")


## Save the data --------------------------------------------------------------

# Create list of results:
results_list <- list("Raw Protein" = tidy_dm, "BioID Results" = results_df)

# Add the mitochondrial proteins that were removed.
df <- raw_prot %>% dplyr::filter(Accession %in% mito_prot) %>% 
	dplyr::select(Accession)
results_list[["Mitochondrial Contaiminants"]] <- df

# Write to file.
message("\nSaving results.")
myfile <- file.path(tabsdir,"S1_WASH_BioID_Results.xlsx")
write_excel(results_list,myfile)

# Save results as rdata for downstream analysis.
dep_results <- results_df
myfile <- file.path(rdatdir,"WASH_BioID_Results.RData")
saveRDS(dep_results,myfile)

# Save final results for R package in root/data.
wash_interactome <- unique(results_df$Accession[sig & up])
myfile <- file.path(datadir, "wash_interactome.rda")
save(wash_interactome,file=myfile,version=2)
