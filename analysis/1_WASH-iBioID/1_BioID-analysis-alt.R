#!/usr/bin/env Rscript

# title: WASH iBioID Proteomics Analysis
# description: Preprocessing and statistical analysis of WASH1 (Washc1) iBioID.
# author: Tyler W A Bradshaw

## User parameters to change:
FDR_alpha = 0.1 # FDR significance threshold for protein enrichment.
enrichment_threshold = log2(3.0) # enrichment threshold.

## Input data in root/data/BioID.zip/
zipfile = "BioID.zip"
datafile = "BioID_raw_protein.csv" 

## Output in root/tables:
# * WASH_BioID_Results.xlsx

## Output in root/data:
# * wash_interactome.rda # The WASH iBioID proteome.

#---------------------------------------------------------------------
## Misc function - getrd().
#---------------------------------------------------------------------

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
	cv_summary <- tmp_res %>% filter(!is.na(CV)) %>% 
		group_by(Group) %>% 
		summarize('Mean(CV)' = mean(CV,na.rm=TRUE),
			  'SD(CV)' = sd(CV, na.rm=TRUE),
			  'n' = length(CV),.groups="drop")
	return(cv_summary)
}

#-------------------------------------------------------------------------------
## Prepare the workspace.
#-------------------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root)

# Imports.
suppressPackageStartupMessages({
	library(dplyr) # For manipulating the data.
	#library(edgeR) # For statitical comparisons.
	suppressWarnings({
	  library(getPPIs) # For mapping gene identifiers.
	})
	library(geneLists) # For a list of mito proteins.
	library(data.table) # For working with data.tables.
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
if (!dir.exists(rdatdir)){ dir.create(rdatdir) }
if (!dir.exists(downdir)){ dir.create(downdir) }


#-------------------------------------------------------------------------------
## Load the raw data.
#-------------------------------------------------------------------------------

# Extract the raw data from zipped file.
myfile <- file.path(datadir,zipfile)
unzip(myfile) # unzip 

# Read into R with data.table::fread.
myfile <- file.path(getwd(),tools::file_path_sans_ext(zipfile),datafile)
raw_prot <- fread(myfile)

# Clean-up.
myfile <- file.path(downdir,tools::file_path_sans_ext(zipfile))
unlink(myfile,recursive=TRUE)
unlink("./BioID", recursive=TRUE)

# Tidy the data.
message("\nLoading raw Swip BioID protein data.")
tidy_prot <- tidyProt(raw_prot,species="Mus musculus",
		      id.vars=c("Accession","Description","Peptides"))

# Insure that keratins have been removed--typically the proteomics core removes
# most of these.
idx <- grepl("Keratin|keratin",tidy_prot$Description)
keratins <- tidy_prot %>% filter(idx) %>% dplyr::select(Accession) %>% 
	unlist() %>% unique()
warning(paste(length(keratins),
	      "Keratin proteins remain, and  will be removed."))
tidy_prot <- tidy_prot %>% filter(Accession %notin% keratins)

# Load mitochondrial protein list from twesleyb/geneLists.
mito_list <- geneLists("mitocarta2")
if (length(mito_list)==0) { 
	stop("Problem loading mitochrondrial contaiminants.")
}
data(list=mito_list)
mito_entrez <- unlist(mitocarta2,use.names=FALSE)
mito_prot <- getIDs(mito_entrez,from="entrez",to="uniprot",species="mouse")

# Status.
nMito <- sum(mito_prot %in% tidy_prot$Accession)
warning(paste(nMito,"mitochondrial proteins will be removed as contaminants."))

# Remove mitochondrial proteins as contaminants.
tidy_prot <- tidy_prot %>% filter(Accession %notin% mito_prot)

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

stopifnot(!any(is.na(entrez)))

# map to symbols
symbol <- getIDs(entrez,from="Entrez",to="Symbol",species="mouse")

stopifnot(!any(is.na(symbol)))

gene_map <- data.table(uniprot, entrez, symbol)


#-------------------------------------------------------------------------------
## Sample loading normalization. 
#-------------------------------------------------------------------------------

# Sample loading normalization:
message("\nPerforming sample loading normalization.")
SL_prot <- normSL(tidy_prot,groupBy="Sample")

# Check, column sums should now be equal.
message("Total intensity sums are equal after sample loading normalization:")
df <- SL_prot %>% group_by(Sample) %>% 
	summarize("Total Intensity"=sum(Intensity,na.rm=TRUE),.groups="drop")
knitr::kable(df)


#-------------------------------------------------------------------------------
## Sample pool normalization.
#-------------------------------------------------------------------------------

# Perform normalization to QC samples.
message("\nPerforming sample pool normalization used pooled QC samples.")
SPN_prot <- normSP(SL_prot,pool="QC")


#-------------------------------------------------------------------------------
## Protein level filtering.
#-------------------------------------------------------------------------------

# Status.
message("\nFiltering proteins...")

# Remove one hit wonders -- proteins identified by a single peptide.
one_hit_wonders <- unique(SPN_prot$Accession[SL_prot$Peptides == 1])
n_ohw <- length(one_hit_wonders)
filt_prot <- SPN_prot %>% filter(Peptides>1)
message(paste("... Number of one-hit-wonders:",n_ohw))

# Remove proteins with unreliable QC (more than 1 missing values).
df <- filt_prot %>% filter(grepl("QC",Sample)) %>% group_by(Accession) %>%
	summarize(N=length(Intensity),
		  n_missing=sum(is.na(Intensity)),.groups="drop")
out <- unique(df$Accession[df$n_missing>0])
filt_prot <- filt_prot %>% filter(Accession %notin% out)
n_out <- length(out)
message(paste("... Number of proteins with missing QC data:",n_out))

# Remove proteins with more than 50% missingness as these cannot be imputed.
df <- filt_prot %>% filter(!grepl("QC",Sample)) %>% group_by(Accession) %>%
	summarize(N=length(Intensity),
		  n_missing=sum(is.na(Intensity)),.groups="drop")
out <- unique(df$Accession[df$n_missing>0.5*df$N])
filt_prot <- filt_prot %>% filter(Accession %notin% out)
n_out <- length(out)
message(paste("... Number of proteins with too many missing values:",n_out))

# Insure that we are still working with a data.table.
filt_prot <- as.data.table(filt_prot)

# Status.
prots <- unique(filt_prot$Accession)
n_prot <- length(prots)
message(paste0("\nFinal number of quantifiable proteins: ",
	      formatC(n_prot,big.mark=","),"."))


#-------------------------------------------------------------------------------
## Impute missing values. 
#-------------------------------------------------------------------------------

# Proteins with missing values are less abundant than those without. 
# This is evidence that missing values are MNAR and can be imputed with
# the KNN algorithm.
message("\nImputing missing protein values using the KNN algorithm (k=10).")
imp_prot <- imputeKNNprot(filt_prot,k=10,rowmax=0.5,colmax=0.8,quiet=FALSE)

# Status.
knitr::kable(check_group_CV(imp_prot))


#-------------------------------------------------------------------------------
## Statistical testing with DEP::limma
#-------------------------------------------------------------------------------

# munge data into format for DEP
prot_df <- imp_prot %>% 
	filter(!grepl("QC",Sample)) %>% # drop QC
	dcast(Accession ~ Sample, value.var = "Intensity") %>%
	as.matrix(rownames="Accession") %>%  # coerce to table with ID col
	as.data.table(keep.rownames="ID")

# protein data.frame should contain unique protein IDs in 'ID' column
# check for duplicates
stopifnot(!any(duplicated(prot_df$ID))) # there should be no duplicate IDs

# there should be no NA, check for NA
stopifnot(!any(is.na(prot_df$ID))) # there should be no NA

# protein data.frame should contain unique protein names in 'name' column
# we will annotate proteins with gene 'symbol's
prot_df$name <- gene_map$symbol[match(prot_df$ID,gene_map$uniprot)]

# check for NA
stopifnot(!any(is.na(prot_df$name))) # there should be no NA

# Multiple Uniprot Accession IDs may be mapped to the same gene Symbol.
# If any duplicated, make names unique base::make.unique.
prot_df$name <- make.unique(prot_df$name,sep="-")

# check for duplicates
stopifnot(!any(duplicated(prot_df$name))) # there should be no duplicate names


## prepare exp_design for DEP --------------------------------------------------
# build experiment design data.frame from SWIP TMT sample data
# NOTE: exp_design should contain colums for 'label', 'condition' 
# and 'replicate' information
# NOTE: 'condition' is important, it must be the contrast of interest
# NOTE: DEP currently does not support more complicated experimental designs.
# Formula Model:
# 	>>> 	~ 0 + Genotype

# create exp design
samples <- unique(tidy_prot$Sample)
condition <- sapply(strsplit(samples," "),"[",2)
replicate <- sapply(strsplit(samples," "),"[",3)
exp_design <- data.frame(label = samples, condition, replicate) %>%
	filter(condition != "QC")

# check for required columns in input
stopifnot(all(c("label","condition","replicate") %in% colnames(exp_design)))


## build SummarizedExperiment (se) object -------------------------------------

# use DEP helper function to build a SE object
# specify the column indices containing the numeric data (idy)
idy <- which(colnames(prot_df) %notin% c("ID","name"))
prot_se <- DEP::make_se(prot_df,columns=idy,exp_design)

# NOTE: you can access the data contained in a sumamrized experiment object with 
# functions from the SummarizedExperiment package, e.g.:
# library(SummarizedExperiment)
# rowDat(se)
# colData(se)
# assay(se)


## Normalization --------------------------------------------------------------

# Normalize the data
# suppress meanSdPlot message --> FIXME: generate plots
message("\nPerforming VSN normalization.")
suppressMessages({ prot_norm <- DEP::normalize_vsn(prot_se) })


## Protein-level analysis of differential abundance ----------------------------
# Differential enrichment analysis  based on linear models and empherical Bayes
# statistics

# define all contrasts of intrafraction groups

# loop to perform tests for every intra-fraction comparison
results <- DEP::test_diff(prot_norm, type = "manual", test = "WASH_vs_Control")

# Denote significant proteins based on user defined cutoffs
# FIXME: what is p.adjust method?
# from limma::topTable() adjust.method = BH 'Benjamini Hochberg'
results <- DEP::add_rejections(results, alpha = FDR_alpha)

# Generate a results table
dep_results <- DEP::get_results(results)


# Clean-up results ------------------------------------------------------------

# clean-up 
dep_results$significant <- NULL
idy <- apply(dep_results,2,function(x) all(is.na(x)))
colnames(dep_results) <- gsub("WASH_vs_Control_","",colnames(dep_results))
dep_results <- dep_results[,!idy]
colnames(dep_results)[colnames(dep_results) == "p.val"] <- "PValue"
colnames(dep_results)[colnames(dep_results) == "p.adj"] <- "FDR"
colnames(dep_results)[colnames(dep_results) == "ratio"] <- "logFC"
colnames(dep_results)[colnames(dep_results) == "name"] <- "Symbol"
colnames(dep_results)[colnames(dep_results) == "ID"] <- "Accession"
dep_results <- dep_results %>% arrange(PValue,logFC)
dep_results$significant <- NULL

# Convert logCPM column to percent control.
pc <- 100*(2^dep_results$logFC)
dep_results <- tibble::add_column(dep_results,`Percent Control (%)`=pc,
				  .after="logFC")

# Summary:
sig <- dep_results$FDR < FDR_alpha
up <- dep_results$logFC > enrichment_threshold
nsig <- sum(sig & up)
message(paste0("\nNumber of significantly enriched proteins ",
	      "(log2FC > ",round(enrichment_threshold,2),
	      "; FDR < ",FDR_alpha,"): "), nsig,".")

# Add identifiers to data table.
idx <- match(dep_results$Accession,gene_map$uniprot)
dep_results <- tibble::add_column(dep_results,
				  "Entrez"=gene_map$entrez[idx],
			          .after="Accession")

## Save the data --------------------------------------------------------------

# Create list of results:
results_list <- list("Raw Protein" = tidy_dm, "BioID Results" = dep_results)

# Add the mitochondrial proteins that were removed.
df <- raw_prot %>% filter(Accession %in% mito_prot) %>% 
	dplyr::select(Accession)
results_list[["Mitochondrial Contaiminants"]] <- df

# Write to file.
message("\nSaving results.")
myfile <- file.path(tabsdir,"S1_WASH_BioID_Results.xlsx")
write_excel(results_list,myfile)

# Save results as rdata for downstream analysis.
myfile <- file.path(rdatdir,"WASH_BioID_Results.RData")
saveRDS(results,myfile)

# Save final results for R package in root/data.
wash_interactome <- results %>% filter(candidate != "no")
myfile <- file.path(datadir, "wash_interactome.rda")
save(wash_interactome,file=myfile,version=2)
