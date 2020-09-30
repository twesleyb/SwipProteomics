#!/usr/bin/env Rscript

# title: SwipProteomics
# description: analysis of Swip TMT spatial proteomics data with MSstatsTMT
# author: Tyler W Bradshaw <twesleyb10@gmail.com>
# os: windows linux subsystem (WSL)

## Input ----------------------------------------------------------------------

# specify the projects root directory:
root = "~/projects/SwipProteomics"

# the PSM data is in root/data:
input_dir = "data/PSM.zip"

# PSM.zip contains:
input_psm = "PSM/5359_PSM_Report.xlsx"
input_samples = "PSM/5359_Sample_Report.xlsx"

## Options --------------------------------------------------------------------
load_rda = TRUE # load saved objects to speed things up

## Output ---------------------------------------------------------------------
# * several intermediate datasets in root/rdata
# * the MSstatsTMT statistical results for intrafraction comparisons saved as
#   an excel workbook in root/tables

## FUNCTIONS -----------------------------------------------------------------

mkdir <- function(...) {
	# create a new directory
	newdir <- file.path(...)
	if (!dir.exists(newdir)) { 
		dir.create(newdir)
		message(paste("created",newdir))
	}
}


fix_colname <- function(df,old_colname,new_colname) {
	# change a column's name in a data.frame
	colnames(df)[which(colnames(df) == old_colname)] <- new_colname
	return(df)
}


munge1 <- function(x) { # x = samples$ConditionFraction
	# function to reformat the data
	# extract sample 'Condition' annotation from ConditionFraction
	paste0("F",as.numeric(sapply(strsplit(x,"Control|Mutant|SPQC"),"[",2)))
}


munge2 <- function(x) { # x = samples$ConditionFraction 
	# function to reformat the data
	# extract sample 'Fraction' annotation from ConditionFraction
	sapply(strsplit(x,"[0-9]{1,2}"),"[",1)
}


munge3 <- function(x) { # x = samples$Experiment
	# coerce Experiment to replicate ID = R#
	paste0("R",as.numeric(as.factor(samples$Experiment))) 
}


reformat_cols <- function(raw_pd) {
	# make columns look like MSstats by replacing special characters with '.'
	# replace special characters in column names with "."
	chars <- c(" ","\\[","\\]","\\:","\\(","\\)","\\/","\\+","\\#","\\-")
	new_names <- gsub(paste(chars,collapse="|"),".",colnames(raw_pd))
	colnames(raw_pd) <- new_names
	# add 'X' if starts column name starts with ".."
	colnames(raw_pd) <- gsub("^\\.\\.","X..",colnames(raw_pd))
	# return the reformatted data
	return(raw_pd)
}


## Prepare the working environment --------------------------------------------

# project directories
datadir <- file.path(root,"data"); mkdir(datadir)
rdatdir <- file.path(root,"rdata"); mkdir(rdatdir)
downdir <- file.path(root,"downloads"); mkdir(downdir)


# Prepare the R environment ---------------------------------------------------

# load renv
renv::load(root, quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	suppressWarnings({ library(getPPIs) }) # FIXME: annoying warnings!
	library(data.table)
	library(MSstatsTMT)
})

# load functions in root/R
suppressPackageStartupMessages({ devtools::load_all() })


## load sample data in root/rdata --------------------------------------
# recieved this excel spreadsheet from GW (exported from PD?)

unzip(file.path(root,input_dir),exdir=downdir) # unzip into downloads

# pass meaningful colnames to read_excel
col_names <- c("Sample","Mixture","MS.Channel","drop",
	       "Channel","Proteomics ID","ConditionFraction","Experiment")
myfile <- file.path(downdir,input_samples)
samples <- readxl::read_excel(myfile,col_names=col_names)

# clean-up
unlink(myfile)


## read PSM data from excel ------------------------------------------------

# read PSM-level data exported from PD as an excel worksheet
# NOTE: this takes several minutes!
message(paste("\nLoading PSM data from ProteomeDiscoverer.",
	      "This will take several minutes."))
myfile <- file.path(downdir,input_psm)
raw_pd <- readxl::read_excel(myfile,progress=FALSE)

unlink(myfile) # rm input_psm
unlink(tools::file_path_sans_ext(basename(input_dir))) # rmdir ./PSM


## re-format sample metadata annotations for MSstats ---------------------------

# remove un-needed col
samples$drop <- NULL

# Munge 'ConditionFraction' column to 'BioFraction' and 'BioCondition'
# NOTE: BioFraction is the subcellular fraction not MSstats 'Fraction'
# NOTE: BioCondition is the treatment 'Condition' not MSstats 'Condition'
samples$BioFraction <- munge1(samples$ConditionFraction)
samples$BioCondition <- munge2(samples$ConditionFraction)

# this is how MSstatsTMT needs 'Condition' for intra-fraction contrasts:
# >>> BioCondition.BioFraction %eg% Control.F10
samples$Condition <- as.character(interaction(samples$BioCondition,
					      samples$BioFraction)) 
samples$Condition[grepl("SPQC",samples$Condition)] <- "Norm" 

# remove un-needed col
samples$ConditionFraction <- NULL

# clean-up 'Mixture' column by replacing F with M
samples$Mixture <- gsub("F","M",samples$Mixture)

# FIXME: how should BioReplicate be defined?
# for intrafraction contrasts, 'BioReplicate' should be ...
# add R# to indicate the replicate %eg% Control.F4.R2
samples$BioReplicate <- paste(samples$Condition,
			      munge3(samples$Experiment),sep=".")


## re-format PSM data for MSstatsTMT ---------------------------------------------

raw_pd <- reformat_cols(raw_pd)


## drop contaminant proteins ---------------------------------------------------

# drop PSM that are mapped to multiple proteins
idx_drop <- grepl(";",raw_pd$Master.Protein.Accessions)
filt_pd <- raw_pd[!idx_drop,]

# keep mouse proteins and drop ig prots
idx_keep <- grepl("OS=Mus musculus",filt_pd$Protein.Descriptions)
ig_prots <- c("P01631","P01646","P01665","P01680","P01746","P01750",
	      "P01786","P01864","P01878","P03975","P06330","P03987")

# drop human keratins, proteins without genes, pig trypsin, a predicted gene
misc_drop <- c("P13647","P13645","P08779", "P04264","P02533","Q04695",
	       "Q7Z794","P01626","P01637","P00761","P10404","P02538","Q3ZRW6",
	       "P04940","Q6KB66","P01723","Q80WG5") 

# filter data
idx_drop1 <- filt_pd$Master.Protein.Accessions %in% ig_prots 
idx_drop2 <- filt_pd$Master.Protein.Accessions %in% misc_drop
filt_pd <- filt_pd[idx_keep & !idx_drop1 & !idx_drop2,]

# collect all uniprot accession ids
uniprot <- unique(filt_pd$Master.Protein.Accessions)

# create gene_map
if (load_rda) {
  myfile <- file.path(root,"data","msstats_gene_map.rda")
  load(myfile)
  message("Loaded saved 'gene_map'.")
} else {
  # map uniprot to entrez
  # NOTE: this may take several minutes
  entrez <- mgi_batch_query(uniprot,quiet=TRUE)
  names(entrez) <- uniprot
  # Map any remaining missing IDs by hand.
  message("Mapping missing IDs by hand.\n")
  missing <- entrez[is.na(entrez)]
  mapped_by_hand <- c(P05214=22144, 
  		    P0CG14=214987, 
  		    P84244=15078,
  		    P05214=22144) #tub3a
  entrez[names(mapped_by_hand)] <- mapped_by_hand
  # Check: Have we successfully mapped all Uniprot IDs to Entrez?
  check <- sum(is.na(entrez)) == 0
  if (!check) { stop("Unable to map all UniprotIDs to Entrez.") }
  # Map entrez ids to gene symbols using twesleyb/getPPIs.
  symbols <- getPPIs::getIDs(entrez,from="entrez",to="symbol",species="mouse")
  # check there should be no missing gene symbols
  if (any(is.na(symbols))) { stop("Unable to map all Entrez to gene Symbols.") }
  # Create gene identifier mapping data.table.
  gene_map <- data.table(uniprot = names(entrez),
                         entrez = entrez,
  	               symbol = symbols)
  gene_map$id <- paste(gene_map$symbol,gene_map$uniprot,sep="|")
  # save the gene map
  myfile <- file.path(root,"data","msstats_gene_map.rda")
  save(gene_map,file=myfile,version=2)
  message(paste("\nSaved",basename(myfile),"in",dirname(myfile)))
}


## map Spectrum.Files to MS.Channels ------------------------------------------
# Basically, its just complicated. There were 3x 16-plex TMT experiments.
# In each, the concatenated TMT mixture was fractionated into 12 fractions in 
# order to increase analytical depth. Therefore, there were 12 * 3 = 36 mass
# spectrometry runs. Each Run is recorded as a 'Spectrum File' by Proteome
# Discover. Note that each MS run cooresponds to the measurment of all 16 
# samples and therefore the total number of TMT channels for which we 
# have made measurements is 16 x 3 Experiments = 48. In other words, a single
# 'Spectrum.File' cooresponds to 12x MS.Runs and 16x Samples. 

all_files <- filt_pd$Spectrum.File

# collect all Spectrum.Files grouped by Experiment
# split 'Spectrum.File' at first "_" to extract experiment identifiers
exp_files <- lapply(split(all_files, sapply(strsplit(all_files,"_"),"[",1)),
		    unique)
files_dt <- data.table("Experiment ID" = rep(names(exp_files),
					     times=sapply(exp_files,length)),
	               "Run" = unlist(exp_files))

# add Fraction annotation
files_dt$Fraction <- unlist({
	sapply(exp_files, function(x) as.numeric(as.factor(x)),simplify=FALSE) })

# collect all MS.Channels, grouped by Experiment
all_channels <- samples$MS.Channel
exp_channels <- split(all_channels, sapply(strsplit(all_channels,"_"),"[",1))

# as a dt
exp_dt <- data.table("Experiment ID" = rep(names(exp_channels),
					times=sapply(exp_channels,length)),
	             "MS.Channel" = unlist(exp_channels))

# save samples
myfile <- file.path(root,"data","msstats_samples.rda")
save(samples,file=myfile,version=2)
message(paste("\nSaved",myfile))


## build annotation file for MSstats -------------------------------------------
# the annotation data.frame passed to MSstats requires the following columns:
# * Run - indicates which channel within a Spectrum.File; this should match
#     Spectrum.File in raw_pd
# * Fraction - a TMT mixture may be fractionated into multiple fractions to
#     increase analytical depth; column 'Fraction' indicates the MS fraction,
#     it may be all 1
# * TechRepMixture - was a TMT mixture analyzed in technical replicate? 
#     all 1 indicates no replicates of a mixture
# * Mixture - concatenation of TMT labeled samples - an MS experiment
# * Channel - the TMT lables/channels used e.g. 126N, 134N
# * BioReplicate - indicates an individual Biological subject
# * Condition - indicates the treatment condition e.g. WT, MUT or SPQC (Norm)

message("\nBuilding sample annotation dataframe for MSstats.") 

# create annotation_dt from Spectrum.Files and MS.Runs
# add additional freatures from samples
annotation_dt <- left_join(files_dt,exp_dt,by="Experiment ID")
idx <- match(annotation_dt$"MS.Channel", samples$MS.Channel)
annotation_dt$BioFraction <- samples$BioFraction[idx]
annotation_dt$TechRepMixture <- rep(1,length(idx))
annotation_dt$Mixture <- samples$Mixture[idx]
annotation_dt$Condition <- samples$Condition[idx]
annotation_dt$Channel <- samples$Channel[idx]
annotation_dt$BioReplicate <- samples$BioReplicate[idx]
annotation_dt$BioReplicate[grepl("Norm",annotation_dt$BioReplicate)] <- "Norm"
# FIXME: how to pass additional covariates to MSstats?
# FIXME: how to handle repeated measures design? 

# Remove un-needed cols
annotation_dt$"MS.Channel" <- NULL

# save to file
myfile <- file.path(rdatdir,"annotation_dt.rda")
save(annotation_dt,file=myfile,version=2)
message(paste("\nSaved",basename(myfile),"in",dirname(myfile)))


## combine annotation_dt and raw_pd to coerce data to MSstats format -----------

if (load_rda) {
	load(file.path(rdatdir,"data_pd.rda"))
	message("Loaded saved 'data_pd'.")
} else {
	message("\nConverting PSM data to MSstatsTMT format.")
	# NOTE: this takes a considerable amount of time
	data_pd <- PDtoMSstatsTMTFormat(filt_pd, annotation_dt, rmProtein_with1Feature = TRUE)
	# save to file
	myfile <- file.path(rdatdir,"data_pd.rda")
	save(data_pd,file=myfile,version=2)
	message(paste("\nSaved",myfile))
}


# Protein summarize and normalization -----------------------------------------
# use MSstats for protein summarization	
# Sample Summary does not look correct, should be 

# NOTE: this takes a considerable amount of time
# FIXME: Joining, by = ("Run", "Channel") # unexpected output
# FIXME: remove extremely long message about normalization between runs
if (load_rda) {
  load(file.path(root,"data","msstats_prot.rda"))
  message("Loaded saved 'msstats_prot'.")
} else {
  message("\nPerforming normalization and protein sumamrization.")
  msstats_prot <- proteinSummarization(data_pd,
				  method="msstats",	
                                  global_norm=TRUE,	
                                  reference_norm=TRUE,	
                                  remove_norm_channel = TRUE)
  # save to file
  myfile <- file.path(root,"rdata","msstats_prot.rda")
  save(msstats_prot,file=myfile,version=2)
  message(paste("\nSaved",myfile))
}


## Remove protein outliers? -----------------------------------------------------
	
nprot <- unique(msstats_prot$Protein) # There are ~6900 proteins

# Subset for speed:
rand_prots <- sample(unique(msstats_prot$Protein),10)

filt_prot <- msstats_prot %>% filter(Protein %in% rand_prots)


## Protein-level statistical testing -------------------------------------------
# Tests for significant changes in protein abundance across conditions based on
# a family of linear mixed-effects models in TMT experiment. Experimental
# design of case-control study (patients are not repeatedly measured) is
# automatically determined based on proper statistical model.	

# test for all the possible pairs of conditions	
# FIXME: parallelize!
all_results <- groupComparisonTMT(filt_prot)	
# Note: for intrafraction comparisons the ANOVA case is utilized. The model is
# just Abundance ~ Condition
#load(file.path(rdatdir,"all_results.rda"))


## clean-up MSstats results ----------------------------------------------------

# collect relevant contrasts
control_groups <- paste("Control",paste0("F",seq(4,10)),sep=".")
mutant_groups <- paste("Mutant",paste0("F",seq(4,10)),sep=".")
contrasts <- paste(control_groups,mutant_groups,sep="-")
filt_results <- all_results %>% filter(Label %in% contrasts)

# map uniprot accession to gene names using twesleyb/getPPIs
uniprot <- unique(filt_results$Protein)
genes <- suppressWarnings({
	getPPIs::getIDs(uniprot,from="uniprot",to="symbol",species="mouse")
})
# FIXME: clean-up/fix getPPIs namespace causing warnings on import
entrez <- suppressWarnings({
	getPPIs::getIDs(uniprot,from="uniprot",to="entrez",species="mouse")
})
gene_map <- data.frame("uniprot" = uniprot, "symbol" = genes, "entrez" = entrez)

# status
not_mapped <- sum(is.na(genes))
warning(paste0("Unable to map UniProt Accession to Gene Symbol for ",
	      not_mapped, " (",
	      round(100*(not_mapped/length(uniprot)),2),"%) ",
	      "proteins."))

# add to the data
idx <- match(filt_results$Protein,gene_map$uniprot)
Symbol <- gene_map$symbol[idx]
Entrez <- gene_map$entrez[idx]
filt_results <- tibble::add_column(filt_results,Entrez,.after="Protein")
filt_results <- tibble::add_column(filt_results,Symbol,.after="Entrez")

# Can we add protein level data to  stats and save as tidy object?
# format should be swip_msstats
swip_msstats <- left_join(filt_results,msstats_prot,by="Protein")

# split into BioFraction's for ease of comprehension
# extract fraction from label
fraction <- regmatches(filt_results$Label,
		       regexpr("\\F[0-9]{1,2}",filt_results$Label))
filt_results <- tibble::add_column(filt_results,
				   "Fraction" = fraction,
				   .before="Protein")
results_list <- split(filt_results,fraction)

#length(results_list)
# [1] 7
#names(results_list)
# [1] "F10" "F4"  "F5"  "F6"  "F7"  "F8"  "F9"

# sort
results_list <- results_list[c("F4","F5","F6","F7","F8","F9","F10")]

# sort each df by pvalue and drop rows with issues or NA pvals
clean_results <- function(x) { # x = results_list[[1]]
	is_issue <- !is.na(x$issue)
	x <- x[!is_issue,] # drop rows with issue
	x$issue <- NULL
	x <- x[order(x$pvalue),]
	return(x)
}
final_results <- lapply(results_list, clean_results)

# save as excel worksheet
myfile <- file.path(root,"tables","SWIP_MSstatsTMT_Results.xlsx")
write_excel(final_results, myfile)
message(paste("\nSaved",myfile))

# save as rda
myfile <- file.path(root,"data","swip_msstats.rda")
save(swip_msstats, file=myfile,version=2)
message(paste("\nSaved",myfile))

# any protein overlap?
sigprot_list <- lapply(results_list,function(x) {
	      x %>% filter(adj.pvalue<0.05) %>% 
		      dplyr::select(Protein) %>% 
		      unlist() %>% as.character() %>% unique() })

top_sigprots <- Reduce(intersect, sigprot_list)
idx <- match(top_sigprots, gene_map$uniprot)
message(c("\nProteins that are differentially abundant in all fractions: ",
	      paste(gene_map$symbol[idx],collapse=", ")))

# Examine the number of significant proteins
message(paste("\nNumber of differentially abundant (FDR < 0.05) proteins:"))
df <- sapply(results_list,function(x) sum(x$adj.pvalue<0.05,na.rm=TRUE))
knitr::kable(t(df))

## FIXME: flip sign of log2FC
## add percent WT
## Change Label to Contrast or Comparison
## add protein level data?
