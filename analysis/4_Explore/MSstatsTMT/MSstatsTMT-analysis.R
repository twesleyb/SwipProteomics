#!/usr/bin/env Rscript

# title: SwipProteomics
# description: analysis of Swip TMT spatial proteomics data with MSstatsTMT
# author: Tyler W Bradshaw <twesleyb10@gmail.com>
# os: windows linux subsystem (WSL)

## Input ----------------------------------------------------------------------

# specity the projects root directory:
root <- "~/projects/SwipProteomics"

# input data in root/data:
input_dir = "data/PSM.zip"

# PSM.zip contains:
input_psm = "PSM/5359_PSM_Report.xlsx"
input_samples = "PSM/5359_Sample_Report.xlsx"

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
renv::load(root,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(dplyr)
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
message("\nLoading raw PSM data.") 
myfile <- file.path(downdir,input_psm)
raw_pd <- readxl::read_excel(myfile,progress=FALSE)

unlink(myfile) 
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


## map Spectrum.Files to MS.Channels ------------------------------------------
# Basically, its just complicated. There were 3x 16-plex TMT experiments.
# In each, the concatenated TMT mixture was fractionated into 12 fractions in 
# order to increase analytical depth. Therefore, there were 12 * 3 = 36 mass
# spectrometry runs. Each Run is recorded as a 'Spectrum File' by Proteome
# Discover. Note that each MS run cooresponds to the measurment of all 16 
# samples and therefore the total number of TMT channels for which we 
# have made measurements is 16 x 3 Experiments = 48. In other words, a single
# 'Spectrum.File' cooresponds to multiple Samples. 

# collect all Spectrum.Files grouped by Experiment
# split 'Spectrum.File' at first "_" to extract experiment identifiers
all_files <- raw_pd$Spectrum.File
exp_files <- lapply(split(all_files, sapply(strsplit(all_files,"_"),"[",1)),
		    unique)
files_dt <- data.table("Experiment ID" = rep(names(exp_files),
					     times=sapply(exp_files,length)),
	               "Run" = unlist(exp_files))

# collect all MS.Channels, grouped by Experiment
all_channels <- samples$MS.Channel
exp_channels <- split(all_channels, sapply(strsplit(all_channels,"_"),"[",1))
exp_dt <- data.table("Experiment ID" = rep(names(exp_channels),
					times=sapply(exp_channels,length)),
	             "MS.Channel" = unlist(exp_channels))

# use exp_channels to create numeric ID for MS Fraction
# NOTE: this is MSstats 'Fraction'
exp_fraction_list <- lapply(exp_channels, function(x) {
	       setNames(as.numeric(as.factor(x)),x)
	   })
x <- unlist(exp_fraction_list)
named_fractions <- setNames(x,sapply(strsplit(names(x),"\\."),"[",2))


# build annotation file for MSstats -------------------------------------------
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
annotation_dt <- left_join(files_dt,exp_dt,by="Experiment ID")

# add additional freatures from samples
idx <- match(annotation_dt$"MS.Channel", samples$MS.Channel)
annotation_dt$Fraction <- named_fractions[annotation_dt$MS.Channel]
annotation_dt$TechRepMixture <- rep(1,length(idx))
annotation_dt$Mixture <- samples$Mixture[idx]
annotation_dt$Condition <- samples$Condition[idx]
annotation_dt$Channel <- samples$Channel[idx]
annotation_dt$BioFraction <- samples$BioFraction[idx]
annotation_dt$BioReplicate <- samples$BioReplicate[idx]
annotation_dt$BioReplicate[grepl("Norm",annotation_dt$BioReplicate)] <- "Norm"

# FIXME: how to pass additional covariates to MSstats?
# FIXME: how to handle repeated measures design? 

# Remove un-needed cols
annotation_dt$"Experiment ID" <- NULL
annotation_dt$"MS.Channel" <- NULL

# check the annotation file
# this basic design is repeated for each experiment:
knitr::kable(annotation_dt[c(1:16),]) # 3x

# save to file
myfile <- file.path(rdatdir,"annotation_dt.rda")
save(annotation_dt,file=myfile,version=2)
message(paste("\nSaved",myfile))


## combine annotation_dt and raw_pd to coerce data to MSstats format -----------

# NOTE: this takes a considerable amount of time
message("\nConverting PSM data to MSstatsTMT format.")
data_pd <- PDtoMSstatsTMTFormat(raw_pd, annotation_dt)
#load(file.path(rdatdir,"data_pd.rda"))

# save to file
myfile <- file.path(rdatdir,"data_pd.rda")
save(data_pd,file=myfile,version=2)
message(paste("\nSaved",myfile))


# Protein summarize and normalization -----------------------------------------
# use MSstats for protein summarization	
# Sample Summary does not look correct, should be 

# NOTE: this takes a considerable amount of time
# FIXME: Joining, by = ("Run", "Channel") # unexpected output
# FIXME: remove extremely long message about normalization between runs
message("\nPerforming normalization and protein sumamrization.")
data_prot <- proteinSummarization(data_pd,
				  method="msstats",	
                                  global_norm=TRUE,	
                                  reference_norm=TRUE,	
                                  remove_norm_channel = TRUE)
#load(file.path(rdatdir,"data_prot.rda"))

# save to file
myfile <- file.path(rdatdir,"data_prot.rda")
save(data_prot,file=myfile,version=2)
message(paste("\nSaved",myfile))
	

## Protein-level statistical testing -------------------------------------------
# Tests for significant changes in protein abundance across conditions based on
# a family of linear mixed-effects models in TMT experiment. Experimental
# design of case-control study (patients are not repeatedly measured) is
# automatically determined based on proper statistical model.	

# test for all the possible pairs of conditions	
all_results <- groupComparisonTMT(data_prot)	
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

# save as excel worksheet
myfile <- file.path(root,"tables","SWIP_MSstatsTMT_Results.xlsx")
write_excel(results_list, myfile)
message(paste("\nSaved",myfile))

# save as rda
myfile <- file.path(root,"rdata","results_list.rda")
save(results_list, file=myfile,version=2)
message(paste("\nSaved",myfile))

# Examine the number of significant proteins
message(paste("\nNumber of differentially abundant (FDR < 0.05) proteins:"))
df <- sapply(results_list,function(x) sum(x$adj.pvalue<0.05,na.rm=TRUE))
knitr::kable(t(df))

quit()

## compare pvalues to DEP and EdgeR pipelines ----------------------------------

# 1. combine EdgeR and DEP stats
# set global plotting settings
ggtheme(); set_font("Arial", font_path = fontdir)

# load the data in root/data
data(swip_tmt) # the preprocessed data (contains edgeR stats)
data(swip_dep) # the DEP stats

# collect edgeR stats for intra-fraction comparisons
tmpdf <-  swip_tmt %>% 
	group_by(Fraction) %>% 
	dplyr::select(Fraction,Accession,PValue) %>% 
	unique()
# combine edgeR and DEP stats
stats_df <- swip_dep %>% 
	dplyr::select(Fraction,Accession,p.val) %>% 
	left_join(tmpdf,by=c("Fraction","Accession"))
# clean-up colnames
stats_df <- fix_colname(stats_df,"PValue","EdgeR")
stats_df <- fix_colname(stats_df,"p.val","DEP")

# 2. add MSstatsTMT stats
temp_list <- lapply(results_list, function(x) {
			    x %>% dplyr::select(Protein,pvalue) })
df <- bind_rows(temp_list,.id="Fraction")
# fix colnames colnames(df) <- c("Fraction","Accession","MSstatsTMT")

# combine with other stats
all_stats <- left_join(df,stats_df,by=c("Accession","Fraction"))

# drop rows in which pvals from all methods is NA
idx <- apply(all_stats, 1, function(x) all(is.na(x[c(3,4,5)])))
filt_stats <- all_stats[!idx,]

# the spearman rank correlation is very high
rho1 <- cor(x=filt_stats$MSstatsTMT, y=filt_stats$DEP, method="spearman", use='pairwise.complete.obs')
rho2 <- cor(x=filt_stats$MSstatsTMT, y=filt_stats$EdgeR, method="spearman", use='pairwise.complete.obs')
rho3 <- cor(x=filt_stats$DEP, y=filt_stats$EdgeR, method="spearman", use='pairwise.complete.obs')
df <- data.frame(rho1,rho2,rho3)
knitr::kable(df)

## ----------------------------------------------------------------------------

# generate a plot examining coorelation between edgeR and DEP pvalues
plot <- ggplot(data=stats_df,aes(x=EdgeR,y=DEP))
plot <- plot + geom_point()
plot <- plot + xlab("P-Value (edgeR)")
plot <- plot + ylab("P-Value (DEP)")
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(panel.border=element_rect(colour="black",fill=NA,size=1))

# save as pdf
myfile <- file.path(figsdir,"PValue_correlation_scatterplot.pdf")
ggsave(myfile,plot,height=5,width=5)
message(paste("\nSaved",myfile))

## ----------------------------------------------------------------------------

# tidy the data
df <- reshape2::melt(stats_df,id=c("Fraction","Accession"),
		     value.name = "PValue",
		     variable.name="Method")

# loop generate p-value histograms for every intra-fraction comparison
plots <- list()
for (fraction in unique(df$Fraction)) {
	plot <- ggplot(df %>% filter(Fraction == fraction),
		       aes(x=PValue,color=Method))
	plot <- plot + geom_histogram(bins=100)
	plot <- plot + theme(panel.background = element_blank())
	plot <- plot + theme(panel.border=element_rect(colour="black",
						       fill=NA,size=1))
	plot <- plot + ggtitle(paste("Fraction:",fraction))
	plots[[fraction]] <- plot
}

# save as a single pdf
myfile <- file.path(figsdir,"PValue_Histograms.pdf")
ggsavePDF(plots,myfile)
message(paste("\nSaved",myfile))

if (file.exists("Rplots.pdf")) { unlink("Rplots.pdf") }

# DONE!
