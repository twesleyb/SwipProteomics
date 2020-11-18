#!/usr/bin/env Rscript

# title: SwipProteomics
# description: 
# author: twab <twesleyb10@gmail.com>
# os: windows linux subsystem (WSL)

## INPUT ----------------------------------------------------------------------
# specify project's root directory
#root <- "~/projects/SwipProteomics"
root <- "~/projects/SynaptopathyProteomics" # need DEP... cant get it installed

## OPTIONS --------------------------------------------------------------------
##FDR_alpha <- 0.05 # BH significance threshold for protein DA

## FUNCTIONS -----------------------------------------------------------------

fix_colname <- function(df,old_colname,new_colname) {
	# change a column's name in a data.frame
	colnames(df)[which(colnames(df) == old_colname)] <- new_colname
	return(df)
}


mkdir <- function(...) {
	# create a new directory
	newdir <- file.path(...)
	if (!dir.exists(newdir)) { 
		dir.create(newdir)
		message(paste("created",newdir))
	}
}


## Prepare the working environment ---------------------------------------------

# directory for output tables:
tabsdir <- file.path(root,"tables"); mkdir(tabsdir)

# directory for output data:
datadir <- file.path(root,"data"); mkdir(datadir)

# directory for temporary data:
rdatdir <- file.path(root,"rdata"); mkdir(rdatdir)

# directory for figures:
figsdir <- file.path(root,"figs","DEP"); mkdir(figsdir)


## Prepare the R environment ---------------------------------------------------

# load renv
renv::load(root,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(DEP) # twesleyb/DEP 
	library(dplyr) # for manipulating data
	library(data.table) # for working with data.tables
})

# load functions in root/R
proot <- "~/projects/SwipProteomics"
devtools::load_all(project)

# load the data in root/data
data(gene_map)
data(msstats_prot)

# load the raw data
myfile <- file.path(proot,"rdata","msstats_psm.rda")
load(myfile) # msstats_psm


## prepare raw protein data for DEP -------------------------------------------
# NOTE: do not log transform the data

# cast tidy data into a data.table
dm <- msstats_prot %>% filter(Abundance > 0) %>%
	group_by(Protein,Mixture,Channel,Condition) %>% 
	reshape2::dcast(Protein ~ Mixture + Channel + Condition, 
			value.var = "Abundance") %>%
	as.data.table() %>%  na.omit() %>%
	as.matrix(rownames="Protein")
norm_df <- as.data.table(2^dm,keep.rownames="ID") # coerce back to dt with "ID" column


raw_df <- msstats_psm %>%  mutate(Condition = as.character(Condition)) %>%
	filter(Condition != "Norm") %>%  # drop QC
	group_by(ProteinName,Mixture,Channel,Condition) %>% 
	summarize("Abundance" = sum(Intensity),.groups="drop") %>%
	reshape2::dcast(ProteinName ~ Mixture + Channel + Condition, 
			value.var = "Abundance") %>% 
	as.data.table() %>% 
	as.matrix(rownames="ProteinName") %>%
	as.data.table(keep.rownames="ID") # coerce back to dt with "ID" column

# protein data.frame should contain unique protein IDs in 'ID' column
# check for duplicates
#stopifnot(!any(duplicated(prot_df$ID))) # there should be no duplicate IDs

# check for NA
#stopifnot(!any(is.na(prot_df$ID))) # there should be no NA

# protein data.frame should contain unique protein names in 'name' column
# we will annotate proteins with gene 'symbol's
raw_df$name <- gene_map$symbol[match(raw_df$ID,gene_map$uniprot)]
norm_df$name <- gene_map$symbol[match(norm_df$ID,gene_map$uniprot)]


# Multiple Uniprot Accession IDs may be mapped to the same gene Symbol.
# If any duplicated, make names unique base::make.unique.
raw_df$name <- make.unique(raw_df$name,sep="-")
norm_df$name <- make.unique(norm_df$name,sep="-")

# check for duplicates
stopifnot(!any(duplicated(norm_df$name))) # there should be no duplicate names


## create exp_design -------------------------------------------

# create experimental design
namen <- colnames(norm_df)[grep("M",colnames(norm_df))]
exp_design <- data.table(label = namen,
			 Mixture = sapply(strsplit(namen,"_"),"[",1),
			 Channel = sapply(strsplit(namen,"_"),"[",2),
			 condition = sapply(strsplit(namen,"_"),"[",3))
exp_design$treatment <- sapply(strsplit(exp_design$condition,"\\."),"[",1)
exp_design$BioFraction <- sapply(strsplit(exp_design$condition,"\\."),"[",2)
exp_design$replicate <-  interaction(exp_design$Mixture,exp_design$treatment)


# check for required columns in input
stopifnot(all(c("label","condition","replicate") %in% colnames(exp_design)))


## build SummarizedExperiment (se) object -------------------------------------

# use DEP helper function to build a SE object
# specify the column indices containing the numeric data (idy)
idy <- grep("M",colnames(raw_df))
raw_se <- DEP::make_se(raw_df,columns=idy,exp_design)

# THE DATA SHOULD NOT BE LOGGED!
idy <- grep("M",colnames(norm_df))
norm_se <- DEP::make_se(norm_df,columns=idy,exp_design)


## plots ----------------------------------------------------------------------

# extract normalized protein data from se_list
# NOTE: the se in se_list differ only in their rowData

# save as pdf
figsdir <- file.path("~/projects/SwipProteomics","figs")
myfile <- file.path(figsdir,"DEP_Plots.pdf")

# DEP preprocessing plots
pdf(file=myfile,onefile=TRUE)

#print(plot_frequency(raw_se))

#plot_numbers(raw_se)

#plot_coverage(raw_se)

plot_normalization(norm_se)

#plot_missval(raw_se)

plot_detect(raw_se)

print(plot_pca(raw_se))

print(plot_pca(norm_se))

dev.off()

message("\nSaved",myfile)

if (file.exists("Rplots.pdf")) { unlink("Rplots.pdf") }

#ggsavePDF(plots,file=myfile)
