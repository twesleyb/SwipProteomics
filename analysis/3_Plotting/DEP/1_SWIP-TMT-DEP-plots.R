#!/usr/bin/env Rscript

# title: SwipProteomics
# description: 
# author: twab <twesleyb10@gmail.com>
# os: windows linux subsystem (WSL)

## INPUT ----------------------------------------------------------------------
# specify project's root directory
ROOT <- "~/projects/SwipProteomics"

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
tabsdir <- file.path(ROOT,"tables"); mkdir(tabsdir)

# directory for output data:
datadir <- file.path(ROOT,"data"); mkdir(datadir)

# directory for temporary data:
rdatdir <- file.path(ROOT,"rdata"); mkdir(rdatdir)

# directory for figures:
figsdir <- file.path(ROOT,"figs","DEP"); mkdir(figsdir)


## Prepare the R environment ---------------------------------------------------

# load renv
renv::load(ROOT,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(DEP) # twesleyb/DEP 
	library(dplyr) # for manipulating data
	library(data.table) # for working with data.tables
})

# load functions in root/R
devtools::load_all(ROOT)

# load the data in root/data
data(msstats_prot) # the preprocessed data
data(gene_map)

#  load the data in root/rdata
myfile <- file.path(rdatdir,"tidy_peptide.rda")
load(myfile) # tidy_peptide


## prepare raw protein data for DEP -------------------------------------------

# summarize raw protein data
tidy_prot <- summarize_prot(tidy_peptide)

# cast tidy data into a data.table
# NOTE: do not log transform the data
prot_df <- tidy_prot %>% as.data.table() %>%
	dcast(Accession ~ Sample,value.var = "Intensity") %>% 
	as.matrix(rownames="Accession") %>% # coerce to matrix
	as.data.table(keep.rownames="ID") # coerce back to dt with "ID" column

# protein data.frame should contain unique protein IDs in 'ID' column
# check for duplicates
stopifnot(!any(duplicated(prot_df$ID))) # there should be no duplicate IDs

# check for NA
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


## create exp_design -------------------------------------------

exp_design <- data.frame(
			 label = samples$Sample,
			 condition = interaction(samples$Fraction,
						 samples$Treatment),
			 replicate = interaction(samples$Treatment,
						 samples$Experiment),
			 experiment = samples$Experiment,
			 fraction = samples$Fraction,
			 treatment = samples$Treatment
			 )

# check for required columns in input
stopifnot(all(c("label","condition","replicate") %in% colnames(exp_design)))


## build SummarizedExperiment (se) object -------------------------------------

# use DEP helper function to build a SE object
# specify the column indices containing the numeric data (idy)
idy <- grep("Abundance",colnames(prot_df))
raw_se <- DEP::make_se(prot_df,columns=idy,exp_design)


## plots ----------------------------------------------------------------------

# extract normalized protein data from se_list
# NOTE: the se in se_list differ only in their rowData
norm_se <- se_list[[1]]

# save as pdf
myfile <- file.path(figsdir,"DEP_Plots.pdf")

# DEP preprocessing plots
pdf(file=myfile,onefile=TRUE)
print(plot_frequency(raw_se))
print(plot_numbers(raw_se))
print(plot_coverage(raw_se))
print(plot_normalization(raw_se,norm_se))
print(plot_missval(raw_se))
print(plot_detect(raw_se))
print(plot_pca(raw_se))
print(plot_pca(norm_se))
dev.off()

message("\nSaved",myfile)

if (file.exists("Rplots.pdf")) { unlink("Rplots.pdf") }

#ggsavePDF(plots,file=myfile)
