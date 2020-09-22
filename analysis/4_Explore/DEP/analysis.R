#!/usr/bin/env Rscript

# title: SwipProteomics
# description: protein-level statiscal analysis with DEP
# author: Tyler W Bradshaw <twesleyb10@gmail.com>

# Prepare the R environment ---------------------------------------------------

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(DEP)
	library(dplyr)
	library(data.table)
})

# load functions in root/R and data in root/data
devtools::load_all()

# load the preprocessed protein data as a tidy data.table
data(swip_tmt)
data(samples)
data(gene_map)

# fix missing gene name in gene_map
idx <- which(is.na(gene_map$symbol))
if (gene_map$uniprot[idx] == "Q80WG5") {
	gene_map$symbol[idx] <- "Lrrc8a"
}

# load the DEP pipeline
#load("UbiLength_ExpDesign.rda")
#load("data_unique.rda")
#load("data_se.rda")

# cast tidy data into a data.table
# NOTE: do not log transform the data
prot_df <- swip_tmt %>% 
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

# multiple Uniprot Accession IDs may be mapped to the same gene Symbol.
# if any duplicated, make names unique by annotating with #
prot_df$name <- make.unique(prot_df$name,sep="-")

# check for duplicates
stopifnot(!any(duplicated(prot_df$name))) # there should be no duplicate names

# build experiment design data.frame from SWIP TMT sample data
# NOTE: exp_design should contain colums for 'label', 'condition' 
# and 'replicate' information
exp_design <- data.frame(
			 label = samples$Sample,
			 condition = interaction(samples$Treatment,
						 samples$Fraction),
			 replicate = interaction(samples$Treatment,
						 samples$Experiment)
			 )

# check for required columns in input
stopifnot(all(c("label","condition","replicate") %in% colnames(exp_design)))

# use DEP helper function to build a SE object
# specify the column indices containing the numeric data (idy)
idy <- grep("Abundance",colnames(prot_df))
prot_se <- DEP::make_se(prot_df,columns=idy,exp_design)

#plot_frequency(prot_se)
#prot_filt <- filter_missval(prot_se, thr = 0) # ERROR
#Error in Summary.factor(c(1L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L,  : 
#  ‘max’ not meaningful for factors
#Calls: filter_missval -> Summary.factor
#Execution halted

# Normalize the data
# FIXME: need this be done? does it hurt if it doesnt?
#prot_norm <- normalize_vsn(prot_se)
save(prot_norm,file="prot_norm.rda",version=2)

# Protein-level analysis of differential abundance ----------------------------
# Differential enrichment analysis  based on linear models and empherical Bayes
# statistics

# Test every sample versus control
#data_diff <- test_diff(prot_norm, type = "control", control = "Control")

# Test all possible comparisons of samples
#data_diff_all_contrasts <- test_diff(prot_norm, type = "all")

# Test manually defined comparisons
data_diff <- test_diff(prot_se, type = "manual",
                              test = c("Control.F5", "Mutant.F5"))

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

# Plot the first and second principal components
# FIXME: WARNING Message
# Use of `pca_df[[indicate[2]]]` is discouraged. 
# Use `.data[[indicate[2]]]` instead.
#plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

# Plot a frequency plot of significant proteins for the different conditions
#plot_cond(dep)

## results_table -----------------------------------------------------------

# Generate a results table
data_results <- get_results(dep)

colnames(data_results)

sum(data_results$significant)

quit()

#####################################################################




data_diff <- test_diff(prot_se, type = "manual",
                              test = c("Control.F5", "Mutant.F5"))
