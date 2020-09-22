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
			 condition = interaction(samples$Fraction,
						 samples$Treatment),
			 replicate = interaction(samples$Treatment,
						 samples$Experiment),
			 # you can include additional covariates:
			 experiment = samples$Experiment,
			 fraction = samples$Fraction,
			 treatment = samples$Treatment
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
prot_norm <- normalize_vsn(prot_se)

# Protein-level analysis of differential abundance ----------------------------
# Differential enrichment analysis  based on linear models and empherical Bayes
# statistics

# Test manually defined comparisons
all_contrasts <- sapply(unique(exp_design$fraction),
		    paste,unique(exp_design$treatment),simplify=FALSE,sep=".")

# loop to perform tests for every intra-fraction comparison
results <- lapply(all_contrasts, function(comparison) {
			  test_diff(prot_se, type = "manual",
				    test = comparison) })
		    

# Denote significant proteins based on user defined cutoffs
dep_results <- lapply(results,add_rejections, alpha = 0.05)

# Plot the first and second principal components
# FIXME: WARNING Message
# Use of `pca_df[[indicate[2]]]` is discouraged. 
# Use `.data[[indicate[2]]]` instead.
#plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

# Plot a frequency plot of significant proteins for the different conditions
#plot_cond(dep)

## results_table -----------------------------------------------------------

# Generate a results table
data_results <- lapply(dep_results,get_results)

# Save as excel
write_excel(data_results,file.path(root,tables,"DEP_results.xlsx"))

# Significant proteins
knitr::kable(sapply(data_results,function(x) sum(x$significant)))

quit()

#####################################################################

data_diff <- test_diff(prot_se, type = "manual",
                              test = c("Control.F5", "Mutant.F5"))
