#!/usr/bin/env Rscript

# title:
# description:
# author: twab <twesleyb10@gmail.com>
# os: windows linux subsystem (WSL)

## Prepare the working environment ----------------------------------------------

# project root directory:
root <- "~/projects/DEP"

# directory for output tables:
tabsdir <- file.path(root,"tables")
if (!dir.exists(tabsdir)) { 
	dir.create(tabsdir)
	message(paste("mkdir",figsdir))
}

## Prepare the R environment ---------------------------------------------------

# load renv
renv::load(root,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	# dont load DEP, use devtools::load_all()
	library(dplyr)
	library(data.table)
	library(TBmiscR) # for write_excel
})

# load functions in root/R and data in root/data
suppressPackageStartupMessages({
	devtools::load_all()
})

# load the data
data(swip_tmt)
data(samples)
data(gene_map)

# fix missing gene name in gene_map
idx <- which(is.na(gene_map$symbol))
if (gene_map$uniprot[idx] == "Q80WG5") {
	gene_map$symbol[idx] <- "Lrrc8a"
}

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
# if any duplicated, make names unique by annotating with # using make.unique
prot_df$name <- make.unique(prot_df$name,sep="-")

# check for duplicates
stopifnot(!any(duplicated(prot_df$name))) # there should be no duplicate names

# build experiment design data.frame from SWIP TMT sample data
# NOTE: exp_design should contain colums for 'label', 'condition' 
# and 'replicate' information
# NOTE: condition is important, it must be the contrast of interest
exp_design <- data.frame(
			 label = samples$Sample,
			 condition = samples$Treatment,
			 replicate = interaction(samples$Treatment,
						 samples$Experiment,
						 samples$Fraction), 
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

# Normalize the data
# FIXME: need this be done? does it hurt if it doesnt?
norm_prot <- DEP::normalize_vsn(prot_se)
save(norm_prot,file=file.path(root,"data","norm_prot.rda"),version=2)

quit()

# Protein-level analysis of differential abundance ----------------------------
# Differential enrichment analysis  based on linear models and empherical Bayes
# statistics

# FIXME: how to specify more complex design?
#dep <- DEP::test_diff(prot_norm, type = "all",
#		      design_formula=formula("~ fraction + condition"))

#limma::lmFit(dm, design, block="fraction")
		    
# args for testing
se = prot_norm
type = "control"
control = "Control"
design_formula=formula("~ 0 + fraction + condition")
quiet=TRUE

# Denote significant proteins based on user defined cutoffs
# FIXME: what is p.adjust method?
# from limma::topTable() adjust.method = BH 'Benjamini Hochberg'
dep_results <- lapply(results, DEP::add_rejections, alpha = 0.05)

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
TBmiscR::write_excel(data_results,file.path(tabsdir,"DEP_results.xlsx"))

# Summarize Significant proteins
knitr::kable(sapply(data_results,function(x) sum(x$significant)))
