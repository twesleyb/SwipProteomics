#!/usr/bin/env Rscript

# title:
# description:
# author: twab <twesleyb10@gmail.com>
# os: windows linux subsystem (WSL)

## INPUT ----------------------------------------------------------------------
# specify projject's root directory
ROOT <- "~/projects/SwipProteomics"

## OPTIONS --------------------------------------------------------------------
FDR_alpha <- 0.05

## FUNCTIONS -----------------------------------------------------------------

fix_colname <- function(df,old_colname,new_colname) {
	# change a column's name in a data.frame
	colnames(df)[which(colnames(df) == old_colname)] <- new_colname
	return(df)
}

## Prepare the working environment ----------------------------------------------

# directory for output tables:
tabsdir <- file.path(ROOT,"tables")
if (!dir.exists(tabsdir)) { 
	dir.create(tabsdir)
	message(paste("created",figsdir))
}

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
data(swip_tmt)
data(samples)
data(gene_map)


## fix gene_map ----------------------------------------------------------------

# fix missing gene name in gene_map
idx <- which(is.na(gene_map$symbol))
if (gene_map$uniprot[idx] == "Q80WG5") {
	gene_map$symbol[idx] <- "Lrrc8a"
}


## prepare data for DEP --------------------------------------------------------

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


## prepare exp_design for DEP --------------------------------------------------
# build experiment design data.frame from SWIP TMT sample data
# NOTE: exp_design should contain colums for 'label', 'condition' 
# and 'replicate' information
# NOTE: 'condition' is important, it must be the contrast of interest
# NOTE: DEP currently does not support more complicated experimental designs.
# Formula Model:
# 	>>> 	~ 0 + condition (Fraction.Genotype)

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
suppressMessages({ prot_norm <- normalize_vsn(prot_se) })


## Protein-level analysis of differential abundance ----------------------------
# Differential enrichment analysis  based on linear models and empherical Bayes
# statistics

# define all contrasts of intrafraction groups
groups <- sapply(unique(exp_design$fraction),
		    paste,unique(exp_design$treatment),simplify=FALSE,sep=".")
all_contrasts <- sapply(groups, paste, collapse="_vs_")

# loop to perform tests for every intra-fraction comparison
results <- lapply(all_contrasts, function(comparison) {
			  DEP::test_diff(prot_norm, type = "manual",
				    test = comparison) })

# Denote significant proteins based on user defined cutoffs
# FIXME: what is p.adjust method?
# from limma::topTable() adjust.method = BH 'Benjamini Hochberg'
dep_results <- lapply(results, DEP::add_rejections, alpha = FDR_alpha)


## collect results_table -----------------------------------------------------------

# Generate a results table
data_results <- lapply(dep_results,get_results)

# Summarize Significant proteins
message(paste0("\nSummary of significant proteins (FDR < ", FDR_alpha,") ",
	      "for intra-fraction comparisons:"))
df <- data.frame(fraction = names(data_results),
		 n_sig =sapply(data_results, function(x) sum(x$significant)),
		 sig_prots=sapply(data_results,function(x) {
					  paste(x$name[x$significant],collapse=", ")}))
knitr::kable(df,row.names=FALSE)


# Clean-up results ------------------------------------------------------------

# change 'name' and 'ID' columns to 'Symbol' and 'Accession'
data_results <- lapply(data_results, fix_colname,"name","Symbol")
data_results <- lapply(data_results, fix_colname,"ID","Accession")

# FIXME: centered calculation is not working, all NA
drop_cols <- function(df) { return(df[,!grepl("_centered",colnames(df))]) }
data_results <- lapply(data_results,drop_cols)


## save any output ------------------------------------------------------------

# Save DEP results as excel workbook
myfile <- file.path(tabsdir,"DEP_results.xlsx")
write_excel(data_results,myfile)
message(paste("\nSaved",myfile))

# DONE
