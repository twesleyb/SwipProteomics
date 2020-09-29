#!/usr/bin/env Rscript

# title: SwipProteomics
# description: analysis of intra-fraction differential abundance with DEP
# author: twab <twesleyb10@gmail.com>
# os: windows linux subsystem (WSL)

## INPUT ----------------------------------------------------------------------
# specify project's root directory
ROOT <- "~/projects/SwipProteomics"

## OPTIONS --------------------------------------------------------------------
FDR_alpha <- 0.05 # BH significance threshold for protein DA

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


## Prepare the working environment ----------------------------------------------

# directory for output tables:
tabsdir <- file.path(ROOT,"tables"); mkdir(tabsdir)
datadir <- file.path(ROOT,"data"); mkdir(datadir)


## Prepare the R environment --------------------------------------------------

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
#data(gene_map)

# Set Fraction levels
levels(samples$Fraction) <- c("F4","F5","F6","F7","F8","F9","F10")
samples %>% select(Sample,Experiment,Channel,Treatment,Fraction) %>% 
	arrange(Experiment,Treatment) %>% knitr::kable()

# load the MSstats processed data
myfile <- file.path(ROOT,"rdata","data_prot.rda")
load(myfile) # data_prot

# pass MSstats protein data to DEP --> visualization and imputing


## create gene map ------------------------------------------------------------

uniprot <- unique(data_prot$Protein)
symbols <- getPPIs::getIDs(uniprot,from="uniprot",to="symbol",species="mouse")
gene_map <- data.frame(uniprot,symbols)

# FIXME: Need to go upstream in workflow and map uniprot to genes/drop bad prots
# For now, drop any missing gene symbols
is_missing <- gene_map$uniprot[is.na(gene_map$symbol)]
data_prot <- data_prot %>% filter(Protein %notin% is_missing)

## prepare data for DEP -------------------------------------------------------

# cast tidy data into a data.table
# NOTE: do not log transform the data
#prot_df <- swip_tmt %>% 
#	dcast(Accession ~ Sample,value.var = "Intensity") %>% 
#	as.matrix(rownames="Accession") %>% # coerce to matrix
#	as.data.table(keep.rownames="ID") # coerce back to dt with "ID" column
prot_df <- data_prot %>% as.data.table() %>% 
	dcast(Protein ~ Mixture + Channel + Condition, value.var = "Abundance") %>% 
	as.matrix(rownames="Protein") %>% # coerce to matrix
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


## prepare exp_design for DEP --------------------------------------------------
# build experiment design data.frame from SWIP TMT sample data
# NOTE: exp_design should contain colums for 'label', 'condition' 
# and 'replicate' information
# NOTE: 'condition' is important, it must be the contrast of interest
# NOTE: DEP currently does not support more complicated experimental designs.
# Formula Model:
# 	>>> 	~ 0 + condition (Fraction.Genotype)

samples$label <- as.character(interaction(paste(gsub("Exp","M",samples$Experiment),samples$Channel,samples$Treatment,sep="_"),samples$Fraction))

exp_design <- data.frame(
			 label = samples$label,
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
idy <- grep("M[1,2,3]",colnames(prot_df))
prot_se <- DEP::make_se(prot_df,columns=idy,exp_design)

# NOTE: you can access the data contained in a sumamrized experiment object with 
# functions from the SummarizedExperiment package, e.g.:
# library(SummarizedExperiment)
# rowDat(se)
# colData(se)
# assay(se)


## impute missing values -------------------------------------------------------

# Data can be missing at random (MAR) or missing not at random (MNAR).
# MAR means that values are randomly missing from all samples.
# In the case of MNAR, values are missing in specific samples and/or for specific proteins.
# For example, certain proteins might not be quantified in specific conditions, 
# because they are below the detection limit in these specific samples.  
# 
# To mimick these two types of missing values, 
# we introduce missing values randomly over all data points (MAR)
# and we introduce missing values in the control samples of
# 100 differentially expressed proteins (MNAR).
# The variables used to introduce missing values are depicted below.

# Filter proteins based on missing values

#A first consideration with missing values is whether or not to filter out
#proteins with too many missing values.

## Visualize the extend of missing values
#The number of proteins quantified over the samples can be visualized 
#to investigate the extend of missing values. 

# Plot a barplot of the protein quantification overlap between samples
plot_frequency(se)

# Many proteins are quantified in all six samples and 
# only a small subset of proteins were detected in less than half of the samples.

## Filter options

# We can choose to not filter out any proteins at all,
# filter for only the proteins without missing values,
# filter for proteins with a certain fraction of quantified samples, and
# for proteins that are quantified in all replicates of at least one condition.

no_filter <- se

# Filter for proteins that are quantified in all replicates of at least one condition
condition_filter <- filter_proteins(se, "condition", thr = 0)

# Filter for proteins that have no missing values
complete_cases <- filter_proteins(se, "complete")

# Filter for proteins that are quantified in at least 2/3 of the samples.
frac_filtered <- filter_proteins(se, "fraction", min = 0.66)

# To check the consequences of filtering, we calculate the number of background 
# and DE proteins left in the dataset.

# Mean versus Sd plot
meanSdPlot(no_filter)

# Data imputation of missing data

# A second important consideration with missing values is
# whether or not to impute the missing values.

# MAR and MNAR (see [Introduce missing values](#introduce-missing-values)) 
# require different imputation methods.
# See the `r Biocpkg("MSnbase") ` vignette and more specifically 
# the _impute_ function descriptions for detailed information.  
 
# ## Explore the pattern of missing values
 
# To explore the pattern of missing values in the data, 
# a heatmap can be plotted indicating whether values are missing (0) or not (1).
# Only proteins with at least one missing value are visualized.

# Plot a heatmap of proteins with missing values
plot_missval(no_filter)

# The missing values seem to be randomly distributed across the samples (MAR).
# However, we do note a block of values that are missing in all control samples
# (bottom left side of the heatmap).
# These proteins might have missing values not at random (MNAR).  
 
# To check whether missing values are biased to lower intense proteins, 
# the densities and cumulative fractions are plotted for proteins with 
# and without missing values.   

# Plot intensity distributions and cumulative fraction of proteins 
# with and without missing values
plot_detect(no_filter)

# In our example data, there is no clear difference between the two distributions.  

## Imputation options

# DEP borrows the imputation functions from `r Biocpkg("MSnbase") `.
# See the `r Biocpkg("MSnbase") ` vignette and more specifically the _impute_ 
# function description for more information on the imputation methods.  

# All possible imputation methods are printed in an error, if an invalid function name is given.
impute(no_filter, fun = "")


# No imputation
no_imputation <- no_filter

# Impute missing data using random draws from a 
# Gaussian distribution centered around a minimal value (for MNAR)
MinProb_imputation <- impute(no_filter, fun = "MinProb", q = 0.01)

# Impute missing data using random draws from a 
# manually defined left-shifted Gaussian distribution (for MNAR)
manual_imputation <- impute(no_filter, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
knn_imputation <- impute(no_filter, fun = "knn", rowmax = 0.9)

# The effect of data imputation on the distributions can be visualized.

# Plot intensity distributions before and after imputation
plot_imputation(no_filter, MinProb_imputation, 
  manual_imputation, knn_imputation)
