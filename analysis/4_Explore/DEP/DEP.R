#!/usr/bin/env rscript

## functions ------------------------------------------------------------------

getrd <- function(here=getwd(), dpat= ".git") {
	# get the repository's root directory
	in_root <- function(h=here, dir=dpat) {
		check <- any(grepl(dir,list.dirs(h,recursive=false)))
		return(check)
	}
	# loop to find root.
	while (!in_root(here)) {
		here <- dirname(here)
	}
	root <- here
	return(root)
}

## ----required packages, echo = false, warning=false, results="hide"-----------

root <- getrd()
renv::load(root)

suppresspackagestartupmessages({
  library("dep")
  library("dplyr")
  library("biocstyle")
})

## ----install, eval = false----------------------------------------------------
#
#  if (!requirenamespace("biocmanager", quietly=true))
#      install.packages("biocmanager")
#  biocmanager::install("dep")
#
#  library("dep")
#

## ----run_app, eval = false----------------------------------------------------
#  # for lfq analysis
#  run_app("lfq")
#
#  # for tmt analysis
#  run_app("tmt")

## ----load data----------------------------------------------------------------

# the data is provided with the package
data <- ubilength

# we filter for contaminant proteins and decoy database hits, which are indicated by "+" in the columns "potential.contaminants" and "reverse", respectively.
data <- filter(data, reverse != "+", potential.contaminant != "+")

## ----dimension----------------------------------------------------------------
dim(data)

## ----colnames-----------------------------------------------------------------
colnames(data)

## ----unique-------------------------------------------------------------------
# are there any duplicated gene names?
data$gene.names %>% duplicated() %>% any()

# make a table of duplicated gene names
data %>% group_by(gene.names) %>% summarize(frequency = n()) %>%
  arrange(desc(frequency)) %>% filter(frequency > 1)

## ----unique_names-------------------------------------------------------------
# make unique names using the annotation in the "gene.names" column as primary names and the annotation in "protein.ids" as name for those that do not have an gene name.
data_unique <- make_unique(data, "gene.names", "protein.ids", delim = ";")

# are there any duplicated names?
any_duplicates <- data$name %>% duplicated() %>% any()
if (any_duplicates) { stop("there are duplicated gene names!") }

## ----expdesign, echo = false--------------------------------------------------
# display experimental design
knitr::kable(ubilength_expdesign)

## ----to_exprset---------------------------------------------------------------
# generate a summarizedexperiment object using an experimental design
lfq_columns <- grep("lfq.", colnames(data_unique)) # get lfq column numbers
experimental_design <- ubilength_expdesign
data_se <- make_se(data_unique, lfq_columns, experimental_design)

# generate a summarizedexperiment object by parsing condition information from
# the column names
lfq_columns <- grep("lfq.", colnames(data_unique)) # get lfq column numbers
data_se_parsed <- make_se_parse(data_unique, lfq_columns)

# let's have a look at the summarizedexperiment object
data_se


## ----plot_data_nofilt, fig.width = 4, fig.height = 4--------------------------
# plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

## ----filter_missval-----------------------------------------------------------
# filter for proteins that are identified in all replicates of at least one
# condition
data_filt <- filter_missval(data_se, thr = 0)

# Less stringent filtering:
# Filter for proteins that are identified in 2 out of 3 replicates of at least
# one condition
data_filt2 <- filter_missval(data_se, thr = 1)

## ----plot_data, fig.width = 4, fig.height = 4---------------------------------
# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

## ----plot_data2, fig.width = 3, fig.height = 4--------------------------------
# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

## ----normalize----------------------------------------------------------------
# Normalize the data
data_norm <- normalize_vsn(data_filt)

## ----plot_norm, fig.width = 4, fig.height = 5---------------------------------
# Visualize normalization by boxplots for all samples before and after
# normalization
plot_normalization(data_filt, data_norm)

## ----plot_missval, fig.height = 4, fig.width = 3------------------------------
# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

## ----plot_detect, fig.height = 4, fig.width = 4-------------------------------
# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)

## ----impute, results = "hide", message = FALSE, warning = FALSE, error = TRUE----
# All possible imputation methods are printed in an error, if an invalid function name is given.
impute(data_norm, fun = "")

# Impute missing data using random draws from a Gaussian distribution centered
# around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Impute missing data using random draws from a manually defined left-shifted
# Gaussian distribution (for MNAR)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)


## ----plot_imp, fig.width = 4, fig.height = 4----------------------------------
# Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)

## ----statistics---------------------------------------------------------------
# Differential enrichment analysis  based on linear models and empherical Bayes
# statistics

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "Ctrl")

# Test all possible comparisons of samples
data_diff_all_contrasts <- test_diff(data_imp, type = "all")

# Test manually defined comparisons
data_diff_manual <- test_diff(data_imp, type = "manual",
                              test = c("Ubi4_vs_Ctrl", "Ubi6_vs_Ctrl"))

## ----add_reject---------------------------------------------------------------
# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

## ----pca, fig.height = 3, fig.width = 4---------------------------------------
# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

## ----corr, fig.height = 3, fig.width = 5--------------------------------------
# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

## ----heatmap, fig.height = 5, fig.width = 3-----------------------------------
# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE,
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))

## ----heatmap2, fig.height = 5, fig.width = 3----------------------------------
# Plot a heatmap of all significant proteins (rows) and the tested contrasts
# (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE,
             k = 6, col_limit = 10, show_row_names = FALSE)

## ----volcano, fig.height = 5, fig.width = 5-----------------------------------
# Plot a volcano plot for the contrast "Ubi6 vs Ctrl""
plot_volcano(dep, contrast = "Ubi6_vs_Ctrl", label_size = 2, add_names = TRUE)

## ----bar, fig.height = 3.5, fig.width = 3.5-----------------------------------
# Plot a barplot for USP15 and IKBKG
plot_single(dep, proteins = c("USP15", "IKBKG"))

# Plot a barplot for the protein USP15 with the data centered
plot_single(dep, proteins = "USP15", type = "centered")

## ----overlap, fig.height = 4, fig.width = 6-----------------------------------
# Plot a frequency plot of significant proteins for the different conditions
plot_cond(dep)

## ----results_table------------------------------------------------------------
# Generate a results table
data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()

## ----results_table2-----------------------------------------------------------
# Column names of the results table
colnames(data_results)

## ----get_df-------------------------------------------------------------------
# Generate a wide data.frame
df_wide <- get_df_wide(dep)
# Generate a long data.frame
df_long <- get_df_long(dep)

## ----save_load, eval = FALSE--------------------------------------------------
#  # Save analyzed data
#  save(data_se, data_norm, data_imp, data_diff, dep, file = "data.RData")
#  # These data can be loaded in future R sessions using this command
#  load("data.RData")

## ----LFQ, message = FALSE-----------------------------------------------------
# The data is provided with the package
data <- UbiLength
experimental_design <- UbiLength_ExpDesign

# The wrapper function performs the full analysis
data_results <- LFQ(data, experimental_design, fun = "MinProb",
      type = "control", control = "Ctrl", alpha = 0.05, lfc = 1)

## ----report2, eval = FALSE----------------------------------------------------
#  # Make a markdown report and save the results
#  report(data_results)

## ----report1------------------------------------------------------------------
# See all objects saved within the results object
names(data_results)

## ----LFQ_results3-------------------------------------------------------------
# Extract the results table
results_table <- data_results$results

# Number of significant proteins
results_table %>% filter(significant) %>% nrow()

## ----LFQ_results4, fig.height = 5, fig.width = 3------------------------------
# Extract the sign object
full_data <- data_results$dep

# Use the full data to generate a heatmap
plot_heatmap(full_data, type = "contrast", kmeans = TRUE,
             k = 6, col_limit = 4, show_row_names = FALSE)

## ----TMT, eval = FALSE--------------------------------------------------------
#  # Need example data
#  TMTdata <- example_data
#  Exp_Design <- example_Exp_Design
#
#  # The wrapper function performs the full analysis
#  TMTdata_results <- TMT(TMTdata, expdesign = Exp_Design, fun = "MinProb",
#      type = "control", control = "Control", alpha = 0.05, lfc = 1)

## ----session_info, echo = FALSE-----------------------------------------------
sessionInfo()
