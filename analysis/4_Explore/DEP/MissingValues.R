## ----required packages, echo = FALSE, warning=FALSE, results="hide"-----------
suppressPackageStartupMessages({
  library("BiocStyle")
  library("DEP")
  library("dplyr")
  library("tidyr")
  library("purrr")
  library("ggplot2")
  library("SummarizedExperiment")
})

## ----variables----------------------------------------------------------------
# Variables
replicates = 3
bg_proteins = 3000
DE_proteins = 300
log2_mean_bg = 27
log2_sd_bg = 2
log2_mean_DE_control = 25
log2_mean_DE_treatment = 30
log2_sd_DE = 2

## ----simulate_data------------------------------------------------------------
# Loading DEP and packages required for data handling
library("DEP")
library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("SummarizedExperiment")

# Background data (sampled from null distribution)
sim_null <- data_frame(
  name = paste0("bg_", rep(1:bg_proteins, rep((2*replicates), bg_proteins))),
  ID = rep(1:bg_proteins, rep((2*replicates), bg_proteins)),
  var = rep(c("control_1", "control_2", "control_3", 
    "treatment_1", "treatment_2", "treatment_3"), bg_proteins), 
  val = 2^rnorm(2*replicates*bg_proteins, mean = log2_mean_bg, sd = log2_sd_bg))

# Data for DE proteins (sampled from alternative distribution)
sim_diff <- rbind(
  data_frame(
    name = paste0("DE_", rep(1:DE_proteins, rep(replicates, DE_proteins))),
    ID = rep((bg_proteins+1):(bg_proteins+DE_proteins), 
      rep(replicates, DE_proteins)),
    var = rep(c("control_1", "control_2", "control_3"), DE_proteins), 
    val = 2^rnorm(replicates*DE_proteins, 
      mean = log2_mean_DE_control, sd = log2_sd_DE)),
  data_frame(
    name = paste0("DE_", rep(1:DE_proteins, rep(replicates, DE_proteins))),
    ID = rep((bg_proteins+1):(bg_proteins+DE_proteins), 
      rep(replicates, DE_proteins)),
    var = rep(c("treatment_1", "treatment_2", "treatment_3"), DE_proteins),
    val = 2^rnorm(replicates*DE_proteins, 
      mean = log2_mean_DE_treatment, sd = log2_sd_DE)))

# Combine null and DE data
sim <- rbind(sim_null, sim_diff) %>% 
  spread(var, val) %>% 
  arrange(ID)

# Generate experimental design
experimental_design <- data_frame(
  label = colnames(sim)[!colnames(sim) %in% c("name", "ID")],
  condition = c(rep("control", replicates), 
    rep("treatment", replicates)),
  replicate = rep(1:replicates, 2))

## ----missing_values_variables-------------------------------------------------
# Variables
MAR_fraction = 0.05
MNAR_proteins = 100

## ----add_missing_values-------------------------------------------------------
# Generate a MAR matrix
MAR_matrix <- matrix(
  data = sample(c(TRUE, FALSE), 
    size = 2*replicates*(bg_proteins+DE_proteins), 
    replace = TRUE, 
    prob = c(MAR_fraction, 1-MAR_fraction)), 
  nrow = bg_proteins+DE_proteins, 
  ncol = 2*replicates)

# Introduce missing values at random (MAR)
controls <- grep("control", colnames(sim))
treatments <- grep("treatment", colnames(sim))
sim[, c(controls, treatments)][MAR_matrix] <- 0
sim$MAR <- apply(MAR_matrix, 1, any)

# Introduce missing values not at random (MNAR)
DE_protein_IDs <- grep("DE", sim$name)
sim[DE_protein_IDs[1:MNAR_proteins], controls] <- 0
sim$MNAR <- FALSE
sim$MNAR[DE_protein_IDs[1:MNAR_proteins]] <- TRUE

## ----generate_SE--------------------------------------------------------------
# Generate a SummarizedExperiment object
sim_unique_names <- make_unique(sim, "name", "ID", delim = ";")
se <- make_se(sim_unique_names, c(controls, treatments), experimental_design)

## ----plot_data_noFilt, fig.width = 4, fig.height = 4--------------------------
# Plot a barplot of the protein quantification overlap between samples
plot_frequency(se)

## ----filter_proteins----------------------------------------------------------
# No filtering
no_filter <- se

# Filter for proteins that are quantified in all replicates of at least one condition
condition_filter <- filter_proteins(se, "condition", thr = 0)

# Filter for proteins that have no missing values
complete_cases <- filter_proteins(se, "complete")

# Filter for proteins that are quantified in at least 2/3 of the samples.
frac_filtered <- filter_proteins(se, "fraction", min = 0.66)

## ----filter_consequences------------------------------------------------------
# Function to extract number of proteins
number_prots <- function(se) {
  names <- rownames(get(se))
  data_frame(Dataset = se,
    bg_proteins = sum(grepl("bg", names)),
    DE_proteins = sum(grepl("DE", names)))
}

# Number of bg and DE proteins still included
objects <- c("no_filter", 
  "condition_filter",
  "complete_cases",
  "frac_filtered")

map_df(objects, number_prots)

## ----scale_vst, message = FALSE-----------------------------------------------
# Scale and variance stabilize
no_filter <- normalize_vsn(se)
condition_filter <- normalize_vsn(condition_filter)
complete_cases <- normalize_vsn(complete_cases)
frac_filtered <- normalize_vsn(frac_filtered)

## ----plot_norm, fig.width = 4, fig.height = 5---------------------------------
# Mean versus Sd plot
meanSdPlot(no_filter)

## ----plot_missval, fig.height = 4, fig.width = 3------------------------------
# Plot a heatmap of proteins with missing values
plot_missval(no_filter)

## ----plot_detect, fig.height = 4, fig.width = 4-------------------------------
# Plot intensity distributions and cumulative fraction of proteins 
# with and without missing values
plot_detect(no_filter)

## ----imputation_methods, error = TRUE-----------------------------------------
# All possible imputation methods are printed in an error, if an invalid function name is given.
impute(no_filter, fun = "")

## ----impute, results = "hide", message = FALSE, warning = FALSE---------------
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

## ----plot_imp, fig.width = 5, fig.height = 7----------------------------------
# Plot intensity distributions before and after imputation
plot_imputation(no_filter, MinProb_imputation, 
  manual_imputation, knn_imputation)

## ----mixed_impuation, results = "hide", message = FALSE, warning = FALSE------
# Extract protein names with missing values 
# in all replicates of at least one condition
proteins_MNAR <- get_df_long(no_filter) %>%
  group_by(name, condition) %>%
  summarize(NAs = all(is.na(intensity))) %>% 
  filter(NAs) %>% 
  pull(name) %>% 
  unique()

# Get a logical vector
MNAR <- names(no_filter) %in% proteins_MNAR

# Perform a mixed imputation
mixed_imputation <- impute(
  no_filter, 
  fun = "mixed",
  randna = !MNAR, # we have to define MAR which is the opposite of MNAR
  mar = "knn", # imputation function for MAR
  mnar = "zero") # imputation function for MNAR

## ----sample_specific_impuation, results = "hide", message = FALSE, warning = FALSE----
# SummarizedExperiment to MSnSet object conversion
sample_specific_imputation <- no_filter
MSnSet <- as(sample_specific_imputation, "MSnSet")

# Impute differently for two sets of samples
MSnSet_imputed1 <- MSnbase::impute(MSnSet[, 1:3], method = "MinProb")
MSnSet_imputed2 <- MSnbase::impute(MSnSet[, 4:6], method = "knn")

# Combine into the SummarizedExperiment object
assay(sample_specific_imputation) <- cbind(
  MSnbase::exprs(MSnSet_imputed1), 
  MSnbase::exprs(MSnSet_imputed2))

## ----plot_imp2, fig.width = 5, fig.height = 5---------------------------------
# Plot intensity distributions before and after imputation
plot_imputation(no_filter, mixed_imputation, sample_specific_imputation)

## ----DE_analysis, message = FALSE, warning = FALSE----------------------------
# Function that wraps around test_diff, add_rejections and get_results functions
DE_analysis <- function(se) {
  se %>% 
    test_diff(., type = "control", control = "control") %>%
    add_rejections(., alpha = 0.1, lfc = 0) %>% 
    get_results()
}

# DE analysis on no, knn, MinProb and mixed imputation
no_imputation_results <- DE_analysis(no_imputation)
knn_imputation_results <- DE_analysis(knn_imputation)
MinProb_imputation_results <- DE_analysis(MinProb_imputation)
mixed_imputation_results <- DE_analysis(mixed_imputation)

## ----rejections---------------------------------------------------------------
# Function to extract number of DE proteins
DE_prots <- function(results) {
  data_frame(Dataset = gsub("_results", "", results),
    significant_proteins = get(results) %>% 
      filter(significant) %>% 
      nrow())
}

# Number of significant proteins
objects <- c("no_imputation_results", 
  "knn_imputation_results",
  "MinProb_imputation_results",
  "mixed_imputation_results")

map_df(objects, DE_prots)

## ----plot_performance---------------------------------------------------------
# Function to obtain ROC data
get_ROC_df <- function(results) {
  get(results) %>% 
    select(name, treatment_vs_control_p.val, significant) %>% 
    mutate(
      DE = grepl("DE", name),
      BG = grepl("bg", name)) %>% 
    arrange(treatment_vs_control_p.val) %>% 
    mutate(
      TPR = cumsum(as.numeric(DE)) / 300,
      FPR = cumsum(as.numeric(BG)) / 3000,
      method = results)
}

# Get ROC data for no, knn, MinProb and mixed imputation
ROC_df <- map_df(objects, get_ROC_df)

# Plot ROC curves
ggplot(ROC_df, aes(FPR, TPR, col = method)) +
  geom_line() +
  theme_DEP1() +
  ggtitle("ROC-curve")

# Plot ROC curves zoom
ggplot(ROC_df, aes(FPR, TPR, col = method)) +
  geom_line() +
  theme_DEP1() +
  xlim(0, 0.1) +
  ggtitle("ROC-curve zoom")

## ----differences_results------------------------------------------------------
# Function to obtain summary data
get_rejected_proteins <- function(results) {
  get(results) %>% 
    filter(significant) %>% 
    left_join(., select(sim, name, MAR, MNAR), by = "name") %>% 
    mutate(
      DE = grepl("DE", name),
      BG = grepl("bg", name),
      method = results)
}

# Get summary data for no, knn, MinProb and mixed imputation
objects <- c("no_imputation_results", 
  "knn_imputation_results",
  "MinProb_imputation_results",
  "mixed_imputation_results")

summary_df <- map_df(objects, get_rejected_proteins)

# Plot number of DE proteins (True and False)
summary_df %>% 
  group_by(method) %>% 
  summarize(TP = sum(DE), FP = sum(BG)) %>% 
  gather(category, number, -method) %>% 
  mutate(method = gsub("_results", "", method)) %>% 
  ggplot(aes(method, number, fill = category)) +
  geom_col(position = position_dodge()) +
  theme_DEP2() +
  labs(title = "True and False Hits",
    x = "",
    y = "Number of DE proteins",
    fill = "False or True")

## ----plot_missval_histogram---------------------------------------------------
# Plot number of DE proteins with missing values
summary_df %>% 
  group_by(method) %>% 
  summarize(MNAR = sum(MNAR), MAR = sum(MAR)) %>% 
  gather(category, number, -method) %>% 
  mutate(method = gsub("_results", "", method)) %>% 
  ggplot(aes(method, number, fill = category)) +
  geom_col(position = position_dodge()) +
  theme_DEP2() +
  labs(title = "Category of Missing Values",
    x = "",
    y = "Number of DE proteins",
    fill = "")

## ----session_info, echo = FALSE-----------------------------------------------
sessionInfo()

