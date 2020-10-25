#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: fit protein-wise lmer models and perform statistical inferences
# for defined contrasts

## formulae to be fit:
# [1] fx0: Abundance ~ 0 + Condition + (1|Mixture)
# [2] fx1: Aundance ~ 0 + Genotype:BioFraction + (1|Subject) + (1|Mixture)

# To assess 'Mutant-Control' comparisons we should fit formula 2.
# However, 'Subject' and 'Mixture' are confounded.
# Therefore can either choose to account for the effect of 'Subject' or
# 'Mixture', but not both. 
# As 'Mixture' contributes more to the overall variance, it makes sense to
# account for 'Mixture' and not the random effect inherint in the repeated
# mesaures of each 'Subject'.
# Removing 'Subject' from the model, we have:
# [3] fx1: Aundance ~ 0 + Genotype:BioFraction + (1|Mixture)

# this is equivalent to:
# [1] fx0: Abundance ~ 0 + Condition + (1|Mixture)
# When Condition is interaction(Genotype,BioFraction)

# Thus, we can actually perform the contrast of interest using the MSstatsTMT,
# provided the correct contrast matrix.


## prepare the env ------------------------------------------------------------

## input
results_file <- "BioFraction_lmerTestProtein_results.xlsx" # saved in root/rdata
FDR_alpha <- 0.05 # threshold for significance

## load renv
root <- "~/projects/SwipProteomics"
renv::load(root)

devtools::load_all(root)

## load data
data(swip)
data(gene_map)
data(msstats_prot)

## other imports
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
  # require(lme4)
  # require(knitr)
  # require(tibble)
  # require(lmerTest)
  # requre(data.table)
})


## functions ------------------------------------------------------------------

getContrast <- function(fm, negative_index, positive_index) {
  # build a contrast matrix:
  contrast_matrix <- lme4::fixef(fm)
  contrast_matrix[] <- 0
  contrast_matrix[positive_index] <- +1
  contrast_matrix[negative_index] <- -1
  return(contrast_matrix)
} # EOF


expandGroups <- function(conditions, biofractions) {
  # munge to create contrast matrix for fm1
  groups <- apply(expand.grid(conditions, biofractions), 1, paste, collapse = ".")
  idx <- rep(c(1:length(biofraction)), each = length(condition))
  contrast_list <- split(groups, idx)
  return(contrast_list)
} # EOF


## check Swip's fit -----------------------------------------------------------

## formula to be fit:
fx0 <- formula("Abundance ~ 0 + Condition + (1|Mixture)")

message("\nfit: ", 
	paste(as.character(fx0)[2], as.character(fx0)[3], sep = " ~ "))

# save model formula
myfile <- file.path(root, "data", "fx0.rda")
save(fx0, file = myfile, version = 2)


## fit model 0 and get Satterthwaite degrees of freedon
fm0 <- lmerTest::lmer(fx0, msstats_prot %>% filter(Protein == swip))

print(summary(fm0, ddf = "Satterthwaite"))

# save fit model
myfile <- file.path(root, "data", "fm0.rda")
save(fm0, file = myfile, version = 2)

## evaluate gof using a function ported directly from the MuMin package
r2_nakagawa <- r.squaredGLMM.merMod(fm0)

knitr::kable(rbind(c("marginal/fixef", "conditional/total"), r2_nakagawa))

## FIXME: add some more info about calculation

## munge to create contrast matrices for intrafraction comparisons:
condition <- c("ConditionControl", "ConditionMutant")
biofraction <- c("F4", "F5", "F6", "F7", "F8", "F9", "F10")
contrasts <- lapply(expandGroups(condition, biofraction), function(x) {
  getContrast(fm0, x[1], x[2])
})
names(contrasts) <- sapply(contrasts, function(x) {
  paste(names(x)[x == +1], names(x)[x == -1], sep = "-")
})


# save contrasts
cm0 <- contrasts
myfile <- file.path(root, "data", "cm0.rda")
save(cm0, file = myfile, version = 2)


# check the results for swip
message("\nResults for WASHC4:")
results_df <- lmerTestProtein(swip, fx0, msstats_prot, contrasts)

results_df %>% knitr::kable()

# you can also directly test a contrast:
#lmerTestContrast(fm0,contrasts[[1]])


## loop to fit all proteins ----------------------------------------------------

prots <- unique(as.character(msstats_prot$Protein))
message("\nAnalyzing ", formatC(length(prots),big.mark=","), " proteins.")

n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(cores = n_cores)

t0 <- Sys.time() # start

results_list <- foreach(protein = prots) %dopar% {
  suppressMessages({
    try(lmerTestProtein(protein, fx0, msstats_prot, contrasts), silent = TRUE)
  })
} # EOL

t1 <- Sys.time() # stop

# elapsed time
message("\nTime to analyze ", length(prots), " proteins:")
difftime(t1, t0)


## process results ------------------------------------------------------------

# collect results
idx <- unlist(sapply(results_list, class)) != "try-error"
filt_list <- results_list[which(idx)]
results_df <- do.call(rbind, filt_list)
rownames(results_df) <- NULL

# drop singular
results_df <- results_df %>% filter(!isSingular)
results_df$isSingular <- NULL

## annotate with gene symbols
idx <- match(results_df$Protein, gene_map$uniprot)
results_df <- tibble::add_column(results_df,
  Symbol = gene_map$symbol[idx],
  .after = "Protein"
)

## adjust pvals  for each contrast
results_df <- results_df %>%
  group_by(Contrast) %>%
  mutate("Padjust" = p.adjust(Pvalue, "BH"))


# sort cols
results_df <- results_df %>%
  select(
    Protein, Symbol, Contrast, log2FC,
    percentControl, Pvalue, Padjust, SE, DF
  )

# sort rows
results_df <- results_df %>% arrange(Pvalue)

# status
message(
  "\nTotal instances of significant change: ",
  sum(results_df$Padjust < FDR_alpha)
)

# total (unique) sig prots
sig_prots <- results_df %>%
  ungroup() %>%
  filter(Padjust < FDR_alpha) %>%
  select(Protein) %>%
  unique() %>%
  unlist()
message(
  "\nTotal number of significant proteins: ",
  length(sig_prots)
)


## save results ----------------------------------------------------------------

# results split by BioFraction

results_list <- results_df %>%
  group_by(Contrast) %>%
  group_split()

# munge to get list names
namen <- sapply(strsplit(
  sapply(results_list, function(x) unique(x$Contrast)),
  "\\."
), "[", 3) # get third element after split, F#
names(results_list) <- namen

# sort the results by biofraction
results_list <- results_list[biofraction]

# FIXME: where is strange class from?
# class(results_list)
# [1] "vctrs_list_of" "vctrs_vctr"    "list"
class(results_list) <- "list"

# save as excel
myfile <- file.path(root, "tables", results_file)
write_excel(results_list, myfile)

# status
message("\nSummary of significant proteins for intrafraction comparisons:")
df <- data.frame(nsig = sapply(results_list, function(x) {
  sum(x$Padjust < FDR_alpha)
}))
knitr::kable(df)

# proteins
prot_list <- sapply(results_list, function(x) x$Protein[x$Padjust < FDR_alpha])
common_prots <- Reduce(intersect, prot_list)
if (length(common_prots) > 0) {
  message(
    "\nCommonly significant proteins: ",
    paste(common_prots, collapse = ", ")
  )
}

# save as rda
fit0_results <- results_df
myfile <- file.path(root, "data", "fitBioFraction_results.rda")
save(fit0_results, file = myfile, version = 2)

########################################################################################


## prepare the env ------------------------------------------------------------

## options
save_rda = FALSE
results_file = "lmerTest_Control-Mutant.csv" # saved in root/tables
FDR_alpha <- 0.05 # threshold for significance

## load renv
root <- "~/projects/SwipProteomics"
renv::load(root)

# load project functions
devtools::load_all(root)

## load input data
data(swip)
data(gene_map)
data(msstats_prot)

## other imports
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
  ## other packages used:
  # require(lme4)
  # require(knitr)
  # require(tibble)
  # require(lmerTest)
  # requre(data.table)
})


## functions ------------------------------------------------------------------

getContrast <- function(fm, negative_index, positive_index) {
  # build a contrast matrix:
  contrast_matrix <- lme4::fixef(fm)
  contrast_matrix[] <- 0
  contrast_matrix[positive_index] <- +1
  contrast_matrix[negative_index] <- -1
  return(contrast_matrix)
} # EOF


expandGroups <- function(conditions, biofractions) {
  # munge to create contrast matrix for fm1
  groups <- apply(expand.grid(conditions, biofractions), 1, paste, collapse = ".")
  idx <- rep(c(1:length(biofraction)), each = length(condition))
  contrast_list <- split(groups, idx)
  return(contrast_list)
} # EOF


## check Swip's fit -----------------------------------------------------------


## formula to be fit:
fx1 <- formula("Abundance ~ 0 + Condition + (1|Mixture)")

# status
gene <- gene_map$symbol[which(gene_map$uniprot == swip)]
message(
  "\nlmer: ", as.character(fx1)[2],
  "(", gene, ") ~ ", as.character(fx1)[3]
)

myfile <- file.path(root, "data", "fx1.rda")
save(fx1, file = myfile, version = 2)

## fit model 1
# NOTE: the underlying model will be the same!
fm1 <- lmerTest::lmer(fx1, msstats_prot %>% filter(Protein == swip))

# model summary with Satterthwaite degrees of freedom:
print(summary(fm1, ddf = "Satterthwaite"))

myfile <- file.path(root, "data", "fm1.rda")
save(fm1, file = myfile, version = 2)

# create contrast vector
# NOTE: here's the appropriate contrast matrix:
alt_contrast <- lme4::fixef(fm1)
alt_contrast[] <- 0

idx <- which(grepl("Control",names(alt_contrast)))
alt_contrast[idx] <- -1/length(idx)

idx <- which(grepl("Mutant",names(alt_contrast)))
alt_contrast[idx] <- +1/length(idx)

cm1 <- alt_contrast

if (save_rda) {
  myfile <- file.path(root, "data", "cm1.rda")
  save(cm1, file = myfile, version = 2)
}

# check the results for swip
results <- lmerTestContrast(fm1, cm1)

results %>% knitr::kable()

# goodness-of-fit
r.squaredGLMM.merMod(fm1) %>% knitr::kable()


## loop to fit all proteins ----------------------------------------------------

prots <- unique(as.character(msstats_prot$Protein))

n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(cores = n_cores)

t0 <- Sys.time() # start

results_list <- foreach(protein = prots) %dopar% {
  suppressMessages({
    try(lmerTestProtein(
      protein, fx1,
      msstats_prot, cm1
    ), silent = TRUE)
  })
} # EOL

t1 <- Sys.time() # stop

message("\nTime to analyze ", length(prots), " proteins:")
difftime(t1, t0)


## process results ------------------------------------------------------------

# collect results
idx <- unlist(sapply(results_list, class)) != "try-error"
filt_list <- results_list[which(idx)]
results_df <- do.call(rbind, filt_list)

# drop singular
results_df <- results_df %>% filter(!isSingular)
results_df$isSingular <- NULL

## annotate with gene symbols
idx <- match(results_df$Protein, gene_map$uniprot)
results_df <- tibble::add_column(results_df,
  Symbol = gene_map$symbol[idx],
  .after = "Protein"
)


## Remove redundant contrasts
results_df <- results_df %>% mutate(Contrast = "Mutant-Control") %>% unique()


## adjust pvals
results_df <- tibble::add_column(results_df,
  FDR = p.adjust(results_df$Pvalue, "BH"),
  .after = "Pvalue"
)

results_df <- tibble::add_column(results_df,
  Padjust = p.adjust(results_df$Pvalue, "bonferroni"),
  .after = "FDR"
)

## sort
results_df <- results_df %>% arrange(Pvalue)

# examine top results
results_df %>%
  head() %>%
  knitr::kable()

# status
message(
  "Total number of significant proteins: ",
  sum(results_df$FDR < FDR_alpha), 
  " (FDR < " , FDR_alpha,")."
)

# status
message(
  "Total number of significant proteins: ",
  sum(results_df$Padjust < FDR_alpha),
  " (Bonferroni < " , FDR_alpha,")."
)


## save results ----------------------------------------------------------------

# save as excel
myfile <- file.path(root, "tables", results_file)
write_excel(results_df, myfile)

# save as rda
fit1_results <- results_df
myfile <- file.path(root, "data", "fitGenotype_results.rda")
save(fit1_results, file = myfile, version = 2)
