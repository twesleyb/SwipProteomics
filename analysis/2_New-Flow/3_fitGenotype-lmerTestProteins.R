#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: fit protein-wise lmer models and perform statistical inferences
# for given contrasts

## formulae to be fit:
# [1] fx0: Abundance ~ 0 + Condition + (1|Mixture)

# [?] fx1: Aundance ~ 0 + Genotype:BioFraction + (1|Subject) + (1|Mixture)

# This is the model we really want to fit. But Subject and Mixture are
# confounded. We can either choose to account for the effect of Subject or Mixture. 
# Mixture contributes more to variance, thus it makes sense to use Mixture.
# The model becomres:

# [2] fx1: Aundance ~ 0 + Genotype:BioFraction + (1|Mixture)

# this is equivalent to :
# [1] fx0: Abundance ~ 0 + Condition + (1|Mixture)
# When Condition is interaction(Genotype,BioFraction)

# Thus we can actually perform the contrast of interest using the MSstatsTMT
# given the correct contrast matrix. We define this below.


## prepare the env ------------------------------------------------------------

## options
save_rda = FALSE
results_file = "fitGenotype_lmerTestProtein_results.xlsx" # saved in root/tables
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
#formula("Abundance ~ 0 + Genotype + BioFraction + (1|Subject)")

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

## adjust pvals
results_df <- tibble::add_column(results_df,
  Padjust = p.adjust(results_df$Pvalue, "BH"),
  .after = "Pvalue"
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
  sum(results_df$Padjust < FDR_alpha)
)


## save results ----------------------------------------------------------------

# save as excel
myfile <- file.path(root, "tables", results_file)
write_excel(results_df, myfile)

# save as rda
fit1_results <- results_df
myfile <- file.path(root, "data", "fitGenotype_results.rda")
save(fit1_results, file = myfile, version = 2)
