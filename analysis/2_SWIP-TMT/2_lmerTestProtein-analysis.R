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
save_rda = TRUE
results_file <- "lmerTestProtein_Results.xlsx" # saved in root/tables
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


compileTests <- function(results_list, moderated=FALSE) {
  # FIXME: where is strange class from?
  # class(results_list)
  # [1] "vctrs_list_of" "vctrs_vctr"    "list"
  # clean-up the results list return by lmerTestProtein loop
  class(results_list) <- "list"
  idx <- unlist(sapply(results_list, class)) != "try-error"
  filt_list <- results_list[which(idx)]
  results_df <- do.call(rbind, filt_list)
  # drop singular
  results_df <- results_df %>% filter(!isSingular)
  results_df$isSingular <- NULL
  # annotate with gene symbols
  idx <- match(results_df$Protein, gene_map$uniprot)
  results_df$Symbol <- gene_map$symbol[idx]
  # Remove redundant contrasts
  results_df <- results_df %>% mutate(Contrast = "Mutant-Control") %>% unique()
  # adjust pvals
  # adjust pvals  for each contrast
  #results_df <- results_df %>% group_by(Contrast) %>% 
  #	  mutate("Padjust" = p.adjust(Pvalue, "BH"))
  results_df$FDR <- p.adjust(results_df$Pvalue, "BH")
  results_df$Padjust <- p.adjust(results_df$Pvalue, "bonferroni")
  # sort
  results_df <- results_df %>% arrange(Pvalue)
  return(results_df)
} #EOF


## check Swip's fit -----------------------------------------------------------

## formula to be fit:
fx <- formula("Abundance ~ 0 + Condition + (1|Mixture)")

# status
gene <- gene_map$symbol[which(gene_map$uniprot == swip)]
message(
  "\nlmer: ", as.character(fx)[2],
  "(", gene, ") ~ ", as.character(fx)[3]
)


## fit model 1
fm <- lmerTest::lmer(fx, msstats_prot %>% filter(Protein == swip))

# model summary with Satterthwaite degrees of freedom:
print(summary(fm1, ddf = "Satterthwaite"))


# create contrast
alt_contrast <- lme4::fixef(fm1)
alt_contrast[] <- 0

idx <- which(grepl("Control",names(alt_contrast)))
alt_contrast[idx] <- -1/length(idx)

idx <- which(grepl("Mutant",names(alt_contrast)))
alt_contrast[idx] <- +1/length(idx)
cm1 <- alt_contrast


# examine the results for swip
results <- lmerTestContrast(fm1, cm1)

results %>% mutate(Contrast = "Mutant-Control") %>% unique() %>% knitr::kable()

# goodness-of-fit
r2_nakagawa <- r.squaredGLMM.merMod(fm)

knitr::kable(rbind(c("marginal/fixef", "conditional/total"), r2_nakagawa))


## munge to create contrast matrices for intrafraction comparisons:
condition <- c("ConditionControl", "ConditionMutant")
biofraction <- c("F4", "F5", "F6", "F7", "F8", "F9", "F10")
contrasts <- lapply(expandGroups(condition, biofraction), function(x) {
  getContrast(fm, x[1], x[2])
})
names(contrasts) <- sapply(contrasts, function(x) {
  paste(names(x)[x == +1], names(x)[x == -1], sep = "-")
})
cm0 <- contrasts


# check the results for swip
message("\nResults for WASHC4:")
results_df <- lmerTestProtein(swip, fx, msstats_prot, contrasts)

results_df %>% knitr::kable()


## assess intrafraction comparisons -------------------------------------------

prots <- unique(as.character(msstats_prot$Protein))
message("\nAnalyzing ", formatC(length(prots),big.mark=","), " proteins.")

n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(cores = n_cores)

t0 <- Sys.time() # start

results_list <- foreach(protein = prots) %dopar% {
  suppressMessages({
    try(lmerTestProtein(protein, fx, msstats_prot, contrasts), silent = TRUE)
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

## perform p-value moderation
eb_fit <- limma::squeezeVar(results_df$SE, results_df$DF)

#eb_fit$df.prior

#eb_fit$var.prior

## annotate with gene symbols
idx <- match(results_df$Protein, gene_map$uniprot)
results_df <- tibble::add_column(results_df,
  Symbol = gene_map$symbol[idx],
  .after = "Protein"
)

# sort cols
results_df <- results_df %>%
  select(
    Protein, Symbol, Contrast, log2FC,
    percentControl, Pvalue, Padjust, SE, DF
  )

# sort rows
results_df <- results_df %>% arrange(Pvalue)


# status
message("\nSummary of significant proteins for intrafraction comparisons:")
results_df %>% group_by(Contrast) %>% 
	summarize(nsig = sum(Padjust < FDR_alpha)) %>% knitr::kable()


## loop to assess Mutant-Control contrast -------------------------------------

prots <- unique(as.character(msstats_prot$Protein))

n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(cores = n_cores)

t0 <- Sys.time() # start

results_list <- foreach(protein = prots) %dopar% {
  suppressMessages({
	  input = list(protein,fx,msstats_prot,cm1)
	  try(do.call(lmerTestProtein,input), silent = TRUE)
  })
} # EOL

t1 <- Sys.time() # stop

message("\nTime to analyze ", length(prots), " proteins:")
difftime(t1, t0)


## save results ----------------------------------------------------------------

# save results as excel document
myfile <- file.path(root, "tables", results_file)
write_excel(results_list, myfile)

if (save_rda) {

  myfile <- file.path(root, "data", "fm1.rda")
  save(fm1, file = myfile, version = 2)

  myfile <- file.path(root, "data", "cm0.rda")
  save(cm0, file = myfile, version = 2)

  myfile <- file.path(root, "data", "fx1.rda")
  save(fx1, file = myfile, version = 2)

  myfile <- file.path(root, "data", "cm1.rda")
  save(cm1, file = myfile, version = 2)
}
