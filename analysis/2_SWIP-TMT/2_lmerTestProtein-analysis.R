#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: fit protein-wise lmer models and perform statistical inference

## formulae to be fit:
# [1] fx0: Abundance ~ 0 + Condition + (1|Mixture)
# [2] fx1: Aundance ~ 0 + Genotype:BioFraction + (1|Subject) + (1|Mixture)

# To assess 'Mutant-Control' comparisons we should fit formula 2.
# However, 'Subject' and 'Mixture' are partially confounded.
# Therefore, can either choose to account for the effect of 'Subject' or
# 'Mixture', but not both. 

# As 'Mixture' contributes more to the overall variance, it makes sense to
# account for 'Mixture' and not the random effect inherint in the repeated
# measures of each 'Subject'.

# Removing 'Subject' from the model, we have:
# [3] fx1: Aundance ~ 0 + Genotype:BioFraction + (1|Mixture)

# this is equivalent to:
# [1] fx0: Abundance ~ 0 + Condition + (1|Mixture)
# When Condition is interaction(Genotype,BioFraction)

# Thus, we can actually perform the contrast of interest using the MSstatsTMT,
# provided the correct contrast matrix.


## prepare the env ------------------------------------------------------------

## input
save_rda <- TRUE
FDR_alpha <- 0.05 # threshold for significance

## load renv
root <- "~/projects/SwipProteomics"
renv::load(root)

## load SwipProteomics
devtools::load_all(root)
data(swip)
data(gene_map)
data(msstats_prot)

## other imports
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
  ## other requirements:
  # require(lme4)
  # require(knitr)
  # require(tibble)
  # require(lmerTest)
  # requre(data.table)
})


## functions ------------------------------------------------------------------

getContrast <- function(fm, negative_index, positive_index) {
  # build a contrast matrix from model coefficients
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


compileTests <- function(results_list) {
  # clean-up the results list return by lmerTestProtein loop
  # remove errors
  idx <- unlist(sapply(results_list, class)) != "try-error"
  filt_list <- results_list[which(idx)]
  results_df <- do.call(rbind, filt_list)

  # drop singular
  results_df <- results_df %>% filter(!isSingular)
  results_df$isSingular <- NULL

  # annotate with gene symbols
  idx <- match(results_df$Protein, gene_map$uniprot)
  results_df$Symbol <- gene_map$symbol[idx]

  # Remove any redundant results
  results_df <- results_df %>% unique()

  # for each contrast perform moderation??
  #y = results_df %>% group_by(Contrast) %>% group_split()
  #x = y[[1]]
  #mod_stats = limma::squeezeVar(x$SE^2,x$DF)
  #df_prior = mod_stats$df.prior
  #s2_prior = mod_stats$var.prior
  #lmerTestContrast(fm,kkk

  # adjust pvals for each contrast
  results_df <- results_df %>% group_by(Contrast) %>% 
  	  mutate("FDR" = p.adjust(Pvalue, "BH"),
		 "Padjust" = p.adjust(Pvalue, "bonferroni"))
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

## fit the model
fm <- lmerTest::lmer(fx,msstats_prot %>% filter(Protein == swip))

## evaluate goodness-of-fit
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

## create contrast for Mutant-Control comparison
alt_contrast <- lme4::fixef(fm)
alt_contrast[] <- 0
idx <- which(grepl("Control",names(alt_contrast)))
alt_contrast[idx] <- -1/length(idx)
idx <- which(grepl("Mutant",names(alt_contrast)))
alt_contrast[idx] <- +1/length(idx)

# combine contrats
all_contrasts <- c(contrasts,"Mutant-Control" = list(alt_contrast))

# examine the results for swip
df <- lmerTestProtein(swip,fx,msstats_prot,all_contrasts)
df %>% unique() %>% arrange(Pvalue) %>% knitr::kable()


## perform tests for all proteins ----------------------------------------------

prots <- unique(as.character(msstats_prot$Protein))
message("\nAnalyzing ", formatC(length(prots),big.mark=","), " proteins.")

n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(cores = n_cores)

t0 <- Sys.time() # start

results_list <- foreach(protein = prots) %dopar% {
  suppressMessages({
	  input <- list(protein, fx, msstats_prot, all_contrasts)
	  try(do.call(lmerTestProtein,input), silent = TRUE)
  })
} # EOL

t1 <- Sys.time() # stop

# elapsed time
message("\nTime to analyze ", length(prots), " proteins:")
difftime(t1, t0)


## process results ------------------------------------------------------------

# collect results as a data.frame
lmerTest_results <- compileTests(results_list)

# split into list of results for each contrast
all_results <- lmerTest_results %>% group_by(Contrast) %>% group_split()
names(all_results) <- sapply(all_results,function(x) unique(x$Contrast))

# rename
namen <- gsub("ConditionMutant.F[0-9]{1,2}-ConditionControl.","",
	      names(all_results))
names(all_results) <- namen

# sort
all_results <- all_results[c(biofraction,"Mutant-Control")]

# add protein data
all_results <- c("Normalized Protein" = list(msstats_prot), all_results)


## save results ----------------------------------------------------------------




## save results ----------------------------------------------------------------

# save as excel
myfile <- file.path(root, "tables", "S3_lmerTest_Results.xlsx")
write_excel(all_results, myfile)

# save other key R objects
if (save_rda) {

  #myfile <- file.path(root, "data", "fm1.rda")
  #save(fm1, file = myfile, version = 2)

  #myfile <- file.path(root, "data", "cm0.rda")
  #save(cm0, file = myfile, version = 2)

  #myfile <- file.path(root, "data", "fx1.rda")
  #save(fx1, file = myfile, version = 2)

  #myfile <- file.path(root, "data", "cm1.rda")
  #save(cm1, file = myfile, version = 2)

  myfile <- file.path(root, "data", "lmerTest_results.rda")
  save(lmerTest_results, file = myfile, version = 2)

}
