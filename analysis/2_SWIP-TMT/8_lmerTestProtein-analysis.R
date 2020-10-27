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
FDR_alpha <- 0.05 # threshold for significance

## prepare the environment
root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

## load SwipProteomics data
data(swip)
data(gene_map)
data(msstats_prot)
data(alt_contrast) # MSstatsTMT contrast matrix for 'Mutant-Control' comparison
data(msstats_contrasts)

msstats_alt_contrast <- alt_contrast

## other imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(MSstatsTMT)
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
  # adjust pvals for each contrast
  results_df <- results_df %>%
    group_by(Contrast) %>%
    mutate(
      "FDR" = p.adjust(Pvalue, "BH"),
      "Padjust" = p.adjust(Pvalue, "bonferroni")
    )
  # sort
  results_df <- results_df %>% arrange(Pvalue)
  return(results_df)
} # EOF


## illustrate analysis with lmerTest ------------------------------------------

# fit a model
fm <- lmerTest::lmer(fx1,data=msstats_prot %>% filter(Protein == swip))

# given the model, test a single contrast:
lmerTestContrast(fm, cm1) %>% knitr::kable()

## or assess multiple contrasts with lmerTestProtein
# FIXME: need to add capability for multiple proteins 
# FIXME: need to add moderated capability
lmerTestProtein(swip, fx0, msstats_prot, cm0) %>% knitr::kable() 

# NOTE: cm0 is a list of numeric vectors indicating positive and negative coeff
# defining a comparison to be tested

# we can also pass msstats_contrasts (a matrix) to lmerTestProtein:
lmerTestProtein(swip, fx0, msstats_prot, msstats_contrasts) %>% knitr::kable() 


## [1] fx0 / fitBioFraction / intra-fraction comparisons -----------------------

# store results in a list
results <- list()

## analysis with lmerTest
input = list(swip, fx0, msstats_prot, cm0) 
results[["lmerTest::IntraFraction"]] <- do.call(lmerTestProtein, input)

## analysis with MSstats 
input = list("data" = msstats_prot %>% filter(Protein==swip),
	     "contrast.matrix" = msstats_contrasts,
	     "moderated"=FALSE)
results[["MSstatsTMT::IntraFraction"]] <- suppressMessages({
	do.call(groupComparisonTMT, input) })

# FIXME: Remove message printing formula


# [2] fx1 -- fitGenotype - Mutant vs Control ----------------------------------

## analysis with lmerTest
input = list(swip, fx1, msstats_prot, cm1) 
results[["lmerTest::Mutant-Control"]] <- do.call(lmerTestProtein, input)

## analysis with MSstats 

# create a contrast for assessing difference between Control and Mutant
alt_contrast <- matrix(c(-1/7,-1/7,-1/7,-1/7,-1/7,-1/7,-1/7,
		     1/7,1/7,1/7,1/7,1/7,1/7,1/7), nrow=1)
row.names(alt_contrast) <- "Mutant-Control"
colnames(alt_contrast)<- levels(msstats_prot$Condition)


# do MSstatsTMT groupComparisons
input = list("data" = msstats_prot %>% filter(Protein==swip),
	     "contrast.matrix" = alt_contrast,
	     "moderated"=FALSE)
results[["MSstatsTMT::Mutant-Control"]] <- suppressMessages({
	do.call(groupComparisonTMT, input) })


###############################################################################

# NOTE: the above approaches differ in the way the model and contrasts are
# specified. Will lmerTest approach return the same results?

###############################################################################


## Control-Mutant comparisons with alternative contrast and lmerTest ----------
# the result is not the same. Something is different here.
# somehow the result is different, even though we expected the same result

## analysis with lmerTest
input = list(swip, fx0, msstats_prot, alt_contrast) 
results[["lmerTest::alt-Mutant-Control"]] <- do.call(lmerTestProtein, input)

# ^this result matches MSstats

## Examine results:
lapply(results,knitr::kable)


###############################################################################

# which is correct way to specify the model?

# as ther are no sig prots with lmerTestContrast(fm1), it seems like something
# may be wrong

###############################################################################

## timed comparison -----------------------------------------------------------

# do MSstatsTMT groupComparisons
MSstatsTMT <- function() { 
	input = list("data" = msstats_prot %>% filter(Protein==swip),
		     "contrast.matrix" = msstats_contrasts,"moderated"=FALSE)
        results <- suppressMessages({ do.call(groupComparisonTMT, input) })
}

# do lmerTestProtein
lmerTest <- function() {
	input = list(swip, fx0, msstats_prot, msstats_contrasts) 
	results <- do.call(lmerTestProtein, input)
}

timed_res <- microbenchmark(MSstatsTMT(), lmerTest(), times=100L)   

print(timed_res)

# it appears that lmerTest approach is significantly faster ~1.5x
df <- as.data.table(do.call(cbind,timed_res)) 

# time to evaluate 10,000 proteins
fold_diff <- df %>% group_by(expr) %>% 
	summarize(mean=mean(time),.groups="drop") %>%
	mutate("Duration" = 10000*(mean*10^-9)/60)

# savings of ~7 min
lmerTestProtein(swip,fx0,msstats_prot,contrast)

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
fm <- lmerTest::lmer(fx, msstats_prot %>% filter(Protein == swip))

summary(fm,ddf="Satterthwaite")

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
idx <- which(grepl("Control", names(alt_contrast)))
alt_contrast[idx] <- -1 / length(idx)
idx <- which(grepl("Mutant", names(alt_contrast)))
alt_contrast[idx] <- +1 / length(idx)

# combine contrasts
all_contrasts <- c(contrasts, "Mutant-Control" = list(alt_contrast))

# examine the results for swip
df <- lmerTestProtein(swip, fx, msstats_prot, all_contrasts)
df %>%
  unique() %>%
  knitr::kable()

## Check MSstatsTMT results for Swip ------------------------------------------

# Msstats generates the same results:
message("\nMSstatsTMT intra-fraction results:")
MSstatsTMT::groupComparisonTMT(msstats_prot %>% filter(Protein == swip),
			       msstats_contrasts,moderated=FALSE) %>% 
  knitr::kable()

message("\nMSstatsTMT intra-fraction results:")
MSstatsTMT::groupComparisonTMT(msstats_prot %>% filter(Protein == swip),
			       msstats_alt_contrast,moderated=FALSE) %>% 
  knitr::kable()

# NOTE: MSstatsTMT could do both contrasts simultaneously if we defined the 
# contrast matrix correctly


## perform tests for all proteins ----------------------------------------------

# a real world test of lmerTestProtein.

prots <- unique(as.character(msstats_prot$Protein))
message("\nAnalyzing ", formatC(length(prots), big.mark = ","), " proteins.")

n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(cores = n_cores)

t0 <- Sys.time() # start

results_list <- foreach(protein = prots) %dopar% {
  suppressMessages({
    input <- list(protein, fx, msstats_prot, all_contrasts)
    try(do.call(lmerTestProtein, input), silent = TRUE)
  })
} # EOL

t1 <- Sys.time() # stop

# elapsed time
message("\nTime to analyze ", length(prots), " proteins:")
difftime(t1, t0)
message("...utilizing ", n_cores, " parallel processors.")


## process results ------------------------------------------------------------

# collect results as a data.frame
lmerTest_results <- compileTests(results_list)

# split into list of results for each contrast
all_results <- lmerTest_results %>%
  group_by(Contrast) %>%
  group_split()
names(all_results) <- sapply(all_results, function(x) unique(x$Contrast))

# rename
namen <- gsub(
  "ConditionMutant.F[0-9]{1,2}-ConditionControl.", "",
  names(all_results)
)
names(all_results) <- namen

# sort
all_results <- all_results[c(biofraction, "Mutant-Control")]

# check number of sig prots
sapply(all_results,function(x) sum(x$FDR<0.05)) %>% t() %>% knitr::kable()

# add protein data
all_results <- c("Normalized Protein" = list(msstats_prot), all_results)



