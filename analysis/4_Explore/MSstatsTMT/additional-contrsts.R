#!/usr/bin/env Rscript

# 1. formula("Abundance ~ 0 + Condition + (1 | Mixture)")
# 2. formula("Abundance ~ 0 + Genotype + BioFraction + (1 | Subject)")
# 3. formula("Abundance ~ 0 + Genotype + ... + (1|Protein)")

## prepare the environment
root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)


## load project data
data(swip)
data(gene_map)
data(msstats_prot)
data(fx0); data(cm0)
data(fx1); data(cm1)
data(msstats_contrasts)


## imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(MSstatsTMT)
})


## illustrate analysis with lmerTest

# fit a model
fm <- lmerTest::lmer(fx0,data=msstats_prot %>% filter(Protein == swip))

# the appropriate contrast for Control-Mutant comparison:
contrast <- lme4::fixef(fm)
contrast[] <- 0
idx <- grepl("Control",names(contrast))
idy <- grepl("Mutant",names(contrast))
contrast[idx] <- -1/length(idx)
contrast[idy] <- +1/length(idy)

# given the model, test a single contrast:
lmerTestContrast(fm, contrast) %>% knitr::kable()

## or assess multiple contrasts with lmerTestProtein
# FIXME: need to add capability for multiple proteins 
# FIXME: need to add moderated capability
lmerTestProtein(swip, fx0, msstats_prot, contrast) %>% knitr::kable() 

# NOTE: cm0 is a list of numeric vectors indicating positive and negative coeff
# defining a comparison to be tested
# FIXME: can we report CV?
# that way, we can say reduced to 59.11 % +/- CV WT.

# With this contrast it is clearer how we might assess other interesting
# comparisons, like BioFraction vs all else.
contrast <- lme4::fixef(fm)
contrast[] <- 0
idx <- grepl("F4",names(contrast))
contrast[idx] <- -1/sum(idx)
idy <- !grepl("F4",names(contrast))
contrast[idy] <- +1/sum(idy)

prots = unique(as.character(msstats_prot$Protein))
prots = sample(prots,10)

n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(cores = n_cores)

library(doParallel)

results_list <- foreach(protein = prots) %dopar% {
  suppressMessages({
    try(lmerTestProtein(protein, fx0, msstats_prot, contrast), silent = TRUE)
  })
} # EOL

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
results_df <- results_df %>% mutate(Contrast = "F4") %>% unique()


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



quit()

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
