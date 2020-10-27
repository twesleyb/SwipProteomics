#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: fit protein-wise lmer models and perform statistical inference
# with lmerTest.

## formulae to be fit:
# [1] fx1: Abundance ~ 0 + Genotype:BioFraction + (1|Mixture)
# [2] fx2: Aundance ~ 0 + Genotype:BioFraction + (1|Subject) + (1|Mixture)

# To assess intra-BioFraction comparisons we should fit formula 1.
# This is equivalent to MSstatsTMT formula:
# [0] fx0: Abundance ~ 0 + Condition + (1|Mixture)

# To assess 'Mutant-Control' comparisons we should fit formula 2.
# However, 'Subject' and 'Mixture' are partially confounded.
# Therefore, can either choose to account for the effect of 'Subject' or
# 'Mixture', but not both.

# Removing 'Subject' from the model, we have:
# [3] fx3: Aundance ~ 0 + Genotype:BioFraction + (1|Mixture)

# this is equivalent to the model used by MSstatsTMT:
# [0] fx0: Abundance ~ 0 + Condition + (1|Mixture)
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
data(alt_contrast)
data(msstats_contrasts)

# contrast matrix for 'Mutant-Control' comparison
msstats_alt_contrast <- alt_contrast

## other imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(MSstatsTMT)
  library(doParallel)
  library(microbenchmark)
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


## illustrate analysis with lmerTest ------------------------------------------

## formula to be fit:
fx0 <- formula("Abundance ~ 0 + Condition + (1|Mixture)")

# status
gene <- gene_map$symbol[which(gene_map$uniprot == swip)]
message(
  "\nlmer: ", as.character(fx0)[2],
  "(", gene, ") ~ ", as.character(fx0)[3]
)

## fit the model
idx <- msstats_prot$Protein == swip
fm <- lmerTest::lmer(fx0, msstats_prot[idx,])

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

## assess multiple contrasts with lmerTestProtein
df <- lmerTestProtein(swip, fx0, msstats_prot, all_contrasts)
df %>% as.data.table() %>% unique() %>% knitr::kable()


## Check MSstatsTMT results for Swip ------------------------------------------

# Msstats generates the same results:
message("\nMSstatsTMT intra-fraction results:")
idx <- msstats_prot$Protein == swip
MSstatsTMT::groupComparisonTMT(msstats_prot[idx,], msstats_contrasts, moderated=FALSE) %>% 
  knitr::kable()

message("\nMSstatsTMT intra-fraction results:")
MSstatsTMT::groupComparisonTMT(msstats_prot[idx,], 
			       msstats_alt_contrast, moderated=FALSE) %>% 
  knitr::kable()

# NOTE: MSstatsTMT could do both contrasts simultaneously if we defined the 
# contrast matrix correctly
