#!/usr/bin/env Rscript

# prepare the env
root = "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

# load the data
data(fx0) # protein model
data(fx1) # module model
data(swip)
data(gene_map)
data(partition)
data(msstats_prot)
data(pd_annotation)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(variancePartition)
})


## prepare the data -----------------------------------------------------------

# anno msstats_prot with Module membership
msstats_prot <- msstats_prot %>% filter(Protein %in% names(partition)) %>% 
	mutate(Module = paste0("M",partition[Protein]))

# munge to split Condition (Geno.BioFrac) into geno annotation
samples <- pd_annotation
samples$Genotype <- sapply(strsplit(samples$Condition,"\\."),"[",1)

# munge to create mouse Subject identifier
samples$Subject <- as.numeric(interaction(samples$Genotype,samples$Experiment))

# munge to annotate msstats_prot with Subject (1-6)
msstats_prot$Subject <- as.numeric(interaction(msstats_prot$Genotype,msstats_prot$Mixture))

# cast the data into a matrix.
fx <- formula(Protein ~ Mixture + Genotype + BioFraction)
dm <- msstats_prot %>%
	reshape2::dcast(fx, value.var= "Abundance") %>% 
	as.data.table() %>% na.omit() %>% as.data.table() %>% as.matrix(rownames="Protein")

# munge to create sample info from the dcast fx
info <- as.data.table(do.call(rbind,strsplit(colnames(dm),"_")))
colnames(info) <- strsplit(as.character(fx)[3]," \\+ ")[[1]]

# calculate variance explained by major covariates
form <- formula(~ (1|Mixture) + (1|Genotype) + (1|BioFraction))
prot_varpart <- fitExtractVarPartModel(dm, form, info)

# collect results
df <- as.data.frame(prot_varpart)
varpart_df <- as.data.table(df,keep.rownames="Protein")

# annotate with gene ids
idx <- match(varpart_df$Protein,gene_map$uniprot)
varpart_df$Symbol <- gene_map$symbol[idx]
varpart_df$Entrez <- gene_map$entrez[idx]
# sort cols
varpart_df <- varpart_df %>%
	select(Protein,Symbol,Entrez,Mixture,Genotype,BioFraction,Residuals) %>%
	arrange(desc(Genotype))

## fit all protein models -----------------------------------------------------

## ----------------------------------------------------------------------------


# Combine varpart_df with R2.mermod for each model -- save as goodness of fit in
# module results

# examine the overall goodness-of-fit in the sense of how much variation we can
# attribute to each of the major covariates in the data in  the full model:
fx1 <- Abundance ~ (1|Mixture) + (1|Genotype) + (1|BioFraction) + (1|Module) + (1|Protein)
fm1 <- lmerTest::lmer(fx1, msstats_prot %>% filter(Module != "M0"))

rho <- calcVarPart(fm1)
knitr::kable(round(rho,3))


## save ------------------------------------------------------------------------

myfile <- file.path(root,"data","prot_varpart.rda")
save(prot_varpart, file=myfile, version=2)
#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: fit protein-wise lmer models and perform statistical inference
# with lmerTest.

## formulae to be fit:
# [1] fx0: Abundance ~ 0 + Genotype:BioFraction + (1|Mixture)
# [2] fx1: Aundance ~ 0 + Genotype:BioFraction + (1|Mixture) + (1|Protein)

# To assess intra-BioFraction comparisons we should fit formula 1.
# This is equivalent to MSstatsTMT formula:
# [0] fx0: Abundance ~ 0 + Condition + (1|Mixture)
# When Condition is interaction(Genotype,BioFraction)

# To assess 'Mutant-Control' comparisons we should fit formula 1 as well.
# As Subject and Mixture are confounded we choose to model mixture and not
# subject. Thus, we can actually perform the contrast of interest using 
# MSstatsTMT, provided the correct contrast matrix.


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

model_summary <- summary(fm,ddf="Satterthwaite")

df <- model_summary$coefficients
df %>% as.data.table(keep.rownames="Coefficient") %>% mutate(Coeffi


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
