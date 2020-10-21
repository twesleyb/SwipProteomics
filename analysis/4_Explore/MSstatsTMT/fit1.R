#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: fit lmer model for comparisons between Control and Mutant mice
# adjusted for differences in subcellular fraction.

## prepare the env
root <- "~/projects/SwipProteomics"
renv::load(root); devtools::load_all(root)

## inputs:
data(msstats_prot)
data(gene_map)
data(swip)

## other imports
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
})


## check Swip's fit -----------------------------------------------------------

# fit a model
fx1 <- formula("Abundance ~ 0 + Genotype + BioFraction + (1|Subject)")
fm1 <- lmerTest::lmer(fx1, msstats_prot %>% filter(Protein == swip))

model_summary <- summary(fm1, ddf = "Satterthwaite")
model_summary

# FIXME: replace with better code
knitr::kable(r.squaredGLMM.merMod(fm1))

# generate contrast
cm1 <- lme4::fixef(fm1)
cm1[] <- 0
cm1["GenotypeControl"] <- -1
cm1["GenotypeMutant"] <- +1

# test the model
model1 <- lmerTestProtein(swip, fx1, msstats_prot, cm1)
model1$stats %>% knitr::kable()


## loop to fit all proteins ----------------------------------------------------

n_cores <- parallel::detectCores() - 1
BiocParallel::register(BiocParallel::SnowParam(n_cores))

prots = unique(as.character(msstats_prot$Protein))

results_list <- foreach(protein = prots) %dopar% {
	suppressMessages({
	  try(lmerTestProtein(protein, fx1, msstats_prot, cm1),silent=TRUE)
	})
} # EOL


## collect results ------------------------------------------------------------

idx <- unlist(sapply(results_list,class)) != "try-error"
filt_list <- results_list[which(idx)]
results_df <- bind_rows(sapply(filt_list,"[[","stats"))

# drop singular
results_df <- results_df %>% filter(!isSingular)
results_df$isSingular <- NULL

## annotate with gene symbols
idx <- match(results_df$protein,gene_map$uniprot)
results_df <- tibble::add_column(results_df,
  				 symbol=gene_map$symbol[idx],
  				 .after="protein")

## adjust pvals 
results_df <- tibble::add_column(results_df, 
			 Padjust=p.adjust(results_df$Pvalue,"BH"),
			 .after="Pvalue")

## sort
results_df <- results_df %>% arrange(Pvalue)

# examine top results
results_df %>% head() %>% knitr::kable()

# status:
message("Total number of significant proteins: ",
	sum(results_df$Padjust < 0.05))

# save
myfile <- file.path(root,"rdata","fit1_results.csv")
fwrite(results_df, myfile)
