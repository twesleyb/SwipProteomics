#!/usr/bin/env Rscript

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root)

# load project
devtools::load_all(root)

# load project's data
data(swip)
data(gene_map)
data(partition)
data(msstats_prot)

# other imports
suppressPackageStartupMessages({
  library(dplyr)
  library(lmerTest)
  library(data.table)
})


## remove any unclustered proteins
msstats_filt <- msstats_prot %>% filter(Protein %in% names(partition))


## annotate data with module membership
msstats_filt$Module <- paste0("M",partition[msstats_filt$Protein])

# wash prots
washc <- gene_map$uniprot[grepl("Washc*",gene_map$symbol)]

# the model to be fit:
fx <- formula(paste(c("Abundance ~ 0 + Genotype:BioFraction + (1|Mixture)",
		      " + (1|Protein)"), collapse=" "))

fm <- lmerTest::lmer(fx, data=msstats_filt %>% filter(Protein %in% washc))

summary(fm,ddf="Satterthwaite")

# goodness of fit
r.squaredGLMM.merMod(fm) %>% knitr::kable()

## build a contrast
contrast = lme4::fixef(fm)
contrast[] <- 0
idx <- which(grepl("Control",names(contrast)))
contrast[idx] <- -1/length(idx)
idy <- which(grepl("Mutant",names(contrast)))
contrast[idy] <- +1/length(idy)

# asses contrast:
lmerTestContrast(fm, contrast) %>% 
	mutate(Contrast = "Mutant-Control") %>% unique() %>% knitr::kable() 


## loop through all modules
results_list <- list()
for (module in names(modules)){
  input <- list(fx, data=msstats_filt %>% filter(Module==module))
  fm <- suppressMessages({ # about boundary fits
	  try(do.call(lmerTest::lmer,input), silent=TRUE)
  })
  df <- lmerTestContrast(fm, contrast) %>% 
	mutate(Contrast = "Mutant-Control") %>% unique()
  results_list[[module]] <- df
}

## collect results
results_df <- do.call(rbind,results_list) 

warning(sum(results_df$isSingular),
	" modules with singular fit will be removed.")

results_df <- results_df %>% filter(!isSingular) %>% 
	arrange(Pvalue) %>% 
	as.data.table(keep.rownames="Module") %>% 
	mutate(FDR = p.adjust(Pvalue,method="BH"),
	       PAdjust = p.adjust(Pvalue,method="bonferroni"))

## examine top results
knitr::kable(head(results_df))
