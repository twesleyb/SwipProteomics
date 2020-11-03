#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: generate some gof statistics for protein and module-level models

save_rda = TRUE 

# prepare the env
root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

# load the data
data(fx1) # module-level model
data(swip)
data(gene_map)
data(sigprots)
data(partition)
data(poor_prots) # proteins with poor fm0 lmer fit (R2<thresh)
data(msstats_prot)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(doParallel)
  library(variancePartition)
})



#------------------------------------------------------------------------------
# above we calculated the R2 for fixed and total affects, but we also are
# curious to know how much variation is explained by our clustering

# how much variation is explained by partition?

prot_df <- msstats_prot %>% 
	filter(Protein %in% names(partition)) %>% 
	mutate(Module = paste0("M",partition[Protein]))

# fit to protein
form0 <- formula(Abundance ~ (1|Protein)) 
fm0 <- lmer(form0, data=prot_df %>% filter(Module != "M0"))
rho1 <- calcVarPart(fm1)

# fit to partition
form1 <- formula(Abundance ~ (1|Module)) 
fm1 <- lmer(form1, data=prot_df %>% filter(Module != "M0"))
rho0 <- calcVarPart(fm0)

t(rho0) %>% knitr::kable() # [1]

t(rho1) %>% knitr::kable() # [2] 

# shuffle partition - refit Module
prot_df <- prot_df %>% mutate(Module = sample(Module))
fm2 <- lmer(form1, data=prot_df %>% filter(Module != "M0"))
rho2 <- calcVarPart(fm2)

t(rho2) %>% knitr::kable() # bad

# so we can do alot better than random but we are not approaching the ideal case
# of partition == Protein

# partition == protein - refit Module
prot_df <- msstats_prot %>% 
	filter(Protein %in% names(partition)) %>% 
	mutate(Module = as.factor(Protein))
fm3 <- lmer(form1, prot_df %>% filter(Protein %notin% names(which(partition==0))))
rho3 <- calcVarPart(fm3)
t(rho3) %>% knitr::kable()  # same as [1]
