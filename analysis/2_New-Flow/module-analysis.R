#!/usr/bin/env Rscript

# load renv
root = "~/projects/SwipProteomics"
renv::load(root)

# load project
devtools::load_all(root)

# load data
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

# the function to be fit:
#fx0 <- formula("Abundance ~ 0 + Genotype")
#fx1 <- formula("Abundance ~ 0 + Genotype + BioFraction")
#fx2 <- formula("Abundance ~ 0 + Genotype + BioFraction + (1|Protein)")
#fx3 <- formula("Abundance ~ 0 + Genotype + BioFraction + (1|Mixture) + (1|Protein)")
#anova(fm2,fm3) # justify inclusion of Mixture

module = paste0("M",partition[swip])

fx <- formula("Abundance ~ 0 + Genotype + BioFraction + (1|Mixture) + (1|Protein)")
fm <- lmerTest::lmer(fx, data=msstats_filt %>% filter(Module==module))

summary(fm, ddf="Satterthwaite")


## examine gof:
r.squaredGLMM.merMod(fm) %>% knitr::kable()


## build a contrast
contrast = lme4::fixef(fm)
contrast[] <- 0
contrast["GenotypeMutant"] <- +1
contrast["GenotypeControl"] <- -1

#contrast = lme4::fixef(fm)
#contrast[] <- -(1/5)
#contrast["BioFractionF6"] <- +1
#contrast["GenotypeMutant"] <- 0
#contrast["GenotypeControl"] <- 0


# asses contrast:
lmerTestContrast(fm, contrast) %>% knitr::kable()


## loop to fit modules
fit_list <- list()
modules <- unique(msstats_filt$Module)

for (module in modules) {
	fm <- lmerTest::lmer(fx, data=msstats_filt %>% filter(Module==module))
	fit_list[[module]] <- fm
}


## loop to assess contrast
results <- list()
for (i in seq(fit_list)) {
	fm <- fit_list[[i]]
	results[[i]] <- lmerTestContrast(fm, contrast)
}
names(results) <- names(fit_list)	


## collect results
df <- do.call(rbind,results) %>% arrange(Pvalue) %>% 
	as.data.table(keep.rownames="Module")


## examine top results
knitr::kable(head(df))


# wash prots
washc = gene_map$uniprot[grepl("Washc*",gene_map$symbol)]

fm <- lmerTest::lmer(fx, data=msstats_filt %>% filter(Protein %in% washc))

summary(fm,ddf="Satterthwaite")

# goodness of fit
r.squaredGLMM.merMod(fm) %>% knitr::kable()

# asses contrast:
lmerTestContrast(fm, contrast) %>% knitr::kable()


fx0 <- formula("Abundance ~ 0 + Genotype + BioFraction + (1|Mixture) + (1|Protein)")
fx1 <- formula("Abundance ~ 0 + BioFraction + (1|Mixture) + (1|Protein)")

fm0 <- lmerTest::lmer(fx0, data=msstats_filt %>% filter(Protein %in% washc))
fm1 <- lmerTest::lmer(fx1, data=msstats_filt %>% filter(Protein %in% washc))

anova(fm0,fm1)

fx0 <- formula("Abundance ~ 0 + Genotype + BioFraction + (1|Mixture) + (1|Protein)")
fx1 <- formula("Abundance ~ 0 + BioFraction + (1|Mixture) + (1|Protein)")

prots <- names(which(partition==19))
fm0 <- lmerTest::lmer(fx0, data=msstats_filt %>% filter(Protein %in% prots))
fm1 <- lmerTest::lmer(fx1, data=msstats_filt %>% filter(Protein %in% prots))

av <- anova(fm0,fm1)
summary(av)
