#!/usr/bin/env Rscript

root = "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

data(swip)
data(partition)
data(msstats_prot)

suppressPackageStartupMessages({
  library(dplyr)
  library(lmerTest)
  library(data.table)
})


# remove any unclustered proteins
msstats_filt <- msstats_prot %>% filter(Protein %in% names(partition))

# annotate data with module membership
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

# examine gof:
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
lmerTestContrast(fm, contrast)


################################################################################

# loop through all modules
modules <- unique(msstats_filt$Module)

# the function to be fit:
fx <- formula("Abundance ~ 0 + Genotype + BioFraction + (1|Protein)")

fit_list <- list()
for (module in modules) {
	fm <- lmerTest::lmer(fx, data=msstats_filt %>% filter(Module==module))
	fit_list[[module]] <- fm
}

# extract a fit
fm = fit_list[["M1"]]
summary(fm,ddf="Satterthwaite")

# examine gof:
gof <- setNames(as.numeric(r.squaredGLMM.merMod(fm)),nm=c("R2fixef","R2total"))
gof["R2mixef"] <- gof["R2total"] - gof["R2fixef"]
gof <- gof[c("R2mixef","R2fixef","R2total")]
t(gof) %>% knitr::kable()

# build a contrast
contrast = lme4::fixef(fm)
contrast[] <- 0
contrast["GenotypeMutant"] <- +1
contrast["GenotypeControl"] <- -1

# asses contrast:
lmerTestContrast(fm, contrast)

# p-values do not seem correct -- too small

# loop to do 
results <- list()
for (i in seq(fit_list)) {
	fm <- fit_list[[i]]
	results[[i]] <- lmerTestContrast(fm, contrast)
}

names(results) <- names(fit_list)	

df = do.call(rbind,results) %>% arrange(Pvalue) %>% 
	as.data.table(keep.rownames="Module")

knitr::kable(head(df))

################################################################################

data(gene_map)

#washc = gene_map$symbol[grepl("Washc*",gene_map$symbol)]
washc = gene_map$uniprot[grepl("Washc*",gene_map$symbol)]

fx <- formula("Abundance ~ 0 + Genotype + BioFraction + (1|Mixture) + (1|Protein)")
fm <- lmerTest::lmer(fx, data=msstats_filt %>% filter(Protein %in% washc))

summary(fm,ddf="Satterthwaite")

r.squaredGLMM.merMod(fm) %>% knitr::kable()


# build a contrast
contrast = lme4::fixef(fm)
contrast[] <- 0
contrast["GenotypeMutant"] <- +1
contrast["GenotypeControl"] <- -1

# asses contrast:
lmerTestContrast(fm, contrast)

prots <- names(partition[partition==partition[swip]])
fm <- lmerTest::lmer(fx, data=msstats_filt %>% filter(Protein %in% prots))
summary(fm,ddf="Satterthwaite")
r.squaredGLMM.merMod(fm) %>% knitr::kable()



