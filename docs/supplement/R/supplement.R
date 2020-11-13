#!/usr/bin/env Rscript

# title: supplment.R
# description: code associated with supplement.Rnw
# author: twab

# NOTE: code chunks are delimited with '## ---- LABEL' where LABEL cooresponds
# to its label in the main *.Rnw file. Declare code chunk options in the main
# R noweave file. 

# NOTE: it seems best if your approach each chunk as a stand-alone bit of code.

# prepare the renv
root <- "~/projects/swipproteomics"
renv::load(root)
devtools::load_all(root)

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})


## ---- fit0

## fit the protein-level model to WASHC4

# load dependencies
library(dplyr)
library(lmerTest)

# load SwipProteomics
data(swip)
data(msstats_prot)

# LMM formula
fx0 <- 'Abundance ~ 0 + Genotype:BioFraction + (1|Mixture)'

# fit the model
fm0 <- lmer(fx0, data = msstats_prot %>% subset(Protein == swip))

# examine the model's summary
summary(fm0, ddf = "Satterthwaite")

## ---- tab0

library(grid)
library(gtable)
library(gridExtra)
library(ggplot2)
library(cowplot)

df <- summary(fm0, ddf="Satterthwaite")[["coefficients"]] %>%
	as.data.table(keep.rownames="Coefficient")
df$Coefficient <- gsub("Genotype|BioFraction","",df$Coefficient)
colnames(df)[colnames(df)=="Pr(>|t|)"] <- "p value"
colnames(df)[colnames(df)=="Std. Error"] <- "SE"
colnames(df)[colnames(df)=="df"] <- "DF"
df$"p value" <- formatC(df$"p value",digits=3)
df$Estimate <- round(df$Estimate,2)
df$SE <- round(df$SE,3)
df$DF <- round(df$DF,2)
df$"t value" <- round(df$"t value",2)
tt <- gridExtra::ttheme_default(
      base_size = 11,
      core = list(bg_params = list(fill = "white"))
    )
tab <- tableGrob(df, rows = NULL, theme = tt)
# add outside border
g <- gtable_add_grob(tab, grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)), 
		     t = 1, b = nrow(tab), l = 1, r = ncol(tab))
# add rows -- i have no idea how gtable works
for (i in c(1:nrow(tab))) {
	g <- gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)), 
		     t = i, b = i, l = 1, r = ncol(tab))
}
plot_grid(g)

## ---- alt0

# fit the model
fx0 <- 'Abundance ~ 0 + Genotype + BioFraction + (1|Mixture)'
fm0 <- lmer(fx0, data = msstats_prot %>% subset(Protein == swip))

# examine the model's summary
summary(fm0, ddf = "Satterthwaite")


## ---- contrast7 

## compare 'Mutant:F7' and 'Control:F7' conditions

# create a contrast
coeff <- lme4::fixef(fm0)
contrast7 <- setNames(rep(0,length(coeff)), nm = names(coeff))
contrast7["GenotypeMutant:BioFractionF7"] <- +1 # positive coeff
contrast7["GenotypeControl:BioFractionF7"] <- -1 # negative coeff

# evaluate contrast
lmerTestContrast(fm0, contrast7)


## ---- contrast8 

# create a contrast to compare 'Mutant' versus 'Control'
contrast8 <- getContrast(fm0, "Mutant","Control")

# evaluate contrast
lmerTestContrast(fm0, contrast8)


## ---- fit1

# the module-level formula to be fit:
fx1 <- 'Abundance ~ 0 + Condition + (1|Mixture) + (1|Protein)'

# load WASH Complex proteins
data(washc_prots)

fm1 <- lmer(fx1, data=msstats_prot %>% subset(Protein %in% washc_prots))

# assess 'Mutant-Control' comparison
lmerTestContrast(fm1, contrast8)


## ---- flexibility

# the module-level formula to be fit:
fx2 <- 'Abundance ~ 0 + Genotype + BioFraction + (1|Mixture) + (1|Protein)'
fm2 <- lmer(fx2, data = msstats_prot %>% subset(Protein %in% washc_prots))
lT <- getContrast(fm2,"Mutant","Control")
lmerTestContrast(fm2, lT)


## ---- nakagawa

# assess gof with Nakagawa coefficient of determination
r.squaredGLMM.merMod(fm0)

r.squaredGLMM.merMod(fm1)


## ---- variancePartition

# load variancePartition
suppressPackageStartupMessages({
	library(variancePartition)
})

# calculate partitioned variance
form <- "Abundance ~ (1|Genotype) + (1|BioFraction) + (1|Mixture) + (1|Protein)"
fit <- lmer(form, data = msstats_prot %>% filter(Protein %in% washc_prots))

calcVarPart(fit)
