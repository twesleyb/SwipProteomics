#!/usr/bin/env Rscript

## ---- input
root = "~/projects/SwipProteomics"

## ---- renv
if (!dir.exists(file.path(root,"renv"))){
  msg <- c("This executable expects to be run in the SwipProteomics project tree.",
	 "This executable expectes that R dependencies are managed by renv in root/renv.")
  stop(mgs)
} else {
  renv::load(root,quiet=TRUE)
  devtools::load_all(root,quiet=TRUE)
}

## ---- load data

data(swip)
data(tmt_prot)

## ---- fit the linear model

fx = Abundance ~ 0 + Condition

fm = lm(fx, tmt_prot %>% subset(Protein == swip))


## ---- assess a contrast

LT = getContrast(fm,"MUT","WT")

lmTestContrast(fm,LT)
