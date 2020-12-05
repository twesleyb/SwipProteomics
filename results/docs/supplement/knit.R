#!/usr/bin/env Rscript

root <- "~/projects/SwipProteomics"
renv::load(root)

knitr::knit("supplement.Rnw",quiet=TRUE)
