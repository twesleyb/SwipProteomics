#!/usr/bin/env Rscript

root <- "~/projects/SwipProteomics"
renv::load(root)

knitr::knit("response.Rnw",quiet=TRUE)
