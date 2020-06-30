#!/usr/bin/env Rscript

load_renv <- function(root,...){
	renv::load(root,...)
	devtools::load_all()
}

root <- dirname(dirname(getwd()))
load_renv(root,quiet=TRUE)

# Compare sig prots at level 1 and level 2. 
data(tmt_protein)

data(sig_proteins)
