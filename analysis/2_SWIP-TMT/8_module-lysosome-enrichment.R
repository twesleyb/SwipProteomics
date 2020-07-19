#!/usr/bin/env Rscript

#' ---
#' title: Swip TMT Proteomics
#' description: analysis of modules for enrichment of WASH proteome
#' authors: Tyler W Bradshaw
#' ---

## Optional parameters:
BH_alpha = 0.10 # FDR significance threshold.

#---------------------------------------------------------------------
## Misc function - getrd().
#---------------------------------------------------------------------

# Get the repository's root directory.
getrd <- function(here=getwd(), dpat= ".git") {
	in_root <- function(h=here, dir=dpat) { 
		check <- any(grepl(dir,list.dirs(h,recursive=FALSE))) 
		return(check)
	}
	# Loop to find root.
	while (!in_root(here)) { 
		here <- dirname(here) 
	}
	root <- here
	return(root)
}

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
	library(geneLists) # for lysosome gene list
})

# Functions.
suppressWarnings({ devtools::load_all() })

# Directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Load the data from root/data.
data(hlgd) # human and mouse Lysosome genes from http://lysosome.unipg.it/ 
data(gene_map) # gene mapping data
data(partition) # graph partition
data(tmt_protein) # the proteomics data

#--------------------------------------------------------------------
## Do work.
#--------------------------------------------------------------------

# Load lysosome entrez ids.
lyso_prots <- unique(unlist(hlgd))

knitr::kable(data.table("Lysosome Proteins" = length(lyso_prots)))

# FIXME: GO BACK AND MAKE SUURE THAT LYSO GENES ARE MOUSE!

# Map module uniprot to entrez.
idx <- match(names(partition),gene_map$uniprot)
modules <- split(gene_map$entrez[idx],partition)[-1]

# Background is union of WASH BioID and lysosome proteins in network.
background <- unique(c(unlist(modules),lyso_prots))

# Loop to perform hypergeometric test for enrichment.
results_list <- list()
message(paste("\nPerforming hypergeometric tests to assess",
	     "modules for lysosome protein enrichment."))

for (i in c(1:length(modules))) {
	       results_list[[i]] <- hyperTest(lyso_prots,modules[[i]],background)
}
names(results_list) <- paste0("M",names(modules))

# Collect results.
hyper_dt <- as.data.table(do.call(rbind,results_list),keep.rownames="Module")

# Adjust p-values.
hyper_dt$FDR <- p.adjust(hyper_dt$"P-value",method="BH")
hyper_dt$P.adjust <- p.adjust(hyper_dt$"P-value",method="bonferroni")

# Calculate n lyso proteins per module.
hyper_dt <- tibble::add_column(hyper_dt, N=sapply(modules,length),
				.after="Module")
n <- sapply(modules,function(x) sum(x %in% lyso_prots))
hyper_dt <- tibble::add_column(hyper_dt, "n Lysosomal Proteins"=n,.after="N")

# Pretty print:
#message("Modules that are enriched for Lysosome genes:")
knitr::kable(head(hyper_dt %>% arrange(`P-value`)))

fwrite(hyper_dt,file=file.path(rdatdir,"lyso_gse.csv"))
