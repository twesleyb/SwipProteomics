#!/usr/bin/env Rscript

#' ---
#' title: Swip TMT proteomics
#' description: Analysis of modules for GO enrichment
#' authors: Tyler W A Bradshaw
#' ---

## User parameters to change:
FDR_alpha = 0.05 # significance threshold for gse enrichment.

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

#---------------------------------------------------------------------
## Set-up the workspace.
#---------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global options and imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
})

# Load additional functions in root/R.
suppressWarnings({ devtools::load_all() })

# Project directories:
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

#---------------------------------------------------------------------
## Load the data.
#---------------------------------------------------------------------

# Load the partition and tmt data.
data(partition)
data(tmt_protein)

#---------------------------------------------------------------------
## Module GO enrichment.
#---------------------------------------------------------------------

# List of modules.
module_list <- split(partition,partition)
names(module_list) <- paste0("M",names(module_list))

# Collect all sig prots.
sig_prots <- tmt_protein %>% filter(FDR < 0.1) %>% 
	select(Accession) %>% unlist() %>% unique()

# collect list of entrez ids and gene symbols.
gene_list <- lapply(module_list,function(x) {
			    tmt_protein$Entrez[match(names(x),tmt_protein$Accession)]
			    })

gene_symbols <- lapply(module_list,function(x) {
			    tmt_protein$Gene[match(names(x),tmt_protein$Accession)]
			    })

# Build a GO reference collection:
message("\nBuilding a mouse GO reference collection with anRichment.")
gocollection <- suppressPackageStartupMessages({
	anRichment::buildGOcollection(organism="mouse")
})

# perform module go enrichment.
# NOTE: background defaults to all genes in gene_list.
go_gse <- gse(gene_list,gocollection)

# collect significant modules.
idx <- sapply(go_gse,function(x) any(x$FDR < FDR_alpha))
nsig_go <- sum(idx)
message(paste("\nNumber of modules with significant",
	      "GO term enrichment:",nsig_go))

# top go term for each module.
fe <- sapply(go_gse,function(x) round(x$enrichmentRatio[1],2))
fdr <- sapply(go_gse,function(x) round(x$FDR[1],2))
m <- sapply(go_gse, function(x) x$shortDataSetName[1])
top_go <- data.table(module=names(go_gse),
		       term = m,
		       fe=fe,fdr=fdr)

# summary:
message("\nModules with significant go enrichment:")
knitr::kable(top_go[idx,])

# save significant results.
results_list <- list("GO Results" = bind_rows(go_gse[idx],.id="Module"))
write_excel(results_list,file.path(tabsdir,"Swip_TMT_Module_GO_Results.xlsx"))

# Save as rda.
module_GO <- results_list[["GO Results"]]
myfile <- file.path(root,"data","module_GO.rda")
save(module_GO,file=myfile,version=2)
