#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

## Optional parameters:
FDR_alpha = 0.10 # FDR significance threshold.

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

start <- Sys.time()
message(paste("Starting analysis at:",start))

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
})

# Functions.
suppressWarnings({ devtools::load_all() })

# Directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

#--------------------------------------------------------------------
## Load the proteomics data.
#--------------------------------------------------------------------

# Load the data from root/data.
data(tmt_protein)

# Load the graph partition:
data(partition)

# Load WASH iBioID interactome.
data(wash_interactome)
wash_prots <- unique(wash_interactome$Accession) # Get uniprot accession

# List of all modules. 
module_list <- split(names(partition),partition)[-1] # drop M0

# Background is union of WASH BioID and proteins in network.
background <- unique(c(unlist(module_list),wash_prots))

# Loop to perform hypergeometric test for enrichment.
message("\nPerforming hypergeometric tests to assess WASH BioID enrichment.")
results_list <- lapply(c(1:length(module_list)), function(x) {
	       hyperTest(wash_prots,module_list[[x]],background)
}) 
names(results_list) <- paste0("M",names(module_list))

# Collect results.
hyper_dt <- as.data.table(do.call(rbind,results_list),keep.rownames="Module")

# Adjust p-values.
hyper_dt$FDR <- p.adjust(hyper_dt$"P-value",method="BH")
hyper_dt$P.adjust <- p.adjust(hyper_dt$"P-value",method="bonferroni")

# Calculate n BioID proteins per module.
hyper_dt <- tibble::add_column(hyper_dt, N=sapply(module_list,length),
				.after="Module")
n <- sapply(module_list,function(x) sum(x %in% wash_prots))
hyper_dt <- tibble::add_column(hyper_dt, "n BioID Proteins"=n,.after="N")

# Pretty print:
message("Modules that are enriched for WASH iBioID proteins:")
knitr::kable(hyper_dt %>% filter(P.adjust < FDR_alpha))

# Done!
end <- Sys.time()
message(paste("\nCompleted analysis at:",end))
message(paste("Elapsed time:",
	      round(difftime(end,start,units="mins"),2),"minutes."))
