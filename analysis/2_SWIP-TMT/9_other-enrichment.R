#!/usr/bin/env Rscript

#' ---
#' title: Swip TMT Proteomics
#' description: analysis of modules for enrichment of WASH proteome
#' authors: Tyler W Bradshaw
#' ---

## Optional parameters:
FDR_alpha = 0.10 # FDR significance threshold.

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
	library(geneLists) # for bioid gene lists
})

# Load functions in root/R and data in root/data.
suppressWarnings({ devtools::load_all() })

# Directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Load the geneList data from geneLists.
data(list="iPSD") # iPSD$Arhgef9, Gphn, InSyn1, iPSD
data(list="ciPSD") # ciPSD$ciPSD
data(list="spence2019") # spence2019$"wrp_P5-interactome"
data(list="gao2020") # gao2020$PV, SST, CamkII
data(list="dube2020") # dube2020$"syp-preesynapse"

# Load the data from root/data.
data(wash_interactome) # Courtland et al., 2020 - WASH1 BioID.
data(gene_map) # gene mapping data
data(partition) # graph partition
data(tmt_protein) # the proteomics data

#--------------------------------------------------------------------
## Do work.
#--------------------------------------------------------------------

# Wash genes.
wash_prots <- unique(wash_interactome$Accession)
wash_genes <- na.omit(gene_map$entrez[match(wash_prots,gene_map$uniprot)])

# Collect list of entrez ids.
gene_lists <- list(
		   "Arhgef9"=iPSD$Arhgef9,
		   "Gphn" = iPSD$Gphn,
		   "Insyn1" = iPSD$InSyn1,
		   "Wrp-P5" = spence2019$"wrp_P5-interactome",
		   "PV-Dlg4" = gao2020$PV,
		   "SST-Dlg4" = gao2020$SST,
		   "CamkII-Dlg4" = gao2020$CamkII,
		   "Synaptosome-Syp" = dube2020$"syp-presynapse",
		   "Wash1" = unique(wash_genes)
		   )

# collect list of modules.
modules <- split(names(partition),partition)[-1]
module_entrez <- lapply(modules,function(x) gene_map$entrez[match(x,gene_map$uniprot)])

# Loop to perform GSE for each pathway.
results <- list()
for (experiment in names(gene_lists)){
	# Get pathway specific genes.
	pathway_genes <- gene_lists[[experiment]]
	# Background is union of WASH BioID and lysosome proteins in network.
	background <- unique(c(unlist(module_entrez),pathway_genes))
	# Loop to perform hypergeometric test for enrichment.
	results_list <- list()
	for (i in c(1:length(module_entrez))) {
		results_list[[i]] <- hyperTest(pathway_genes,module_entrez[[i]],background)
	}
	names(results_list) <- paste0("M",names(modules))
	# Collect results in a data.table.
	hyper_dt <- as.data.table(do.call(rbind,results_list),keep.rownames="Module")
	# Adjust p-values.
	hyper_dt$FDR <- p.adjust(hyper_dt$"P-value",method="BH")
	hyper_dt$P.adjust <- p.adjust(hyper_dt$"P-value",method="bonferroni")
	# Add module size annotation.
	sizes <- sapply(module_entrez,length)
	hyper_dt <- tibble::add_column(hyper_dt,"Module Size"=sizes,.after="Module")
	# Add pathway annotation.
	hyper_dt <- tibble::add_column(hyper_dt,Pathway=experiment,.after="Module Size")
	# Add total number of pathway genes.
	hyper_dt <- tibble::add_column(hyper_dt,"N Pathway Genes"=length(pathway_genes),.after="Pathway")
	n <- sapply(module_entrez,function(x) sum(pathway_genes %in% x))
	hyper_dt <- tibble::add_column(hyper_dt,"n Pathway Genes in Module"=n,.after="N Pathway Genes")
	# Return the results.
	results[[experiment]] <- hyper_dt
}

# Collect the results in a single data.table.
dt <- bind_rows(results)
sig_modules <- dt %>% filter(FDR < 0.05)

# Pretty print.
knitr::kable(sig_modules)
