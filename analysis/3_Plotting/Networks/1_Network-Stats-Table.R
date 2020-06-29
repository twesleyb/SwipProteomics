#!/usr/bin/env Rscript

#' ---
#' title: Swip Proteomics Plotting
#' description: Plot an overview of the network.
#' authors: Tyler W Bradshaw
#' ---

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
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. 

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Imports.
suppressPackageStartupMessages({
	library(gridExtra)
	library(data.table)
	library(ggplot2)
	library(grid)
})

suppressWarnings({ devtools::load_all() })

# Set plotting theme.
set_font("Arial",font_path=file.path(root,"fonts"))

# Load the data.
data(partition)
data(tmt_protein)
data(module_stats)

# Summarize key module statistics.
n_nodes <- formatC(length(unique(tmt_protein$Accession)),big.mark=",")
n_modules <- length(unique(partition))-1 # M0 is not a module
p_clustered <- round(100*sum(partition!=0)/length(partition),2)
median_pve <- median(module_stats$PVE)
median_size <- median(module_stats$Nodes)

# Create table.
df <- data.table("N Nodes" = n_nodes,
		 "N Modules" = n_modules,
		 "Percent Clustered" = p_clustered,
		 "Median Size" = median_size)
gtab <- tableGrob(df, rows=NULL, theme=ttheme_default())

# Save.
myfile <- file.path(root,"figs","Networks","Network_Summary.pdf")
ggsaveTable(gtab,myfile)
