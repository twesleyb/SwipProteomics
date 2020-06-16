#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

## Analysis options:

## R Options:
options(renv.config.synchronized.check = FALSE) # skip renv::check(repo).
options(renv.settings.snapshot.type = "simple") # use simple renv::snapshot.

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
  library(RCy3)
  library(dplyr)
  library(igraph)
  library(data.table)
})

# Functions.
suppressWarnings({ devtools::load_all() })

# Project directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

## Output directory for cytoscape networks.
netwdir <- file.path(root,"networks")
if (!dir.exists(netwdir)) {
	dir.create(netwdir)
}

#--------------------------------------------------------------------
## Load the data.
#--------------------------------------------------------------------

# Load the data from root/data.
data(tmt_protein)

# Load the graph partition:
data(partition)

# Load wash interactome.
data(wash_interactome)
wash_prots <- unique(wash_interactome$Accession) # Get uniprot accession

# Load networks.
data(ne_adjm)
data(ppi_adjm)

# Load gene map
data(gene_map)

#--------------------------------------------------------------------
## Create igraph graph objects.
#--------------------------------------------------------------------

# List of all modules. 
module_list <- split(names(partition),partition)[-1] # drop M0
names(module_list) <- paste0("M",names(module_list))

# Insure that matrices are in matching order.
check <- all(colnames(ne_adjm) == colnames(ppi_adjm))
if (!check) { stop() }

# Create igraph graph objects.
netw_g <- graph_from_adjacency_matrix(ne_adjm,mode="undirected",diag=FALSE,
				      weighted=TRUE)
ppi_g <- graph_from_adjacency_matrix(ppi_adjm,mode="undirected",diag=FALSE,
				     weighted=TRUE)

# Annotate graph's with gene symols.
symbols <- gene_map$symbol[match(names(V(netw_g)),gene_map$uniprot)]
netw_g <- set_vertex_attr(netw_g,"symbol",value = symbols)
symbols <- gene_map$symbol[match(names(V(ppi_g)),gene_map$uniprot)]
ppi_g <- set_vertex_attr(ppi_g,"symbol",value = symbols)

#--------------------------------------------------------------------
## Annotate graphs with additional meta data.
#--------------------------------------------------------------------

idx <- match(names(V(netw_g)),tmt_protein$Accession)


#--------------------------------------------------------------------
## Create Cytoscape graphs.
#--------------------------------------------------------------------

## FIXME: remove existing images.
imgsdir <- file.path(netwdir,"Modules")

if (!dir.exists(imgsdir)) {
	# Create the directory.
	dir.create(imgsdir)
} else {
	# Remove any existing figures.
	invisible({ file.remove(list.files(imgsdir,full.name=TRUE)) })
}

# Loop to create graphs:
for (module_name in names(module_list)){
	nodes = module_list[[module_name]]
	createCytoscapeGraph(netw_g, ppi_g, nodes, module_name, 
			     netwdir=netwdir,imgsdir=imgsdir)
}

# When done, save cytoscape session.
myfile <- file.path(netwdir,paste0("Modules.cys"))
winfile <- gsub("/mnt/d/","D:/",myfile)
saveSession(winfile)

# Convert to svg to pdf.

## Combine modules as a single pdf.

# Directory for grouped pdfs.
output_dir <- file.path(figsdir,"Final-Markers")
dir.create(output_dir,recursive=TRUE)

message("Combining markers and saving as single pdf...")

for (module in names(marker_modules)){
# Get subset of proteins.
prots <- paste0(gsub("\\|","_",marker_modules[[module]]),".pdf")
check <- all(prots %in% names(all_plots))
if (!check) { stop("We are missing some plots!") } 
# Create directory for grouped plots.
to_dir <- file.path(output_dir,gsub(" ","_",module))
dir.create(to_dir)
# Remove any existing plots.
unlink(list.files(to_dir,full.names=TRUE))
# Copy to new directory.
namen <- names(all_plots[prots])
response <- file.copy(from=all_plots[prots],
		      to=file.path(to_dir,namen))
if (!all(response)) { stop("Problem saving pdfs.") }

# Combined into a single pdf with ghostscript cli utility.
# Create gs command.
myfile <- file.path(to_dir,paste0(gsub(" ","_",module),".pdf"))
cmd <- c("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite",
	 "-sOutputFile=OUTPUT.pdf", "INPUT.pdf")
cmd <- gsub("OUTPUT.pdf",myfile,cmd)
cmd <- gsub("INPUT.pdf",file.path(to_dir,"*.pdf"),cmd)
cmd <- paste(cmd,collapse=" ")

# Execute system command
response <- system(cmd,intern=TRUE)

# Remove input pdfs.
unlink(file.path(to_dir,namen))

# Done!
end <- Sys.time()
message(paste("\nCompleted analysis at:",end))
message(paste("Elapsed time:",
	      round(difftime(end,start,units="mins"),2),"minutes."))
