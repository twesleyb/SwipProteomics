#!/usr/bin/env Rscript

# Analysing Swip (WASH) BioID proteomics.

## User parameters to change:
FDR_alpha = 0.1
enrichment_threshold = log2(3.5)

## Input data in root/rdata
data_file = "WASH_BioID_Results.RData"

## Output in root/tables:
# * WASH_Network_PPIs.xlsx

#-------------------------------------------------------------------------------
## Prepare the workspace.
#-------------------------------------------------------------------------------

start <- Sys.time()
message(paste("Starting analysis at:",start))

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Imports
suppressPackageStartupMessages({
	library(dplyr) # For manipulating data.
	library(igraph) # For networks.
	#library(TBmiscr) # For misc. utilities.
	library(getPPIs) # For mapping gene names.
	library(data.table) # For data.tables.
	library(openxlsx)
})

# Directories.
datadir <- file.path(root,"data")
rdatdir <- file.path(root,"rdata")
tabsdir <- file.path(root,"tables")

# Load any additional functions in root/R.
devtools::load_all()

# Load BioID results.
myfile <- list.files(rdatdir,pattern=data_file,full.names=TRUE)[1]
results <- readRDS(myfile)

#---------------------------------------------------------------------
## Create PPI graph.
#---------------------------------------------------------------------

# Protein identifiers are Gene|Uniprot.
ids <- paste(results$Gene,results$Accession,sep="|")
names(ids) <- results$Entrez

# Get entrez identifiers for proteins of interest.
sig <- results$FDR < FDR_alpha
up <- results$logFC > enrichment_threshold
entrez <- names(ids)[sig & up]

# Load PPIs from getPPIs.
data(musInteractome)

# Subset dataset--get PPIs among our proteins of entrest.
ppis <- musInteractome %>% 
	filter(osEntrezA %in% entrez & osEntrezB %in% entrez) %>% 
	select(osEntrezA,osEntrezB,Interactor_A_Taxonomy,
	       Interaction_detection_methods,
	       Publications,Interaction_type,Confidence_score,
	       Source_database,Methods)

# Create edges data.frame -- aka SIF or simple interaction file.
sif <- ppis 
colnames(sif)[c(1,2,3)] <- c("EntrezA","EntrezB","Taxonomy")
sif$Weight <- 1

# Map Entrez to protein ids.
sif <- tibble::add_column(sif,ProteinA=ids[as.character(sif$EntrezA)],
			  .before="EntrezA")
sif <- tibble::add_column(sif,ProteinB=ids[as.character(sif$EntrezB)],
			  .after="ProteinA")

# Node attributes -- first column must be matching protein identifiers.
noa <- results %>% filter(Entrez %in% entrez)
noa <- tibble::add_column(noa,Vertex=ids[as.character(noa$Entrez)],.before=1)

# Create graph.
message("\nCreating PPI graph from WASH iBioID interactome.")
g <- graph_from_data_frame(sif,vertices=noa,directed=FALSE)

# Summarize number of nodes and edges.
message(paste("\nNumber of Nodes:",length(V(g))))
message(paste("\nNumber of Edges:",length(E(g))))

# data in root/tables:
# Save sif as a sheet in WASH_BioID_Results.xlsx
myfile <- file.path(root,"tables","WASH_BioID_Results.xlsx")
wb <- loadWorkbook(file = myfile)
addWorksheet(wb, sheetName = "PPIs")
writeData(wb,sheet=2,sif,rowNames=TRUE,colNames=TRUE)
saveWorkbook(wb, file=myfile, overwrite=TRUE)

# Done!
end <- Sys.time()
message(paste("\nCompleted analysis at:",end))
