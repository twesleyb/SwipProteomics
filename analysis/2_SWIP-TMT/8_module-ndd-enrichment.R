#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

## Optional parameters:
FDR_alpha = 0.05 # Significance threshold.
single_group = FALSE

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(geneLists)
	library(data.table)
})

# Project functions.
suppressWarnings({ devtools::load_all() })

# Directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

#--------------------------------------------------------------------
## Load the data.
#--------------------------------------------------------------------

# Load the data from root/data.
data(tmt_protein)

# Load the graph partition:
data(partition)

# Load gene map.
data(gene_map)

#--------------------------------------------------------------------
## Build a NDD-gene reference collection.
#--------------------------------------------------------------------

# Use the geneLists() function from my geneLists() package to collect
# DBD-associated gene lists.
dataset_list <- list("ALS" = geneLists("als"),
		     "PD" = geneLists("pd"),
		     "ALZ" = geneLists("alz"),
		     "MS" = geneLists("ms"),
		     "NDD"= geneLists("NDD"))

# Coerce to a named vector.
datasets <- unlist(dataset_list,use.names=FALSE)
datasets <- setNames(rep(names(dataset_list),times=sapply(dataset_list,length)),
		     nm=datasets)

# Load the data from twesleyb/geneLists.
data(list=names(datasets))

# Define a function that converts a given geneList into a dataframe.
geneList_2df <- function(geneList){
	gene_list <- eval(parse(text=geneList)) # dynamically set gene_list
	max_list_length <- max(sapply(gene_list,length))
	buffer_NA <- lapply(gene_list,function(x) {
				    rep(NA, each = max_list_length - length(x))
		     })
	df <- as.data.frame(mapply(c,gene_list,buffer_NA))
	df$dataset <- geneList # Name of geneList
	# Tidy-up.
	geneList_df <- reshape2::melt(df, id.var="dataset",value.name="entrez",
				      variable.name="disorder_association")
	return(geneList_df %>% filter(!is.na(entrez))) # Drop NA.
}

# Coerce geneLists to dataframes.
suppressWarnings({ data_list <- lapply(names(datasets),geneList_2df) })
names(data_list) <- names(datasets)

# Collect as a single data.table. 
# Ignore warnings about coercing variable to character.
dt <- suppressWarnings({ bind_rows(data_list) })
dt$disorder_association <- as.character(dt$disorder_association)

# Clean-up disorder associations.
idx <- grepl("amyotrophic lateral sclerosis",dt$disorder_association)
if (sum(idx) > 0) { dt$disorder_association[idx] <- "ALS" }

idx <- grepl("multiple sclerosis",dt$disorder_association)
if (sum(idx) > 0) { dt$disorder_association[idx] <- "MS" }

idx <- grepl("alzheimer",dt$disorder_association)
if (sum(idx) > 0) { dt$disorder_association[idx] <- "ALZ" }

idx <- grepl("parkinson disease",dt$disorder_association)
if (sum(idx) > 0) { dt$disorder_association[idx] <- "PD" }

idx <- grepl("huntington disease",dt$disorder_association)
if (sum(idx) > 0) { dt$disorder_association[idx] <- "HD" }

# Summarize number of unique genes per category.
# Summary.
message(paste("Compiled",
	formatC(length(unique(dt$entrez)),big.mark=","),
	"genes associated with neurodegenerative disorders."))
dt %>% group_by(disorder_association) %>% 
	summarize(nGenes=formatC(length(unique(entrez)),big.mark=","),
		  .groups="drop") %>% 
	knitr::kable()

# Single group: 
if (single_group) { 
	message("Combined NDD genes into single group!")
	gene_list <- list("NDD" = unique(dt$entrez))
} else {
	# Multiple groups:
	gene_list <- split(dt$entrez,dt$disorder_association)
}

# Create gene sets -- utilizes anRichment library.
geneSets <- lapply(names(gene_list),function(x) {
			    createGeneSet(gene_list[[x]],x,
			 description="NDD genes compiled from 6 databases",
			 data_source="this script")})

# Combine gene sets into single collection:
PLgroup <- newGroup(name="Compiled NDD Genes",
		    description="NDD genes compiled from 6 databases.",
		    source = "this script")
NDDcollection <- newCollection(dataSets=geneSets,groups=list(PLgroup))

#---------------------------------------------------------------------
## Protein NDD annotations.
#---------------------------------------------------------------------

# Collect protein NDD annotations for proteins in the data
dt$uniprot <- gene_map$uniprot[match(dt$entrez,gene_map$entrez)]
dt$symbol <- gene_map$symbol[match(dt$entrez,gene_map$entrez)]
tmp_dt <- dt %>% filter(!is.na(uniprot)) %>% group_by(uniprot) %>% 
	summarize(NDD = paste(unique(disorder_association),collapse="|"))
NDD_proteins <- tmp_dt$NDD
names(NDD_proteins) <- tmp_dt$uniprot

# Save as rda.
save(NDD_proteins,file=file.path(datadir,"NDD_proteins.rda"),version=2)

#---------------------------------------------------------------------
## Module NDD enrichment analysis.
#---------------------------------------------------------------------

# Perform enrichment analysis
# NOTE: Background is all genes in partition.
message("\nPerforming NDD enrichment analysis for every module...")
NDD_results <- moduleGOenrichment(partition, gene_map, 
				  GOcollection = NDDcollection,
				  partition.ids = "uniprot")

# Number of significant terms per module.
nSig <- sapply(NDD_results, function(x) sum(x$FDR < FDR_alpha))

# Modules with any sig terms:
sigModules <- names(which(nSig > 0))
nSig_Modules <- length(sigModules)

# Status.
message(paste("\nNumber of modules with any significant DBD-gene enrichment:",
	    nSig_Modules))

# Get sig results.
temp_list <- lapply(NDD_results[sigModules],function(x) {
	       as.data.table(x) %>% filter(FDR < FDR_alpha) %>%
		       dplyr::select(class,dataSetName,nCommonGenes,
			      nCommonGenes,pValue,FDR,enrichmentRatio)
				  })

# Significant modules.
message("\nSignificant modules:")
lapply(temp_list,knitr::kable)

#--------------------------------------------------------------------
## Save results.
#--------------------------------------------------------------------

# Save as excel workbook.
message("\nSaving data...")
myfile <- file.path(tabsdir,"Swip_TMT_Module_NDD_Results.xlsx")
write_excel(NDD_results,myfile)
