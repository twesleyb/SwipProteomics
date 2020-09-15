#!/usr/bin/env Rscript

#' ---
#' title: Swip TMT Proteomics
#' description: generate protein networks
#' authors: Tyler W A Bradshaw
#' ---

## INPUT:
# tmt_protein in root/data

## OPTIONS:
os_keep = as.character(c(9606,10116,1090)) # keep ppis from human, rat, and mus.

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

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Imports.
suppressPackageStartupMessages({
	library(dplyr) # for manipulating data
	library(WGCNA) # for bicor function
	library(neten) # for network enhancement
	library(getPPIs) # for PPI data
	library(data.table) # for working with tables
})

# Project imports.
devtools::load_all()

#--------------------------------------------------------------------
## Create protein covariation network.
#--------------------------------------------------------------------

# Load the proteomics data.
data(tmt_protein)

# Cast to a data.matrix.
dm <- tmt_protein %>% as.data.table() %>%
	dcast(Sample ~ Accession, value.var="Intensity") %>%
	as.matrix(rownames=TRUE) %>% log2()

# Create correlation (adjacency) matrix.
message("\nGenerating protein co-variation matrix using bicor().")
adjm <- WGCNA::bicor(dm)

# Enhanced network.
message("\nPerforming network enhancement with to denoise network.")
ne_adjm <- neten::neten(adjm)

#--------------------------------------------------------------------
## Create PPI network.
#--------------------------------------------------------------------

# Load mouse PPIs.
message("\nCreating protein-protein interaction network.")
data(musInteractome)

# Collect all entrez cooresponding to proteins in our network.
proteins <- colnames(adjm)
entrez <- tmt_protein$Entrez[match(proteins,tmt_protein$Accession)]
names(proteins) <- entrez

# Collect PPIs among all proteins.
ppi_data <- musInteractome %>%
	filter(Interactor_B_Taxonomy %in% os_keep) %>%
	filter(Interactor_B_Taxonomy %in% os_keep) %>%
	filter(osEntrezA %in% entrez) %>% 
	filter(osEntrezB %in% entrez)

# Save to excel.
myfile <- file.path(root,"tables","Swip_TMT_Network_PPIs.xlsx")
write_excel(list("Network PPIs" = ppi_data),file=myfile)

# Create simple edge list (sif) and matrix with node attributes (noa).
sif <- ppi_data %>% select(osEntrezA, osEntrezB)
sif$uniprotA <- proteins[as.character(sif$osEntrezA)]
sif$uniprotB <- proteins[as.character(sif$osEntrezB)]

# Create igraph object from sif.
g <- graph_from_data_frame(sif[,c("uniprotA","uniprotB")], directed = FALSE)
g <- simplify(g)

# Extract as adjm.
ppi_adjm <- as.matrix(as_adjacency_matrix(g))

# Fill matrix.
all_proteins <- colnames(adjm)
missing <- all_proteins[all_proteins %notin% colnames(ppi_adjm)]
x <- matrix(nrow=dim(ppi_adjm)[1],ncol=length(missing))
colnames(x) <- missing
ppi_adjm <- cbind(ppi_adjm,x)
y <- matrix(nrow=length(missing),ncol=dim(ppi_adjm)[2])
rownames(y) <- missing
ppi_adjm <- rbind(ppi_adjm,y)
ppi_adjm <- ppi_adjm[all_proteins,all_proteins]
ppi_adjm[is.na(ppi_adjm)] <- 0

# Check, all should be the same.
c1 <- all(colnames(adjm) == colnames(ne_adjm)) 
c2 <- all(colnames(ne_adjm) == colnames(ppi_adjm))
if (!(c1 & c2)){ stop() }

# Number of edges and nodes.
n_edges <- sum(ppi_adjm[upper.tri(ppi_adjm)])
n_nodes <- ncol(ppi_adjm)

#--------------------------------------------------------------------
## Save the data.
#--------------------------------------------------------------------

message("\nSaving the data.")

# Save adjacency matrices as rda objects.
# To reduce the file size of a large matrix saved as a csv file, 
# the matrix is saved as an edge list, after removing the diagonal and 
# lower half of the matrix are removed. 
# This dataframe is saved as an rda object. It can be cast
# back into a N x N matrix with the convert_to_adjm function.

myfile <- file.path(root,"data","adjm.rda")
save_adjm_as_rda(round(adjm,5), myfile)

myfile <- file.path(root,"data","ne_adjm.rda")
save_adjm_as_rda(ne_adjm, myfile)

myfile <- file.path(root,"data","ppi_adjm.rda")
save_adjm_as_rda(ppi_adjm, myfile)

# Save adjm as csv.
adjm %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(file.path(root,"rdata","adjm.csv"))

# Save enhanced adjm as csv.
ne_adjm %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(file.path(root,"rdata","ne_adjm.csv"))

# Save ppi network as csv.
ppi_adjm %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(file.path(root,"rdata","ppi_adjm.csv"))

# Save norm_protein as matrix. 
norm_protein <- tmt_protein %>% as.data.table() %>%
	dcast(Accession ~ Sample, value.var = "Intensity") %>%
	fwrite(file.path(root,"rdata","norm_protein.csv"))
