#!/usr/bin/env Rscript

# title: Swip TMT Proteomics
# description: generate protein networks
# authors: Tyler W A Bradshaw


## INPUTs ----------------------------------------------------------------------
root = "~/projects/SwipProteomics"

# * msstats_prot
# * msstats_gene_map
# * musInteractome from twesleyb/getPPIs
# * 

## OPTIONs ---------------------------------------------------------------------
os_keep = as.character(c(9606,10116,1090)) # keep ppis from human, rat, and mus.


## Prepare the workspace ------------------------------------------------------

# Load renv
root <- getrd()
renv::load(root,quiet=TRUE)

# Imports
suppressPackageStartupMessages({
	library(dplyr) # for manipulating data
	library(WGCNA) # for bicor function
	library(neten) # for network enhancement
	suppressWarnings({ library(getPPIs) }) # for PPI data 
	library(data.table) # for working with tables
})

# Project imports
devtools::load_all()

# Load the proteomics data
data(msstats_prot)
data(msstats_gene_map)
data(musInteractome)


## Create protein covariation network -----------------------------------------

# Cast protein data into a data.matrix. No need to log2 transform.
# MSstats has done this already. Rownames are unique for each Sample.
dm <- msstats_prot %>% as.data.table() %>%
	dcast(interaction(Mixture,BioFraction,Genotype) ~ Protein, value.var="Abundance") %>%
	as.matrix(rownames=TRUE)

# drop rows with any missing values
# NOTE: the data was transposed for bicor
out <- apply(dm,2,function(x) any(is.na(x)))
subdm <- dm[,!out]
warning(formatC(sum(out),big.mark=",")," proteins with any missing values were removed.")

# Create correlation (adjacency) matrix
message("\nGenerating protein co-variation network.")
adjm <- WGCNA::bicor(subdm)

# Enhanced network
# NOTE: this can take a couple minutes
message("\nPerforming network enhancement.")
ne_adjm <- neten::neten(adjm)


## Create PPI network ---------------------------------------------------------
message("\nCreating protein-protein interaction network.")

# Collect all entrez cooresponding to proteins in our network.
proteins <- colnames(adjm)
idx <- match(proteins,gene_map$uniprot)
entrez <- gene_map$entrez[idx]
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
# If a protein was unconnected it gets dropped by igraph. Add it back so that
# all the matrices look the same.

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

message("\nPPI graph:")
data.table("Edges" = formatC(n_edges, format = "d", big.mark=","), 
	   "Nodes" = formatC(n_nodes, format = "d", big.mark=",")) %>% knitr::kable()


## Save the data --------------------------------------------------------------
# save data in root/rdata

myfile <- file.path(root,"rdata","adjm.rda")
save(adjm,file=myfile,version=2)

myfile <- file.path(root,"rdata","ne_adjm.rda")
save(ne_adjm,file=myfile,version=2)

myfile <- file.path(root,"rdata","ppi_adjm.rda")
save(ppi_adjm,file=myfile,version=2)

# Save adjm as csv.
adjm %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(file.path(root,"rdata","adjm.csv"))

# Save enhanced adjm as csv.
ne_adjm %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(file.path(root,"rdata","ne_adjm.csv"))

# Save ppi network as csv.
ppi_adjm %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(file.path(root,"rdata","ppi_adjm.csv"))
