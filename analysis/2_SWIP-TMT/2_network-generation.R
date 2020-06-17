#!/usr/bin/env Rscript

## INPUT:
# tmt_protein in root/data

## ANALYSIS OPTIONS:
os_keep = as.character(c(9606,10116,1090)) # keep ppis from human, rat, and mus.

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------

start <- Sys.time()
message(paste("Starting analysis at:",start))

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(WGCNA)
	library(neten)
	library(getPPIs)
	library(data.table)
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

# To shrink the size of the adjm, melt it to an edge list, remove diag, and
# lower.tri, coerce to a matrix and save as rda.

# Define a function that casts data back into a matrix.
# This is packaged with the data.
convert_to_adjm <- function(edge_df) {
	suppressPackageStartupMessages({ library(data.table) })
	message("Casting edge list into adjacency matrix.")
	dm <- edge_df %>% as.data.table() %>% 
		dcast.data.table(Var1 ~ Var2, value.var = "value") %>% 
		as.matrix(rownames="Var1")
	return(dm)
}

# Define a function that saves adjm as an edge list 
# data.frame as an rda object.
save_adjm_as_rda <- function(adjm,file) {
	diag(adjm) <- 0 # ZERO is smaller than NA
	adjm[lower.tri(adjm)] <- 0
	edges <- reshape2::melt(adjm)
	save(edges,file=file,version=2)
}

# Save adjmatrices.
save_adjm_as_rda(ne_adjm,file="ne_adjm.rda")

# load adjm.
load("ne_adjm.rda")
adjm <- convert_to_adjm(edges)

adjm <- convert(adjm)


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

#--------------------------------------------------------------------
## Save the data.
#--------------------------------------------------------------------

message("\nSaving the data.")

# Save adjm as csv and rda.
adjm %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(file.path(root,"rdata","adjm.csv"))
save(adjm,file=file.path(root,"rdata","adjm.rda"), version=2)

# Save enhanced adjm as csv and rda.
ne_adjm %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(file.path(root,"rdata","ne_adjm.csv"))
save(ne_adjm,file=file.path(root,"rdata","ne_adjm.rda"), version=2)

# Save ppi network as csv and rda.
ppi_adjm %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(file.path(root,"rdata","ppi_adjm.csv"))
save(ppi_adjm,file=file.path(root,"rdata","ppi_adjm.rda"), version=2)

# Save norm_protein as matrix. 
# FIXME: do we need this?
norm_protein <- tmt_protein %>% as.data.table() %>%
	dcast(Accession ~ Sample, value.var = "Intensity") %>%
	fwrite(file.path(root,"rdata","norm_protein.csv"))

# Done!
end <- Sys.time()
message(paste("\nCompleted analysis at:",end))
