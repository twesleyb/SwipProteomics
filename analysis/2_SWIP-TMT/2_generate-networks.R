#!/usr/bin/env Rscript

# title: Swip TMT Proteomics
# description: generate protein networks
# authors: Tyler W A Bradshaw


## INPUTs ----------------------------------------------------------------------
root = "~/projects/SwipProteomics"

# * data(msstats_prot)
# * data(gene_map)
# * data(musInteractome) from twesleyb/getPPIs

## OPTIONs ---------------------------------------------------------------------
os_keep = as.character(c(9606,10116,1090)) # keep ppis from human, rat, and mus.

## OUTPUTs in root/rdata:

# * adjm.csv
# * ne_adjm.csv --> leidenalg clustering
# * ppi_adjm.csv --> leidenalg clustering

# * adjm.rda
# * ne_adjm.rda
# * ppi_adjm.rda

# * norm_prot in root/data

## functions ------------------------------------------------------------------

mkdir <- function(...,warn=TRUE,report=FALSE) {
	# create a new directory
	newdir <- file.path(...)
	if (warn & dir.exists(newdir)) {
		warning("dir exists")
	} else if (!dir.exists(newdir)) { 
		dir.create(newdir)
		if (report) {
		message(paste("created",newdir))
		}
	}
}


## Prepare the workspace ------------------------------------------------------

# Load renv
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
data(gene_map)
data(msstats_prot)
data(musInteractome)

# data for output tables
tabsdir <- file.path(root,"tables")
if (!dir.exists(tabsdir)) { mkdir(tabsdir) }


## Create protein covariation network -----------------------------------------

# Cast protein data into a data.matrix. No need to log2 transform.
# MSstats has done this already. Rownames are unique for each Sample.
dm <- msstats_prot %>% as.data.table() %>%
	dcast(interaction(Mixture,BioFraction,Genotype) ~ Protein, 
	      value.var="Abundance") %>%
	as.matrix(rownames=TRUE)

# normalize such that sum is 1? normalize each protein to its max!
norm_dm = dm
#norm_dm = apply(dm,2,function(x) x/max(x,na.rm=TRUE))

# drop rows with any missing values
# NOTE: the data was transposed for bicor
out <- apply(norm_dm,2,function(x) any(is.na(x)))
subdm <- norm_dm[,!out]
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


## cast data into a matrix for permutation testing -----------------------------------

# load the data and save as norm protein for permutation testing
norm_prot <- msstats_prot %>% as.data.table() %>%
	dcast(Protein ~ interaction(Mixture,Channel,Genotype),value.var="Abundance") %>%
	na.omit() %>% as.matrix(rownames="Protein")


## Save the data --------------------------------------------------------------
# all data is saved in root/rdata

# adjm
myfile <- file.path(root,"rdata","adjm.rda")
save(adjm,file=myfile,version=2)
message("\nSaved ", basename(myfile)," in ", dirname(myfile))

# ne_adjm
myfile <- file.path(root,"rdata","ne_adjm.rda")
save(ne_adjm,file=myfile,version=2)
message("\nSaved ", basename(myfile)," in ", dirname(myfile))

# ppi_adjm
myfile <- file.path(root,"rdata","ppi_adjm.rda")
save(ppi_adjm,file=myfile,version=2)
message("\nSaved ", basename(myfile)," in ", dirname(myfile))

# adjm
adjm %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(file.path(root,"rdata","adjm.csv"))
message("\nSaved ", basename(myfile)," in ", dirname(myfile))

# ne_adjm
ne_adjm %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(file.path(root,"rdata","ne_adjm.csv"))
message("\nSaved ", basename(myfile)," in ", dirname(myfile))

# ppi_adjm
ppi_adjm %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(file.path(root,"rdata","ppi_adjm.csv"))
message("\nSaved ", basename(myfile)," in ", dirname(myfile))

# norm_prot
myfile <- file.path(root,"data","norm_prot.rda")
save(norm_prot,file=myfile,version=2)
message("\nSaved ", basename(myfile)," in ", dirname(myfile))
