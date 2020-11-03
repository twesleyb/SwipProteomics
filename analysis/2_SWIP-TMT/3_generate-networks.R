#!/usr/bin/env Rscript

# title: Swip TMT Proteomics
# description: generate protein networks
# authors: Tyler W A Bradshaw


## INPUTs ----------------------------------------------------------------------
root <- "~/projects/SwipProteomics"

# * data(gene_map)
# * data(msstats_prot)
# * data(musInteractome) from twesleyb/getPPIs

## HELLO ## 
# this script seems to be a hub. it generates adjm, ne_adjm, and norm_prot --
# all used in the network building process

## OPTIONs ---------------------------------------------------------------------
# species to keep when building the ppi graph:
os_keep <- as.character(c(9606, 10116, 1090)) # human, rat, and mus.

## OUTPUT

## saved in root/rdata:
# * adjm.csv
# * ne_adjm.csv --> leidenalg clustering
# * ppi_adjm.csv

## saved in root/rdata:
# * adjm.rda
# * ne_adjm.rda
# * ppi_adjm.rda

## saved in root/data
# * norm_prot in root/data --> module preservation


## functions ------------------------------------------------------------------

mkdir <- function(..., warn = TRUE, report = FALSE) {
  # create a new directory
  newdir <- file.path(...)
  if (warn & dir.exists(newdir)) {
    warning("dir exists")
  } else if (!dir.exists(newdir)) {
    dir.create(newdir)
    if (report) {
      message(paste("created", newdir))
    }
  }
}


## Prepare the workspace ------------------------------------------------------

# Load renv
renv::load(root, quiet = TRUE)

# Imports
suppressPackageStartupMessages({
  library(dplyr) # for manipulating data
  library(neten) # for network enhancement
  library(getPPIs) # for PPI data
  library(data.table) # for working with tables
})

# Project imports
devtools::load_all()

# Load the proteomics data
data(gene_map)
data(poor_prots)
data(msstats_prot)
data(musInteractome)

# NOTE: We will remove proteins with poor fit before creating network

# data for output tables
tabsdir <- file.path(root, "tables")
if (!dir.exists(tabsdir)) {
  mkdir(tabsdir)
}


## Create protein covariation network -----------------------------------------

# Cast protein data into a data.matrix. No need to log2 transform.
# MSstats has done this already. Rownames are unique for each Sample.
# Use norm_Abundance (the adjusted protein-level data)

# NOTE: we also drop proteins with poor fit
# this is the data used to create adjm, ne_adjm, and ppi_adjm

message("\nProteins with poor fit are not used to build network.")

dm <- msstats_prot %>% 
  # rm poor fits
  filter(Protein %notin% poor_prots) %>%
  as.data.table() %>%
  dcast(interaction(Mixture, BioFraction, Genotype) ~ Protein,
    value.var = "norm_Abundance"
  ) %>%
  as.matrix(rownames = TRUE)

# drop rows with any missing values
out <- apply(dm, 2, function(x) any(is.na(x)))
subdm <- dm[, !out]
warning(formatC(sum(out), big.mark = ","), " proteins with any missing values were removed.")

# status
knitr::kable(cbind(samples = dim(subdm)[1], proteins = dim(subdm)[2]))

# Create correlation (adjacency) matrix
message("\nGenerating protein co-variation network.")
adjm <- cor(subdm, method = "pearson")

# Enhanced network
# NOTE: this can take a couple minutes
# NOTE: we really need to generate a visualization...
message("\nPerforming network enhancement.")
ne_adjm <- neten::neten(adjm,alpha=0.9,diffusion=1) 


## Create PPI network ---------------------------------------------------------
message("\nCreating protein-protein interaction network.")

# NOTE: the PPI network is not used in the clustering of the data.
# But it can be.

# Collect all entrez cooresponding to proteins in our network.
proteins <- colnames(adjm)
idx <- match(proteins, gene_map$uniprot)
entrez <- gene_map$entrez[idx]
names(proteins) <- entrez

# Collect PPIs among all proteins.
ppi_data <- musInteractome %>%
  filter(Interactor_B_Taxonomy %in% os_keep) %>%
  filter(Interactor_B_Taxonomy %in% os_keep) %>%
  filter(osEntrezA %in% entrez) %>%
  filter(osEntrezB %in% entrez)


# Create simple edge list (sif) and matrix with node attributes (noa)
sif <- ppi_data %>% select(osEntrezA, osEntrezB)
sif$uniprotA <- proteins[as.character(sif$osEntrezA)]
sif$uniprotB <- proteins[as.character(sif$osEntrezB)]

# Create igraph object from sif
g <- graph_from_data_frame(sif[, c("uniprotA", "uniprotB")], directed = FALSE)
g <- simplify(g)

# extract as adjm
ppi_adjm <- as.matrix(as_adjacency_matrix(g))

## munge to fill matrix
# If a protein was unconnected it gets dropped by igraph. Add it back so that
# all the matrices look the same.
all_proteins <- colnames(adjm)
missing <- all_proteins[all_proteins %notin% colnames(ppi_adjm)]
x <- matrix(nrow = dim(ppi_adjm)[1], ncol = length(missing))
colnames(x) <- missing
ppi_adjm <- cbind(ppi_adjm, x)
y <- matrix(nrow = length(missing), ncol = dim(ppi_adjm)[2])
rownames(y) <- missing
ppi_adjm <- rbind(ppi_adjm, y)
ppi_adjm <- ppi_adjm[all_proteins, all_proteins]
ppi_adjm[is.na(ppi_adjm)] <- 0

# Check, all should be the same.
c1 <- all(colnames(adjm) == colnames(ne_adjm))
c2 <- all(colnames(ne_adjm) == colnames(ppi_adjm))
if (!(c1 & c2)) {
  stop()
}

# Number of edges and nodes.
n_edges <- sum(ppi_adjm[upper.tri(ppi_adjm)])
n_nodes <- ncol(ppi_adjm)

message("\nPPI graph:")
data.table(
  "Edges" = formatC(n_edges, format = "d", big.mark = ","),
  "Nodes" = formatC(n_nodes, format = "d", big.mark = ",")
) %>% knitr::kable()


## cast data into a matrix for permutation testing -----------------------------------

# load the data and save as norm protein for permutation testing
# NOTE: we use normalized Abundance!
# NOTE: we rm poor fits
# NOTE: prots with missing vals are dropped!
norm_prot <- msstats_prot %>% as.data.table() %>% 
  filter(Protein %notin% poor_prots) %>%
  dcast(Protein ~ interaction(Mixture, Channel, Genotype), value.var = "norm_Abundance") %>%
  na.omit() %>%
  as.matrix(rownames = "Protein")


## Save the data --------------------------------------------------------------
# data > 50 MB are saved in root/rdata

# adjm
myfile <- file.path(root, "rdata", "adjm.rda")
save(adjm, file = myfile, version = 2)
message("\nSaved ", basename(myfile), " in ", dirname(myfile))

# ne_adjm
myfile <- file.path(root, "rdata", "ne_adjm.rda")
save(ne_adjm, file = myfile, version = 2)
message("\nSaved ", basename(myfile), " in ", dirname(myfile))

# ppi_adjm
myfile <- file.path(root, "rdata", "ppi_adjm.rda")
save(ppi_adjm, file = myfile, version = 2)
message("\nSaved ", basename(myfile), " in ", dirname(myfile))

# adjm
adjm %>%
  as.data.table(keep.rownames = "Accession") %>%
  fwrite(file.path(root, "rdata", "adjm.csv"))
message("\nSaved ", basename(myfile), " in ", dirname(myfile))

# ne_adjm
ne_adjm %>%
  as.data.table(keep.rownames = "Accession") %>%
  fwrite(file.path(root, "rdata", "ne_adjm.csv"))
message("\nSaved ", basename(myfile), " in ", dirname(myfile))

# ppi_adjm
ppi_adjm %>%
  as.data.table(keep.rownames = "Accession") %>%
  fwrite(file.path(root, "rdata", "ppi_adjm.csv"))
message("\nSaved ", basename(myfile), " in ", dirname(myfile))

# norm_prot
myfile <- file.path(root, "data", "norm_prot.rda")
save(norm_prot, file = myfile, version = 2)
message("\nSaved ", basename(myfile), " in ", dirname(myfile))
