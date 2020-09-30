#!/usr/bin/env Rscript

# title: SwipProteomics
# description: 
# author: twab <twesleyb10@gmail.com>
# os: windows linux subsystem (WSL)

## INPUT ----------------------------------------------------------------------
# specify project's root directory
ROOT = "~/projects/SwipProteomics"

## OPTIONS --------------------------------------------------------------------
# keep PPIs from human, rat, and mus
os_keep = as.character(c(9606,10116,1090)) 

## FUNCTIONS -----------------------------------------------------------------

mkdir <- function(...) {
	# create a new directory
	newdir <- file.path(...)
	if (!dir.exists(newdir)) { 
		dir.create(newdir)
		message(paste("created",newdir))
	}
}


## Prepare the working environment ----------------------------------------------

# directory for output tables:
#tabsdir <- file.path(ROOT,"tables"); mkdir(tabsdir)
#datadir <- file.path(ROOT,"data"); mkdir(datadir)


## Prepare the R environment --------------------------------------------------

# load renv
renv::load(ROOT,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(dplyr) # for manipulating data
	library(WGCNA) # for bicor function
	library(neten) # for network enhancement
	suppressWarnings( library(getPPIs) )# for PPI data
	library(data.table) # for working with tables
})

# load functions in root/R
devtools::load_all(ROOT)

# load the data in root/data
data(samples) # sample metadata
data(msstats_prot) # load the MSstats processed data
data(msstats_gene_map) # big gene_map 


## Create protein covariation network -----------------------------------------

data_prot <- msstats_prot

# Cast protein data into a matrix
dm <- data_prot %>% as.data.table() %>%
	dcast(Protein ~ Mixture + Channel + Condition, value.var="Abundance") %>% 
	as.matrix(rownames="Protein") # coerce to matrix

# use complete cases; drop rows with na
idx <- apply(dm,1,function(x) any(is.na(x)))
dm <- dm[!idx,]

# create correlation (adjacency) matrix
message("\nGenerating protein co-variation matrix using bicor().")
adjm <- WGCNA::bicor(t(dm))

# perform network enhancement
message("\nPerforming network enhancement with to denoise network.")
ne_adjm <- neten::neten(adjm)


## Create PPI network ---------------------------------------------------------

# Load mouse PPIs
message("\nCreating protein-protein interaction network.")
data(musInteractome) # from twesleyb/getPPIs

# collect all entrez cooresponding to proteins in our network
proteins <- colnames(adjm)
entrez <- gene_map$entrez[match(proteins,gene_map$uniprot)]
names(proteins) <- entrez

# collect PPIs among entrez genes cooresponding to all proteins in adjm
ppi_data <- musInteractome %>%
	filter(Interactor_B_Taxonomy %in% os_keep) %>%
	filter(Interactor_B_Taxonomy %in% os_keep) %>%
	filter(osEntrezA %in% entrez) %>% 
	filter(osEntrezB %in% entrez)

# save to excel workbook
#myfile <- file.path(root,"tables","Swip_TMT_Network_PPIs.xlsx")
#write_excel(list("Network PPIs" = ppi_data),file=myfile)

# Create simple edge list (sif) and matrix with node attributes (noa).
sif <- ppi_data %>% dplyr::select(osEntrezA, osEntrezB)
sif$uniprotA <- proteins[as.character(sif$osEntrezA)]
sif$uniprotB <- proteins[as.character(sif$osEntrezB)]

# create igraph object from sif
g <- graph_from_data_frame(sif[,c("uniprotA","uniprotB")], directed = FALSE)
g <- simplify(g)

# extract as adjm
ppi_adjm <- as.matrix(as_adjacency_matrix(g))

# fill matrix -- add unconnected proteins
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

# check, all should be the same
c1 <- all(colnames(adjm) == colnames(ne_adjm)) 
c2 <- all(colnames(ne_adjm) == colnames(ppi_adjm))
if (!(c1 & c2)){ stop() }

# number of edges and nodes:
n_edges <- sum(ppi_adjm[upper.tri(ppi_adjm)])
n_nodes <- ncol(ppi_adjm)


## Save the data --------------------------------------------------------------

message("\nSaving the data.")

# Save adjacency matrices as rda objects.
# To reduce the file size of a large matrix saved as a csv file, 
# the matrix is saved as an edge list, after removing the diagonal and 
# lower half of the matrix are removed. 
# This dataframe is saved as an rda object. It can be cast
# back into a N x N matrix with the convert_to_adjm function.

## NOTE: FILES ARE TOO BIG, SAVE IN RDATA

myfile <- file.path(ROOT,"rdata","adjm.rda")
save(adjm,file=myfile,version=2)

myfile <- file.path(ROOT,"rdata","ne_adjm.rda")
save(ne_adjm,file=myfile,version=2)

myfile <- file.path(ROOT,"rdata","ppi_adjm.rda")
save(ppi_adjm,file=myfile,version=2)

# Save adjm as cvv
adjm %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(file.path(ROOT,"rdata","adjm.csv"))

# Save enhanced adjm as csv.
ne_adjm %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(file.path(ROOT,"rdata","ne_adjm.csv"))

# Save ppi network as csv.
ppi_adjm %>% as.data.table(keep.rownames="Accession") %>%
	fwrite(file.path(ROOT,"rdata","ppi_adjm.csv"))

# map colnames of dm (normalize protein data) to Sample names 
samples$label <- paste(paste(gsub("Exp","M",samples$Experiment),
			     samples$Channel,sep="_"),
		       paste(samples$Treatment,samples$Fraction,sep="."),
		       sep="_")
idx <- match(colnames(dm),samples$label)
colnames(dm) <- samples$Sample[idx]

# Save norm_protein as matrix. 
norm_protein <- dm %>% as.data.table(keep.rownames="Accession")
myfile <- file.path(ROOT,"rdata","norm_protein.csv")
fwrite(norm_protein, myfile)
