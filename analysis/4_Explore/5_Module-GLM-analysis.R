#!/usr/bin/env Rscript

#' ---
#' title: Swip TMT Proteomics
#' description: module level statistical analysis
#' authors: Tyler W A Bradshaw
#' ---

## Options:
BF_alpha = 0.05 # Significance threshold.

#---------------------------------------------------------------------
## Misc function - getrd().
#---------------------------------------------------------------------

getrd <- function(here=getwd(), dpat= ".git") {
	# Get the repository's root directory.
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

# Load renv -- use renv::load NOT activate!
root <- getrd()
renv::load(root,quiet=TRUE)

# Load required packages and functions.
suppressPackageStartupMessages({
	library(dplyr) # For manipulating data.
	library(edgeR) # For statistical analysis.
	library(data.table) # For working with tables.
})

# Load additional functions.
suppressWarnings({ devtools::load_all() })

# Project directories:
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")
suppdir <- file.path(root, "supplement")
figsdir <- file.path(root, "figs","Modules")

# If necessary, create dir for figs.
if (!dir.exists(figsdir)){ dir.create(figsdir, recursive = TRUE) }

# Load the data.
data(tmt_protein)

# Load the gene map.
data(gene_map)

# Load WASH BioID results.
data(wash_interactome)
wash_prots <- unique(wash_interactome$Accession)

# Load the graph partition.
data(partition)

# Load the network.
data(ne_adjm)
ne_adjm <- convert_to_adjm(edges)

#---------------------------------------------------------------------
## Perform module level statistical analysis.
#---------------------------------------------------------------------

# Calculate the number of proteins per module.
module_sizes <- sapply(split(partition,partition), length)

# Annotate data with module membership.
tmt_protein$Module <- partition[tmt_protein$Accession]

# Cast data into a dm, summarize (sum) all proteins in a module.
dm <- tmt_protein %>% group_by(Module, Genotype, Fraction) %>% 
	dplyr::summarize(Sum.Intensity=sum(Intensity),.groups="drop") %>% 
	as.data.table() %>%
	dcast(Module ~ Fraction + Genotype,value.var="Sum.Intensity") %>%
	as.matrix(rownames=TRUE)

# Create dge object.
dge <- DGEList(counts=dm)

# Sample to group mapping.
groups <- rownames(dge$samples)
dge$samples$fraction <- as.factor(sapply(strsplit(groups,"_"),"[",1))
dge$samples$genotype <- as.factor(sapply(strsplit(groups,"_"),"[",2))
levels(dge$samples$genotype) <- c("WT","MUT")

# Create design matrix.
design <- model.matrix( ~ fraction + genotype, data=dge$samples)

# Flip sign of comparison.
design[,"genotypeMUT"] <- abs(design[,"genotypeMUT"]-1)

# Estimate dispersion.
dge <- estimateDisp(dge, design, robust=TRUE)

# Fit a model.
fit <- glmQLFit(dge,design,robust=TRUE)

# Default comparison is last contrast: genotype.
qlf <- glmQLFTest(fit)

x = log2(qlf$fitted.values["19",])
x/max(x)

# Collect results.
glm_results <- topTags(qlf,n=Inf,sort.by="p.value")$table %>%
	as.data.table(keep.rownames="Module")

# Add deviance -- a measure of the goodness of fit of the model.
glm_results$Deviance <- qlf$deviance

# Drop M0.
glm_results <- glm_results %>% filter(Module != 0)

# Adjust p-values for n module comparisons.
glm_results$PAdjust <- p.adjust(glm_results$PValue,method="bonferroni")

# Fix logCPM column -- > convert to percentWT.
idy <- which(colnames(glm_results)=="logCPM")
colnames(glm_results)[idy] <- "PercentWT"
glm_results$PercentWT <- 2^glm_results$logFC

# Number of nodes per module:
n <- module_sizes[as.character(glm_results$Module)]
glm_results <- tibble::add_column(glm_results,Nodes=n,.after="Module")

# Number of significant modules.
message(paste("\nTotal number of modules:",length(unique(partition)) -1))
nsig <- sum(glm_results$PAdjust < BF_alpha)
message(paste0("\nNumber of significant ",
	      "(p-adjust < ", BF_alpha,") ",
	      "modules: ", nsig,"."))

# Pretty print sig results:
glm_results %>% filter(PAdjust < BF_alpha) %>%
	knitr::kable()

