#!/usr/bin/env Rscript

#' ---
#' title: Swip Proteomics
#' description:
#' authors: Tyler W Bradshaw
#' ---

## Options:
alpha = 0.1 # Significance threshold.

## R options.
options(renv.config.synchronized.check = FALSE) # skip renv::check(repo).
options(renv.settings.snapshot.type = "simple") # use simple renv::snapshot.

#FIXME: ADD ADDITIONAL MODULE METRICS

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. 

start <- Sys.time()
message(paste("Starting analysis at:",start))

# Load renv -- use renv::load NOT activate!
root <- getrd()
renv::load(root,quiet=TRUE) # NOTE: getrd is a f(x) in .Rprofile.

# Load required packages and functions.
suppressPackageStartupMessages({
	library(dplyr) # For manipulating data.
#	library(ggplot2) # For making plots.
	library(edgeR) # For statistical analysis.
	library(data.table) # For working with tables.
})

# Load additional functions.
suppressWarnings({ devtools::load_all() })

# Project directories:
datadir <- file.path(root, "data")
fontdir <- file.path(root, "fonts")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")
downdir <- file.path(root, "downloads")
figsdir <- file.path(root, "figs","Modules")

# Global plotting settings.
ggtheme()
set_font("Arial", font_path = fontdir)

# Load the data.
data(tmt_protein)

# Load the graph partition.
data(partition)

# Calculate the number of proteins per module.
module_sizes <- sapply(split(partition,partition), length)

# Annotate data with module membership.
tmt_protein$Module <- partition[tmt_protein$Accession]

# Cast data into a dm, summaraze all proteins in a module.
dm <- tmt_protein %>% group_by(Module, Genotype, Fraction) %>% 
	summarize(Sum.Intensity=sum(Intensity),.groups="drop") %>% 
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

# Default is last contrast?
qlf <- glmQLFTest(fit)

# Collect results.
glm_results <- topTags(qlf,n=Inf,sort.by="p.value")$table %>%
	as.data.table(keep.rownames="Module")
glm_results$PAdjust <- p.adjust(glm_results$PValue,method="bonferroni")

# Fix logCPM column.
idy <- which(colnames(glm_results)=="logCPM")
colnames(glm_results)[idy] <- "PercentWT"
glm_results$PercentWT <- 2^glm_results$logFC

# Number of nodes per module:
n <- module_sizes[as.character(glm_results$Module)]
glm_results <- tibble::add_column(glm_results,Nodes=n,.after="Module")

# Number of significant modules.
nsig <- sum(glm_results$PAdjust < alpha)
message(paste0("\nNumber of significant ",
	      "(p-adjust < ", alpha,") ",
	      "modules: ", nsig,"."))

# Sig results:
glm_results %>% filter(PAdjust < alpha) %>%
	knitr::kable()

# Save results.
myfile <- file.path(tabsdir,"Module_GLM_results.csv")
fwrite(glm_results,myfile)

# Done!
end <- Sys.time()
message(paste("\nCompleted analysis at:",end))
message(paste("Elapsed time:",
	      round(difftime(end,start,units="mins"),2),"minutes."))
