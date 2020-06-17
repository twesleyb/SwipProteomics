#!/usr/bin/env Rscript

#' ---
#' title: Swip Proteomics
#' description:
#' authors: Tyler W Bradshaw
#' ---

## Options:
FDR_alpha = 0.1 # Significance threshold.

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

# Load WASH BioID results.
data(wash_interactome)
wash_prots <- unique(wash_interactome$Accession)

# Load the graph partition.
data(partition)

# Load the networks.
data(adjm)
adjm <- convert_to_adjm(edges)

data(ne_adjm)
ne_adjm <- convert_to_adjm(edges)

# Calculate the number of proteins per module.
module_sizes <- sapply(split(partition,partition), length)

# Annotate data with module membership.
tmt_protein$Module <- partition[tmt_protein$Accession]

# Cast data into a dm, summarize all proteins in a module.
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

# Collect results.
glm_results <- topTags(qlf,n=Inf,sort.by="p.value")$table %>%
	as.data.table(keep.rownames="Module")

# Adjuste p-values.
glm_results$PAdjust <- p.adjust(glm_results$PValue,method="bonferroni")

# Drop M0.
glm_results <- glm_results %>% filter(Module != 0)

# Fix logCPM column -- > convert to percentWT.
idy <- which(colnames(glm_results)=="logCPM")
colnames(glm_results)[idy] <- "PercentWT"
glm_results$PercentWT <- 2^glm_results$logFC

# Number of nodes per module:
n <- module_sizes[as.character(glm_results$Module)]
glm_results <- tibble::add_column(glm_results,Nodes=n,.after="Module")

# Number of significant modules.
nsig <- sum(glm_results$PAdjust < FDR_alpha)
message(paste0("\nNumber of significant ",
	      "(p-adjust < ", FDR_alpha,") ",
	      "modules: ", nsig,"."))

# Pretty print sig results:
glm_results %>% filter(PAdjust < FDR_alpha) %>%
	knitr::kable()

#--------------------------------------------------------------------
## Other module properties.
#--------------------------------------------------------------------

#  Calculate module eigengenes.
dm <- tmt_protein %>% as.data.table() %>% 
	dcast(Sample ~ Accession, value.var = "Intensity") %>% 
	as.matrix(rownames="Sample") %>% log2()
ME_data <- WGCNA::moduleEigengenes(dm, colors = partition, 
				   excludeGrey = TRUE, softPower = 1 ,
				   impute = FALSE)

# Extract PVE.
pve <- as.numeric(ME_data$varExplained)
names(pve) <- gsub("X","M",names(ME_data$varExplained))

# Add to glm_results.
glm_results <- tibble::add_column(glm_results,
				  PVE=pve[paste0("M",glm_results$Module)],
				  .after="Nodes")

#--------------------------------------------------------------------
# Add additional meta data.
#--------------------------------------------------------------------

# Annotate with module proteins.
data(gene_map)
#
symbols <- unique(tmt_protein$Symbol)
names(symbols) <- gene_map$uniprot[match(symbols,gene_map$symbol)]
#
module_list <- split(names(partition),partition)
names(module_list) <- paste0("M",names(module_list))
named_module_list <- lapply(module_list,function(x) symbols[x])
#
module_prots <- lapply(named_module_list,function(x){ 
			       paste(x,names(x),collapse="|",sep=";") })
glm_results$Proteins <-  module_prots[paste0("M",glm_results$Module)]
#
n_wash_prots <- sapply(module_list,function(x) sum(x %in% wash_prots))
glm_results$nWASH <- n_wash_prots[paste0("M",glm_results$Module)]

# Sig 85 - Significant intrafraction comparisons.
sig_prots <- tmt_protein %>% filter(FDR < 0.1) %>% 
	select(Accession) %>% unlist() %>% unique()
sig85 <- sig_prots
n_sig85 <- sapply(module_list,function(x) sum(x %in% sig_prots))
glm_results$nSig85 <- n_sig85[paste0("M",glm_results$Module)]

# Sig 62 -- sig WT v KO with > +/- 20% percent change.
sig62 <- tmt_protein %>% filter(Adjusted.FDR < 0.1) %>% 
	filter(Adjusted.logFC > log2(1.2) | Adjusted.logFC < log2(0.8)) %>%
	select(Accession) %>% unlist() %>% unique()
glm_results$nSig62 <- sapply(paste0("M",glm_results$Module),function(x){
				     sum(x %in% sig62) })

# Sig 968 --  sig WT v KO - NO log2FC threshold.
sig968 <- tmt_protein %>% filter(Adjusted.FDR < 0.1) %>% 
	select(Accession) %>% unlist() %>% unique()
glm_results$nSig968 <- sapply(paste0("M",glm_results$Module),function(x){
				      sum(x %in% sig968) })

# Combine as single list.
sig_proteins <- list()

# Save as rda object.

#--------------------------------------------------------------------
# Save results.
#--------------------------------------------------------------------

# Save as excel table
myfile <- file.path(tabsdir,"Module_GLM_results.xlsx")
results <- list("results" = glm_results)
write_excel(results,file=myfile)

# Save as rda object.
module_stats <- glm_results
myfile <- file.path(datadir,"module_stats.rda")
save(module_stats,file=myfile,version=2)

# Done!
end <- Sys.time()
message(paste("\nCompleted analysis at:",end))
