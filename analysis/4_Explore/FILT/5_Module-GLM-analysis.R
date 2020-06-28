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

# Save sig modules.
myfile <- file.path(root,"data","sig_modules.rda")
sig_modules <- paste0("M",glm_results$Module[glm_results$PAdjust < BF_alpha])
save(sig_modules,file=myfile,version=2)

#--------------------------------------------------------------------
## Calculate Module PVE
#--------------------------------------------------------------------

# Calculate module eigengenes.
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
# Determine module hubs.
#--------------------------------------------------------------------

# All modules.
module_list <- split(partition,partition)
names(module_list) <- paste0("M",names(module_list))

# Get top three proteins.
# Sorted by node weighted degree.
module_hubs <- lapply(module_list, function(x) {
	       prots <- names(x)
	       subadjm <- ne_adjm[prots,prots]
	       node_degree <- apply(subadjm,2,sum)
	       node_degree <- node_degree[order(node_degree,decreasing=TRUE)]
	       hubs <- names(head(node_degree,3))
	       symbols <- gene_map$symbol[match(hubs,gene_map$uniprot)]
	       ids <- paste(symbols,hubs,sep="|")
	       return(ids)
				  } )

# We will add hubs to the data in the code below.
				  
#--------------------------------------------------------------------
## Add module level data.
#--------------------------------------------------------------------

df <- tmt_protein %>% group_by(Module,Genotype,Fraction) %>% 
	dplyr::summarize(Intensity = sum(Adjusted.Intensity),.groups="drop")
df$.groups <- NULL
dm <- df %>% as.data.table() %>%
	dcast(Module ~ Fraction + Genotype, value.var = "Intensity")
dm$Module <- as.character(dm$Module)

# Sort columns.
idy <- grepl("F[0-9]{1,2}",colnames(dm))
lvls <- c("F4","F5","F6","F7","F8","F9","F10")
f <- factor(sapply(strsplit(colnames(dm)[idy],"_"),"[",1),levels=lvls)
g <- factor(sapply(strsplit(colnames(dm)[idy],"_"),"[",2),levels=c("WT","MUT"))
idy <- c("Module",gsub("\\.","_",as.character(levels(interaction(f,g)))))
dm <- dm %>% select(all_of(idy))

# Add to the data.
glm_results <- left_join(glm_results,dm,by="Module")

# Annotate with hubs.
hubs_list <- module_hubs[paste0("M",glm_results$Module)]
glm_results$Hubs <- sapply(hubs_list,paste,collapse="; ")

# Annotate with module proteins.
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

#--------------------------------------------------------------------
## Calculate mean and SEM of groups.
#--------------------------------------------------------------------

df <- glm_results
cols <- c(grep("F[0-9]{1,2}_WT",colnames(df)),grep("MUT",colnames(df)))
dm <- df %>% dplyr::select(Module,all_of(cols)) %>% 
	as.data.table() %>% as.matrix(rownames="Module")
idy <- grepl("WT",colnames(dm))
WT_means <- apply(dm,1,function(x) log2(mean(x[idy])))
WT_SEM <- apply(dm,1,function(x) sd(log2(x[idy])))/WT_means
idy <- grepl("MUT",colnames(dm))
MUT_means <- apply(dm,1,function(x) log2(mean(x[idy])))
MUT_SEM <- apply(dm,1,function(x) sd(log2(x[idy])))/MUT_means
df <- tibble::add_column(df,"WT Mean" = WT_means, .after="PVE")
df <- tibble::add_column(df,"WT SEM" = WT_SEM, .after="WT Mean")
df <- tibble::add_column(df,"MUT Mean" = MUT_means, .after="WT SEM")
df <- tibble::add_column(df,"MUT SEM" = MUT_SEM, .after="MUT Mean")
glm_results <- df

#--------------------------------------------------------------------
# Save results.
#--------------------------------------------------------------------

# FIXME: Save sig proteins.
#sig_proteins <- ""
#myfile <- file.path(root,"data","sig_proteins.rda")
#save(sig_proteins,file=myfile,version=2)

# Data.table describing graph partition.
Uniprot <- names(partition)
idx <- match(Uniprot,gene_map$uniprot)
Entrez <- gene_map$entrez[idx]
Symbol <- gene_map$symbol[idx]
part_dt <- data.table(Uniprot,Entrez,Symbol,Module=partition)

# Save as excel table.
myfile <- file.path(tabsdir,"Swip_TMT_Module_GLM_Results.xlsx")
results <- list("Network Partition" = part_dt,
		"Module GLM Results" = glm_results)
write_excel(results,file=myfile)

# Save a copy in root/supplment.
myfile <- file.path(suppdir,"S3_Swip_TMT_Module_GLM_Results.xlsx")
write_excel(results,file=myfile)

# Save as rda object.
module_stats <- glm_results
myfile <- file.path(datadir,"module_stats.rda")
save(module_stats,file=myfile,version=2)
