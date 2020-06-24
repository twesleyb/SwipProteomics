#!/usr/bin/env Rscript

#' ---
#' title: Swip Proteomics
#' description:
#' authors: Tyler W Bradshaw
#' ---

## Options:
BF_alpha = 0.1 # P.adjust significance threshold for module DA.

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. 

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

# If necessary, create dir for figs.
if (!dir.exists(figsdir)){ dir.create(figsdir, recursive = TRUE) }

# Load the data.
data(tmt_protein)

# Load the gene map.
data(gene_map)

# Load the graph partition.
data(partition)

# Load the networks.
data(adjm)
adjm <- convert_to_adjm(edges)
data(ne_adjm)
ne_adjm <- convert_to_adjm(edges)

#---------------------------------------------------------------------
## Perform module level statistical analysis with EdgeR GLM.
#---------------------------------------------------------------------
# Modules are summarized as the sum of all 
# protein Adjusted.Intensities in the module.

# Calculate the number of proteins per module.
module_sizes <- sapply(split(partition,partition), length)
names(module_sizes) <- paste0("M",names(module_sizes))

# Annotate data with module membership.
tmt_protein$Module <- paste0("M",partition[tmt_protein$Accession])

# Cast tmt_data into a dm, summarizing module Intensities grouped 
# by fraction.
dm <- tmt_protein %>% group_by(Module, Genotype, Fraction) %>% 
	dplyr::summarize(Sum.Intensity=sum(Intensity),.groups="drop") %>% 
	as.data.table() %>%
	dcast(Module ~ Fraction + Genotype,value.var="Sum.Intensity") %>%
	as.matrix(rownames=TRUE)

# Create dge object with the module-level data.
dge <- DGEList(counts=dm)

# Sample to group mapping.
groups <- rownames(dge$samples)
dge$samples$fraction <- as.factor(sapply(strsplit(groups,"_"),"[",1))
dge$samples$genotype <- as.factor(sapply(strsplit(groups,"_"),"[",2))
levels(dge$samples$genotype) <- c("WT","MUT")

# Create design matrix.
design <- model.matrix( ~ fraction + genotype, data=dge$samples)

# Flip the sign of the comparison.
# Report logFC MUT/WT.
design[,"genotypeMUT"] <- abs(design[,"genotypeMUT"]-1)

# Estimate dispersion.
dge <- estimateDisp(dge, design, robust=TRUE)

# Fit a glm model.
fit <- glmQLFit(dge,design,robust=TRUE)

# Test for DA modules:
# Default comparison is last contrast: genotype.
qlf <- glmQLFTest(fit)

# Collect results.
glm_results <- topTags(qlf,n=Inf,sort.by="p.value")$table %>%
	as.data.table(keep.rownames="Module")

# Drop M0.
glm_results <- glm_results %>% filter(Module != "M0")

# Adjust p-values for n module comparisons.
glm_results$`PAdjust (Bonferroni)` <- p.adjust(glm_results$PValue,method="bonferroni")

# Drop FDR column.
glm_results$FDR <- NULL

# Fix logCPM column -- > convert to percentWT.
# NOTE: these values are small, but reflect the average affect of GENOTYPE for
# all proteins in the module.
idy <- which(colnames(glm_results)=="logCPM")
colnames(glm_results)[idy] <- "PercentWT"
glm_results$PercentWT <- 2^glm_results$logFC

# Annotate with the number of nodes per module.
glm_results <- tibble::add_column(glm_results,
				  N=module_sizes[glm_results$Module],
				  .after="Module")

# Number of significant modules.
message(paste("\nTotal number of modules:",length(unique(partition)) -1))
nsig <- sum(glm_results$`PAdjust (Bonferroni)` < BF_alpha)
message(paste0("\nNumber of significant ",
	      "(Bonferroni p-adjust < ", BF_alpha,") ",
	      "modules: ", nsig,"."))

# Save sig modules.
sig_modules <- glm_results %>% 
	filter(`PAdjust (Bonferroni)` < BF_alpha) %>% 
	select(Module) %>% unlist() %>% unique()
myfile <- file.path(datadir,"sig_modules.rda")
save(sig_modules,file=myfile,version=2)

# Pretty print sig results:
glm_results %>% filter(`PAdjust (Bonferroni)` < BF_alpha) %>%
	knitr::kable()

#--------------------------------------------------------------------
## Add module level data to Module stats.
#--------------------------------------------------------------------

# Cast the module-level data into a matrix.
dm <- tmt_protein %>% group_by(Module,Genotype,Fraction) %>% 
	summarize(Intensity = sum(Adjusted.Intensity)) %>% 
	as.data.table() %>%
	dcast(Module ~ Fraction + Genotype, value.var = "Intensity") %>%
	as.matrix(rownames="Module")

# Sort the data columns.
f <- sapply(strsplit(colnames(dm),"_"),"[",1)
f <- as.factor(f)
levels(f) <- c("F4","F5","F6","F7","F8","F9","F10")
g <- sapply(strsplit(colnames(dm),"_"),"[",2)
g <- as.factor(g)
levels(g) <- c("WT","MUT")
col_order <- gsub("\\.","_",as.character(levels(interaction(f,g))))
dm_sorted <- dm[,col_order]

dt <- as.data.table(dm_sorted,keep.rownames="Module")
glm_results <- left_join(glm_results,dt,by="Module")

#--------------------------------------------------------------------
## These results will be saved in an excel workbook. In addition to 
## the statistical results, lets include a sheet with module summary
## statistics. We will call this object module_noa for module node 
## attributes.
#--------------------------------------------------------------------

module_list <- split(names(partition),partition)
names(module_list) <- paste0("M",names(module_list))

module_noa <- data.table("Module" = names(module_list))

module_noa$`Module Size` <- module_sizes[module_noa$Module]

#--------------------------------------------------------------------
## Calculate Module PVE
#--------------------------------------------------------------------

# PVE is the perccent of the variance explained by the modules first
# principle component.

#  Calculate module eigengenes using WGCNA::mooduleEigengenes().
dm <- tmt_protein %>% as.data.table() %>% 
	dcast(Sample ~ Accession, value.var = "Intensity") %>% 
	as.matrix(rownames="Sample") %>% log2()
ME_data <- WGCNA::moduleEigengenes(dm, colors = partition, 
				   excludeGrey = TRUE, softPower = 1 ,
				   impute = FALSE)

# Calculate bicor coorelation matrix.
ME_adjm <- WGCNA::bicor(t(ME_data$eigengenes))
ME_ne_adjm <- neten::neten(ME_adjm)

# Save.
myfile <- file.path(rdatdir,"ME_adjm.csv")
ME_adjm %>% as.data.table(keep.rownames=TRUE) %>% fwrite(myfile)
myfile <- file.path(rdatdir,"ME_ne_adjm.csv")
ME_ne_adjm %>% as.data.table(keep.rownames=TRUE) %>% fwrite(myfile)

# Extract PVE.
pve <- as.numeric(ME_data$varExplained)
names(pve) <- gsub("X","M",names(ME_data$varExplained))

# Add to module noa.
module_noa$`Module PC1 PVE` <- pve[module_noa$Module]

#--------------------------------------------------------------------
## Determine module hubs.
#--------------------------------------------------------------------

# Hubbiness statistic = node weighted degree
# Node weighted degree is the sum of a nodes edges.

# For each module, calculate node (protein) hubbiness, get the top 
# k proteins in the module and map them to GENE|ACCESSION. 

# Arbitrarily choose k = 3 HUBS.
k = 3
module_hubs <- lapply(module_list, function(prots) {
	       subadjm <- ne_adjm[prots,prots]
	       node_degree <- apply(subadjm,2,sum)
	       node_degree <- node_degree[order(node_degree,decreasing=TRUE)]
	       hubs <- names(head(node_degree,k)) # Get top three proteins.
	       symbols <- gene_map$symbol[match(hubs,gene_map$uniprot)]
	       ids <- paste(symbols,hubs,sep="|")
	       return(ids)
				  } )
idx <- module_noa$Module
module_noa$`Top 3 Hubs` <- sapply(module_hubs,paste,collapse=";")[idx]

#--------------------------------------------------------------------
## Annotate with number of sigprots at various thresholds.
#--------------------------------------------------------------------

# Sig 85: Significant intrafraction comparisons.
sig85 <- tmt_protein %>% filter(FDR < 0.1) %>% 
	select(Accession) %>% unlist() %>% unique()

# Sig 62: Sig WT v KO with > +/- 20% percent change.
sig62 <- tmt_protein %>% filter(Adjusted.FDR < 0.1) %>% 
	filter(Adjusted.logFC > log2(1.2) | Adjusted.logFC < log2(0.8)) %>%
	select(Accession) %>% unlist() %>% unique()

# Sig 968: Sig WT v KO - NO log2FC threshold.
sig968 <- tmt_protein %>% filter(Adjusted.FDR < 0.1) %>% 
	select(Accession) %>% unlist() %>% unique()

# Combine as single list.
sig_proteins <- list(sig85=sig85,sig62=sig62,sig968=sig968)

# Save as rda object.
myfile <- file.path(datadir,"sig_proteins.rda")
save(sig_proteins, file=myfile,version=2)

#--------------------------------------------------------------------
## Annotate modules with protein identifiers.
#--------------------------------------------------------------------

# List of uniprot accession.
protein_list <- module_list

# Collect gene symbols for the uniprot ids in every module. Paste these
# together as GENE|UNIPROT and collapse this seperated  by ;.
module_proteins <- sapply(protein_list, function(x) {
				  paste(paste(gene_map$symbol[match(x,gene_map$uniprot)],x,sep="|"),collapse=";")})

module_noa$Proteins <- module_proteins[module_noa$Module]

# Drop "M0"
module_noa <- module_noa %>% filter(Module != "M0")

# Save as rda.
myfile <- file.path(root,"data","module_noa.rda")
save(module_noa,file=myfile,version=2)

# Save module stats.
module_stats <- glm_results
myfile <- file.path(datadir,"module_stats.rda")
save(module_stats,file=myfile,version=2)

#--------------------------------------------------------------------
## Save results as excel document.
#--------------------------------------------------------------------

# Add to final results.
results <- list("Module Summary Stats" = module_noa, 
		"Module GLM Results" = glm_results)

# Save as excel table
myfile <- file.path(tabsdir,"Swip_TMT_Module_GLM_Results.xlsx")
write_excel(results,file=myfile)
