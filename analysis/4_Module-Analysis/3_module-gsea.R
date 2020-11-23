#!/usr/bin/env Rscript

# title: Swip TMT Proteomics
# description: analysis of modules for GSE
# author: Tyler W Bradshaw

## Optional parameters:
BF_alpha <- 0.05 # Significance threshold for enrichment
Padjust_alpha <- 0.05 


## Set-up the workspace -------------------------------------------------------

# Load renv
root <- "~/projects/SwipProteomics"
renv::load(root, quiet = TRUE)

# Global imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(geneLists) # for gene lists (pathways)
})

# Load functions in root/R and data in root/data
devtools::load_all()

# Project Directories
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Load the gene lists from geneLists.
# FIXME: geneLists need Pubmed IDs!
data(list = "corum") # CORUM protein complexes. [1]
data(list = "lopitDCpredictions") # Predicted subcellular locations from Geledaki. [2]
data(list = "takamori2006SV") # Presynaptic proteome from Takamori et al. [3]
data(list = "ePSD") # Uezu et al., 2016 [4]
data(list = "iPSD") # Uezu et al., 2016 [5]

# Add retriever complex
# Retriever complex from McNally et al., 2017. [6]
retriever <- c("Vps35l", "Vps26c", "Vps29")

# Load the data from root/data.
data(gene_map) # gene mapping data
data(partition) # the graph partition
data(msstats_prot) # the proteomics data
data(msstats_results) # protein stats
data(module_results) # module-level statistics
data(wash_interactome) # WASH1 BioID from this study, Courtland et al., 2020. [7]


## add sigprots
sig_prots <- list("SigProts" = msstats_results %>% ungroup() %>%
	filter(Contrast == "Mutant-Control", FDR < 0.05) %>% 
	select(Entrez) %>% unlist() %>% as.character() %>% unique())


## Do work --------------------------------------------------------------------

# get entrez ids for retriever prots
names(retriever) <- gene_map$entrez[match(retriever, gene_map$symbol)]

# Collect list of modules, map Uniprot accession to Entrez.
modules <- split(names(partition), partition)[-1]
module_entrez <- lapply(modules, function(x) {
  gene_map$entrez[match(x, gene_map$uniprot)]
})
names(module_entrez) <- paste0("M", names(module_entrez))
all_entrez <- unlist(module_entrez, use.names = FALSE)

# Collect WASH BioID genes.
wash_prots <- wash_interactome
wash_genes <- na.omit(gene_map$entrez[match(wash_prots, gene_map$uniprot)])

# Clean up lopit dc predictions.
names(lopitDCpredictions) <- paste("LopitDC:", names(lopitDCpredictions))

# Clean-up corum names:
names(corum) <- paste("CORUM:", names(corum))

# Clean-up takamori names.
names(takamori2006SV) <- paste("Takamori et al., 2006:", names(takamori2006SV))

# Clean-up iPSD names.
names(iPSD) <- paste("Uezu et al., 2016:", names(iPSD))

# Clean-up ePSD names.
names(ePSD) <- paste("Uezu et al., 2016:", names(ePSD))

# Collect list of entrez ids for pathways of interest.
gene_lists <- c(
  list("WASH-iBioID" = wash_genes), # 1
  corum, # 2
  lopitDCpredictions, # 3 and 4_
  list("McNally et al., 2017: Retriever Complex" = names(retriever)),
  takamori2006SV, # 5
  iPSD, # 6
  ePSD, # 7
  sig_prots
)

# Remove lists with less than 3 proteins.
part_entrez <- setNames(partition,
			nm=mapID(names(partition),"uniprot","entrez"))
idx <- which(sapply(gene_lists, function(x) sum(x %in% names(part_entrez))<3))
gene_lists <- gene_lists[-idx]

# Loop to perform GSE for each pathway.
message("\nPerforming GSE analysis for all modules:")
results <- list()
pbar <- txtProgressBar(max = length(gene_lists), style = 3)
for (experiment in names(gene_lists)) {
  # Get pathway specific genes.
  pathway_genes <- gene_lists[[experiment]]
  # Background is union of all network genes and pathway genes.
  background <- unique(c(all_entrez, pathway_genes))
  # Loop to perform hypergeometric test for enrichment.
  results_list <- list()
  for (i in c(1:length(module_entrez))) {
    results_list[[i]] <- hyperTest(
      pathway_genes,
      module_entrez[[i]],
      background
    )
  }
  names(results_list) <- paste0("M", names(modules))
  # Collect results in a data.table.
  hyper_dt <- as.data.table(do.call(rbind, results_list),
    keep.rownames = "Module"
  )
  # Adjust p-values.
  hyper_dt$FDR <- p.adjust(hyper_dt$"P-value", method = "BH")
  hyper_dt$Padjust <- p.adjust(hyper_dt$"P-value", method = "bonferroni")
  # Add module size annotation.
  sizes <- sapply(module_entrez, length)
  hyper_dt <- tibble::add_column(hyper_dt,
    "Module Size" = sizes,
    .after = "Module"
  )
  # Add pathway annotation.
  hyper_dt <- tibble::add_column(hyper_dt,
    Pathway = experiment,
    .after = "Module Size"
  )
  # Add total number of pathway genes.
  hyper_dt <- tibble::add_column(hyper_dt,
    "Total Pathway Genes" = length(pathway_genes),
    .after = "Pathway"
  )
  hyper_dt <- tibble::add_column(hyper_dt,
    "N Pathway Genes" = sum(pathway_genes %in% all_entrez),
    .after = "Total Pathway Genes"
  )
  # Number of pathway genes in each module.
  n <- sapply(module_entrez, function(x) sum(pathway_genes %in% x))
  hyper_dt <- tibble::add_column(hyper_dt,
    "n Pathway Genes in Module" = n,
    .after = "N Pathway Genes"
  )
  # Sort by fold enrichment.
  hyper_dt <- hyper_dt %>% arrange(desc(`Fold enrichment`))
  Entrez <- sapply(module_entrez[hyper_dt$Module], function(x) x[x %in% as.numeric(gene_lists[[experiment]])])
  hyper_dt$Genes <- lapply(Entrez, function(x) {
    paste(gene_map$symbol[match(x, gene_map$entrez)], collapse = "; ")
  })
  # Return the results.
  results[[experiment]] <- hyper_dt
  setTxtProgressBar(pbar, value = match(experiment, names(gene_lists)))
}
close(pbar)

# Collect the results in a single data.table.
dt <- bind_rows(results)

# only sig + enriched results:
sig_dt <- dt %>% filter(Padjust < Padjust_alpha) %>% 
	filter(`Fold enrichment` > 1) 

# Status:
n_mods <- length(unique(sig_dt$Module))
message(paste(
  "\nNumber of modules with something interesting going on:",
  n_mods
))


# Save the data.
myfile <- file.path(rdatdir, "Module_GSEA_Results.csv")
fwrite(sig_dt, myfile)

# save as rda
module_gsea <- dt
myfile <- file.path(root, "data", "module_gsea.rda")
save(module_gsea, file = myfile, version = 2)

# summarize top pathway for sig modules
sig_modules <- module_results %>% 
	filter(Padjust < 0.05) %>% 
	select(Module) %>% unlist() %>% unique()

message("\nSignificant Modules with significant gse:")
message("(",sum(sig_modules %in% sig_dt$Module)," of ",
	length(sig_modules)," significant modules.)")

# summary of sig modules with sig enrichment:
sig_dt %>% 
	group_by(Module) %>% 
	filter(Module %in% sig_modules) %>% 
	arrange(Padjust) %>% 
	summarize(TopPathway = head(Pathway,1),
		  moduleSize = head(`Module Size`,1),
		  n = head(`n Pathway Genes in Module`,1),
		  N = head(`Total Pathway Genes`,1),
		  FE = head(`Fold enrichment`,1),
		  Padjust = head(Padjust,1),
		  .groups="drop") %>% arrange(desc(FE)) %>%
	knitr::kable()

## summarize top LopitDC predictions:
message("\nModules with significant LopitDC gse:")
sig_dt[grepl("LopitDC",sig_dt$Pathway),] %>% 
	group_by(Pathway) %>%  arrange(Padjust) %>%
	summarize(TopModule = head(Module,1), 
		  FE = head(`Fold enrichment`,1),
		  Padjust = head(Padjust,1),
		  .groups = "drop") %>%
	knitr::kable()

# Save as excel
tmp_df <- data.table(Pathway=names(gene_lists),
		Entrez = sapply(gene_lists,paste,collapse=";"))
tmp_list <- list("Pathways" = tmp_df,"Module GSEA" = sig_dt)
myfile <- file.path(root,"tables","S4_Module_GSEA_Results.xlsx")
write_excel(tmp_list,myfile)

# examine NS modules with sig GSE enrichment
dt %>% filter(`Fold enrichment` > 1) %>% filter(Padjust < 0.05) %>%
	filter(Module %notin% sig_modules) %>% 
	group_by(Module) %>% 
	summarize(TopPathway = head(Pathway,1),
		  moduleSize = head(`Module Size`,1),
		  n = head(`n Pathway Genes in Module`,1),
		  N = head(`Total Pathway Genes`,1),
		  FE = head(`Fold enrichment`,1),
		  Padjust = head(Padjust,1),
		  .groups="drop") %>%
	arrange(as.numeric(gsub("M","",Module))) %>% knitr::kable()


sig_mods <- sig_dt %>% filter(Pathway == "SigProts") %>% select(Module) %>% unlist() %>% unique()
myfile <- file.path(root,"data","sig_mods.rda")
save(sig_mods,file=myfile,version=2)
