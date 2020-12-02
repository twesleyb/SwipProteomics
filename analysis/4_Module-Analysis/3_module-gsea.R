#!/usr/bin/env Rscript

# title: Swip TMT Proteomics
# description: analysis of modules for GSE
# author: Tyler W Bradshaw

## ---- Input parameters

BF_alpha <- 0.05 
save_results <- TRUE
input_part <- "ne_surprise_partition"

## ---- Set-up the workspace 

# Load renv
root <- "~/projects/SwipProteomics"
renv::load(root, quiet = TRUE)

# Global imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(geneLists) # for gene lists (pathways) and hyperTest function
})

# Load functions in root/R and data in root/data
devtools::load_all(root, quiet = TRUE)

# Load the data from root/data
data(list=input_part)

data(gene_map)
data(sig_prots)
data(sig_modules)
data(msstats_prot)
data(msstats_results)
data(wash_interactome) 

# Project Directories
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Load the gene lists from twesleyb/geneLists
# FIXME: geneLists need Pubmed IDs!
data(list = "corum") # CORUM protein complexes [1]
data(list = "lopitDCpredictions") # protein predicted subcellular loc [2]
data(list = "takamori2006SV") # Presynaptic proteome from Takamori et al. [3]
data(list = "ePSD") # Uezu et al., 2016 [4]
data(list = "iPSD") # Uezu et al., 2016 [5]

# Add retriever complex
# Retriever complex from McNally et al., 2017. [6]
retriever <- c("Vps35l", "Vps26c", "Vps29")


## ---- Do work 

# get entrez ids for retriever prots
names(retriever) <- gene_map$entrez[match(retriever, gene_map$symbol)]

# Collect list of modules, map Uniprot accession to Entrez.
modules <- split(names(partition), partition)[-1]
module_entrez <- lapply(modules, function(x) {
  gene_map$entrez[match(x, gene_map$uniprot)]
})
names(module_entrez) <- paste0("M", names(module_entrez))
all_entrez <- unlist(module_entrez, use.names = FALSE)

# Collect WASH BioID genes
wash_prots <- wash_interactome
wash_genes <- na.omit(gene_map$entrez[match(wash_prots, gene_map$uniprot)])

# Clean-up lopit dc predictions
names(lopitDCpredictions) <- paste("LopitDC:", names(lopitDCpredictions))

# Clean-up corum names
names(corum) <- paste("CORUM:", names(corum))

# Clean-up Takamori et al. pathway names
names(takamori2006SV) <- paste("Takamori et al., 2006:", names(takamori2006SV))

# Clean-up iPSD names
names(iPSD) <- paste("Uezu et al., 2016:", names(iPSD))

# Clean-up ePSD names
names(ePSD) <- paste("Uezu et al., 2016:", names(ePSD))

# Collect list of entrez ids for pathways of interest
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

# Remove lists with less than 3 proteins
part_entrez <- setNames(partition,
			nm=mapID(names(partition),"uniprot","entrez"))
idx <- which(sapply(gene_lists, function(x) sum(x %in% names(part_entrez))<3))
gene_lists <- gene_lists[-idx]

# Loop to perform GSE for each pathway
message("\nPerforming GSE analysis for all modules:")
results <- list()
pbar <- txtProgressBar(max = length(gene_lists), style = 3)
for (experiment in names(gene_lists)) {
  # Get pathway specific genes
  pathway_genes <- gene_lists[[experiment]]
  # Background is union of all network genes and pathway genes
  background <- unique(c(all_entrez, pathway_genes))
  # Loop to perform hypergeometric test for enrichment
  results_list <- list()
  for (i in c(1:length(module_entrez))) {
    results_list[[i]] <- hyperTest(
      pathway_genes,
      module_entrez[[i]],
      background
    )
  }
  names(results_list) <- paste0("M", names(modules))
  # Collect results in a data.table
  hyper_dt <- as.data.table(do.call(rbind, results_list),
    keep.rownames = "Module"
  )
  # Adjust p-values
  hyper_dt$FDR <- p.adjust(hyper_dt$"P-value", method = "BH")
  hyper_dt$Padjust <- p.adjust(hyper_dt$"P-value", method = "bonferroni")
  # Add module size annotation
  sizes <- sapply(module_entrez, length)
  hyper_dt <- tibble::add_column(hyper_dt,
    "Module Size" = sizes,
    .after = "Module"
  )
  # Add pathway annotation
  hyper_dt <- tibble::add_column(hyper_dt,
    Pathway = experiment,
    .after = "Module Size"
  )
  # Add total number of pathway genes
  hyper_dt <- tibble::add_column(hyper_dt,
    "Total Pathway Genes" = length(pathway_genes),
    .after = "Pathway"
  )
  hyper_dt <- tibble::add_column(hyper_dt,
    "N Pathway Genes" = sum(pathway_genes %in% all_entrez),
    .after = "Total Pathway Genes"
  )
  # Number of pathway genes in each module
  n <- sapply(module_entrez, function(x) sum(pathway_genes %in% x))
  hyper_dt <- tibble::add_column(hyper_dt,
    "n Pathway Genes in Module" = n,
    .after = "N Pathway Genes"
  )
  # Sort by fold enrichment
  hyper_dt <- hyper_dt %>% arrange(desc(`Fold enrichment`))
  Entrez <- sapply(module_entrez[hyper_dt$Module], function(x) x[x %in% as.numeric(gene_lists[[experiment]])])
  hyper_dt$Genes <- lapply(Entrez, function(x) {
    paste(gene_map$symbol[match(x, gene_map$entrez)], collapse = "; ")
  })
  # Return the results
  results[[experiment]] <- hyper_dt
  setTxtProgressBar(pbar, value = match(experiment, names(gene_lists)))
}
close(pbar)


## ---- Collect the results in a single data.table

dt <- bind_rows(results)

# only sig + enriched results:
sig_dt <- dt %>% filter(Padjust < BF_alpha) %>% 
	filter(`Fold enrichment` > 1) 


## ---- status

m <- length(unique(sig_dt$Module))
M <- length(modules)
message(m, " of ", M, " modules exhibit some significant GSE.")

# modules that are enriched for sig_prots
sig_gsea <- sig_dt %>% filter(Pathway == "SigProts") %>% 
	select(Module) %>% unlist() %>% unique()

# modules with sig lopitDC enrichment
message("Modules enriched for LopitDC subcellular compartments:")
idx <- grepl("LopitDC", sig_dt$Pathway)
sig_dt %>% filter(idx) %>% 
	select(Module, Pathway, Padjust, `Fold enrichment`) %>% 
	knitr::kable()


## ---- save results

if (save_results) {

  # save sig_gsea as rda
  namen <- gsub("partition","sig_gsea.rda",input_part)
  myfile <- file.path(root,"data", namen)
  save(sig_gsea,file=myfile,version=2)

  # save as rda
  module_gsea <- sig_dt
  namen <- gsub("partition","module_gsea.rda",input_part)
  myfile <- file.path(root, "data", namen)
  save(module_gsea, file = myfile, version = 2)

  # save as excel
  idx <- order(as.numeric(gsub("M","",sig_dt$Module)))
  sig_dt <- sig_dt[idx,] # sort by module
  tmp_df <- data.table(Pathway=names(gene_lists),
 		  Entrez = sapply(gene_lists,paste,collapse=";"))
  tmp_list <- list("Module GSEA" = sig_dt,"Pathways" = tmp_df)
  namen <- gsub("partition","S6_SWIP-TMT_Module_GSEA.xlsx",input_part)
  myfile <- file.path(root,"tables", namen)
  write_excel(tmp_list,myfile)

} # EIS
