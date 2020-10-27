#!/usr/bin/env Rscript

#' title: Swip TMT Proteomics
#' description: analysis of modules for enrichment of WASH proteome
#' authors: Tyler W Bradshaw

## Optional parameters:
BF_alpha <- 0.05 # Significance threshold for enrichment

## Misc function - getrd() ----------------------------------------------------

# Get the repository's root directory.
getrd <- function(here = getwd(), dpat = ".git") {
  in_root <- function(h = here, dir = dpat) {
    check <- any(grepl(dir, list.dirs(h, recursive = FALSE)))
    return(check)
  }
  # Loop to find root.
  while (!in_root(here)) {
    here <- dirname(here)
  }
  root <- here
  return(root)
}

## Set-up the workspace -------------------------------------------------------

# Load renv
root <- getrd()
renv::load(root, quiet = TRUE)

# Global imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(geneLists) # for gene lists (pathways)
})

# Load functions in root/R and data in root/data
suppressWarnings({
  devtools::load_all()
})

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
data(module_results) # module-level statistics
data(wash_interactome) # WASH1 BioID from this study, Courtland et al., 2020. [7]


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

# Clean-up names of pre- and post- synaptic go terms in SynGO.
# idx1 <- which(grepl("GO:0098794",names(synGO))) # postsynapse
# idx2 <- which(grepl("GO:0098793",names(synGO))) # presynapse
# names(synGO)[idx1] <- paste0("SYNGO:Postsynapse (",names(synGO)[idx1],")")
# names(synGO)[idx2] <- paste0("SYNGO:Presynapse (",names(synGO)[idx2],")")

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
  ePSD # 7
)

# Remove lists with less than 3 proteins.
drop <- which(sapply(gene_lists,length) < 3)
gene_lists <- gene_lists[-drop]

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
  hyper_dt$P.adjust <- p.adjust(hyper_dt$"P-value", method = "bonferroni")
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
  hyper_dt$Proteins <- lapply(Entrez, function(x) {
    paste(gene_map$symbol[match(x, gene_map$entrez)], collapse = "; ")
  })
  # Return the results.
  results[[experiment]] <- hyper_dt
  setTxtProgressBar(pbar, value = match(experiment, names(gene_lists)))
}
close(pbar)

# Collect the results in a single data.table.
dt <- bind_rows(results)

# only sig results:
sig_dt <- dt %>% filter(P.adjust < BF_alpha)

# NOTE: No retriever enrichment
#unique(sapply(strsplit(dt$Pathway,":"),"[",1)) %notin% unique(sapply(strsplit(sig_dt$Pathway,":"),"[",1))
# NOTE: all lopitDC compartments are represented
#all(names(lopitDCpredictions) %in% sig_dt$Pathway)

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

# summarize top pathway for all modules (sig)
sig_dt %>% 
	group_by(Module) %>% 
	arrange(P.adjust) %>% 
	summarize(TopPathway = head(Pathway,1),.groups="drop") %>%
	knitr::kable()

## summarize top LopitDC predictions:
sig_dt[grepl("LopitDC",sig_dt$Pathway),] %>% 
	group_by(Pathway) %>%  arrange(P.adjust) %>%
	summarize(TopModule = head(Module,1), 
		  FE = head(`Fold enrichment`,1),
		  P.adjust = head(P.adjust,1)) %>%
	knitr::kable()

# DONE!
