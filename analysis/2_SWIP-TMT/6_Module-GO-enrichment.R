#!/usr/bin/env Rscript

#' ---
#' title: Module Analysis
#' description: Analysis of co-expression modules
#' authors: Tyler W Bradshaw
#' ---

## User parameters to change:
FDR_alpha = 0.05 # significance threshold for GO enrichment.

# NOTE: GO enrichment is performed with the Langfelder Lab's 
# anRichment::enrichmentAnalysis function which utilizes 
# Fisher's exact test in order to evaluate the statistical enrichment
# of GO BP MF and CC terms in a module. The background for the test 
# is all 'given' genes. All 5,84x protein in the dataset.

#---------------------------------------------------------------------
## Set-up the workspace.
#---------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global options and imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
})

# Load additional functions in root/R.
suppressWarnings({ devtools::load_all() })

# Project directories:
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

#---------------------------------------------------------------------
## Load the data.
#---------------------------------------------------------------------

# Load the partition and tmt data.
data(partition)
data(tmt_protein)

#---------------------------------------------------------------------
## Module GO enrichment.
#---------------------------------------------------------------------

# List of modules.
module_list <- split(partition,partition)
names(module_list) <- paste0("M",names(module_list))

# Drop M0.
module_list <- module_list[-which(names(module_list)=="M0")]

# collect list of entrez ids for each module.
entrez_list <- lapply(module_list,function(x) {
			    idx <- match(names(x),tmt_protein$Accession)
			    entrez <- tmt_protein$Entrez[idx]
			    return(entrez)
			    })

# Build a GO reference collection:
message("\nBuilding a mouse GO reference collection with anRichment.")
gocollection <- suppressPackageStartupMessages({
	anRichment::buildGOcollection(organism="mouse")
})

# perform module go enrichment.
# NOTE: background defaults to all genes in gene_list.
go_gse <- gse(entrez_list, gocollection)

# collect significant modules.
idx <- sapply(go_gse,function(x) any(x$FDR < FDR_alpha))
nsig_go <- sum(idx)
message(paste("\nTotal number of modules:",length(module_list)))
message(paste("\nNumber of modules with significant",
	      "GO term enrichment:",nsig_go))

# top sig go term for each module.
fe <- sapply(go_gse,function(x) round(x$enrichmentRatio[1],2))
fdr <- sapply(go_gse,function(x) round(x$FDR[1],2))
m <- sapply(go_gse, function(x) x$shortDataSetName[1])
top_go <- data.table(module=names(go_gse),
		       term = m,
		       fe=fe,fdr=fdr)

# summary:
message("\nModules with significant go enrichment:")
knitr::kable(top_go[idx,])

# Combine results, keep sig results.
go_results <- bind_rows(go_gse,.id="Module") %>% filter(FDR < FDR_alpha)

# Clean-up.
rename_column <- function(col,new_name) {
	# A function that renames a column in a data.frame.
	idx <- which(colnames(go_results)==col)
	colnames(go_results)[idx] <- new_name
	return(go_results)
}

# Remove some un-needed columns.
go_results$class <- NULL
go_results$fracOfEffectiveClassSize <- NULL
go_results$expectedFracOfEffectiveClassSize <- NULL
go_results$effectiveClassSize <- NULL
go_results$fracOfEffectiveSetSize <- NULL
go_results$effectiveSetSize <- NULL
go_results$rank <- NULL

# Rename some columns.
go_results <- rename_column("dataSetID","ID")
go_results <- rename_column("dataSetName","Name")
go_results <- rename_column("inGroups","Ontology")
go_results <- rename_column("pValue","PValue")
go_results <- rename_column("nCommonGenes","n Genes")
go_results <- rename_column("enrichmentRatio","Fold Enrichment")
go_results <- rename_column("classSize","Module Size")
go_results <- rename_column("shortDataSetName","Short Name")
go_results <- rename_column("overlapGenes","Entrez")
go_results$Entrez <- gsub("\\|",";",go_results$Entrez)

# Reorganize.
go_results <- go_results %>% dplyr::select(Module,`Short Name`,Name, ID, 
					   Ontology, `Module Size`, `n Genes`, 
					`Fold Enrichment`, PValue, FDR, Entrez)

#--------------------------------------------------------------------
## Save GO results for significant modules in a seperate sheet.
#--------------------------------------------------------------------

data(sig_modules)
sub_results <- go_results %>% filter(Module %in% sig_modules)

#--------------------------------------------------------------------
## Save the data.
#--------------------------------------------------------------------

# save significant results.
write_excel(list("Module GO Results" = go_results, 
		 "Sig Module GO Results"=sub_results),
	    file.path(tabsdir,"Swip_TMT_Module_GO_Results.xlsx"))

# Save as rda.
module_go <- go_results
myfile <- file.path(root,"data","module_GO.rda")
save(module_go,file=myfile,version=2)
