#!/usr/bin/env Rscript

#' ---
#' title: 
#' description:
#' authors: Tyler W Bradshaw
#' ---

#---------------------------------------------------------------------
## ARGUMENTS
#---------------------------------------------------------------------

## INPUT

## OPTION

## OUTPUT

#---------------------------------------------------------------------
## FUNCTIONS
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
## IMPORTS
#---------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root, quiet=TRUE)

# Global Imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
})

# Local Imports.
suppressMessages({ devtools::load_all() })

#---------------------------------------------------------------------
## LOAD DATA 
#---------------------------------------------------------------------

# Let's add our own data.
data(list=c("partition","gene_map"))
data(sig_modules)

# Spatial proteomics datasets:
data(itzhak2016)
data(itzhak2016markers)
data(itzhak2016predictions)

data(itzhak2017)
data(itzhak2017markers)
data(itzhak2017predictions)

data(lopitDCmarkers)
data(lopitDCpredictions)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

# Map names(parition) from Uniprot to Entrez.
# Drop unclustered proteins from partition.
data(partition)
partition <- partition[!partition == 0]
idx <- match(names(partition),gene_map$uniprot)
entrez <- gene_map$entrez[idx]
names(partition) <- entrez

# List of modules -- this is the same format as the gene lists.
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))
modules <- modules[!names(modules)=="M0"]

# Combine all spatial proteomics datasets:
gene_lists <- list(tmtp = modules,
		   i16m = itzhak2016markers,
		   i16p = itzhak2016predictions,
		   i17m = itzhak2017markers,
		   i17p = itzhak2017predictions,
		   lpm = lopitDCmarkers,
		   lpp = lopitDCpredictions)

# All entrez gene identifiers.
# Complicated bc we need to unnest the list but avoid changing the names.
all_entrez = unlist(gene_lists,recursive = FALSE)
x <- unlist(all_entrez,use.names=FALSE)
names(x) <- rep(names(all_entrez),times=sapply(all_entrez,length))
all_entrez <- x

# Collect the data from all studies. 
# The variable class is the 'known'/marker or predicted subcellular localization.
dt <- data.table(study = sapply(strsplit(names(all_entrez),"\\."),"[",1),
		 class = sapply(strsplit(names(all_entrez),"\\."), "[",2),
		 entrez = all_entrez)

# Number of duplicates per study (they can be duplicated multiple times).
tmp_dt <- dt %>% group_by(study) %>% summarize(n=sum(duplicated(entrez)))

# Drop duplicated entrez ids.
dt <- unique(dt,by=c("study","entrez"))

# We can only compare the proteins that are present in all of the studies.
common_entrez <- Reduce(intersect,split(dt$entrez,dt$study))
message(paste("\nNumber of genes in all studies:",length(common_entrez)))

# Subset the data.
subdt <- dt %>% filter(entrez %in% common_entrez) %>% as.data.table()

# check, each group (study) should have the same number of proteins.
nprot <- unique(sapply(split(subdt$entrez,subdt$study),length)) 
if (!nprot == length(common_entrez)) { stop() } 

# Work it.
tmtp <- subdt %>% filter(study == "tmtp")
subdt <- subdt %>% filter(study != "tmtp")
subdt$class <- toupper(subdt$class)
subdt <- subdt %>% filter(class != "UNKNOWN")
# Combine ER.
subdt$class[grepl("ER_HIGH_CURVATURE",subdt$class)] <- "ER"
# Combine nucleus.
subdt$class[grepl("NUCLEUS/CHROMATIN",subdt$class)] <- "NUCLEUS"
subdt$class[grepl("NUCLEAR PORE COMPLEX",subdt$class)] <- "NUCLEUS"
subdt$class[grepl("NUCLEAR PORE COMPLEX/NUCLEAR",subdt$class)] <- "NUCLEUS"
# Combine golgi.
subdt$class[grepl("^GA$",subdt$class)] <- "GOLGI"

# Keep proteins that are consistent across all studies -- prediction is
# conserved.
df <- subdt %>% group_by(entrez) %>% summarize(n_studies = length(unique(study)),
					     prediction = unique(class)[1],
					     isConserved=length(unique(class))==1,
					     .groups="drop")
subdf <- df %>% filter(isConserved)

#  collect thes predictions
predictions <- subdf$prediction
names(predictions) <- subdf$entrez

# subset partition, keep genes with predictions
idx <- names(partition) %in% names(predictions)
part <- partition[idx]

# The result.
df <- data.table("Entrez" = names(part),
		 "Module" = part,
		 "Prediction" = predictions[names(part)])
df <- df %>% arrange(Module)

# summarize genes
subdf <- df %>% group_by(Module) %>% 
	summarize('nGenes'=length(Entrez),
		  'nPred'=length(unique(Prediction)),
		  'cPred'=paste(unique(Prediction),collapse="; "),
		  .groups = "drop")

# annotate with sig modules.
subdf$isSig <- subdf$Module %in% as.numeric(gsub("M","",sig_modules))
subdf %>% filter(grepl("LYSOSOME",subdf$cPred)) %>% knitr::kable()
