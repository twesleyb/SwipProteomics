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
# data(partition)
# data(tmt_protein)
# data(sig_modules)

## OPTION
# * None

## OUTPUT
# * enhanced ME adjmatrix of Significant (37) modules.

#---------------------------------------------------------------------
## FUNCTIONS
#---------------------------------------------------------------------

# Get the repository's root directory.
getrd <- function(here=getwd(), dpat= ".git") {
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
	library(igraph)
	library(data.table)
})

# Project Imports.
suppressMessages({ devtools::load_all() })

#---------------------------------------------------------------------
## LOAD DATA
#---------------------------------------------------------------------

data(partition)
data(tmt_protein)
data(sig_modules)

#---------------------------------------------------------------------
## DO WORK
#---------------------------------------------------------------------

# Collect the data and cast into a matrix.
tmt_protein$Module <- paste0("M",partition[tmt_protein$Accession])
subdt <- tmt_protein #%>% filter(Module %in% sig_modules)
subpart <- partition[names(partition) %in% subdt$Accession]
dm <- subdt %>% dcast(Sample ~ Accession, value.var="Intensity") %>%
	as.matrix(rownames="Sample") %>% log2()

# Calculate module eigengenes and create adj matrix.
ME_data <- WGCNA::moduleEigengenes(dm,colors=subpart)
ME_dm <- ME_data$eigengenes
colnames(ME_dm) <- gsub("ME","M",colnames(ME_dm))
ME_adjm <- WGCNA::bicor(ME_dm)
ME_ne_adjm <- neten::neten(ME_adjm)
adjm <- ME_ne_adjm[sig_modules,sig_modules] %>% 
	as.data.table(keep.rownames="Module") 

#---------------------------------------------------------------------
## SAVE
#---------------------------------------------------------------------

# Save to file for LA clustering.
myfile <- file.path(root,"rdata","ne_me_adjm.csv")
fwrite(adjm,myfile)
