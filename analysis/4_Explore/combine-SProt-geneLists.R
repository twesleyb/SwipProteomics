#!/usr/bin/env Rscript

#' ---
#' title: 
#' description:
#' authors: Tyler W Bradshaw
#' ---

# Trying to evaluate our model.
# Does being in the same 'organelle' predict a proteins module membership?

# Q? How well does protein A and protein B being in the same subcellular compartment predict whether or not they are in the same module. Very poorly.

#---------------------------------------------------------------------
## ARGUMENTS
#---------------------------------------------------------------------

## INPUT
# Spatial proteomics data: predicted subcellular localization of proteins from
# several studies.

## DEFAULT

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

vec2edges <- function(m) {
	dm <- 1 - diag(length(m))
	dm[lower.tri(dm,diag=TRUE)] <- NA
	colnames(dm) <- rownames(dm) <- as.character(m)
	df <- reshape2::melt(dm,na.rm=TRUE)
	colnames(df) <- c("EntrezA","EntrezB","weight")
	return(df) # edge df
}

#--------------------------------------------------------------------
## PATH
#--------------------------------------------------------------------

root <- getrd()

#---------------------------------------------------------------------
## IMPORTS
#---------------------------------------------------------------------

# Load renv.
renv::load(root, quiet=TRUE)

# Global Imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(geneLists)
	library(data.table)
})

# Local Imports.
#suppressMessages({ devtools::load_all() })

#---------------------------------------------------------------------
## LOAD DATA 
#---------------------------------------------------------------------

# Let's add our own data.
data(list=c("partition","gene_map"))

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
## Map partition UniprotIDs to Entrez.
#---------------------------------------------------------------------

# Drop unclustered proteins from partition.
data(partition)
partition <- partition[!partition == 0]
idx <- match(names(partition),gene_map$uniprot)
entrez <- gene_map$entrez[idx]
names(partition) <- entrez
head(partition)

# List of modules -- this is the same format as the gene lists.
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))
modules <- modules[!names(modules)=="M0"]

# All datasets:
gene_lists <- list(tmtp = modules,
		   i16m = itzhak2016markers,
		   i16p = itzhak2016predictions,
		   i17m = itzhak2017markers,
		   i17p = itzhak2017predictions,
		   lpm = lopitDCmarkers,
		   lpp = lopitDCpredictions)

# Test our model.
# E(ProteinA,ProteinB) ~ [SP]{0,1} | [PPI]{0,1} | [NE(bicor)]

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
knitr::kable(tmp_dt)

# Drop duplicated entrez ids.
dt <- unique(dt,by=c("study","entrez"))

# We can only compare the proteins that are present in all of the studies.
common_entrez <- Reduce(intersect,split(dt$entrez,dt$study))
message(paste("\nNumber of genes:",length(common_entrez)))

# Subset the data.
subdt <- dt %>% filter(entrez %in% common_entrez) %>% as.data.table()

# check, each group (study) should have the same number of proteins.
nprot <- unique(sapply(split(subdt$entrez,subdt$study),length))

if (!nprot == length(common_entrez)) { stop() }

# If two proteins are in the same class, then the edge between them is 1.
# Lets break up the data into groups.
data_list <- subdt %>% group_by(study) %>% group_split()
names(data_list) <- sapply(data_list,function(x) unique(x$study))

# Loop to collect adjm for each study.
adjm_list <- list()
for (study in names(data_list)) {
	## There was probably a simpler way to do this...
	tmp_dt <- data_list[[study]]
	module_list <- split(tmp_dt$entrez,tmp_dt$class)
	edge_list <- lapply(module_list,vec2edges)
	edge_dt <- bind_rows(edge_list)
	edge_dm <- edge_dt %>% select(EntrezA,EntrezB) %>% as.matrix()
	edge_dm <- apply(edge_dm,2,as.character)
	g <- igraph::graph_from_edgelist(edge_dm)
	adjm <- as.matrix(igraph::as_adjacency_matrix(g))
	# If any module contained only a single node, then this protein
	# become unconnected to any other and is droped from the adjm.
	# add these back as columns/rows == 0
	is_missing <- tmp_dt$entrez %notin% colnames(adjm)
	missing_entrez <- tmp_dt$entrez[is_missing]
	tmp_cols <- matrix(0,ncol=length(missing_entrez),nrow=nrow(adjm))
	colnames(tmp_cols) <- missing_entrez
	rownames(tmp_cols) <- rownames(adjm)
	adjm <- cbind(adjm,tmp_cols)
	tmp_rows <- matrix(0,nrow=length(missing_entrez),ncol=ncol(adjm))
	rownames(tmp_rows) <- missing_entrez
	adjm <- rbind(adjm,tmp_rows)
	adjm_list[[study]] <- adjm
}

# Insure the matrices are sorted in the same order. 
prots <- common_entrez
adjm_list <- lapply(adjm_list,function(x) x[prots,prots])
dt_list <- lapply(adjm_list,reshape2::melt)

# Evaluate the correlation between each.
contrasts <- combn(names(dt_list),2,simplify=FALSE)

myfun <- function(contrast) { 
	x <- dt_list[[contrast[1]]]
	y <- dt_list[[contrast[2]]]
	r <- cor(x$value,y$value)
	return(r)
}

result <- data.table(studyA=sapply(contrasts,"[",1),
		     studyB=sapply(contrasts,"[",2),
		     cor = sapply(contrasts,myfun))

head(result)

# Group if they contain tmtp (our data).
result$group <- as.factor(result$studyA == "tmtp" | result$studyB == "tmtp")

# Rank sum test.
wilcox.test(cor ~ group, data = result, alternative="greater")

