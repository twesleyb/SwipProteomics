#!/usr/bin/env Rscript

# title: Module Preservation
# description: evaluate module preservation by permutation testing
# authors: Tyler W A Bradshaw

## ---- INPUT data

root <- "~/projects/SwipProteomics"

adjm_file <- file.path(root, "rdata", "adjm.rda")
netw_file <- file.path(root, "rdata", "ne_adjm.rda")
data_file <- file.path(root, "data", "msstats_prot.rda")
part_file <- file.path(root, "data", "ne_surprise_surprise_partition.rda")


## ---- OUTPUT saved in root/data
# Partition of the network with preservation enforced.
# Indices of modules that are not preserved are set to 0.

# * partition.rda


## ---- OPTIONS

# Statistics by which module preservation is enforced
pres_stats <- c(1, 6, 7) 

## ---- Description of NetRep Permutation Statistics:

# 1. avg.weight (average edge weight) - Calculated from network. Assumes edge
#    weights are positive.
# 2. coherence - Calculated from the data. Quantifies the percent variance
#    explained by a modules summary vector.
# 3. cor.cor (concordance of correlation structure) - Calculated from
#    correlation matrix. Sensitive to small changes in correlation coefficients.
#    DONT USE--sensitve to small changes in edge weight- may result in false
#    positives.
# 4. cor.degree (concordance of weighted degree) - Calculated from network.
#    Assumes edge weights are positive. Sensitive to small changes in
#    weighted degree. DONT USE. Sensitive to small changes- may result in flase
#    positives.
# 5. cor.contrib (concordance of node contribution) - Calculated from the data.
#    Sensitve to small changes in node contribution. DONT USE! Same as 3,4.
# 6. avg.cor (density of correlation structure) - Calculated from correlation
#    matrix.
# 7. avg.contrib (average node contribution) - Quantifies how similar nodes are
#    to summary profile.


## ---- functions

inputDataMunge <- function(data) {
	# data = msstats_tmt
	# the function cannot tolerate missing vals. So we summarize the median
	# of the three replicates and remove the small (2) rows with any
	# remaining missing vals
	dm <- data %>% group_by(Protein, Condition) %>% 
		summarize(med_Abundance = median(Abundance),.groups="drop") %>% 
   	        reshape2::dcast(Protein ~ Condition, value.var = "med_Abundance") %>% 
		as.data.table() %>% 
		as.matrix(rownames="Protein") %>% t()
	idy <- apply(dm, 2, function(x) any(is.na(x)))
	warning("Removing ", sum(idy), " proteins with missing values",
		" before permutation testing.")
	return(dm[,!idy])
}


checkStat <- function(data) {
	data %>% mutate("p.adjust"=p.adjust(p.value,method="bonferroni")) %>% 
		# check if p.adjust < threshold and if observed > null
		mutate(isSig = `p.adjust` < 0.05 & observed > null)
}


## ---- Set-up the workspace

# Load renv
root <- "~/projects/SwipProteomics"
renv::load(root, quiet = TRUE)

devtools::load_all(quiet = TRUE)

# imports
suppressPackageStartupMessages({
  library(dplyr) 
  library(NetRep) # for permutation testing
  library(data.table) 
})

# Directories:
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")

# Number of threads for parallel processing
nThreads <- parallel::detectCores() - 1


## ---- Collect input for permutation testing

## Prepare input for NetRep

# 1. Load expression data
load(data_file)
data <- inputDataMunge(msstats_prot) 
# a matrix in which rows are samples and columns are proteins

# 2. Load correlation matrix
load(adjm_file)
# adjm

# 3. Load interaction/adjacency network
load(netw_file)
netw <- ne_adjm 

stopifnot(!any(netw<0)) # there should be no negative vals in netw

# 4. Load network partition
load(part_file)

# partition is from python which uses 0 based index. 
# Add 1 such that min = 1.
part <- partition + 1

# only keep proteins that are in all datasets
tmp_list <- list(colnames(adjm), colnames(netw), colnames(data), names(part))
clust_prot <- Reduce(intersect, tmp_list)

# NetRep expects lists as input arguments
netw_list <- list("discovery" = netw[clust_prot, clust_prot])
adjm_list <- list("discovery" = adjm[clust_prot, clust_prot])
data_list <- list("discovery" = data[, clust_prot]) # cols = prots
part_list <- list("discovery" = part[clust_prot])


## ---- Permutation testing with NetRep

# Perform permutation test for module preservation
suppressWarnings({ # Suppress warnings about missing values
  results_list <- NetRep::modulePreservation(
    network = netw_list,
    data = data_list,
    correlation = adjm_list,
    moduleAssignments = part_list,
    modules = NULL, # analyze all modules
    backgroundLabel = 0,
    discovery = "discovery",
    test = "discovery",
    selfPreservation = TRUE,
    nThreads = nThreads,
    nPerm = NULL, # determined automatically by the function
    null = "overlap",
    alternative = 'greater', # should be 'greater' for self-preservation
    simplify = TRUE,
    verbose = TRUE
  )
})


## ---- parse the response

# collect permutation pvalues
dt_pval <- results_list[["p.values"]] %>%
  as.data.table(keep.rownames = "Module") %>%
  data.table::melt(id.vars = "Module", variable.name = "Statistic", value.name = "p.value")

# collect observed statistics
dt_obs <- results_list[["observed"]] %>%
  as.data.table(keep.rownames = "Module") %>%
  data.table::melt(id.vars = "Module", variable.name = "Statistic", value.name = "observed")

# collect permutation null values
dt_nulls <- results_list[["nulls"]] %>% 
	apply(2,function(x) apply(x,1,function(y) mean(y))) %>% 
	as.data.table(keep.rownames = "Module") %>%
	data.table::melt(id.vars = "Module", variable.name = "Statistic", value.name = "null")

# combine as dt_stats
dt_tmp <- left_join(dt_pval, dt_obs, by = c("Module", "Statistic"))
dt_stats <- left_join(dt_tmp, dt_nulls, by = c("Module", "Statistic"))


# collect a list of results for each statistic
stats_list <- dt_stats %>% group_by(Statistic) %>% group_split()

# for each statistic, check if sig and obs > null
stats_list <- lapply(stats_list, checkStat)

# collect modules with preservation (all stats sig)
pres_df <- bind_rows(stats_list) %>% 
	group_by(Module) %>% 
	summarize(isPres = all(isSig), .groups="drop")

# if NA then FALSE
pres_df$isPres[is.na(pres_df$isPres)] <- FALSE

# set ns modules to 0 in final partition
final_part <- part_list[["discovery"]]
not_preserved <- pres_df$Module[!pres_df$isPres]
final_part[final_part %in% as.numeric(not_preserved)] <- 0

# status
message(sum(pres_df$isPres)," of ", nrow(pres_df), " modules are self-preserved.")


## ---- Save partition as rda

partition <- final_part
myfile <- file.path(datadir, "partition.rda")
save(partition, file = myfile, version = 2)
