#!/usr/bin/env Rscript

# title: testPreservation
# author: twab
# description: evaluate module preservation/divergence by permutation testing
# using NetRep


## NOTE: all input data is in root/rdata

root = "~/projects/SwiProteomics"

rownames_col = "Protein" 

adjmA = "wt_adjm.csv"
adjmB = "mut_adjm.csv"

netwA = "wt_netw.csv"
netwB = "mut_netw.csv"

partA = "wt_partition.csv"
partB = "mut_partition.csv"

dataA = "wt_protein.csv"
dataB = "mut_protein.csv"

# number of threads for parallel processing
nThreads = parallel::detectCores() - 1
message("utilizing: ", nThreads, " cores for parallel processing")


## ---- description of NetRep permutation statistics:

# 1. avg.weight (average edge weight) - Calculated from network. Assumes edge
#    weights are positive.
# 2. coherence - Calculated from the data. Quantifies the percent variance
#    explained by a modules summary vector.
# 3. cor.cor (concordance of correlation structure) - Calculated from
#    correlation matrix. Sensitive to small changes in correlation coefficients.
#    This statistic is sensitve to small changes in edge weight and may result 
#    in false positives.
# 4. cor.degree (concordance of weighted degree) - Calculated from network.
#    Assumes edge weights are positive. Sensitive to small changes in
#    weighted degree. This statistic is sensitive to small changes in edge
#    weights and may result in false positives.
# 5. cor.contrib (concordance of node contribution) - Calculated from the data.
#    Sensitve to small changes in node contribution. DONT USE! Same as 3,4.
# 6. avg.cor (density of correlation structure) - Calculated from correlation
#    matrix.
# 7. avg.contrib (average node contribution) - Quantifies how similar nodes are
#    to summary profile.


## ---- functions

getrd <- function(here = getwd(), dpat = ".git") {
  # get the repository's root directory
  in_root <- function(h = here, dir = dpat) {
    check <- any(grepl(dir, list.dirs(h, recursive = FALSE)))
    return(check)
  }
  while (!in_root(here)) {
    here <- dirname(here)
  }
  root <- here
  return(root)
}


rmSmall <- function(part, min_size=5) {
	# set modules with size < mize size to 0
	too_small <- as.numeric(names(which(table(part)< min_size)))
	nodes <- names(part[part %in% too_small])
	part[nodes] <- 0
	return(part)
}


compilePermResults <- function(res) {
  # a function that collects permutation statistics from NetRep
  dt_pval <- res$p.values %>% as.data.table(keep.rownames="Module") %>% 
  	melt(id.vars="Module",variable.name="Statistic",value.name="p.value")
  dt_obs <- res$observed %>% as.data.table(keep.rownames="Module") %>% 
  	melt(id.vars="Module",variable.name="Statistic",value.name="observed")
  # combine pvals and observed module statistics
  dt_stats <- left_join(dt_pval,dt_obs,by=c("Module","Statistic"))
  # compute median of null distribution for each module and each statistic
  dt_nulls <- apply(res$nulls, 2, function(x) apply(x, 1, median)) %>%
  	as.data.table(keep.rownames="Module") 
  # combine observed and null stats
  df <- left_join(dt_stats, dt_nulls, by = "Module") %>% group_by(Statistic) %>% 
  	mutate(p.adjust = p.adjust(p.value, "BH")) %>%  ungroup() %>%
  	reshape2::melt(id.vars = c("Module","Statistic", "observed",
  			           "p.value","p.adjust"),
  		       value.vars = colnames(res$observed), 
  		       value.name="null.median") %>% 
  	filter(Statistic == variable) %>%
  	dplyr::select(Module,Statistic,observed,null.median,p.value,p.adjust)
  return(df)
}


## ---- set-up the R env

suppressPackageStartupMessages({
	library(dplyr) 
	library(data.table)
	library(NetRep)
})


## ---- parse inputs

# all data are in root/rdata

root = getrd()

# full paths to input data

input_adjm = c(discovery=file.path(root,"rdata",adjmA),
	       test=file.path(root,"rdata", adjmB))

input_netw = c(discovery=file.path(root,"rdata", netwA),
	       test=file.path(root,"rdata", netwB))

input_part = c(discovery=file.path(root,"rdata", partA),
	       test=file.path(root,"rdata", partB))

input_data = c(discovery=file.path(root,"rdata", dataA),
	       test=file.path(root,"rdata", dataB))


## ---- prepare input for NetRep


# 1. Load expression data (dm)
stopifnot(all(file.exists(as.character(input_data))))
data_list <- lapply(input_data, function(x) {
			    fread(x, nThread = nThreads) %>% 
				    as.matrix(rownames=rownames_col)
})

# 2. load correlation matrices (adjm)
stopifnot(all(file.exists(as.character(input_adjm))))
adjm_list <- lapply(input_adjm, function(x) {
			    fread(x, nThread = nThreads) %>% 
				    as.matrix(rownames=rownames_col)
})

# check, adjm rownames == colnames?
invisible(lapply(adjm_list, function(x) {
			 stopifnot(all(colnames(x) == rownames(x)))
  }))


# 3. load interaction networks (netw)
netw_list <- lapply(input_netw, function(x) {
			    fread(x, nThread = nThreads) %>% 
				    as.matrix(rownames=rownames_col)
})

# FIXME: why not utilize multiple threads?

# check, netw is symmetric, rownames == colnames?
invisible(lapply(netw_list, function(x) {
			 stopifnot(all(colnames(x) == rownames(x)))
  }))

# check, network must contain positive edge weights (no negative vals)
invisible(lapply(netw_list, function(x) {
	       stopifnot(!any(x < 0))
  }))

# 4. load network partitions
# NOTE: this assumes that partition is a 2 x n protein + 1 matrix
# NOTE: the first column stores information about network resolution
# NOTE: this assumes that the min membership index is 0 from python

# load partition
part_list <- lapply(input_part, function(x) unlist(fread(x, drop=1) + 1))

# enforce minimum module size, these modules are set to 0
part_list <- lapply(part_list, rmSmall)

# Insure names are equal
A <- names(part_list[[1]])
B <- names(part_list[[2]])

stopifnot(length(A) == length(B))

stopifnot(all(A %in% B & B %in% A))

# sort the input

nodes <- A # == B

data_list <- lapply(data_list, function(x) t(x[nodes, ]))

adjm_list <- lapply(adjm_list, function(x) x[nodes, nodes])

netw_list <- lapply(netw_list, function(x) x[nodes, nodes])

part_list <- lapply(part_list, function(x) x[nodes])


## ---- perform permutation testing

# Perform permutation test for module preservation

# the chunk below will test for:
# * test modules identified in A for preservation in network B
# * test modules identified in B for preservation in network A
# * alternative = two.sided 

#  evidence of preservation if 
# * observed > NULL 
#  evidence of divergence if 
# * observed < NULL 

res_list <- NetRep::modulePreservation(
    network = netw_list,
    data = data_list,
    correlation = adjm_list,
    moduleAssignments = part_list,
    modules = NULL,
    backgroundLabel = 0,
    discovery = c("discovery","test"),
    test = c("test","discovery"),
    selfPreservation = FALSE,
    nThreads = nThreads,
    nPerm = NULL, 
    null = "overlap",
    alternative = "two.sided", 
    simplify = TRUE,
    verbose = TRUE
  )

# can analyze WT modules in KO network or MUT modules in 
results <- lapply(res_list, compilePermResults)
res1 <- results[[1]]
res2 <- results[[2]]

  
# check preservation and divergence for res1
pres1 <- res1 %>% group_by(Module) %>% 
	summarize(nSig = sum(p.adjust < 0.10 & observed > null.median)) %>%
	summarize(strong = sum(nSig ==7), weak = sum(nSig > 1))
div1 <- res1 %>% group_by(Module) %>% 
	summarize(nSig = sum(p.adjust < 0.10 & observed < null.median)) %>%
	summarize(strong = sum(nSig ==7), weak = sum(nSig > 1))

message("preservation of module from network A in network B:")
data.table('preservation'=pres1, 'divergence'=div1) %>% knitr::kable()

# check preservation and divergence for res2
pres2 <- res2 %>% group_by(Module) %>% 
	summarize(nSig = sum(p.adjust < 0.10 & observed > null.median)) %>%
	summarize(strong = sum(nSig ==7), weak = sum(nSig > 1))
div2 <- res2 %>% group_by(Module) %>% 
	summarize(nSig = sum(p.adjust < 0.05 & observed < null.median)) %>%
	summarize(strong = sum(nSig ==7), weak = sum(nSig > 1))

message("preservation of module from network B in network A:")
data.table('preservation'=pres2, 'divergence'=div2) %>% knitr::kable()
