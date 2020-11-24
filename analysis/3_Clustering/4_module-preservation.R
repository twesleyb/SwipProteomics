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

stats <- c(1, 6, 7) # Statistics by which module preservation is enforced
strength <- "strong" # Criterion for pres; weak = any(sig), strong = all(sig)
negative_edges <- "zero" # How will negative edges be replacedR? Abs val or zero
replace_zero_index <- TRUE # If min(module index)==0, +1 such that all indices >0


## ORGANIZATION of the permutation test:
# Names of discovery_data and test_data can be anything.
# These names are used to generate the output filename.
# If performing self-preservation, then discovery should be == test.
discovery_data <- "Swip"
test_data <- "Swip"

## Other NetRep DEFAULTS:
verbose <- FALSE
nPerm <- NULL
null <- "overlap"
backgroundLabel <- 0 # Modules with backgroundLabel will be ignored. *See NOTE.
alternative <- "greater" # Greater or less, for preservation use 'greater'.


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

cast2dm <- function(data, transpose = TRUE) {
	if (!transpose) { 
		stop("The data should be transposed such that ",
		     "rows = Samples, cols = Proteins.")
	}
	data %>% 
		reshape2::dcast(Protein ~ Condition + Mixture, 
				value.var = "Abundance") %>% 
	as.data.table() %>% 
	as.matrix(rownames="Protein") %>% t()
}


reset_index <- function(partition) {
  # reset partition indices
  # min index is 0
  if (min(x) == 0) {
    v <- c(0:length(unique(x)))
  } else {
    # min index is 1
    v <- c(1:length(unique(x)))
  }
  # map old indices to new ones
  namen <- unique(x)
  names(v) <- namen[order(namen)]
  y <- as.numeric(v[as.character(x)])
  names(y) <- names(x)
  return(y)
}


# Function to coerce partition matrix into named vector.
part_matrix_to_vec <- function(part_matrix, reset = TRUE) {
  colNames <- colnames(part_matrix)
  part_list <- split(part_matrix, rownames(part_matrix))
  part_list <- lapply(part_list, function(x) {
    part <- as.numeric(x)
    names(part) <- colNames
    return(part)
  })
  # Reset module numbers.
  if (reset) {
    part_list <- lapply(part_list, reset_index)
  }
  return(part_list[[1]])
}

# Function to relace negative edges in network list.
replace_negative <- function(dm, replace = c("abs", "zero")) {
  # Enforce that network (edges) are positive by replacing negative edges.
  if (replace == "abs") {
    # Replace negative edges as absolute value.
    dm_new <- abs(dm)
  } else if (replace == "zero") {
    # Replace negative edges with zero.
    dm_new <- dm
    dm_new[dm < 0] <- 0
  } else {
    stop("Please specify 'abs' or 'zero'.")
  }
  return(dm_new)
}

# Function that sets the order of the data to be the same.
set_matching_node_order <- function(adjm, netw, data, part) {
  #
  if (!is.matrix(adjm)) {
    stop("Input 'adjm' should be a matrix!")
  }
  #
  if (dim(adjm)[1] != dim(adjm)[2]) {
    stop("Input 'adjm' is not symmetric!")
  }
  #
  if (!all(rownames(adjm) %in% colnames(adjm))) {
    stop("Names of input 'adjm' are not symmetric!")
  }
  #
  if (is.null(names(part))) {
    stop("Input 'part(ition)' should be a named vector.")
  }
  #
  new_names <- names(part)
  new_adjm <- adjm[new_names, new_names]
  new_netw <- netw[new_names, new_names]
  new_data <- data[, new_names]
  #
  return(list(
    "data" = new_data, "adjm" = new_adjm,
    "part" = part, "netw" = new_netw
  ))
}


# Function to parse preservation result, set NS modules to 0.
check_preservation <- function(results, part_list, stats,
                               discovery_data, test_data) {
  # Wrapper around check_modules()
  partition <- part_list[["discovery"]]
  preservedParts <- check_modules(results, strength, stats)
  nModules <- length(unique(partition[which(partition != 0)]))
  out <- names(preservedParts)[preservedParts == "ns"]
  partition[partition %in% out] <- 0
  nPreserved <- nModules - length(out)
  message(paste(
    "...", nPreserved, "of", nModules, discovery_data,
    "modules are", c(greater = "preserved", less = "divergent")[alternative],
    "in the", test_data, "network.\n"
  ))
  # Reset module numbering.
  clean_partition <- reset_index(partition)
  # Return partition with NS modules set to 0.
  return(clean_partition)
}


check_modules <- function(x, strength = "strong", stats = c(1:7), alpha = 0.05) {
  #' check_modules
  #'
  #' Check modules for evidence of preservation.
  #'
  #' @param selfPreservation result returned by NetRep.
  #'
  #' @return none
  #'
  #' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
  #'
  #' @references none
  #'
  #' @keywords none
  #'
  #' check_modules(selfPreservation)
  # Collect observed values, nulls, and p.values -> p.adj.
  obs <- x$observed[, stats]
  nulls <- apply(x$nulls, 2, function(x) apply(x, 1, mean))[, stats]
  q <- apply(x$p.values, 2, function(x) p.adjust(x, "bonferroni"))[, stats]
  q[is.na(q)] <- 1
  # If testing more than one statistic, consider strong or weak preservation.
  fx <- c("strong" = "all", "weak" = "any")[strength]
  if (length(stats) > 1) {
    sig <- apply(q < alpha, 1, eval(fx))
    greater <- apply(obs > nulls, 1, eval(fx))
    less <- apply(obs < nulls, 1, eval(fx))
  } else {
    # If testing a single statistic...
    sig <- q < alpha
    greater <- obs > nulls
    less <- obs < nulls
  }
  # Define preserved, divergent, and ns modules.
  nModules <- length(x$nVarsPresent)
  v <- rep("ns", nModules)
  v[greater & sig] <- "preserved"
  v[less & sig] <- "divergent"
  names(v) <- names(x$nVarsPresent)
  return(v)
} # EOF


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

# load data
data(partition)
data(msstats_prot)

## Prepare input for NetRep.

# 1. Load expression data
#data(msstats_prot)
data <- msstats_prot %>% cast2dm()

# 2. Load correlation matrix
load(adjm_file)
# adjm

# 3. Load interaction/adjacency network
load(netw_file)
netw <- ne_adjm 

stopifnot(!any(netw<0)) # there should be no negative vals in netw

# 4. Load network partition
# data(partition)

# Add 1 if minimum partition index is 0.
if (replace_zero_index) {
  part <- partition + 1
}

# enforce minimum module size
stopifnot(!any(table(part) < min_size))

# only keep proteins that are in all datasets
clust_prot <- Reduce(intersect,list(colnames(adjm),colnames(netw),colnames(data),names(part)))

netw_list[["discovery"]] <- netw[clust_prot, clust_prot]
adjm_list[["discovery"]] <- adjm[clust_prot, clust_prot]
data_list[["discovery"]] <- data[, clust_prot] # cols = prots
part_list[["discovery"]] <- part[clust_prot]


## ---- status report

# Module preservation stats.
module_stats <- paste(c(
  "avg.weight", "coherence", "cor.cor", "cor.degree",
  "cor.contrib", "avg.cor", "avg.contrib"
)[stats], collapse = ", ")

# Status report:
message(paste0(
  "\nModule statistic(s) used to evaluate module ",
  c(greater = " preservation", less = " divergence")[alternative], ":\n",
  module_stats
), ".")
message(paste0(
  "\nCriterion for module",
  c(greater = " preservation", less = " divergence")[alternative], ": ",
  strength, ".", "\n"
))
message(paste(
  "\nEvaluating", c(greater = " preservation", less = " divergence")[alternative],
  "of", discovery_data,
  "modules in the", test_data, "network..."
))


## ---- Permutation testing

# Check, if performing self-preservation test, then discovery == test.
if (discovery_data == test_data) {
  self_preservation <- TRUE
  discovery <- "discovery"
  test <- "discovery"
} else {
  self_preservation <- FALSE
  discovery <- "discovery"
  test <- "test"
}


# Perform permutation test for module preservation
suppressWarnings({ # Suppress warnings about missing values
  results <- NetRep::modulePreservation(
    network = netw_list,
    data = data_list,
    correlation = adjm_list,
    moduleAssignments = part_list,
    modules = NULL,
    backgroundLabel = backgroundLabel,
    discovery = discovery,
    test = test,
    selfPreservation = self_preservation,
    nThreads = nThreads,
    nPerm = nPerm, # If null then determined automatically by the function
    null = null,
    alternative = alternative, # should be 'greater' for self-preservation
    simplify = TRUE,
    verbose = verbose
  )
})

# Collect permutation stats
dt_pval <- results$p.values %>%
  as.data.table(keep.rownames = "Module") %>%
  melt(id.vars = "Module", variable.name = "Statistic", value.name = "p.value")
dt_obs <- results$observed %>%
  as.data.table(keep.rownames = "Module") %>%
  melt(id.vars = "Module", variable.name = "Statistic", value.name = "observed")
dt_stats <- left_join(dt_pval, dt_obs, by = c("Module", "Statistic"))
fwrite(dt_stats, file.path(rdatdir, "Permutation_Statistics.csv"))

# Check how many modules are preserved.
partition <- check_preservation(
  results, part_list,
  stats, discovery_data, test_data
)

# Generate filename and save as RData or csv.
file_ext <- c(csv = ".csv", "rdata" = ".RData")
if (self_preservation) {
  output_name <- paste0(
    discovery_data,
    "_Module_Self_Preservation", file_ext[save_as]
  )
} else {
  output_name <- paste0(
    discovery_data, "_", test_data,
    "_Module_Preservation", file_ext[save_as]
  )
}
if (save_as == "rdata") {
  # Save as RData.
  saveRDS(partition, file.path(rdatdir, output_name))
} else if (save_as == "csv") {
  # Save as csv.
  # Saves partition in same format as output from LA clustering script.
  fwrite(as.data.table(t(partition)),
    file.path(rdatdir, output_name),
    row.names = TRUE
  )
} else if (save_as == "rda") {
  # Save as rda.
  myfile <- file.path(datadir, paste0(output_prefix, "partition.rda"))
  save(partition, file = myfile, version = 2)
} else {
  stop("ut oh")
}
