#!/usr/bin/env Rscript

#' ---
#' title: Module Preservation
#' description: Evaluate module preservation by permutation testing.
#' authors: Tyler W A Bradshaw
#' ---

## OPTIONS:
stats = c(1,6,7) # Statistics by which module preservation is enforced. 
strength = "strong" # Criterion for pres; weak = any(sig), strong = all(sig).
negative_edges = "zero" # How will negative edges be replacedR? Abs val or zero.
replace_zero_index = TRUE # If min(module index)==0, +1 such that all indices >0.
log_data = FALSE # Should input data be log2 transformed?
min_size = 5 # Minimum size of a module.

## ORGANIZATION of the permutation test:
# Names of discovery_data and test_data can be anything.
# These names are used to generate the output filename.
# If performing self-preservation, then discovery should be == test.
discovery_data = "Swip"
test_data = "Swip"

## Other NetRep DEFAULTS:
verbose = FALSE
nPerm = NULL
null = "overlap"
backgroundLabel = 0 # Modules with backgroundLabel will be ignored. *See NOTE.
alternative = "greater" # Greater or less, for preservation use 'greater'.

## INPUT data:
# All data are in root/rdata.
adjm = c(discovery="adjm.csv", test=NULL)
netw = c(discovery="ne_adjm.csv", test=NULL)
part = c(discovery="Modularity_Surprise_partition.csv", test=NULL)
data = c(discovery="norm_protein.csv", test=NULL)

## OUTPUT saved in root/rdata:
# Partition of the network with preservation enforced. 
# NOTE: Indices of modules that are not preserved are set to 0.
save_as = "rda" # Output format for partition: can be RData or csv
output_name = "MS"

#---------------------------------------------------------------------
## Description of NetRep Permutation Statistics:
#---------------------------------------------------------------------

# 1. avg.weight (average edge weight) - Calculated from network. Assumes edge
#    weights are positive.
# 2. coherence - Calculated from the data. Quantifies the percent variance
#    explained by a modules summary vector.
# 3. cor.cor (concordance of correlation structure) - Calculated from
#    correlation matrix. Sensitive to small changes in correlation coefficients.
#    DONT USE -- sensitve to small changes in edge weight--results in false positives.
# 4. cor.degree (concordance of weighted degree) - Calculated from network.
#    Assumes edge weights are positive. Sensitive to small changes in
#    weighted degree. DONT USE. Sensitive to small changes--results in false positives.
# 5. cor.contrib (concordance of node contribution) - Calculated from the data.
#    Sensitve to small changes in node contribution. DONT USE! Same as 3,4.
# 6. avg.cor (density of correlation structure) - Calculated from correlation
#    matrix.
# 7. avg.contrib (average node contribution) - Quantifies how similar nodes are
#    to summary profile.

#---------------------------------------------------------------------
## Misc function - getrd().
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

## Additional misc functions. 

# Function to add source directory to file path.
append_paths <- function(source_dir,input_file){
	# Like file.path, but using sapply to keep names.
	output_file <- sapply(input_file,function(x) {
				       file.path(source_dir,x) })
	return(output_file)
}

# Function to coerce partition matrix into named vector.
part_matrix_to_vec <- function(part_matrix,reset=TRUE){
	colNames <- colnames(part_matrix)
	part_list <- split(part_matrix,rownames(part_matrix))
	part_list <- lapply(part_list,function(x) {
				    part <- as.numeric(x)
				    names(part) <- colNames
                                    return(part) })
	# Reset module numbers.
	if (reset) { 
		part_list <- lapply(part_list, reset_index)
	} 
	return(part_list[[1]])
}

# Function to relace negative edges in network list.
replace_negative <- function(dm,replace=c("abs","zero")){
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
	if (!is.matrix(adjm)) { stop("Input 'adjm' should be a matrix!") }
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
	new_adjm <- adjm[new_names,new_names]
	new_netw <- netw[new_names,new_names]
	new_data <- data[,new_names]
	#
	return(list("data"=new_data,"adjm" = new_adjm, 
		    "part"=part, "netw" = new_netw))
}

# Function to parse preservation result, set NS modules to 0.
check_preservation <- function(results, part_list, stats, 
			       discovery_data,test_data) {
  # Wrapper around check_modules()
  partition <- part_list[["discovery"]]
  preservedParts <- check_modules(results, strength, stats)
  nModules <- length(unique(partition[which(partition != 0)]))
  out <- names(preservedParts)[preservedParts == "ns"]
  partition[partition %in% out] <- 0
  nPreserved <- nModules - length(out)
  message(paste(
    "...", nPreserved, "of", nModules, discovery_data,
    "modules are", c(greater="preserved",less="divergent")[alternative],
    "in the", test_data, "network.\n"
  ))
  # Reset module numbering.
  clean_partition <- reset_index(partition)
  # Return partition with NS modules set to 0.
  return(clean_partition)
}

#---------------------------------------------------------------------
## Set-up the workspace.
#---------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global options and imports.
suppressPackageStartupMessages({
	library(dplyr) # For manipulating the data.
	library(NetRep) # For permutation testing.
	library(data.table) # For working with tables.
})

# Additional functions.
devtools::load_all()

# Directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")

# Number of threads for parallel processing.
nThreads <- parallel::detectCores() - 1

#---------------------------------------------------------------------
## Collect input for permutation testing.
#---------------------------------------------------------------------

## Prepare input for NetRep.
# All data should be in datadir/.

# 1. Load expression data.
data_files <- append_paths(rdatdir,data)
data_list <- lapply(data_files, function(x) fread(x))
data_list <- lapply(data_list,function(x) x 
		    %>% as.matrix(rownames="Accession"))
data_list <- lapply(data_list, function(x) t(x))

# 2. Load correlation matrices.
adjm_files <- append_paths(rdatdir,adjm)
adjm_list <- lapply(adjm_files, function(x) fread(x))
adjm_list <- lapply(adjm_list, function(x) x %>% 
		    as.matrix(rownames="Accession"))
adjm_list <- lapply(adjm_list, function(x) { 
			    rownames(x) <- colnames(x)
			    return(x) })

# 3. Load interaction/adjacency networks.
netw_files <- append_paths(rdatdir,netw)
netw_list <- lapply(netw_files, function(x)  fread(x))
netw_list <- lapply(netw_list, function(x) x %>% 
		    as.matrix(rownames="Accession"))
netw_list <- lapply(netw_list,function(x) {
			    replace_negative(x,replace=negative_edges)})
netw_list <- lapply(netw_list, function(x) { 
			    rownames(x) <- colnames(x) 
			    return(x) })

# 4. Load network partitions.
part_files <- append_paths(rdatdir,part)
part_list <- lapply(part_files, function(x) fread(x,drop=1))
part_list <- lapply(part_list, part_matrix_to_vec)

# Insure that all module indices are > 0.
if (replace_zero_index) {
	# Add 1 if minimum partition index is 0.
	part_list <- lapply(part_list[which(sapply(part_list,min)==0)],function(x) x+1)
}

# Coerce to named vector.
part_list <- lapply(part_list, unlist)

# Enforce minimum module size.
too_small <- lapply(part_list, function(x) { 
			    as.numeric(names(which(sapply(split(x,x),length) < min_size))) 
			    })
# Loop to remove small modules.
message(paste("\nRemoving modules that contain less than",min_size,"nodes."))
for (i in c(1:length(part_list))) {
	part <- part_list[[i]]
	part[part %in% too_small[[1]]] <- 0
	part_list[[i]] <- part
}

# Insure that names of the and data, adjm,  netw, and partitions match.
adjm <- adjm_list[["discovery"]]
netw <- netw_list[["discovery"]]
data <- data_list[["discovery"]]
part <- part_list[["discovery"]]

# Insure names are equal.
new_data <- set_matching_node_order(adjm,netw,data,part)
data_list[["discovery"]] <- new_data$data
adjm_list[["discovery"]] <- new_data$adjm
netw_list[["discovery"]] <- new_data$netw
part_list[["discovery"]] <- new_data$part

# Should data be log transformed?
if (log_data == TRUE) {
	data_list <- lapply(data_list,log2)
	message("\nInput data was log2 transformed.")
}

#---------------------------------------------------------------------
## Permutation testing.
#---------------------------------------------------------------------

# Module preservation stats.
module_stats <- paste(c(
  "avg.weight", "coherence", "cor.cor", "cor.degree",
  "cor.contrib", "avg.cor", "avg.contrib"
)[stats], collapse = ", ")

# Status report:
message(paste0(
  "\nModule statistic(s) used to evaluate module ", 
  c(greater=" preservation",less=" divergence")[alternative],":\n",
  module_stats
), ".")
message(paste0(
  "\nCriterion for module",
  c(greater=" preservation",less=" divergence")[alternative],": ",
  strength, ".", "\n"
))
message(paste(
  "\nEvaluating",c(greater=" preservation",less=" divergence")[alternative],
  "of", discovery_data,
  "modules in the", test_data, "network..."
))

# Check, if performing self-preservation test, then discovery == test.
if (discovery_data == test_data) {
	self_preservation = TRUE
	discovery <- "discovery"
	test <- "discovery"
} else { 
	self_preservation = FALSE
	discovery <- "discovery"
	test <- "test"
}

# Perform permutation test for module preservation.
# Suppress warnings about missing values.
suppressWarnings({
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
    nPerm = nPerm, # If null then determined automatically by the function.
    null = null,
    alternative = alternative, # Greater for self-preservation.
    simplify = TRUE,
    verbose = verbose
  )
})

# Collect permutation stats.
dt_pval <- results$p.values %>% as.data.table(keep.rownames="Module") %>% 
	melt(id.vars="Module",variable.name="Statistic",value.name="p.value")
dt_obs <- results$observed %>% as.data.table(keep.rownames="Module") %>% 
	melt(id.vars="Module",variable.name="Statistic",value.name="observed")
dt_stats <- left_join(dt_pval,dt_obs,by=c("Module","Statistic"))
fwrite(dt_stats,file.path(rdatdir,"Permutation_Statistics.csv"))

# Check how many modules are preserved.
partition <- check_preservation(results, part_list, 
				stats, discovery_data, test_data)

# Generate filename and save as RData or csv.
file_ext <- c(csv = ".csv", "rdata" = ".RData")
if (self_preservation) {
	output_name <- paste0(discovery_data,
			      "_Module_Self_Preservation",file_ext[save_as])
} else {
	output_name <- paste0(discovery_data,"_",test_data,
			      "_Module_Preservation",file_ext[save_as])
}
if (save_as == "rdata") {
	# Save as RData.
	saveRDS(partition, file.path(rdatdir, output_name))
} else if (save_as == "csv") {
	# Save as csv.
	# Saves partition in same format as output from LA clustering script.
	fwrite(as.data.table(t(partition)),
	       file.path(rdatdir,output_name), row.names=TRUE)
} else if (save_as == "rda") {
	# Save as rda.
	myfile <- file.path(datadir,paste0(output_name,"partition.rda",sep="_"))
	save(partition,file=myfile,version=2)
} else {
	stop("ut oh")
}
