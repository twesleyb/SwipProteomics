#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: generate some gof statistics for protein and module-level models

# prepare the env
root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

# load the data
data(partition)
data(msstats_prot)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(doParallel)
  library(variancePartition)
})


## functions ------------------------------------------------------------------


moduleQuality <- function(module,partition,msstats_prot){
  # a function the evaluates gof of modules with variancePartition
  # prepare input arguments to lme4::lmer
  prot_df <- msstats_prot %>% filter(Protein %in% names(partition)) %>% 
	  mutate(Module=paste0("M",partition[Protein]))
  # NOTE: variancePartition expects all factors to be modeled as mixed effects
  form1 <- Abundance ~ (1|Mixture) + (1|Genotype) + (1|BioFraction) + (1|Protein)
  lmer_args <- list(formula=form1, 
		    data = prot_df %>% filter(Module==module), 
		    control = lme4::lmerControl(check.conv.singular = "ignore"))
  # fit the model for variancePartition
  fm1 <- tryCatch(expr = { do.call(lme4::lmer,lmer_args) },
		  # error
		  error = function(e) {}, 
		  # warning
		  warning = function(w) {} ) 
  # if error or warning, return NULL
  if (is.null(fm1)) { 
    return(NULL) 
  }
  # otherwise calculate variance explained by protein
  r2_protein <- as.numeric(variancePartition::calcVarPart(fm1)["Protein"])
  # prepare input arguments to lme4::lmer
  prot_df <- msstats_prot %>% filter(Protein %in% names(partition)) %>% 
	  mutate(Module=paste0("M",partition[Protein]))
  # NOTE: variancePartition expects all factors to be modeled as mixed effects
  form2 <- Abundance ~ 0 + Genotype:BioFraction + (1|Mixture) + (1|Protein)
  lmer_args <- list(formula=form2, 
		    data = prot_df %>% filter(Module==module), 
		    control = lme4::lmerControl(check.conv.singular = "ignore"))
  fm2 <- tryCatch(expr = { do.call(lme4::lmer,lmer_args) },
		  # error
		  error = function(e) {}, 
		  # warning
		  warning = function(w) {} ) 
  r2_fixef <- as.numeric(r.squaredGLMM.merMod(fm2)[,"R2m"])
  q <- r2_fixef/r2_protein
  return(q) 
} #EOF


###############################################################################
## explore module quality at a single resolution
###############################################################################

# register parallel backend
n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(n_cores)

# input partition csv file
part_file <- "cpm_partition.csv"

# load partitions
myfile <- file.path(root,"rdata",part_file)
part_dm <- fread(myfile) %>% as.matrix(rownames="V1") + 1
part_list <- unlist(apply(part_dm,1,function(x) list(x)),recursive=FALSE)


# for a given partition...
i = sample(length(part_list),1)
message("\nExamining resolution: ",i)
partition <- part_list[[i]]

# ... remove modules < 5 nodes
partition[partition %in% which(table(partition) < 5)] <- 0

# for every module in that partition...
modules <- split(names(partition),partition)
k <- sum(names(modules) != "M0")
modules <- setNames(modules,nm=paste0("M",names(modules)))
message("\nNumber of clusters: ", k)

# ... evaluate its quality
q_list <- foreach(module = names(modules)) %dopar% { 
	moduleQuality(module, partition,msstats_prot) }
names(q_list) <- names(modules)

# collect module gof results, exlude NA
# FIXME: what about NULL
idx <- sapply(q_list,function(x) any(is.na(x))) # drop NA
q <- unlist(q_list[!idx])
Q <- sum(q[which(names(q)!="M0")])/k
message("\nQuality: ", round(Q,5))


###############################################################################
## loop to evaluate goodness of fit for all modules defined at each resolution
###############################################################################

# register parallel backend
n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(n_cores)

# loop args
min_size <- 5

# parallelized nested loop: dopar + do
t1 <- Sys.time()
results <- foreach(i = seq(part_list)) %dopar% {
  # for a given partition...
  partition <- part_list[[i]]
  # ... remove modules < min_size nodes
  partition[partition %in% which(table(partition) < min_size)] <- 0
  # for every module in that partition...
  modules <- split(names(partition),partition)
  k <- sum(names(modules) != "M0")
  modules <- setNames(modules,nm=paste0("M",names(modules)))
  # ... evaluate its quality
  q_list <- foreach(module = names(modules)) %do% { 
 	moduleQuality(module, partition,msstats_prot) }
  names(q_list) <- names(modules)
  # collect module gof results, exlude non-numeric results
  idx <- sapply(q_list,function(x) any(!is.numeric(x)))
  if (sum(idx)>0) { warning(sum(idx),"models removed.") } # null or NA
  q <- unlist(q_list[!idx])
  Q <- sum(q[which(names(q)!="M0")])/k
  return(Q)
}
difftime(Sys.time(),t1) 
# expected time to complete 100 resolutions ~ 20 minutes

## does this statistic make sense???
# contribution of protein to a modules variance and maximize the variance
# attributable to the fixed effects (genotype and biofraction).
unlist(results) # i think it might... basically we want to minimize the 
