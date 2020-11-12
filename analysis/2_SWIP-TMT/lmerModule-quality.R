#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: generate some gof statistics for protein and module-level models

# input partition csv file
part_file <- "ne_surprise_surprise_partition.csv" 


# prepare the env
root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

# load the data
#data(partition)
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
  ## does this statistic make sense???
  # contribution of protein to a modules variance and maximize the variance
  # attributable to the fixed effects (genotype and biofraction).
  prot_df <- msstats_prot %>% filter(Protein %in% names(partition)) %>% 
	  mutate(Module=paste0("M",partition[Protein]))
  # NOTE: variancePartition expects all factors to be modeled as mixed effects
  form <- Abundance ~ (1|Mixture) + (1|Genotype) + (1|BioFraction) + (1|Protein)
  lmer_args <- list(formula=form, 
		    data = prot_df %>% filter(Module==module), 
		    control = lme4::lmerControl(check.conv.singular = "ignore"))
  # fit the model for variancePartition
  fm <- tryCatch(expr = { do.call(lme4::lmer,lmer_args) },
		  # error
		  error = function(e) {}, 
		  # warning
		  warning = function(w) {} ) 
  # if error or warning, return NULL
  if (is.null(fm)) { 
    return(NULL) 
  }
  # otherwise calculate variance explained by protein
  vpart <- variancePartition::calcVarPart(fm)
  q <- (vpart["Genotype"] + vpart["BioFraction"]) / vpart["Protein"]
  return(q) 
} #EOF


moduleQ2 <- function(module,partition,msstats_prot) {

	# we only fit the model once and directly compute the things we are
	# intersted in... should be much faster
	prot_df <- msstats_prot %>% filter(Protein %in% names(partition)) %>% 
		mutate(Module=paste0("M",partition[Protein]))
	form1 <- Abundance ~ 0 + Genotype:BioFraction + (1|Mixture) + (1|Protein)

	lmer_args <- list(formula=form1, 
		  data = prot_df %>% filter(Module==module), 
		  control = lme4::lmerControl(check.conv.singular = "ignore"))
	fm1 <- tryCatch(expr = { do.call(lme4::lmer,lmer_args) },
			error = function(e) {}, warning = function(w) {} ) 

	if (is.null(fm1)) { 
		return(NULL) 
	}

	# compute some model statistics
	rho <- list(
		    beta = lme4::fixef(fm1),
		    X = lme4::getME(fm1, "X"),
		    vc = as.data.frame(lme4::VarCorr(fm1))
		    )

	var_protein <- rho$vc$vcov[rho$vc$grp == "Protein"] # the variance attributable to protein
	var_fixef <- stats::var(as.vector(rho$beta %*% t(rho$X))) # the variance attributable 

	q <- var_fixef/var_protein 

	return(q)
}


###############################################################################
## loop to evaluate goodness of fit for all modules defined at each resolution
###############################################################################

# load partitions
myfile <- file.path(root,"rdata",part_file)
part_dm <- fread(myfile) %>% as.matrix(rownames="V1") + 1
part_list <- unlist(apply(part_dm,1,function(x) list(x)),recursive=FALSE)
names(part_list) <- paste0("R",names(part_list))

message("\nAnalyzing partition quality at ", 
	length(part_list), " resolutions.")

# size range
module_sizes <- sapply(part_list,function(x) length(unique(x)))
t(setNames(range(module_sizes),nm=c("min k","max k"))) %>% 
	knitr::kable()


# register parallel backend
n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(n_cores)

# loop args
min_size <- 5

# parallelized nested loop: dopar + do
t1 <- Sys.time()
partition_quality <- foreach(i = seq(part_list)) %dopar% {

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
 	moduleQuality(module, partition, msstats_prot) 
  }
  names(q_list) <- names(modules)

  # collect module gof results, exlude non-numeric results
  idx <- sapply(q_list,function(x) any(!is.numeric(x)))
  if (sum(idx)>0) { warning(sum(idx)," models removed.") } # null or NA
  q <- unlist(q_list[!idx])
  Q <- sum(q[which(names(q)!="M0")])/k

  return(Q)
}
difftime(Sys.time(),t1) 

# NOTE: expected time to complete 100 resolutions ~ 20 minutes

# collect results
Q <- unlist(partition_quality)
names(Q) <- names(part_list)

# best resolution
rbest <- seq(Q)[Q==max(Q)] 
Qbest <- Q[rbest]
pbest <- part_list[[rbest]]
pbest[pbest %in% which(table(pbest) < min_size)] <- 0
mbest <- split(pbest,pbest)
kbest <- sum(names(mbest) != "M0")
cbind(rbest, kbest, Qbest) %>% knitr::kable()

# save results
myfile <- file.path(root,"rdata",
		    paste0("quality_",
			   strsplit(part_file,"\\.")[[1]][1],".rda"))
save(Q,file=myfile,version=2)
