#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: fit protein-wise lmer models and perform statistical inferences
# for a given contrast

# fx0: Abundance ~ 0 + Condition + (1|Mixture)

## prepare the env ------------------------------------------------------------

## input 
results_file = "fit0_results.csv" # saved in root/rdata
FDR_alpha = 0.05 # threshold for significance

## load renv
root <- "~/projects/SwipProteomics"
renv::load(root); devtools::load_all(root)

## load data
data(swip)
data(gene_map)
data(msstats_prot)

## other imports
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
  # require(lme4)
  # require(knitr)
  # require(tibble)
  # require(lmerTest)
  # requre(data.table)
})


## functions ------------------------------------------------------------------

getContrast <- function(fm, negative_index, positive_index) {
  # build a contrast matrix:
  contrast_matrix <- lme4::fixef(fm)
  contrast_matrix[] <- 0
  contrast_matrix[positive_index] <- +1 
  contrast_matrix[negative_index] <- -1
  return(contrast_matrix)
} #EOF


expandGroups <- function(conditions,biofractions) {
	# munge to create contrast matrix for fm1
	groups <- apply(expand.grid(conditions,biofractions),1,paste,collapse=".")
	idx <- rep(c(1:length(biofraction)),each=length(condition))
	contrast_list <- split(groups, idx)
	return(contrast_list)
} #EOF


## check Swip's fit -----------------------------------------------------------

## formulae to be fit:
fx0 <- formula("Abundance ~ 0 + Condition + (1|Mixture)")
#fx1 <- formula("Abundance ~ 0 + Genotype + BioFraction + (1|Subject)")

## fit model 0
fm0 <- lmerTest::lmer(fx0, msstats_prot %>% filter(Protein == swip))
print(summary(fm0, ddf = "Satterthwaite"))

# cm1
#cm1 <- getContrast(fm1,"GenotypeControl","GenotypeMutant")

## munge to create contrast matrices for intrafraction comparisons:
condition <- c("ConditionControl","ConditionMutant")
biofraction <- c("F4","F5","F6","F7","F8","F9","F10")
contrasts <- lapply(expandGroups(condition,biofraction), function(x) {
			  getContrast(fm0,x[1],x[2]) })
names(contrasts) <- sapply(contrasts, function(x) {
				 paste(names(x)[x==+1],names(x)[x==-1],sep="-") })

# check the results
results <- lmerTestProtein(protein,fx0,msstats_prot,contrasts)
results$stats %>% knitr::kable()


## loop to fit all proteins ----------------------------------------------------

n_cores <- parallel::detectCores() - 1
BiocParallel::register(BiocParallel::SnowParam(n_cores))

prots = unique(as.character(msstats_prot$Protein))

results_list <- foreach(protein = prots) %dopar% {
	suppressMessages({
	  try(lmerTestProtein(protein, fx0, msstats_prot, contrasts),silent=TRUE)
	})
} # EOL


## process results ------------------------------------------------------------

idx <- unlist(sapply(results_list,class)) != "try-error"
filt_list <- results_list[which(idx)]
results_df <- bind_rows(sapply(filt_list,"[[","stats"))

# drop singular
results_df <- results_df %>% filter(!isSingular)
results_df$isSingular <- NULL

## annotate with gene symbols
idx <- match(results_df$protein,gene_map$uniprot)
results_df <- tibble::add_column(results_df,
  				 symbol=gene_map$symbol[idx],
  				 .after="protein")

## adjust pvals 
results_df <- tibble::add_column(results_df, 
			 Padjust=p.adjust(results_df$Pvalue,"BH"),
			 .after="Pvalue")

## sort
results_df <- results_df %>% arrange(Pvalue)

# examine top results
results_df %>% head() %>% knitr::kable()

# status
message("Total number of significant proteins: ",
	sum(results_df$Padjust < FDR_alpha))


## save results ----------------------------------------------------------------

# save as csv
myfile <- file.path(root,"rdata",results_file)
fwrite(results_df, myfile)
#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: examine the fits for WASHC4 (aka Swip) and save contrast matrices
# to be used by lmerTestProtein

## prepare the env
root <- "~/projects/SwipProteomics"
renv::load(root); devtools::load_all(root)

## data
data(swip)
data(msstats_prot)

## other imports
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
  # require(lme4)
  # require(knitr)
  # require(tibble)
  # require(lmerTest)
  # requre(data.table)
})


## functions ------------------------------------------------------------------

getContrast <- function(fm, negative_index, positive_index) {
  # build a contrast matrix:
  contrast_matrix <- lme4::fixef(fm)
  contrast_matrix[] <- 0
  contrast_matrix[positive_index] <- +1 
  contrast_matrix[negative_index] <- -1
  return(contrast_matrix)
}


expandGroups <- function(conditions,biofractions) {
	# munge to create contrast matrix for fm1
	groups <- apply(expand.grid(conditions,biofractions),1,paste,collapse=".")
	idx <- rep(c(1:length(biofraction)),each=length(condition))
	contrast_list <- split(groups, idx)
	return(contrast_list)
}

## save results ----------------------------------------------------------------

# save contrast matrices
myfile <- file.path(root,"data","cm_list.rda")
save(cm_list,file=myfile,version=2)
