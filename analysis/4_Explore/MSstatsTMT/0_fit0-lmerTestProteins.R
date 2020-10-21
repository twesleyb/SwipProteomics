#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: fit protein-wise lmer models and perform statistical inferences
# for a given contrast

## formulae to be fit:
# [+] fx0: Abundance ~ 0 + Condition + (1|Mixture)
# [-] fx1: Aundance ~ 0 + Genotype + BioFraction + (1|Subject)


## prepare the env ------------------------------------------------------------

## input 
results_file = "fit0_lmerTestProtein_results.xlsx" # saved in root/rdata
FDR_alpha = 0.05 # threshold for significance

## load renv
root <- "~/projects/SwipProteomics"
renv::load(root)
 
suppressWarnings({ # FIXME: warnings about replacing imports
	devtools::load_all(root)
})

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

## fit model 0
fm0 <- lmerTest::lmer(fx0, msstats_prot %>% filter(Protein == swip))
#print(summary(fm0, ddf = "Satterthwaite"))

## munge to create contrast matrices for intrafraction comparisons:
condition <- c("ConditionControl","ConditionMutant")
biofraction <- c("F4","F5","F6","F7","F8","F9","F10")
contrasts <- lapply(expandGroups(condition,biofraction), function(x) {
			  getContrast(fm0,x[1],x[2]) })
names(contrasts) <- sapply(contrasts, function(x) {
				 paste(names(x)[x==+1],names(x)[x==-1],sep="-") })

# check the results for swip
results <- lmerTestProtein(swip,fx0,msstats_prot,contrasts)
results$stats %>% knitr::kable()


## loop to fit all proteins ----------------------------------------------------

prots = unique(as.character(msstats_prot$Protein))
#prots = sample(prots, 23)

n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(cores=n_cores)

results_list <- foreach(protein = prots) %dopar% {
	suppressMessages({
	  try(lmerTestProtein(protein, fx0, msstats_prot, contrasts),silent=TRUE)
	})
} # EOL


## process results ------------------------------------------------------------

# collect results
idx <- unlist(sapply(results_list,class)) != "try-error"
filt_list <- results_list[which(idx)]
results_df <- bind_rows(sapply(filt_list,"[","stats")) 

# drop singular
results_df <- results_df %>% filter(!isSingular)
results_df$isSingular <- NULL

## annotate with gene symbols
idx <- match(results_df$protein,gene_map$uniprot)
results_df <- tibble::add_column(results_df,
  				 symbol=gene_map$symbol[idx],
  				 .after="protein")

## adjust pvals 
results_df <- results_df %>% 
	group_by(contrast) %>% 
	mutate("Padjust"=p.adjust(Pvalue,"BH"))


# sort cols
results_df <- results_df %>% 
	select(protein,symbol,contrast,log2FC,percentControl,Pvalue,Padjust,SE,DF)

# sort rows
results_df <- results_df %>% arrange(Pvalue)

# status
message("\nTotal number of significant proteins: ",
	sum(results_df$Padjust < FDR_alpha))


## save results ----------------------------------------------------------------

# results split by BioFraction

results_list <- results_df %>% group_by(contrast) %>% group_split()

# munge to get list names
namen <- sapply(strsplit(sapply(results_list, function(x) unique(x$contrast)),
			 "\\."),"[",3) # get third element after split, F#
names(results_list) <- namen 

# sort the results by biofraction
results_list <- results_list[biofraction]

# FIXME: where is strange class from?
# class(results_list)
# [1] "vctrs_list_of" "vctrs_vctr"    "list"    
class(results_list) <- "list"

# save as excel
myfile <- file.path(root,"tables",results_file)
write_excel(results_list,myfile)

# status
message("\nSummary of significant proteins for intrafraction comparisons:")
df <- data.frame(nsig=sapply(results_list, function(x) sum(x$Padjust < FDR_alpha)))
knitr::kable(df)
