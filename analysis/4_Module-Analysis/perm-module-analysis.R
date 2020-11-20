#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: analysis of module-level changes

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root,quiet=TRUE)

# load project
devtools::load_all(root,quiet=TRUE)

# load project's data
data(partition)
data(msstats_prot)


# other imports
suppressPackageStartupMessages({
  library(dplyr)
  library(lmerTest)
  library(data.table)
  library(doParallel)
})


## prepare the data -----------------------------------------------------------

# all modules
names(partition) <- sample(names(partition))
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))

# register parallel backend
n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(n_cores)

# loop through all modules
perm_list <- list()
pbar <- txtProgressBar(max=1000,style=3)
for (perm in seq(1000)) {
  results_list <- foreach(module = names(modules)) %dopar% {
    # fit full model
    input <- list()
    input[["formula"]] <- formula(Abundance ~ 1 + Condition + (1|Mixture) + (1|Protein))
    input[["data"]] <- msstats_prot %>% filter(Protein %in% modules[[module]])
    fm <- lmerFit(input)
    # test the contrast
    LT <- getContrast(fm,"Mutant","Control")
    df <- as.data.table(calcSatterth(fm,LT))
    df$log2FC <- LT %*% lme4::fixef(fm)
    df$Contrast <- "Mutant-Control"
    return(df)
  } #EOL
  names(results_list) <- names(modules)
  df <- dplyr::bind_rows(results_list,.id="Module")
  perm_list[[perm]] <- df$Fstat
  setTxtProgressBar(pbar,perm)
} #EOL
close(pbar)

myfile<- file.path(root,"rdata","perm_list.rda")
save(perm_list,file=myfile,version=2)
