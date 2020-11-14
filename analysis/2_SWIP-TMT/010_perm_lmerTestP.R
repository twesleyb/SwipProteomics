#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: analysis of module-level changes

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root)

# load project
devtools::load_all(root)

# load project's data
data(gene_map)
data(partition)
data(washc_prots)
data(msstats_prot)

# other imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(doParallel)
})


## prepare the data -----------------------------------------------------------

# all modules
names(partition) <- sample(names(partition))
all_modules <- split(names(partition),partition)
names(all_modules) <- paste0("M",names(all_modules))

# modules -- without M0
modules <- all_modules[-1]

# annotate data with module membership
msstats_filt <- msstats_prot %>% filter(Protein %in% names(partition))
msstats_filt$Module <- paste0("M", partition[msstats_filt$Protein])


## fit WASH complex as an example ----------------------------------------------

# the formula to be fit:
fx1 <- "Abundance ~ 1 + Condition + (1|Mixture) + (1|Protein)"

## fit the model:
fit <- lmerTest::lmer(fx1, msstats_filt %>% filter(Protein %in% washc_prots))

lT <- lme4::fixef(fit)
lT[] <- 0 
lT[grepl("Control",names(lT))] <- -1/7
lT[grepl("Mutant",names(lT))] <- +1/7


## module level analysis for all modules --------------------------------------

# register parallel backend
n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(n_cores)

# loop through all modules
results_list <- foreach(module = names(modules)) %dopar% {
  # fit full model
  input <- list(fx1, data = msstats_filt %>% filter(Module == module))
  suppressMessages({ # about boundary fits
    fm <- try(do.call(lmerTest::lmer, input), silent = TRUE)
  })
  # test the contrast
  df <- lmerTestContrast(fm, lT) %>%
	  mutate(Contrast = "Mutant-Control") %>%
	  unique()
  # return the results
  return(df)
}
names(results_list) <- names(modules)


## collect results
results_df <- do.call(rbind, results_list) %>% 
	as.data.table(keep.rownames = "Module")
