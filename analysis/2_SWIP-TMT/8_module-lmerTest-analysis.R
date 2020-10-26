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
data(swip)
data(gene_map)
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

# remove any unclustered proteins
msstats_filt <- msstats_prot %>% filter(Protein %in% names(partition))

# annotate data with module membership
msstats_filt$Module <- paste0("M", partition[msstats_filt$Protein])


## fit WASH complex as an example ----------------------------------------------

# wash prots
washc <- gene_map$uniprot[grepl("Washc*", gene_map$symbol)]

# the formula to be fit:
fx1 <- formula(paste(c(
  "Abundance ~ 0 + Genotype:BioFraction + (1|Mixture)",
  " + (1|Protein)"
), collapse = " "))

message("\nlmer: ", as.character(fx1)[2], " ~ ", as.character(fx1)[3])

# fit the model:
fm1 <- lmerTest::lmer(fx1, data = msstats_filt %>% filter(Protein %in% washc))

summary(fm1, ddf = "Satterthwaite")

# goodness of fit
r.squaredGLMM.merMod(fm1) %>% knitr::kable()

## build a contrast
contrast <- lme4::fixef(fm1)
contrast[] <- 0
idx <- which(grepl("Control", names(contrast)))
contrast[idx] <- -1 / length(idx)
idy <- which(grepl("Mutant", names(contrast)))
contrast[idy] <- +1 / length(idy)

# asses contrast:
lmerTestContrast(fm1, contrast) %>%
  mutate(Contrast = "Mutant-Control") %>%
  mutate(isSingular = NULL) %>%
  mutate("nProteins" = length(washc)) %>%
  unique() %>%
  knitr::kable()


## goodness of fit!
qqnorm(residuals(fm1))


## module level analysis for all modules --------------------------------------

# loop through all modules
modules <- split(names(partition), partition)[-1]
names(modules) <- paste0("M", names(modules))

# register parallel backend
n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(n_cores)

results_list <- foreach(module = names(modules)) %dopar% {
  input <- list(fx1, data = msstats_filt %>% filter(Module == module))
  suppressMessages({ # about boundary fits
    fm <- try(do.call(lmerTest::lmer, input), silent = TRUE)
  })
  df <- lmerTestContrast(fm, contrast) %>%
	  mutate(Contrast = "Mutant-Control") %>%
	  unique()
  return(df)
}
names(results_list) <- names(modules)

## collect results
results_df <- do.call(rbind, results_list)

# singular results will be removed
warning(
  sum(results_df$isSingular),
  " modules with singular fits will be removed."
)

# drop singular and perform p.adjust
results_df <- results_df %>%
  filter(!isSingular) %>%
  arrange(Pvalue) %>%
  as.data.table(keep.rownames = "Module") %>%
  mutate(
    FDR = p.adjust(Pvalue, method = "BH"),
    PAdjust = p.adjust(Pvalue, method = "bonferroni")
  )

# annotate with module size
module_size <- sapply(modules,length)[results_df$Module]
results_df <- tibble::add_column(results_df,Size=module_size,.after="Module")

## examine top results
knitr::kable(head(results_df))

message("Number of significant modules (Bonferroni<0.05): ",
	sum(results_df$PAdjust<0.05))


# annotate with protein identifiers
results_df$Proteins <- sapply(modules[results_df$Module],paste,collapse=";")
results_df$Symbols <- sapply(modules[results_df$Module], function(x) {
	       paste(gene_map$symbol[match(x,gene_map$uniprot)],collapse=";")
	})


# save results as an excel workbook
df <- data.table(Protein = names(partition), Module = partition)
results_list <- list("Partition" = df, "Mutant-Control" = results_df)
myfile <- file.path(root, "tables", "S3_Module_Results.xlsx")
write_excel(results_list, myfile)

# save
module_results <- results_df
myfile <- file.path(root, "data", "module_results.rda")
save(module_results, file = myfile, version = 2)

# save
myfile <- file.path(root, "data", "fx1.rda")
save(fx1, file = myfile, version = 2)
