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
  library(ggplot2)
  library(lmerTest)
  library(data.table)
  library(doParallel)
})


## functions ------------------------------------------------------------------

mapGeneIDs <- function(id_in,identifiers,new_ids) {
  # requires: gene_map
  # mapping between gene identifiers
  # prots, uniprot, symbols
  if (length(id_in)==1) { type = "one" }
  if (length(id_in)>1) { type = "many" }
idx <- switch(type,
	      one = grepl(id_in, gene_map[[identifiers]]),
	      many = match(id_in,gene_map[[identifiers]]))
return(gene_map[[new_ids]][idx])
} #EOF


## prepare the data -----------------------------------------------------------

# all modules
all_modules <- split(names(partition),partition)
names(all_modules) <- paste0("M",names(all_modules))

# modules -- without M0
modules <- all_modules[-1]

# data.table describing partition
part_df <- data.table(Protein = names(partition), Module = paste0("M",partition))
part_df$Size <- sapply(all_modules,length)[part_df$Module]
part_df$Symbol <- mapGeneIDs(part_df$Protein,"uniprot","symbol")
part_df$Entrez <- mapGeneIDs(part_df$Protein,"uniprot","entrez")

# remove any unclustered proteins
msstats_filt <- msstats_prot %>% filter(Protein %in% names(partition))

# annotate data with module membership
msstats_filt$Module <- paste0("M", partition[msstats_filt$Protein])

# summary:
message("\nNumber of modules : ", length(modules))
message("\nPercent clustered : ", round(sum(partition!=0)/length(partition),3))
message("\nMedian module size: ", median(sapply(modules,length)))


## fit WASH complex as an example ----------------------------------------------

# wash complex prots
washc_prots <- mapGeneIDs("Washc*","symbol","uniprot")

# wash complex genes
washc_genes <- sort(mapGeneIDs(washc_prots,"uniprot","symbol"))


## the formula to be fit:
fx1 <- formula(paste(c(
  "Abundance ~ 0 + Genotype:BioFraction + (1|Mixture) + (1|Protein)"
), collapse = " "))

message("\nlmer fit to WASH complex (", paste(washc_genes,collapse=", "),") proteins:\n",
	">>>\t",as.character(fx1)[2], " ~ ", as.character(fx1)[3])


## fit the model:
fm1 <- lmerTest::lmer(fx1, 
		      data = msstats_filt %>% filter(Protein %in% washc_prots))

## examine the model
p = "Pr(>|t|)"
dm <- summary(fm1, ddf = "Satterthwaite")[["coefficients"]]
df <- t(apply(dm,1,round,3)) %>% as.data.table(keep.rownames="Term")
df <- df %>% mutate(Term=gsub("Genotype","",Term))
idx <- order(sapply(strsplit(df$Term,"\\:"),"[",1))
df <- df[idx,]
df[[p]] <- formatC(dm[,which(colnames(dm) == p)])
colnames(df)[colnames(df) == "Std. Error"] <- "SE"
colnames(df)[colnames(df) == "df"] <- "DF"
colnames(df)[colnames(df) == "t value"] <- "Tvalue"
colnames(df)[colnames(df) == p ] <- "Pvalue"
df %>% knitr::kable(format="markdown")


## goodness of fit
r.squaredGLMM.merMod(fm1) %>% knitr::kable(format="markdown")
message("R2m: Marginal; variation explained by fixed effects.")
message("R2c: Conditional; total variation explained by the model.")


## build a contrast for assessing 'Mutant-Control' comparison
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
  mutate("nProteins" = length(washc_prots)) %>%
  unique() %>%
  knitr::kable()


## module level analysis for all modules --------------------------------------

# register parallel backend
n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(n_cores)

message("\nAssessing module-level contrasts with lmerTest.")

t0 <- Sys.time()

# loop through all modules
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

t1 <- Sys.time()
message("\nTime to analyze ", length(results_list)," modules:")
difftime(t1,t0)

## collect results
results_df <- do.call(rbind, results_list) %>% 
	as.data.table(keep.rownames = "Module")

# singular results will be removed
warning(
  sum(results_df$isSingular),
  " modules with singular fits will be removed."
)

# record singular moduels
drop <- results_df$Module[results_df$isSingular]
part_df$isSingular <- part_df$Module %in% drop

# drop singular and perform p.adjust
results_df <- results_df %>%
  filter(!isSingular) %>%
  arrange(Pvalue) %>%
  mutate(
    FDR = p.adjust(Pvalue, method = "BH"),
    Padjust = p.adjust(Pvalue, method = "bonferroni")
  )

results_df$isSingular <- NULL

# save final modules -- the modules we have fitted models for
final_modules <- unique(results_df$Module)
myfile <- file.path(root,"data","final_modules.rda")
save(final_modules,file=myfile,version=2)

# examine the results
k <- unique(results_df$Module)
m <- modules[k]
p <- partition

# summary:
message("\nFinal number of modules : ", length(k))
message("\nFinal percent clustered : ", round(length(unlist(m))/length(p),3))
message("\nFinal Median module size: ", median(sapply(modules,length)))
message("\nWashc4 assigned to module: ", paste0("M",partition[swip]))

# annotate results with module size
module_size <- sapply(modules,length)[results_df$Module]
results_df <- tibble::add_column(results_df,Size=module_size,.after="Module")

## examine top results
results_df %>% filter(Padjust<0.05) %>% knitr::kable()

message("\nModules with greater than 10% change:")
idx <- abs(results_df$log2FC) > log2(1.1)  
results_df[idx,] %>% knitr::kable()

message("Number of significant modules (Bonferroni<0.05): ",
	sum(results_df$Padjust<0.05))


# annotate results with protein identifiers
results_df$Proteins <- sapply(modules[results_df$Module],paste,collapse=";")
results_df$Symbols <- sapply(modules[results_df$Module], function(x) {
	       paste(gene_map$symbol[match(x,gene_map$uniprot)],collapse=";")
	})


# sort columns
results_df <- results_df %>% 
	select(Module,Size,Contrast,log2FC,percentControl,
	       Pvalue,FDR,Padjust,Tstatistic,SE,DF,Symbols)

# save results as an excel workbook
results_list <- list("Partition" = part_df, "Mutant-Control" = results_df)
myfile <- file.path(root, "tables", "S3_SWIP_Module_Results.xlsx")
write_excel(results_list, myfile)

# save module_results as rda
module_results <- results_df
myfile <- file.path(root, "data", "module_results.rda")
save(module_results, file = myfile, version = 2)

# save formula for module-level contrasts as rda
myfile <- file.path(root, "data", "fx1.rda")
save(fx1, file = myfile, version = 2)
