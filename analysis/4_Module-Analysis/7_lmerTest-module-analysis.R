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
data(washc_prots)
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
names(partition) <- sample(names(partition))
all_modules <- split(names(partition),partition)
names(all_modules) <- paste0("M",names(all_modules))

# modules -- without M0
modules <- all_modules[-1]

# data.table describing partition
part_df <- data.table(Protein = names(partition), Module = paste0("M",partition))
part_df$Size <- sapply(all_modules,length)[part_df$Module]
part_df$Symbol <- mapGeneIDs(part_df$Protein,"uniprot","symbol")
part_df$Entrez <- mapGeneIDs(part_df$Protein,"uniprot","entrez")


## examine the washc_module
fx <- Abundance ~ 1 + Condition + (1|Protein) + (1|Mixture)
fm <- lmerTest::lmer(fx,data=msstats_prot%>% subset(Protein %in% washc_prots))

vp <- getVariance(fm)
vp/sum(vp)

# network summary:
nProts <- formatC(length(names(partition)),big.mark=",")
kModules <- length(modules)
pClustered <- round(sum(partition!=0)/length(partition),3)
medSize <- median(sapply(modules,length))
knitr::kable(cbind(nProts,kModules,pClustered,medSize))


## fit WASH complex as an example ----------------------------------------------


# fit the model:
fx <- "Abundance ~ 1 + Condition + (1|Mixture) + (1|Protein)"
fm <- lmerTest::lmer(fx, msstats_prot %>% subset(Protein %in% washc_prots))

## assess contrast
LT <- getContrast(fm,"Mutant","Control")
lmerTestContrast(fm, LT) %>% mutate(Contrast="Mutant-Control") %>% 
	unique() %>% knitr::kable()


calcSatter


## examine the model
p = "Pr(>|t|)"
dm <- summary(fm, ddf = "Satterthwaite")[["coefficients"]]
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
message("R2m: Marginal; variation explained by fixed effects.")
message("R2c: Conditional; total variation explained by the model.")
r.squaredGLMM.merMod(fm) %>% knitr::kable(format="markdown")


## module level analysis for all modules --------------------------------------

# register parallel backend
n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(n_cores)

message("\nAssessing module-level contrasts with lmerTest.")

t0 <- Sys.time()

# loop through all modules
results_list <- foreach(module = names(modules)) %dopar% {
  # fit full model
  input <- list(fx, data = msstats_filt %>% filter(Module == module))
  suppressMessages({ # about boundary fits
    fm <- try(do.call(lmerTest::lmer, input), silent = TRUE)
  })
  ## if singular, fit reduced model
  #if (lme4::isSingular(fm)) {
  #  input <- list(fx0, data = msstats_filt %>% filter(Module == module))
  #  suppressMessages({
  #    fm <- try(do.call(lmerTest::lmer, input), silent = TRUE)
  #  })
  #}
  # test the contrast
  df <- lmerTestContrast(fm, LT) %>%
	  mutate(Contrast = "Mutant-Control") %>%
	  unique()
  # return the results
  return(df)
}

names(results_list) <- names(modules)

# summary
t1 <- Sys.time()
message("\nTime to analyze ", length(results_list)," modules:")
difftime(t1,t0)

## collect results
results_df <- do.call(rbind, results_list) %>% 
	as.data.table(keep.rownames = "Module")

# singular results will be removed
warning(
  sum(results_df$isSingular),
  " modules with singular fits."
)

# drop singular and perform p.adjust
results_df <- results_df %>%
  arrange(Pvalue) %>%
  mutate(
    FDR = p.adjust(Pvalue, method = "BH"),
    Padjust = p.adjust(Pvalue, method = "bonferroni")
  )

# save modules
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))
myfile <- file.path(root,"data","modules.rda")
save(modules,file=myfile,version=2)

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
