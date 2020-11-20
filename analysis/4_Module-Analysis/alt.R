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
data(ne_surprise_surprise_partition)
data(washc_prots)
#data(msstats_prot)
data(swip_tmt)


# other imports
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lmerTest)
  library(data.table)
  library(doParallel)
})


## prepare the data -----------------------------------------------------------

# all modules
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))

# data.table describing partition
part_df <- data.table(Protein = names(partition), Module = paste0("M",partition))
part_df$Size <- sapply(modules,length)[part_df$Module]
part_df$Symbol <- mapID(part_df$Protein,"uniprot","symbol")
part_df$Entrez <- mapID(part_df$Protein,"uniprot","entrez")


## examine the washc_module
swip_tmt$Abundance <- log2(swip_tmt$Intensity)
swip_tmt$Condition <- interaction(swip_tmt$Genotype,swip_tmt$BioFraction)

fx <- Abundance ~ 1 + Condition + (1|Protein)
fm <- lmerTest::lmer(fx,data=swip_tmt%>% subset(Protein %in% washc_prots))

vp <- getVariance(fm)
knitr::kable(t(vp/sum(vp))) # pve by each covariate


# network summary:
nProts <- formatC(length(names(partition)),big.mark=",")
kModules <- length(modules)
pClustered <- round(sum(partition!=0)/length(partition),3)
medSize <- median(sapply(modules,length))
knitr::kable(cbind(nProts,kModules,pClustered,medSize))


## fit WASH complex as an example ----------------------------------------------

# fit the model:
fx <- "Abundance ~ 1 + Condition + (1|Protein)"
fm <- lmer(fx, swip_tmt %>% subset(Protein %in% washc_prots))

## assess contrast
LT <- getContrast(fm,"Mutant","Control")
lmerTestContrast(fm, LT) %>% mutate(Contrast="Mutant-Control") %>% unique()

# its the same result!
Ftest <- calcSatterth(fm, LT) %>% unlist() 
as.data.table(t(Ftest)) %>% mutate(pvalue=formatC(pvalue)) %>% knitr::kable()

## goodness of fit
r.squaredGLMM.merMod(fm) %>% knitr::kable()


## module level analysis for all modules --------------------------------------

# register parallel backend
n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(n_cores)

message("\nAssessing module-level contrasts with lmerTest.")

t0 <- Sys.time()

modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))


too_small = names(which(sapply(modules,length) < 5))
modules = modules[!names(modules) %in% too_small]


# loop through all modules
results_list <- foreach(module = names(modules)) %dopar% {
  # fit full model
  input <- list(formula = fx)
  input[["data"]] <- swip_tmt %>% filter(Protein %in% modules[[module]])
  input[["control"]] <- lme4::lmerControl(check.conv.singular = "ignore",
					  calc.derivs = FALSE, 
					  check.rankX = "stop.deficient")
  fm <- do.call(lmerTest::lmer, input)
  # test the contrast
  df <- lmerTestContrast(fm, LT) %>% mutate(Contrast = "Mutant-Control")
  df <- df %>% unique()
  # add fstatitistic
  df$Fstat <- calcSatterth(fm,LT)[["Fstat"]]
  # return the results
  return(unique(df))
}
names(results_list) <- names(modules)

# summary
t1 <- Sys.time()
message("\nTime to analyze ", length(results_list)," modules:")
difftime(t1,t0)

## collect results
results_df <- do.call(rbind, results_list) %>% 
	as.data.table(keep.rownames = "Module")

# examine the results
k <- unique(results_df$Module)
m <- modules[k]
p <- partition

# summary:
message("\nFinal number of modules : ", length(k))
message("\nFinal Median module size: ", median(sapply(modules,length)))
message("\nWashc4 assigned to module: ", paste0("M",partition[swip]))

# annotate results with module size
module_size <- sapply(modules,length)[results_df$Module]
results_df <- tibble::add_column(results_df,Size=module_size,.after="Module")

results_df <- results_df %>% arrange(Pvalue)

# annotate results with protein identifiers
results_df$Proteins <- sapply(modules[results_df$Module],paste,collapse=";")
results_df$Symbols <- sapply(modules[results_df$Module], function(x) {
	       paste(gene_map$symbol[match(x,gene_map$uniprot)],collapse=";")
	})

results_df$Padjust <- p.adjust(results_df$Pvalue,"bonferroni")
results_df$FDR <- p.adjust(results_df$Pvalue, "BH")

# sort columns
results_df <- results_df %>% 
	select(Module,Size,Contrast,log2FC,percentControl,
	       Pvalue,FDR,Padjust,Tstatistic,SE,DF,Symbols)

# save results as an excel workbook
results_list <- list("Partition" = part_df, "Mutant-Control" = results_df)
myfile <- file.path(root, "tables", "S4_SWIP_Module_Results.xlsx")
write_excel(results_list, myfile)

# save module_results as rda
module_results <- results_df
myfile <- file.path(root, "data", "module_results.rda")
save(module_results, file = myfile, version = 2)
