#!/usr/bin/env Rscript


## ---- prepare the env

root = "~/projects/SwipProteomics"

renv::load(root)
devtools::load_all(root)


## ---- load the data

data(swip)
data(swip_tmt)
data(swip_gene_map)


## ---- the model to be fit


# Simple, Linear model
prot <- swip
fx0 <- Abundance ~ Condition
fm0 <- lm(fx0, data = swip_tmt %>% subset(Protein == prot))
LT <- getContrast(fm0, "Mutant","Control")
lmTestContrast(fm0, LT) %>%
	dplyr::mutate(Contrast='Mutant-Control') %>%
	dplyr::mutate(Pvalue = formatC(Pvalue)) %>% unique() %>% knitr::kable()


# loop for all prots
proteins <- unique(swip_tmt$Protein)

pbar <- txtProgressBar(max=length(proteins),style=3)
results_list <- list()
for (prot in proteins) {
  fm <- lm(fx0, swip_tmt %>% subset(Protein == prot))
  LT <- getContrast(fm, "Mutant","Control")
  res <- lmTestContrast(fm, LT) %>%
	  dplyr::mutate(Contrast='Mutant-Control') %>% unique()
  results_list[[prot]] <- res
  setTxtProgressBar(pbar,match(prot,proteins))
}
close(pbar)

# collect results and calc FDR
df <- dplyr::bind_rows(results_list,.id="Protein") %>%
	dplyr::mutate(FDR = p.adjust(Pvalue, method = "fdr")) %>%
	dplyr::mutate(Symbol = gene_map$symbol[match(Protein, gene_map$uniprot)]) %>%
	dplyr::arrange(Pvalue)


# check nsig
sum(df$FDR<0.05)

# check top sig
head(df$Symbol)
