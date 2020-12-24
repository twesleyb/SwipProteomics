#!/usr/bin/env Rscript


## ---- prepare the env

root = "~/projects/SwipProteomics"

renv::load(root)
devtools::load_all(root)


## ---- load the data

data(swip)
data(msstats_prot)
data(swip_gene_map)

# add subject annot and calc rel intensity
swip_tmt <- msstats_prot %>%
	dplyr::mutate(Subject = as.numeric(interaction(Mixture,Genotype))) %>%
	group_by(Protein) %>%
	dplyr::mutate(rel_Intensity = Intensity/sum(Intensity))


## ---- linear model

# NOTE: swip_tmt is msstats_prot!

prot <- swip

fx0 <- log2(rel_Intensity) ~ 0 + Condition

fm0 <- lm(fx0, data = swip_tmt %>% subset(Protein == prot))
LT <- getContrast(fm0, "Mutant","Control")
lmTestContrast(fm0, LT) %>% dplyr::mutate(Contrast='Mutant-Control') %>% unique() %>% knitr::kable()

# mixed-model with Mixture
fx1 <- log2(rel_Intensity) ~ 0 + Condition + (1|Mixture)
fm1 <- lmerTest::lmer(fx1, data = swip_tmt %>% subset(Protein == prot))
LT <- getContrast(fm1, "Mutant","Control")
lmerTestContrast(fm1, LT) %>% dplyr::mutate(Contrast='Mutant-Control') %>% unique() %>% knitr::kable()

# mixed-model with Mixture and Subject
fx2 <- log2(rel_Intensity) ~ 0 + Condition + (1|Mixture) + (1|Subject)
fm2 <- lmerTest::lmer(fx2, data = swip_tmt %>% subset(Protein == prot))
LT <- getContrast(fm2, "Mutant","Control")
lmerTestContrast(fm2, LT) %>% dplyr::mutate(Contrast='Mutant-Control') %>% unique() %>% knitr::kable()

# mixed-model with Mixture and Subject nested within Mixture
fx3 <- log2(rel_Intensity) ~ 0 + Condition + (1|Mixture) + (1|Mixture:Subject)
fm3 <- lmerTest::lmer(fx3, data = swip_tmt %>% subset(Protein == prot))
LT <- getContrast(fm3, "Mutant","Control")
lmerTestContrast(fm3, LT) %>% dplyr::mutate(Contrast='Mutant-Control') %>% unique() %>% knitr::kable()

# mixed-model with Subject nested within Mixture + Subject
# fails to converge: we should avoid this model
#fx4 <- Abundance ~ 1 + Condition + (1|Mixture) + (1|Mixture:Subject) + (1|Subject)
#fm4 <- lmerTest::lmer(fx4, data = swip_tmt %>% subset(Protein == prot))
#LT <- getContrast(fm4, "Mutant","Control")
#lmerTestContrast(fm4, LT) %>% dplyr::mutate(Contrast='Mutant-Control') %>% unique() %>% knitr::kable()

# mixed-model with only Subject
fx5 <- log2(rel_Intensity) ~ 0 + Condition + (1|Subject)

fm5 <- lmerTest::lmer(fx5, data = swip_tmt %>% subset(Protein == prot))
LT <- getContrast(fm5, "Mutant","Control")
lmerTestContrast(fm5, LT) %>% dplyr::mutate(Contrast='Mutant-Control') %>% unique() %>% knitr::kable()

# mixed-model with Subject nested within mixture
fx6 <- log2(rel_Intensity) ~ 0 + Condition + (1|Mixture:Subject)

fm6 <- lmerTest::lmer(fx6, data = swip_tmt %>% subset(Protein == prot))
LT <- getContrast(fm6, "Mutant","Control")
lmerTestContrast(fm6, LT) %>% dplyr::mutate(Contrast='Mutant-Control') %>% unique() %>% knitr::kable()


# anova analysis
anova(fm6, fm5, fm3, fm2, fm1, fm0)


# loop for all prots
proteins <- unique(swip_tmt$Protein)
pbar <- txtProgressBar(max=length(proteins),style=3)
lmer_control <- lme4::lmerControl(check.conv.singular="ignore")
results_list <- list()
for (prot in proteins) {
  fm <- lmerTest::lmer(fx6, swip_tmt %>% subset(Protein == prot), control = lmer_control)
  LT <- getContrast(fm, "Mutant","Control")
  res <- lmerTestContrast(fm, LT) %>%
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


sum(df$FDR<0.05)

head(df$Symbol)
