#!/usr/bin/env Rscript

# title: SwipProteomics
# description: preprocessing and statistical analysis of WASHC1 BioID proteomics
# experiment performed by JC
# author: twab

## ---- imports

suppressPackageStartupMessages({
	library(dplyr) # for manipulating data
	library(data.table) # for working with tables
	library(doParallel) # for parallel processing
})

## ---- input

library(tidyProt) # soderling-lab/tidyProt for statistical fun
library(geneLists) # soderling-lab/geneLists for gene mapping fun


## ---- load the raw data

# input data from soderling-lab/SwipProteomics:
datadir <- system.file("extdata", "BioID.zip", package = "SwipProteomics")

unzip(datadir) # unzip into cwd

# list.files(tools::file_path_sans_ext(basename(datadir)))

# read TMT.csv data into R with data.table::fread
protein <- fread(list.files("BioID", pattern="protein",full.names=T))


## ---- create sample meta data from columns in protein


## ---- tidy the data

# melt abundance matrix into long format
abundance_cols <- grepl("Abundance",colnames(peptides))
id_variables <- colnames(peptides)[!abundance_cols]

tidy_pep <- peptides %>% 
	# melt the Abundance columns into a tidy df
	reshape2::melt(id.var=id_variables, 
				 variable.name="Sample", 
				 value.name="Intensity") %>%
        # merge with sample data by column 'Sample'
        left_join(samples, by = "Sample")

dm = tidy_pep %>% 
	reshape2::dcast(Sample ~ Accession + Sequence, value.var="Intensity", 
			fun="sum") %>%
	as.data.table() %>% as.matrix(rownames="Sample")
pca(

# NOTE: samples should contain the following cols:
# col Sample - should match Abundance columns
stopifnot(all(colnames(peptides)[abundance_cols] %in% samples$Sample))

# perform SL normalization
sl_peptide <- normSL(tidy_pep)

# protein summarization
tidy_prot <- sumProt(sl_peptide)

# NOTE: additional meta data is lost at this step unless its added back

# perform IRS normalization
norm_prot <- normIRS(tidy_prot, controls = "SPQC")

# now we can filter prots
filt_prot <- norm_prot %>% subset(Peptides > 1)

# drop QC before fitting models
swip_tmt <- filt_prot %>% subset(Condition != "SPQC")

# the linear mixed model to be fit:
fx <- log2(Intensity) ~ 0 + Condition + (1|BioFraction) + (1|Mixture)

# fit the model with lmerTest
fm <- lmerTest::lmer(fx, data = swip_tmt %>% subset(Accession == swip))

# get a contrast
LT <- getContrast(fm, "Mutant","Control")

# assess the statistical comparison
lmerTestContrast(fm, contrast = LT) %>% knitr::kable()


## ---- end

# before performing stats for all proteins, we should remove any with missing
# values... cast the data into a matrix and then check for prots w/ missing vals
dm = swip_tmt %>% dcast(Mixture + Condition + BioFraction ~ Accession, 
			value.var = "Intensity") %>% 
                  mutate(Sample = paste(Mixture,Condition,BioFraction)) %>%
		  select(-Mixture, - Condition, - BioFraction) %>%
		  select(Sample,everything()) %>% as.data.table() %>%
		  as.matrix(rownames="Sample")
missing_vals <- names(which(apply(dm,2,function(x) any(is.na(x)))))


swip_tmt <- swip_tmt %>% filter(!(Accession %in% missing_vals))


## ---- loop to test all prots

prots <- unique(swip_tmt$Accession)
n <- length(prots)

# fit the model with some lmer control
lmer_control <- lme4::lmerControl(check.conv.singular="ignore", 
				  check.conv.grad="ignore")

# register parallel backend
registerDoParallel(parallel::detectCores()-1)

# this is computationally intense...
res_list <- foreach(prot = prots) %dopar% {
	fm <- lmerTest::lmer(fx, control = lmer_control,
			     data = swip_tmt %>% subset(Accession == prot))
	res <- lmerTestContrast(fm, contrast = LT)
	return(res)
}
names(res_list) <- prots

# collect the results
results <- bind_rows(res_list, .id = "Protein") %>%
	# calc padjust and FDR
	mutate(Padjust = p.adjust(Pvalue, method = "bonferroni")) %>%
	mutate(FDR = p.adjust(Pvalue, method = "hochberg")) %>%
	# annot with geneLists::getIDs
       	mutate(Symbol = getIDs(Protein, "uniprot", "symbol", "mouse")) %>%
	# sort
	arrange(FDR) 

# check top results
head(results) %>% knitr::kable()

# check n sig at FDR < 0.05
knitr::kable(data.table(nsig = sum(results$FDR < 0.05)))
