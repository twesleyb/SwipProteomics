#!/usr/bin/env Rscript

# title: Swip TMT Proteomics
# description: preprocessing and statistical analysis of Swip TMT proteomics
# experiment performed by JC
# author: Tyler W A Bradshaw

## ---- input

library(tidyProt) # twesleyb/tidyProt for statistical fun
library(geneLists) # twesleyb/geneLists for gene mapping fun

# input data from twesleyb/SwipProteomics:
system.file("extdata", "TMT.zip", package = "SwipProteomics")

## ---- prepare the R workspace

# load required packages and functions
suppressPackageStartupMessages({
	library(dplyr) # for manipulating data
	library(data.table) # for working with tables
	library(doParallel) # for parallel processing
})

## ---- load the raw data and sample info

# datadir [1] ~/projects/SwipProteomics/data"

# extract the raw TMT data from zipped file
unzip(myfile, exdir=downdir) # unzip into root/downloads/

# read TMT.csv data into R with data.table::fread
myfile <- file.path(downdir, tools::file_path_sans_ext(zip_file), input_data)
peptides <- data.table::fread(myfile)

# THIS IS THE RAW DATA FROM PD! 
# The file size is fairly large once unzipped.
o = object.size(peptides)
message("object.size(peptides):"); print(o, units="auto")

# load sample information
myfile <- file.path(downdir,tools::file_path_sans_ext(zip_file),input_meta)
samples <- data.table::fread(myfile)

# This fragment of code does it all
abundance_cols <- grepl("Abundance",colnames(peptides))
id_variables <- colnames(peptides)[!abundance_cols]
tidy_pep <- peptides %>% 
	# melt the Abundance columns into a tidy df
	reshape2::melt(id.var=id_variables, 
				 variable.name="Sample", 
				 value.name="Intensity") %>%
        # merge with sample data by column 'Sample'
        left_join(samples, by = "Sample")

# NOTE: samples should contain the following cols:
# col Sample - should match Abundance columns
stopifnot(all(colnames(peptides)[abundance_cols] %in% samples$Sample))

# The tidy data is tidy, but LARGE! this is why its stored the way it is in the
# PD export -- its the most efficient
print(object.size(tidy_pep), units="auto")

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

library(tidyProt)

fx <- log2(Intensity) ~ 0 + Condition + (1|BioFraction) + (1|Mixture)

fm <- lmerTest::lmer(fx, data = swip_tmt %>% subset(Accession == swip))

LT <- getContrast(fm, "Mutant","Control")
lmerTestContrast(fm, contrast = LT) %>% knitr::kable()

# why does that feel too easy...

dm = swip_tmt %>% dcast(Mixture + Condition + BioFraction ~ Accession, 
			value.var = "Intensity") %>% 
                  mutate(Sample = paste(Mixture,Condition,BioFraction)) %>%
		  select(-Mixture, - Condition, - BioFraction) %>%
		  select(Sample,everything()) %>% as.data.table() %>%
		  as.matrix(rownames="Sample")

missing_vals <- names(which(apply(dm,2,function(x) any(is.na(x)))))

# need to account for missing values?

swip_tmt <- swip_tmt %>% filter(!(Accession %in% missing_vals))

prots <- unique(swip_tmt$Accession)
n <- length(prots)

lmer_control <- lme4::lmerControl(check.conv.singular="ignore", 
				  check.conv.grad="ignore")
library(doParallel)
registerDoParallel(parallel::detectCores()-1)

# this is computationally intense...
# its like my computer used to be able to handle it but not its more
# diffucult...
res_list <- foreach(prot = prots) %dopar% {
	fm <- lmerTest::lmer(fx, control = lmer_control,
			     data = swip_tmt %>% subset(Accession == prot))
	res <- lmerTestContrast(fm, contrast = LT)
	return(res)
}
names(res_list) <- prots

results <- bind_rows(res_list, .id = "Protein") %>%
	mutate(Padjust = p.adjust(Pvalue, method = "bonferroni")) %>%
	mutate(FDR = p.adjust(Pvalue, method = "hochberg")) %>%
       	mutate(Symbol = getIDs(Protein, "uniprot", "symbol", "mouse")) %>%
	arrange(FDR) 

head(results) %>% knitr::kable()

knitr::kable(data.table(nsig = sum(results$FDR < 0.05)))
