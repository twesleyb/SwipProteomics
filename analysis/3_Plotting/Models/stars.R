#!/usr/bin/env Rscript

root = "~/projects/SwipProteomics"
renv::load(root)

devtools::load_all(quiet=TRUE)

data(swip)
data(msstats_prot)

subdat <- function(protein,required=c("msstats_prot")) {
	# subset the protein-level data for a given protein
	# FIXME: error message about env is not enformative
	stopifnot(all(sapply(required,exists)))
	require(dplyr, quietly=TRUE)
	return(msstats_prot %>% dplyr::filter(Protein == protein))
} #EOF

fm <- lmerTest::lmer("Abundance ~ 0 + (1|BioFraction) + Genotype", subdat(swip))



