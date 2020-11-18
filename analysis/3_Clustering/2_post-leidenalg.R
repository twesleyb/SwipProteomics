#!/usr/bin/env Rscript

root <- "~/projects/SwipProteomics"
renv::load(root)

library(dplyr)
library(data.table)

devtools::load_all(root)

myfile <- file.path(root,"rdata","cpm_partition.csv")
part_dm <- fread(myfile,drop=1) %>% as.matrix()
part_list = unlist(apply(part_dm,1,list),recursive=FALSE)


# module quality:
data(washc_prots)

proteins <- unique(msstats_prot$Protein)
prots <- sample(proteins,5)

part <- part_list[[99]]+1

q <- sum(calcModQuality(part, msstats_prot))
q

calcModQuality <- function(part, msstats_prot) {
	# module function
	fx <- Abundance ~ 0 + BioFraction:Genotype + (1|Mixture) + (1|Protein)
	q <- vector("numeric", length(unique(part)))
	for (m in unique(part)) {
		prots <- names(part[part==m])
		fm <- lmerTest::lmer(fx, msstats_prot %>% subset(Protein %in% prots))
	        vp <- getVariance(fm)
	        q[m] <- as.numeric(vp["Fixed"]/(vp["Mixture"] + vp["Protein"] + vp["Residual"]))
	}
	return(q)
}

# maximize variance attributable to fixef (Condition)
# minimize variance attributable to mixef (Mixture + Protein)
# minimize unexplained variance (minimize Residual)


