#!/usr/bin/env Rscript

## ---- input 

os_keep <- c(9606, 10116, 10090)


## ---- prepare the renv

root <- "~/projects/SwipProteomics"
renv::load(root, quiet=TRUE)

# library(SwipProteomics)
devtools::load_all(root, quiet=TRUE)


## ---- load the data

data(swip)
data(gene_map)
data(washc_prots)
data(wash_interactome)

library(getPPIs)
data(musInteractome)


## ---- collect WASHC4 ppis

washc_entrez <- mapID(washc_prots,"uniprot","entrez")
washc1 <- mapID("Washc1","symbol","entrez")


bioid_entrez <- mapID(wash_interactome, "uniprot","entrez")
bioid_entrez[is.na(bioid_entrez)] <- 68112

# collect interactions between swip and wash_interactome proteins
washc1_ppis <- musInteractome %>% 
	subset(osEntrezA %in% washc_entrez | osEntrezB %in% washc_entrez) %>% 
	subset(osEntrezA %in% bioid_entrez & osEntrezB %in% bioid_entrez)

x = unique(c(washc1_ppis$osEntrezA, washc1_ppis$osEntrezB))
mapID(x,"entrez","symbol") %>% unique()


data(bioid_anno) 

filt_anno <- bioid_anno %>% subset(Protein %in% wash_interactome)

filt_anno %>% 
	group_by(Annotation) %>% 
	summarize(n=length(unique(Protein)))

filt_anno %>% subset(Annotation == "Synapse")


##

washc_ppis <- musInteractome %>% 
	subset(osEntrezA %in% washc_entrez | osEntrezB %in% washc_entrez) %>% 
	subset(osEntrezA %in% bioid_entrez & osEntrezB %in% bioid_enterz)
entrez = unique(c(washc_ppis$osEntrezA, washc_ppis$osEntrezB))
prots <- getPPIs::getIDs(entrez,"entrez","uniprot","mouse")
prots["228767"] <- "A0A571BEE7"
sum(prots %in% wash_interactome)



length(wash_interactome)
