#!/usr/bin/env Rscript

# author: twab
# description: create edge and noa files for Cytoscape
# title: SwipProteomics

## ---- prepare the R env

root <- "~/projects/SwipProteomics"
renv::load(root, quiet=TRUE)
devtools::load_all(root, quiet=TRUE)


## ---- load the data

#library(SwipProteomics)
data(bioid_anno)
data(bioid_results)
data(bioid_gene_map)

library(getPPIs) # twesleyb/getPPIs

data(musInteractome)


## ---- create node annotation (noa) data frame

anno <- dcast(bioid_anno, Protein ~ Annotation, value.var = "PMID")

noa <- bioid_results %>% 
	filter(Protein %in% wash_interactome) %>% 
	left_join(anno, by="Protein") %>% 
	mutate(Symbol = toupper(Symbol))

myfile <- file.path(root,"rdata","noa.csv")
fwrite(noa, myfile)
message("saved: ", myfile)


## ---- collect ppis

# keep ppis from human, mouse, and rat
os_keep <- c(9606, 10116, 10090)

# map wash interactome uniprot IDs to entrez
entrez <- mapID(wash_interactome,"uniprot","entrez")

# collect interactions between swip and wash_interactome proteins
wash_ppis <- musInteractome %>% 
	filter(Interactor_A_Taxonomy %in% os_keep) %>% 
	filter(Interactor_B_Taxonomy %in% os_keep) %>%
	subset(osEntrezA %in% entrez & osEntrezB %in% entrez)
edge_df <- wash_ppis %>% select(osEntrezA, osEntrezB, Publications)

# entrez IDs to UniProt
protA = mapID(edge_df$osEntrezA,'entrez','uniprot')
protB = mapID(edge_df$osEntrezB,'entrez','uniprot')

# add to table
edge_df <- tibble::add_column(edge_df, protA, .before='osEntrezA')
edge_df <- tibble::add_column(edge_df, protB, .after='protA')

# combine as edge data.frame
df1 = edge_df %>% select(protA, protB, Publications)
df2 = data.table(protA=mapID("Washc1"),
	   protB=wash_interactome,
	   Publications='Courtland et al. 2020')
edge_df = rbind(df1,df2)

# save
myfile <- file.path(root,"rdata","edges.csv")
fwrite(edge_df, myfile)
message("saved: ", myfile)
