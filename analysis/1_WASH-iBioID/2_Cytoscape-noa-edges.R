#!/usr/bin/env Rscript

## ---- edges and noa files for Cytoscape

renv

data(bioid_anno)
data(bioid_results)

library(getPPIs)

data(musInteractome)


anno <- dcast(bioid_anno, Protein ~ Annotation, value.var = "PMID")
noa <- bioid_results %>% 
	filter(Protein %in% wash_interactome) %>% 
	left_join(anno, by="Protein") %>% 
	mutate(Symbol = toupper(Symbol))

fwrite(noa, "noa.csv")


os_keep <- c(9606, 10116, 10090)

data(bioid_gene_map)

entrez <- mapID(wash_interactome,"uniprot","entrez")

# collect interactions between swip and wash_interactome proteins
wash_ppis <- musInteractome %>% 
	filter(Interactor_A_Taxonomy %in% os_keep) %>% 
	filter(Interactor_B_Taxonomy %in% os_keep) %>%
	subset(osEntrezA %in% entrez & osEntrezB %in% entrez)
edge_df <- wash_ppis %>% select(osEntrezA, osEntrezB, Publications)

protA = mapID(edge_df$osEntrezA,'entrez','uniprot')
protB = mapID(edge_df$osEntrezB,'entrez','uniprot')

edge_df <- tibble::add_column(edge_df, protA, .before='osEntrezA')
edge_df <- tibble::add_column(edge_df, protB, .after='protA')

df1 = edge_df %>% select(protA, protB, Publications)
df2 = data.table(protA=mapID("Washc1"),
	   protB=wash_interactome,
	   Publications='Courtland et al. 2020')

edge_df = rbind(df1,df2)

fwrite(edge_df, "edges.csv")
