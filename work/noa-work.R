#!/usr/bin/env Rscript

root <- "~/projects/SwipProteomics"
renv::load(root, quiet=TRUE)

data(gene_map)
data(gene_lists)
data(sig_modules)
data(swip_partition)
data(swip_module_gsea)


sig_paths <- module_gsea %>% 
	filter(Module %in% sig_modules) %>% 
	select(Pathway) %>% 
	unique() %>% 
	unlist(use.names=FALSE)

paths <- names(gene_lists)[names(gene_lists) %in% sig_paths]

#lysosome <- gene_lists[paths[grep("lysosome",tolower(paths))]]
#endosome <- gene_lists[paths[grep("endosome",tolower(paths))]]
#wash <- gene_lists[paths[grep("wash",tolower(paths))]]
#mito <- gene_lists[paths[grep("mito",tolower(paths))]]
#er <- gene_lists[paths[grep("endoplasmic reticulum",tolower(paths))]]

# select attributes
# M17
df <- data.table(Protein = names(partition)) %>% 
	mutate(Entrez = mapID(Protein,"uniprot","entrez"))
idx <- df$Entrez %in% gene_lists[["LopitDC: ER"]] 
idy <- df$Entrez %in% gene_lists[["Uniprot: Endoplasmic reticulum membrane"]]
df[["ER"]] <- idx | idy
df[["EXOCYST"]] <- df$Entrez %in% gene_lists[["CORUM: Exocyst complex"]] 
# M36

idx <- df$Entrez %in% gene_lists[["Uniprot: Lysosome"]] 
idy <- df$Entrez %in% gene_lists[["LopitDC: LYSOSOME"]]
df[["LYSOSOME"]] <- idx | idy

"Man2b2" %in% mapID(df$Protein[df$LYSOSOME],"uniprot","symbol")

lyso <- c(gene_lists[["Uniprot: Lysosome"]],
	      gene_lists[["LopitDC: LYSOSOME"]])
lysosome <- mapID(lyso,"entrez","uniprot")
df <- data.table(Protein = names(partition))
df$LYSOSOME  <- df$Protein %in% lysosome
fwrite(df,"noa.csv")

module_gsea %>% filter(Module == "M22") %>% select(Pathway)

module_gsea %>% filter(Module == "M22")

df <- data.table(Protein = names(partition))
df$BioID <- df$Protein %in% mapID(gene_lists[["WASH-iBioID"]],"entrez","uniprot")
df$ePSD <- df$Protein %in% mapID(gene_lists[["Uezu et al., 2016: ePSD"]], "entrez","uniprot")
fwrite(df,"noa.csv")

df <- data.table(Protein = names(partition))
df$RETROMER <- df$Protein %in% mapID(gene_lists[["CORUM: Retromer complex (SNX1, SNX2, VPS35, VPS29, VPS26B)"]],"entrez","uniprot")
df$ERM <- df$Protein %in% mapID(gene_lists[["Uniprot: Endoplasmic reticulum membrane"]],"entrez","uniprot")
fwrite(df,"noa.csv")

df$LYSOSOME  <- df$Protein %in% lysosome
fwrite(df,"noa.csv")

df[["EE"]] <- df$Entrez %in% gene_lists[["Uniprot: Early endosome"]]
df[["RETRIEVER"]] <- df$Entrez %in% gene_lists[["McNally et al., 2017: Retriever Complex"]]
df[["CCC"]] <- df$Entrez %in% gene_lists[["CORUM: CCC complex"]]
df[["CORVET"]] <- df$Entrez %in% gene_lists[["CORUM: CORVET complex"]]
df[["HOPS"]] <- df$Entrez %in% gene_lists[["CORUM: HOPS complex"]]
df[["WASH"]] <- df$Entrez %in% gene_lists[["CORUM: WASH complex"]]

names(gene_lists)[grep("Retromer",names(gene_lists))]

fwrite(df,"noa.csv")

mapID(df$Protein[df$CCC],"uniprot","symbol")

module_gsea %>% filter(Pathway == "CORUM: HOPS complex")
