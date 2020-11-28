renv
data(ne_surprise_partition)
partition==39
partition[partition==39]
partition[partition==39] %>% names()
prots = partition[partition==39] %>% names()
prots
fm
args_list = list()
args_list[["formula"]] = Abundance ~ 0 + Genotype:BioFraction + (1|Protein) + (1|Mixture)
lmerFit(args_list)
args_list[["data"]] = msstats_prot %>% filter(Protein %in% prots))
args_list[["data"]] = msstats_prot %>% filter(Protein %in% prots)
data(msstats_prot)
args_list[["data"]] = msstats_prot %>% filter(Protein %in% prots)
args_list
names(args_list)
lmerFit(args_list)
fm = lmerFit(args_list)
LT = getContrast(fm,"Mutant","Control")
lmerTestContrast(fm,LT)
prots
modules
partition[prots]
partition
table(partition)
data(module_gsea)
module_gsea %>% filter(Module == "M39")
data(module_results)
module_results$Module
data(corum)
library(geneLists)
data(corum)
corum["Fibrinogen complex"]
fibrinogen <- corum["Fibrinogen complex"]
fibrinogen <- corum[["Fibrinogen complex"]]
fibrinogen
mapID(fibrinogen,"entrez","uniprot")
data(gene_map)
mapID(fibrinogen,"entrez","uniprot")
fibrinogen = mapID(fibrinogen,"entrez","uniprot")
save(fibrinogen,file="fibrinogen.rda",version=2)
fm = lmerTest::lmer(fx, msstats_prot %>% subset(Protein %in% firbrinogen))
fx
fx=Abundance ~ 1 + Condition + (1|Mixture) + (1|Protein)
fm = lmerTest::lmer(fx, msstats_prot %>% subset(Protein %in% firbrinogen))
fm = lmerTest::lmer(fx, msstats_prot %>% subset(Protein %in% fibrinogen))
fm
LT = getContrast(fm,"Mutant","Control")
LT
sum(LT)
lmerTestContrast(fm,LT)
savehistory("hist.R")
