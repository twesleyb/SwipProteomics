renv

data(gene_map)
data(module_gsea)
data(sig_modules)
data(ne_surprise_partition)

x = module_gsea %>% filter(Module %in% sig_modules) %>% select(Pathway) %>% unique() %>% unlist(use.names=F)
y = x[grep("CORUM",x)]

library(geneLists)
data(corum)

data(corum)
corum["Exocyst complex"]
x = corum["Exocyst complex"]

mapID(x,"entrez","uniprot")

x = corum[["Exocyst complex"]]

mapID(x,"entrez","uniprot")
exocyst = mapID(x,"entrez","uniprot")


partition[exocyst]

fx = Abundance~1+Condition+(1|Mixture)+(1|Protein)

data(msstats_prot)

fm = lmerTest::lmer(fx,msstats_prot %>% filter(Protein %in% exocyst))

LT = getContrast(fm,"Mutant","Control")

lmerTestContrast(fm,LT)

r.squaredGLMM.merMod(fm)

getVariance(fm)
lmerTestContrast(fm,LT)
