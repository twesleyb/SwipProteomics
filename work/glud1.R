renv
data(partition)
data(msstats_prot)
data(gene_map)
mapID("Glud1","symbol","uniprot")
glud1 <- mapID("Glud1","symbol","uniprot")
glud1
fx = Abundance ~ 1 + Condition + (1|Mixture)
lmerTest::lmer(fx, msstats_prot %>% filter(Protein == glud1)
lmerTest::lmer(fx, msstats_prot %>% filter(Protein == glud1))
fm = lmerTest::lmer(fx, msstats_prot %>% filter(Protein == glud1))
fm
LT = getContrast(fm, "Mutant","Control")
lmerTestContrast(fm,LT)
myfile = file.path(root,"rdata","adjm.rda"))
myfile = file.path(root,"rdata","adjm.rda")
root = getrd()
myfile = file.path(root,"rdata","adjm.rda")
myfile
load(myfile)
dim(adjm)
partition
M1 = partition[partition==1]
M1
partition
data(ne_surprise_partition)
M1 = partition[partition==1]
M1
M1 = names(partition)[partition==1]
M1
subadjm = adjm[M1,M1]
importance = apply(subadjm,2,sum)
importance
max(importance)
hist(importance)
myfile = file.path(root,"rdata","ne_adjm.rda")
load(myfile)
sub_netw = ne_adjm[M1,M1]
imp2 = apply(sub_netw,2,sum)
imp2
max(imp2)
hist(imp2)
importance[glud1]
max(importance)
adjm[glud1,]
subadjm[glud1,]
x = subadjm[glud1,]
mean(X)
mean(x)
hist(x)
dir
savehistory("glud1.R")
