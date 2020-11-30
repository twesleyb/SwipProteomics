renv
data(modularity_partition)
partition
table
class(partition)
length(partition)
x = partition[[1]]
x
table(x)
partition = x
partition = partition + 1
save(partition,file="modularity_partition.rda",version=2_
save(partition,file="modularity_partition.rda",version=2)
data(msstats_prot)
fx = Abundance ~ 0 + Condition + (1|Mixture) + (1|Protein)
prots = partition[partition==1]
prots = namespartition[partition==1])
prots = names(partition[partition==1])
prots
length(prots)
fx
fm = lmerTest::lmer(fx, data = msstats_prot %>% subset(Protein %in% prots))
fm
summary(fm)
fx
getContrast(fm,"Mutant","Control")
LT = getContrast(fm,"Mutant","Control")
lmerTestContrast(fm,LT)
?history
savehistory("modularity.R")
