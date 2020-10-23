#!/usr/bin/env Rscript

root = "~/projects/SwipProteomics"

renv::load(root)
devtools::load_all(root)

data(fx0)
data(fx1)
data(cm0)
data(cm1)
data(swip)
data(msstats_prot)
data(msstats_contrasts)

library(dplyr)

## analysis with MSstats

levels(msstats_prot$Condition)

# Compare condition Control.F10 and Mutant.F10
comparison<-matrix(c(-1,0,0,0,0,0,0,
                     1,0,0,0,0,0,0),nrow=1)

# Set the names of each row
row.names(comparison)<-"Control.F10-Mutant.F10"
colnames(comparison)<- levels(msstats_prot$Condition)

subdat <- msstats_prot %>% filter(Protein == swip)

test_contrast <- MSstatsTMT::groupComparisonTMT(data = subdat,
                                    contrast.matrix = comparison,
                                    moderated = FALSE)

test_contrast %>% knitr::kable()

lmerTestProtein(swip, fx0, msstats_prot, cm0) %>% knitr::kable()

# Compare condition Control and Mutant
comparison<-matrix(c(-1/7,-1/7,-1/7,-1/7,-1/7,-1/7,-1/7,
                     1/7,1/7,1/7,1/7,1/7,1/7,1/7),nrow=1)

# Set the names of each row
row.names(comparison)<-"Mutant-Control"
colnames(comparison)<- levels(msstats_prot$Condition)

test_contrast <- MSstatsTMT::groupComparisonTMT(data = subdat,
                                    contrast.matrix = comparison,
                                    moderated = FALSE)

test_contrast %>% knitr::kable()

lmerTestProtein(swip, fx1, msstats_prot, cm1) %>% knitr::kable()

# the results differ in their stder
