## DEP 
Protein-wise linear models combined with empirical Bayes statistics are used for
the differential enrichment analysis (or differential expression analysis). The
test_diff function introduced here uses r Biocpkg("limma") and automatically
generates the contrasts to be tested. 

## DEqMS
DEqMS builds on top of Limma, a widely-used R package for microarray data
analysis (Smyth G. et al 2004), and improves it with proteomics data specific
properties, accounting for variance dependence on the number of quantified
peptides or PSMs for statistical testing of differential protein expression.

Limma assumes a common prior variance for all proteinss, the function
spectraCounteBayes in DEqMS package estimate prior variance for proteins
quantified by different number of PSMs.

A documentation of all R functions available in DEqMS is detailed in the PDF
reference manual on the DEqMS Bioconductor page.

## MSstatsTMT
mixed linear models

## EdgeR
* exact test or
* glm - flexible
