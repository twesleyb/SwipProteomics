# SwipProteomics
This repository contains the proteomics data and source code for the analysis 
performed by [Courtland _et al._,2020](https://www.biorxiv.org/content/10.1101/2020.08.06.239517v1).

## The WASH complex member SWIP: an intellectual disability gene
__Ropers _et al.___ [[1]](./refs/Ropers_2011.pdf) identified a non-synonymous mutation in the **_WASHC4_** gene of
seven children with intellectual impairment, learning and social deficits, and delayed motor development. 
This genetic mutation causes a single amino acid substitution in __SWIP__,
a <ins>S</ins>trumpellin and <ins>W</ins>ASH <ins>I</ins>nteracting <ins>P</ins>rotein and component of the pentameric __WASH complex__. 

Below is a model of the WASH complex as well as a model of SWIP's predicted 
3D structure with the position of the __P1019R__ amino acid substitution 
highlighted in red.

| __WASH Complex__ | __SWIP<sup>P1019R</sup>__ |
|------------------|---------------------------|
|![](./figs/github/WASH_Complex.png)|![](./figs/github/Swip.gif)|

The P1019R missense mutation results in the exchange of proline (Pro, P), an internal neutral
amino acid with an arginine (Arg, R), a charged amino acid. This mutation is found in the WASHC4 
domain that is thought be be critical for its interaction with Strumpellin, and
is predicted to be structurally disruptive. Accordingly, we found that MUT SWIP disrupts its 
ntraction Strumpellin and WASH1.

<p align="center">
  <img src="./figs/github/DNA_Protein.png"/>
  </p>

## Spatial proteomics
We performed subcellular fractionation and TMT proteomics, or spatial
proteomics, by adapting the protocol established by Geladaki _et al_.[[2]](refs/Geladaki_2019.pdf).
We prepared 7 subceullar fractions from brain, in all quantifying __5,894__ proteins.
By analyzing the correlated abundance profiles of proteins across all
seven fractions we reoncstructed a brain spatial proteomics network in which
modules of highly correlated proteins represent groups of interacting,
colocalizing, functionally related proteins.

<p align="center">
  <img src="./figs/github/TMT_Design.png"/>
  </p>

## Network construction
A symmetric, signed protein covaration matrix was built using the `WGCNA::bicor`
function, a robust alternative to Pearson's coorelation.

To remove noise from the graph, network enhancement [[3]](./refs/Wang_2018.pdf), 
was performed using an R fork [(github)](https://github.com/microbma/neten) of the original Matlab code. 

<p align="center">
  <img src="./figs/github/Network_Enhancement.png" height="250" />
  </p>

## Community Detection
Modules were identified in the enhanced protein covaration graph using 
the Leiden algorithm as implemented by the `leidenalg` package 
[(github)](https://github.com/vtraag/leidenalg) [[4]](refs/Traag_2019.pdf).

Leidenalg utilizes a quality statistic to identify optimal partitions 
of the graph. We utilized the `Surprise` metric to optimize clustering 
[[5]](refs/Traag_2015.pdf).  

Surprise quantifies the probability of observing a module with at least 
is the negative logarithm probablility of drawing m edges without replacement.

<p align="center">
  <img src="./figs/github/Leiden.png" height="250" />
  </p>

## Module Preservation
To enforce module quality, permutation testing was performed to identify 
and remove modules whose underlying toplogy was not different from 10,000 random 
permutations of the graph using `NetRep` 
[(Cran)](https://cran.r-project.org/web/packages/NetRep/vignettes/NetRep.html)
[[6]](refs/Ritchie_2016.pdf).  

<p align="center">
  <img src="./figs/github/Permutation_Histogram.png" height="250" />
  </p>

## Differential Protein and Module Abundance
The `edgeR` package 
[(BioConductor)](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
used to model proteins and modules and test for differential abundance 
between WT and SWIP MUT mice.

## Additional Methods

#### Protein-Protein Interactions (PPIs)
Interactions among proteins were compiled using `getPPIs` 
[(github)](https://github.com/twesleyb/getPPIs),
a compilation of PPIs from the [HitPredict](http://www.hitpredict.org/) database.

#### Module Colors
The python `randomcolor` library [(github)](https://github.com/kevinwuhoo/randomcolor-py)
was used to generate a large number of aestethically appealing colors. Download
`randomcolor` from [conda](https://anaconda.org/conda-forge/randomcolor).

#### Gene Set Enrichment analysis
Gene set enrichment analysis was performed using lists of genes from 
`geneLists` [(github)](https://github.com/twesleyb/geneLists) and the
hypergeometric function in R.

## Reproducibility 
Effort was made to make the analysis as reproducible as possible.   
Please not however, that randomness in the permutation testing proceedure means
that the exact result may not be reproduced if the code were run again.

This repository can be downloaded as an R package using devtools.
NOTE: currently the package contains all of the project's data in `data` and is
quite large. The download is 121.8 Mb and may take several minutes.
```R
devtools::install.packages()
```
The raw data and key datasets can be accessed 
within R using the `data` function.
```R

library(SwipProteomics)

data(tmt_proteomics) # the final normalized TMT data

data(adjm) # the bicor protein covariation matrix

data(ne_adjm) # the enhanced network

data(partition) # the Leidenalg partition of the protein covariation graph
```

It's recommended to try and reproduce the research environment using 
[conda](https://docs.anaconda.com/anaconda/install/) 
and [renv](https://anaconda.org/conda-forge/r-renv).   

Create a conda environment and then install renv with conda:
```
(SwipProteomics) $ conda install -c conda-forge r-renv
```

All additional Python dependences should be installed with conda (e.g. [Leidenalg](https://anaconda.org/conda-forge/leidenalg)). 
All additional R dependencies should be installed in the R environment managed by renv [renv](https://github.com/rstudio/renv).

## References
[1] [Ropers _et al._, 2011](refs/Ropers_2011.pdf)  
[2] [Geladaki _et al._, 2019](refs/Geladaki_2019.pdf)  
[3] [Wang _et al._, 2018](refs/Wang_2018.pdf)  
[4] [Traag _et al._, 2015](refs/Traag_2015.pdf)  
[5] [Traag _et al._, 2019](refs/Traag_2019.pdf)  
[6] [Ritchie _et al._, 2016](refs/Ritchie_2016.pdf)  

## Liscense
This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-shield]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
