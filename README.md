# SwipProteomics
This repository contains the proteomics data and source code for the analysis 
performed by [Courtland _et al._, 2020](https://www.biorxiv.org/content/10.1101/2020.08.06.239517v1).

<<<<<<< HEAD
## The WASH complex member SWIP: an intellectual disability gene
__Ropers _et al.___ [[1]](./refs/Ropers_2011.pdf) identified a non-synonymous mutation in the **_WASHC4_** gene of
seven children with intellectual impairment, learning and social deficits, and delayed motor development. 
This mutation causes a single amino acid substitution in __SWIP__,
a <ins>S</ins>trumpellin and <ins>W</ins>ASH <ins>I</ins>nteracting <ins>P</ins>rotein 
and component of the pentameric __WASH complex__. 

Below is a schematic of the WASH complex as well as a model of SWIP's predicted 
3D structure with the position of the __P1019R__ amino acid substitution 
highlighted in red.

| __WASH Complex__ | __SWIP<sup>P1019R</sup>__ |
|------------------|---------------------------|
|![](./figs/github/WASH_Complex.png)|![](./figs/github/Swip.gif)|

The P1019R missense mutation results in the exchange of proline (P), an internal neutral
amino acid, with an arginine (R), a charged amino acid. This mutation is found in the WASHC4 
domain that is thought be be critical for its interaction with Strumpellin 
[(Jia et al., 2010; pmid 20498093)](https://pubmed.ncbi.nlm.nih.gov/20498093/), and
is predicted to be structurally disruptive [(Missense3d)](http://www.sbg.bio.ic.ac.uk/~missense3d/). 
Accordingly, we found that MUT SWIP mice have decreased levels of WASH1 and 
Strumpellin, and the SWIP mutation disrupts its interaction both WASH1 and Strumpellin 
[add some data].

<p align="center">
  <img src="./figs/github/DNA_Protein.png"/>
  </p>

## Spatial proteomics
We performed subcellular fractionation and TMT proteomics, or spatial
proteomics, by adapting the protocol established by Geladaki _et al_. [[2]](refs/Geladaki_2019.pdf)
for brain tissue. We prepared 7 subceullar fractions, in all quantifying __5,894__ proteins.
By analyzing the correlated abundance profiles of proteins across all
seven fractions we constructed a spatial proteomics network in which
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
The `edgeR` package can be installed from [ioConductor](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
and was used test for differential protein and module abundance between WT and SWIP MUT mice.

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
Modules were analyzed for gene set enrichment using lists of genes from 
[twesleyb/geneLists](https://github.com/twesleyb/geneLists) and the
hypergeometric function in R.

## Reproducibility 
Effort was made to make the analysis as reproducible as possible.   
Please note, that that randomness in the permutation testing procedure used to
enforce module quality means that the exact result may not be reproduced 
if the code were run again.

This repository can be downloaded as an R package containing the raw and
normalized data as well as the statistical results. The data can be found
in the `data/` directory and in total is quite large (121.8 Mb). 
Downloading the package may take several minutes.

```R
# Install the SwipProteomics package.
devtools::install_github("twesleyb/SwipProteomics")

# Access the project's data with the data() command. 

library(SwipProteomics)

data(tmt_proteomics) # the final normalized TMT data

data(adjm) # the bicor protein covariation matrix

data(ne_adjm) # the enhanced network

data(partition) # the Leidenalg partition of the protein covariation graph
```
Alternatively, download the data directly by visiting the 
[data/ directory](./data/).

### Dependencies
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
This work is licensed under a [Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-shield]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
=======
Analysis of SWIP<sup>P1019R</sup> (Washc4) mutant mouse brain spatial proteomics

## Swip Pro1019Arg Mutation

<p align="center">
  <img src="./models/Swip.gif" height="250" />
</p>
<p align="center">The SWIP protein with location of P1019R highlighted in red.<p align="center">

## Background
The SWIP P1019R point mutation was identified in humans with intellecutual and
motor impairments. 

## Methods
Using TMT-proteomics we quantify 5,894 proteins across seven mouse brain subcellular
fractions. 

## Results 
We show that this point mutation disrupts expression of the
WASH complex in the brain and results in empaired endosomal trafficking in
neurons. Mice exhibit impaired cognition and motor deficites as well as cellular 
biomarkers of neurodegeneration. 

## Conclusions


## Explore
To explore the data, download or install this repository as an R package.
```
git clone twesleyb/SwipProteomics
```
```R
# Download as an R package.
devtools::install_github('twesleyb/SwipProteomics')

```

## Reproducibility 
Effort was made to make the analysis as reproducible as possible. Reproduce the
R environment used to perform that analysis using [conda]() and [renv]().
>>>>>>> 94dece43b955711d647aa4147600c8d6dc417597
