# SwipProteomics
This repository contains the proteomics data and source code for the analysis 
performed by [Courtland _et al._, 2020](manuscript/main/SWIP_Manuscript_v33_Combined_Neuron.pdf):

* Analysis of  __WASH iBioID__ Proteomics
* Analysis of __SWIP<sup>P1019R</sup> 16-plex TMT__ Spatial Proteomics

## The WASH complex member SWIP: an intellectual disability gene
__Ropers _et al.___ [[1]](./refs/Ropers_2011.pdf) identified a non-synonymous mutation in the **_WASHC4_** gene of
seven children with intellectual impairment, learning and social deficits, and delayed motor development. 
This genetic mutation causes a single amino acid substitution in __SWIP__,
a <ins>S</ins>trumpellin and <ins>W</ins>ASH <ins>I</ins>nteracting <ins>P</ins>rotein and component of the pentameric __WASH complex__. 

Below is a model of the WASH complex as well as a model of SWIP's predicted 
3D structure with the position of the __P1019R__ amino acid substitution 
highlighted in red.


## METHODS

#### Genetic Dysruption of SWIP
To study SWIP<sup>P1019R</sup>, a point mutation was made in the 
mouse __Washc4__ gene using CRISPR genome editing. 

#### Spatial proteomics
Geladaki _et al_., [[2]](refs/Geladaki_2019.pdf).
Using 16-plex TMT-proteomics we quantified __5,894__ proteins from 
__7__ subcellular fractions prepared from WT control and 
SWIP<sup>P1019R<\sup> MUT mice.

#### Network construction
A symmetric, signed protein covaration matrix was built using the `bicor`
function, a robust alternative to Pearson's coorelation.

To remove noise from the graph, network enhancment [[3]](./refs/Wang_2018.pdf), 
was performed using the `neten` R package 
[(github)](https://github.com/twesleyb/neten).

#### Community Detection
Modules were identified in the enhanced protein covaration graph using 
the Leiden algorithm as implemented by the `leidenalg` package 
[(github)](https://github.com/vtraag/leidenalg) [[4]](refs/Traag_2019.pdf).

Leidenalg utilizes a quality statistic to identify optimal partitions 
of the graph. We utilized the `Surprise` metric to optimize clustering 
[[5]](refs/Traag_2015.pdf).  

Surprise quantifies the probability of observing a module with at least 
is the negative logarithm probablility of drawing m edges without replacement.

#### Module Preservation
To enforce module quality, permutation testing was performed to identify 
modules whose underlying toplogy was not different from 10,000 random 
permutations of the graph using `NetRep` 
[(Cran)](https://cran.r-project.org/web/packages/NetRep/vignettes/NetRep.html)
[[6]](refs/Ritchie_2016.pdf).  

#### Differential Protein and Module Abundance
The `edgeR` package 
[(BioConductor)](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
was used to model proteins and modules and test for differential abundance 
between WT and SWIP MUT mice [[7]](refs/McCarthy_2012.pdf).  

#### Protein-Protein Interactions (PPIs)
Interactions among proteins were compiled using `getPPIs` 
[(github)](https://github.com/twesleyb/getPPIs),
a compilation of PPIs from the [HitPredict](http://www.hitpredict.org/) database.

#### Module Colors
The python `randomcolor` library [(github)](https://github.com/kevinwuhoo/randomcolor-py)
was used to generate a large number of aestethically appealing colors. Download
`randomcolor` from [conda](https://anaconda.org/conda-forge/randomcolor).

#### Gene Set Enrichment analysis
Mitochondrial contaiminants were removed from the iBioID proteome. 
Mitochondrial proteins from [MitoCarta2](url) were downloaded from 
`geneLists` [(github)](https://github.com/twesleyb/geneLists).

## Results 
We show that this point mutation disrupts expression of the
WASH complex in the brain and results in empaired endosomal trafficking in
neurons. Mice exhibit impaired cognition and motor deficites as well as cellular 
biomarkers of neurodegeneration. 

## Explore
To explore the data, download or install this repository as an R package.

Using git:
```Bash
# Clone the repository.
git clone twesleyb/SwipProteomics
```

Using R and `devtools`:
```R
# Download as an R package with devtools.
devtools::install_github('twesleyb/SwipProteomics')
```

## Datasets
Key datasets can be accessed within R using the `data` function.
```R

library(SwipProteomics)

data(tmt_proteomics) # the final normalized TMT data

data(partition) # the partition of the protein covariation graph
```

## Reproducibility 
Effort was made to make the analysis as reproducible as possible.   

Reproduce the R environment used to perform that analysis using [conda](https://docs.anaconda.com/anaconda/install/) 
and [renv](https://anaconda.org/conda-forge/r-renv).   

All additional Python dependencies were installed with conda (e.g. [Leidenalg](https://anaconda.org/conda-forge/leidenalg)). 
All additional R dependencies  were installed within the [renv](https://github.com/rstudio/renv) environment.

## References
[1] [Ropers _et al._, 2011](refs/Ropers_2011.pdf)  
[2] [Geladaki _et al._, 2019](refs/Geladaki_2019.pdf)  
[3] [Wang _et al._, 2018](refs/Wang_2018.pdf)  
[4] [Traag _et al._, 2015](refs/Traag_2015.pdf)  
[5] [Traag _et al._, 2019](refs/Traag_2019.pdf)  
[6] [Ritchie _et al._, 2016](refs/Ritchie_2016.pdf)  
[7] [McCarthy _et al._, 2012](refs/McCarthy_2012.pdf)  
