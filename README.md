# SwipProteomics
This repository contains the proteomics data and source code for the analysis 
performed by [Courtland _et al._, 2020](manuscript/main/SWIP_Manuscript_v33_Combined_Neuron.pdf):

* Analysis of  __WASH iBioID__ Proteomics
* Analysis of __SWIP<sup>P1019R</sup> 16-plex TMT__ Spatial Proteomics

### The WASH complex member SWIP: an intellectual disability gene
__Ropers _et al.___ __[[1]](./refs/Ropers_2011.pdf)__ identified a non-synonymous mutation in the **_WASHC4_** gene of
seven children with intellectual impairment, learning and social deficits, and delayed motor development. 
This genetic mutation causes a single amino acid substitution in __SWIP__,
a <ins>S</ins>trumpellin and <ins>W</ins>ASH <ins>I</ins>nteracting <ins>P</ins>rotein 
and component of the pentameric __WASH complex__. 
Below is a model of the WASH complex as well as a model of SWIP's predicted 3D structure with the 
position of the __P1019R__ amino acid substitution highlighted in red.

| __WASH Complex__ | __SWIP<sup>P1019R</sup>__ |
|------------------|---------------------------|
|![](./figs/github/WASH_Complex.png)|![](./models/Swip.gif)|

### Genetic Disruption of WASHC4
How does the SWIP<sup>P1019R</sup> mutation cause intellectual and motor impairment in
humans? To study SWIP<sup>P1019R</sup>, a point mutation was made in the 
mouse __Washc4__ gene using CRISPR genome editing. 
In the brains of SWIP<sup>P1019R</sup> mutant (__MUT__) mice,
the abundance of WASH complex proteins__Strumpellin__ and __WASH1__ are drastically reduced
compared to control, wild-type (__WT__) mice.
These data suggest that the entire WASH complex may be
destabilized in the presence of SWIP<sup>P1019R</sup>.

<p align="center">
  <img src="./figs/github/DNA_Protein.png"/>
  </p>

### Endo-lysosome Network Dysfunction in SWIP<sup>P1019R</sup> Brain
To study the changes to the brain proteome caused by SWIP P1019R, we 
utilized 16-plex TMT-proteomics we quantified __5,894__ proteins from 
__7__ subcellular fractions prepared from WT control and SWIP<sup>P1019R<\sup> MUT mice.

We observed a significant reduction in the abundance of WASH complex proteins
__SWIP__, __WASH1__, __Fam21__, and __Strumpellin__. These data confirm that the 
entire WASH complex* is dysrupted in the presence of SWIP<sup>P1019R<\sup>.

<p align="center">
  <img src="./figs/github/TMT_Design.png"/>
  </p>

## METHODS
#### Spatial proteomics
Using 16-plex TMT-proteomics we quantified __5,894__ proteins from 
__7__ subcellular fractions prepared from WT control and SWIP<sup>P1019R<\sup> MUT mice.

#### Network construction
A symmetric, signed protein covaration matrix was built using the `bicor`
function, a robust alternative to Pearson's coorelation.

To remove noise from the graph, network enhancment __[[3]](./refs/Wang_2018.pdf)__, 
was performed using the `neten` R package [(github)](https://github.com/twesleyb/neten).

<p align="center">
  <img src="./figs/github/Network_Enhancement.png" height="250" />
  </p>

#### Community Detection
Modules were identified in the enhanced protein covaration graph using the `Leiden algorithm` and
the Python leidenalg library [(github)](https://github.com/vtraag/leidenalg)__[4]__.
Leidenalg utilizes a quality statistic to identify optimal partitions of the graph.  

We utilized the [Surprise](refs/Traag_2015.pdf) metric to optimize clustering. Surpsie quantifies 
the probability of observing a module with at least is the negative logarithm 
probablility of drawing m edges without replacement.

<p align="center">
  <img src="./figs/github/Leiden.png" height="250" />
  </p>

#### Module Preservation
To enforce module quality, permutation testing was performed to identify modules
whose underlying toplogy was not different from 10,000 random permutations of
the graph using `NetRep` [cran](https://cran.r-project.org/web/packages/NetRep/vignettes/NetRep.html) __[5]__.

<p align="center">
  <img src="./figs/github/Permutation_Histogram.png" height="250" />
  </p>

#### Differential Protein and Module Abundance
The `edgeR` package [bioconductor](https://bioconductor.org/packages/release/bioc/html/edgeR.html) was used
to model proteins and modules and test for differential abundance between WT and SWIP MUT mice __[6]__.

## iBioID Analysis
Mitochondrial contaiminants were removed from the iBioID proteome. Mitochondrial
proteins from MitoCarta2 were downloaded from `geneLists` [(github)](https://github.com/twesleyb/geneLists).

## Protein-Protein Interactions (PPIs)
Interactions among proteins were compiled using `getPPIs` [(github)](https://github.com/twesleyb/getPPIs),
a compilation of PPIs from the [HitPredict](http://www.hitpredict.org/) database.

## Module Colors
The python `randomcolor` library [(github)](https://github.com/kevinwuhoo/randomcolor-py)
was used to generate a large number of aestethically appealing colors. Download
`randomcolor` from [conda](https://anaconda.org/conda-forge/randomcolor).

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
[4] [Traag _et al._, 2019](refs/Traag_2019.pdf)  
[5] [Ritchie _et al._, 2016](refs/Ritchie_2016.pdf)  
[6] [McCarth _et al._, 2012](refs/McCarth_2012.pdf)  
