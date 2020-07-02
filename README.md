# SwipProteomics

This repository contains data and source code used in the analysis of the __WASH iBioID__
and the __SWIP<sup>P1019R</sup> TMT__ proteomics data associated with [Courtland _et al._, 2020](manuscript/nat-neuro/main/SWIP_paper_JC_v25.pdf).

### BACKGROUND
#### The WASH complex member SWIP: an intellectual disability gene
__Ropers _et al.___ [1] identified a non-synonymous mutation in the __WASHC4__ gene of
humans with intellectual impairment. This mutation causes a __P1019R__ amino acid substitution in SWIP
a __S__trumpellin and __W__ASH __I__nteracting-__P__rotein and component of the WASH complex.
Below is a model of SWIP's predicted 3D structure with the position of the P1019R amino acid substitution highlighted in red.

<p align="center">
  <img src="./models/Swip.gif" height="250" />
</p>

The WASH complex functions in endosome trafficking.

### QUESTION
#### How does SWIP<sup>P1019R</sup> cause intellectual disability?

### RESULTS
#### WASH iBioID

#### Genetic Disruption of WASHC4

#### Endo-lysosome Network Dysfunction

## METHODS
#### Spatial proteomics
Using 16-plex TMT-proteomics we quantified __5,894__ proteins from __7__ subcellular fractions prepared from control and SWIP<sup>P1019R</sup> MUT mice.
The brain of WT and SWIP mut mice was fractionated by differential
centrifugation at increasing velocity following the protocol from Geladaki _et al._, 2019[1].

#### Network construction
Construct a protein covariation matrix using `bicor` [(rdocumentation)](https://www.rdocumentation.org/packages/WGCNA/versions/1.69/topics/bicor),
a robust alternative to Pearson's coorelation.

Remove noise from the graph using `network enhancment` [(Wang et al., 2018)](refs/Wang_2018.pdf). 
My `neten` code [(github)](https://github.com/twesleyb/neten) was forked from an R port [(github.com/microbma/neten)](https://github.com/microbma/neten) 
of the original matlab code [(github.com/wangboyunze/Network_Enhancement)](https://github.com/wangboyunze/Network_Enhancement).

## Community Detection
Modules were identified in the enhanced protein covaration graph using the [Leiden](https://github.com/vtraag/leidenalg) algorithm
([Traag _et al._, 2019](refs/Traag_2019.pdf)). We utilized the [Surprise](refs/Traag_2015.pdf) quality metric to optimize clustering.

* [NetRep](https://cran.r-project.org/web/packages/NetRep/vignettes/NetRep.html) and original paper [Ritchie _et al._, 2016](refs/Ritchie_2016.pdf). 

## Differential protein and Module Abundance
* [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) used
    to model proteins and modules and test for differential abundance between WT
    and SWIP MUT mice.

## Misc
* Remove mitochondrial contaiminants [geneLists](https://github.com/twesleyb/geneLists).
* Interactions among WASH interactome as well as SWIP TMT proteome compiled using [getPPIs](https://github.com/twesleyb/getPPIs) 
  which currates PPIs from the [HitPredict](http://www.hitpredict.org/) database.
* [Colors](https://github.com/kevinwuhoo/randomcolor-py) installed with [conda](https://anaconda.org/conda-forge/randomcolor).


## Results 
We show that this point mutation disrupts expression of the
WASH complex in the brain and results in empaired endosomal trafficking in
neurons. Mice exhibit impaired cognition and motor deficites as well as cellular 
biomarkers of neurodegeneration. 

## Explore
To explore the data, download or install this repository as an R package.

```Bash
# Clone the repository.
git clone twesleyb/SwipProteomics
```

```R
# Download as an R package.
devtools::install_github('twesleyb/SwipProteomics')
```

```R
# Load the data.
devtools::load_all()
data(tmt_proteomics) # the final normalized TMT data
data(partition) # the partition of the protein covariation graph
```

## Reproducibility 
Effort was made to make the analysis as reproducible as possible. Reproduce the
R environment used to perform that analysis using [conda](https://docs.anaconda.com/anaconda/install/) 
and [renv](https://anaconda.org/conda-forge/r-renv). 
Install all additional Python dependencies with conda (e.g. [Leidenalg](https://anaconda.org/conda-forge/leidenalg)). 
Install all additional R dependencies with [renv](https://github.com/rstudio/renv). 
The conda environment can reproduced from the `SwipProteomics.yml` file.
The renv environment can reproduced from `renv.lock` file.

## References
[1] [Ropers _et al._, 2011](refs/Ropers_2011.pdf)
[2] [Geladaki _et al._, 2019](refs/Geladaki_2019.pdf)
