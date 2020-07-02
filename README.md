# SwipProteomics

This repository contains data and code used in the analysis of WASH iBioID
and __SWIP<sup>P1019R</sup>__ TMT proteomics data associated with [Courtland _et al._, 2020](URL).

## SWIP
The SWIP protein (gene: _Washc4_) is a component of the WASH complex.
[Ropers _et al._, 2011](refs/Ropers_2011.pdf). 
A P1019R point mutation was identified in humans with intellectual and motor impairments. 
SWIP<sup>P1019R</sup>

<p align="center">
  <img src="./models/Swip.gif" height="250" />
</p>
<p align="center">A model of SWIP's predicted tertiary structure. The position of the P1019R amino acid substitution is highlighted in red.<p align="center">

## Methods
Using 16-plex TMT-proteomics we quantified __5,894__ proteins from __7__ subcellular fractions prepared from control and SWIP<sup>P1019R<\sup> MUT mice.
Spatial proteomics protocol from [Geladaki _et al._, 2019](refs/Geladaki_2019.pdf).

* bicor function [WGCNA::bicor](URL)
* [Leiden algorithm](URL) [paper](refs/Traag_2019.pdf)  and [Surprise](refs/Traag_2015.pdf) quality metric.
* [NetRep](URL) and original paper [Ritchie et al., 2016](refs/Ritchie_2016.pdf). 
* [neten](https://github.com/twesleyb/neten) my fork of [R port](https://github.com/microbma/neten) of the original 
* [matlab code](https://github.com/wangboyunze/Network_Enhancement) the orignal [paper](refs/Wang_2018.pdf)
* [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
* [geneLists](https://github.com/twesleyb/geneLists)
* [getPPIs](https://github.com/twesleyb/getPPIs) and [HitPredict](http://www.hitpredict.org/) database.
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
