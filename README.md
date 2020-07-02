# SwipProteomics
_Spatial proteomics analysis of __SWIP<sup>P1019R</sup>__ mutant mouse brain._

## SWIP Pro1019Arg
SWIP (gene: Washc4) is a component of the WASH complex. A P1019R point 
mutation was identified in humans with intellectual and motor impairments. 

<p align="center">
  <img src="./models/Swip.gif" height="250" />
</p>
<p align="center">The SWIP protein. The P1019R amino acid substitution is highlighted in red.<p align="center">

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
# Clone the repository.
git clone twesleyb/SwipProteomics
```
```R
# Download as an R package.
devtools::install_github('twesleyb/SwipProteomics')

```
```
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
