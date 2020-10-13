Hi everyone, 

Ting, thank you for your help and proposed collaboration. 
I've cc'd _Scott Soderling_, my PI, as well as _Jamie Courtland_, 
our [paper's](https://www.biorxiv.org/content/10.1101/2020.08.06.239517v1)
main co-author.  

In short, we have two goals, could you please help me: 
(1) extend MSstatsTMT to handle our repeated measures design, and 
(2) extend MSstatsTMT to perform inference at the level of protein groups--
  as you suggested in your recent paper. 

Below I describe our study's basic biological motivation as well as our 
proteomics experimental design. I look forward to working with you and
your colleagues. 

Kind regards,
Tyler


I shall describe our study and experimental design in brief, so that we can 
establish some common understanding.

We are studying a mouse model with a point mutation in the __WASHC4__ gene,
which encodes the __Swip__ protein, a required component of the WASH complex. 
This mutation causes a neurodevelopmental disorder in humans characterized by
intellectual disability and deficits in basic motor skills.

The functional effect of this mutation is loss-of-function of Swip and the WASH 
protein complex. We propose that the WASH complex functions in the brain, 
as it does in other cells/organisms, in endosomal trafficking. 

In close collaboration with our proteomics shared resource, we have performed 3
16-plex TMT experiments. Our experimental design is aimed as understanding the
affects of the WASCH4 mutation (i.e. loss of Swip) on the neuronal proteome.

We anticipate that MUT mice exhibit changes in endosomal compartments which 
may then manifest as changes in subcellular trafficking and the distribution of 
proteins across the cell. Thus, we prepared 7 subcellular fractions from the
brains of 6 mice, 3 Control (aka WT) and 3 MUT mice, and profiled these
biological fractions by performing 3x 16-plex TMT proteomics experiments.

| Feature     | Number |  Note |
| ----------- | ------ | ------|
| Run         | 3      | Sometimes I call this an `Experiment`. |
| Mixture     | 3      | 16-plex |
| TechRepMixture | 1   | There were no technical replicates of mixture. |
| Fractions   | 12     | Each mixture was subfractionated to increase depth of coverage. |
| Channels   | 16      | 14x `Condition`.`BioFraction` + 2x `SPQC` or `Norm` |
| Condition | 2  | aka `Group`, `Treatment`, or `Genotype` |
| BioFraction | 7  | 7 subcellular fractions prepared by centrifugation |

I recognized that currently MSstatsTMT does not handle a repeated measures
design, but that MSstats did. 

To do __(1)__ above, i.e. to assess protein-level changes accounting for the
random effect of subcellular fraction (`Biofraction` -- to distinguish from MS Fraction),
I fit the following mixed effects model to each protein in the normalized data:
`Abundance ~ 0 + (1|BioFraction) + Genotype`

In which `(1|BioFraction)` indicates the mixed-effect of subcellular fraction on 
overall protein abundance.

I believe this is analogous to the `lmer` formula fit by MSstats: 
`ABUNDANCE ~ GROUP + (1|SUBJECT) # for case !TechRep `

The attached script `fit1.R` performs the analysis for a single protein.
It fits a protein with the proposed formula and the stastical comparison 
between conditioned means of the `Control` and `Mutant` groups. These steps
were just extracted from internal MSstatsTMT functions.

Would you mind looking over this script to things to see it makes sense?
Could you please provide some addition description about the steps done to
calculate the p-value for comparisons between groups?

If were to actually try to run the script you will need the data. It's in the
data directory of the explore branch of our git repository:
https://github.com/twesleyb/SwipProteomics/tree/explore/data

We are also interested in __(2)___
extending this mixed effects model approach to perform inference at the level of groups of proteins.
I may also call these `modules` or `clusters`. 

To do so, I fit the following model to the protein-level data from each module:
`Abundance ~ 0 + (1|BioFraction) + (1|Protein) + Genotype`

Here I did not include the effect of Mixture (aka Batch or Experiment) because
it did not seem to contribute much to the variance, but I am not sure about
this.
The code to perform such a module -level analysis is the same as fit1.R, with
some minor substitutions.Note that the contrast matrix is something
like: setNames(c(-1, 1), nm = c("GenotypeControl", "GenotypeMutant"))

I have done this analysis and the results seem to make sense but I am trying 
to evaluate them.
