---
mainfont: Arial.ttf
sansfont: Arial.ttf
monofont: Arial.ttf 
mathfont: Arial.otf 
---

```R

library(dplyr)
library(lmerTest)
library(SwipProteomics)

data(msstats_prot)

washc_prots <- c("Q8C2E7", "Q6PGL7", "Q3UMB9", "Q9CR27", "Q8VDD8")

# the formula to be fit:
fx <- Abundance ~ 0 + Genotype:BioFraction + (1 | Mixture) + (1 | Protein)

# fit a model to WASH complex proteins:
fm <- lmer(fx, data = msstats_prot %>% filter(Protein %in% washc_prots))

# compute Satterthwaite degrees of freedom
model_summary <- summary(fm,ddf="Satterthwaite")

```

|Term                   | Estimate|    SE|    DF| Tvalue|Pvalue    |
|:----------------------|--------:|-----:|-----:|------:|:---------|
|Control:BioFractionF4  |    6.884| 0.151| 6.909| 45.686|2.776e-09 |
|Control:BioFractionF5  |    7.168| 0.151| 6.909| 47.570|7.845e-10 |
|Control:BioFractionF6  |    7.465| 0.151| 6.909| 49.548|2.183e-09 |
|Control:BioFractionF7  |    7.495| 0.151| 6.909| 49.745|5.939e-10 |
|Control:BioFractionF8  |    7.327| 0.151| 6.909| 48.629|1.922e-09 |
|Control:BioFractionF9  |    7.138| 0.151| 6.909| 47.377|4.486e-10 |
|Control:BioFractionF10 |    7.756| 0.151| 6.909| 51.478|1.839e-09 |
|Mutant:BioFractionF4   |    5.729| 0.151| 6.909| 38.025|4.364e-10 |
|Mutant:BioFractionF5   |    5.933| 0.151| 6.909| 39.377|2.197e-09 |
|Mutant:BioFractionF6   |    6.044| 0.151| 6.909| 40.113|5.103e-10 |
|Mutant:BioFractionF7   |    6.083| 0.151| 6.909| 40.370|2.275e-09 |
|Mutant:BioFractionF8   |    5.927| 0.151| 6.909| 39.339|6.108e-10 |
|Mutant:BioFractionF9   |    5.897| 0.151| 6.909| 39.141|1.898e-09 |
|Mutant:BioFractionF10  |    6.055| 0.151| 6.909| 40.186|3.447e-10 |

```
# define a contrast
contrast <- lme4::fixef(fm)
contrast[] <- 0
contrast[grep("Mutant",names(contrast))] <- +1/7
contrast[grep("*Control*",names(contrast))] <- -1/7

# test contrast with lmerTestContrast
results <- lmerTestContrast(fm,contrast)

results %>% mutate(Contrast = "Mutant-Control") %>% unique() %>% knitr::kable()

```

|Contrast       |    log2FC| percentControl| Pvalue| Tstatistic|        SE|  DF|
|:--------------|---------:|--------------:|------:|----------:|---------:|---:|
|Mutant-Control | -1.366434|      0.3878488|      0|  -36.93673| 0.0369939| 190|
