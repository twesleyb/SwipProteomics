## Loading SwipProteomics

Fitting protein-level mixed-models.

## fit: Abundance ~ Condition + (1|Mixture)
```
Linear mixed model fit by REML. t-tests use Satterthwaite's method [
lmerModLmerTest]
Formula: fx
   Data: msstats_prot %>% subset(Protein == swip)

REML criterion at convergence: 3.7

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-1.5030 -0.6089 -0.1463  0.7474  1.4302 

Random effects:
 Groups   Name        Variance Std.Dev.
 Mixture  (Intercept) 0.009596 0.09796 
 Residual             0.034418 0.18552 
Number of obs: 42, groups:  Mixture, 3

Fixed effects:
                    Estimate Std. Error      df t value Pr(>|t|)    
(Intercept)           7.6187     0.1211 17.3059  62.899  < 2e-16 ***
ConditionControl.F4  -0.9077     0.1515 26.0000  -5.993 2.51e-06 ***
ConditionControl.F5  -0.6731     0.1515 26.0000  -4.444 0.000146 ***
ConditionControl.F6  -0.3786     0.1515 26.0000  -2.499 0.019080 *  
ConditionControl.F7  -0.2976     0.1515 26.0000  -1.965 0.060209 .  
ConditionControl.F8  -0.4891     0.1515 26.0000  -3.229 0.003356 ** 
ConditionControl.F9  -0.6642     0.1515 26.0000  -4.385 0.000170 ***
ConditionMutant.F10  -1.8343     0.1515 26.0000 -12.109 3.43e-12 ***
ConditionMutant.F4   -2.2144     0.1515 26.0000 -14.619 4.70e-14 ***
ConditionMutant.F5   -2.0513     0.1515 26.0000 -13.542 2.75e-13 ***
ConditionMutant.F6   -1.9785     0.1515 26.0000 -13.061 6.26e-13 ***
ConditionMutant.F7   -1.9870     0.1515 26.0000 -13.118 5.68e-13 ***
ConditionMutant.F8   -2.1259     0.1515 26.0000 -14.035 1.21e-13 ***
ConditionMutant.F9   -1.8377     0.1515 26.0000 -12.132 3.29e-12 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

```


## fit: Abundance ~ 0 + Condition + (1|Mixture)
```
Linear mixed model fit by REML. t-tests use Satterthwaite's method [
lmerModLmerTest]
Formula: fx0
   Data: msstats_prot %>% subset(Protein == swip)

REML criterion at convergence: 3.7

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-1.5030 -0.6089 -0.1463  0.7474  1.4302 

Random effects:
 Groups   Name        Variance Std.Dev.
 Mixture  (Intercept) 0.009596 0.09796 
 Residual             0.034418 0.18552 
Number of obs: 42, groups:  Mixture, 3

Fixed effects:
                     Estimate Std. Error      df t value Pr(>|t|)    
ConditionControl.F10   7.6187     0.1211 17.3059   62.90   <2e-16 ***
ConditionControl.F4    6.7110     0.1211 17.3059   55.41   <2e-16 ***
ConditionControl.F5    6.9456     0.1211 17.3059   57.34   <2e-16 ***
ConditionControl.F6    7.2401     0.1211 17.3059   59.77   <2e-16 ***
ConditionControl.F7    7.3211     0.1211 17.3059   60.44   <2e-16 ***
ConditionControl.F8    7.1296     0.1211 17.3059   58.86   <2e-16 ***
ConditionControl.F9    6.9545     0.1211 17.3059   57.41   <2e-16 ***
ConditionMutant.F10    5.7844     0.1211 17.3059   47.76   <2e-16 ***
ConditionMutant.F4     5.4043     0.1211 17.3059   44.62   <2e-16 ***
ConditionMutant.F5     5.5674     0.1211 17.3059   45.96   <2e-16 ***
ConditionMutant.F6     5.6402     0.1211 17.3059   46.56   <2e-16 ***
ConditionMutant.F7     5.6317     0.1211 17.3059   46.49   <2e-16 ***
ConditionMutant.F8     5.4928     0.1211 17.3059   45.35   <2e-16 ***
ConditionMutant.F9     5.7810     0.1211 17.3059   47.73   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```


## fit: Abundance ~ 1 + Condition + (1|Mixture)
```
Linear mixed model fit by REML. t-tests use Satterthwaite's method [
lmerModLmerTest]
Formula: fx1
   Data: msstats_prot %>% subset(Protein == swip)

REML criterion at convergence: 3.7

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-1.5030 -0.6089 -0.1463  0.7474  1.4302 

Random effects:
 Groups   Name        Variance Std.Dev.
 Mixture  (Intercept) 0.009596 0.09796 
 Residual             0.034418 0.18552 
Number of obs: 42, groups:  Mixture, 3

Fixed effects:
                    Estimate Std. Error      df t value Pr(>|t|)    
(Intercept)           7.6187     0.1211 17.3059  62.899  < 2e-16 ***
ConditionControl.F4  -0.9077     0.1515 26.0000  -5.993 2.51e-06 ***
ConditionControl.F5  -0.6731     0.1515 26.0000  -4.444 0.000146 ***
ConditionControl.F6  -0.3786     0.1515 26.0000  -2.499 0.019080 *  
ConditionControl.F7  -0.2976     0.1515 26.0000  -1.965 0.060209 .  
ConditionControl.F8  -0.4891     0.1515 26.0000  -3.229 0.003356 ** 
ConditionControl.F9  -0.6642     0.1515 26.0000  -4.385 0.000170 ***
ConditionMutant.F10  -1.8343     0.1515 26.0000 -12.109 3.43e-12 ***
ConditionMutant.F4   -2.2144     0.1515 26.0000 -14.619 4.70e-14 ***
ConditionMutant.F5   -2.0513     0.1515 26.0000 -13.542 2.75e-13 ***
ConditionMutant.F6   -1.9785     0.1515 26.0000 -13.061 6.26e-13 ***
ConditionMutant.F7   -1.9870     0.1515 26.0000 -13.118 5.68e-13 ***
ConditionMutant.F8   -2.1259     0.1515 26.0000 -14.035 1.21e-13 ***
ConditionMutant.F9   -1.8377     0.1515 26.0000 -12.132 3.29e-12 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

## Assessing protein-level intra-BioFraction comparisons.

### F7

fit: Abundance ~ Condition + (1|Mixture)
```
|Contrast |    log2FC| percentControl|        SE| Tstatistic| Pvalue| DF|isSingular |
|:--------|---------:|--------------:|---------:|----------:|------:|--:|:----------|
|F7       | -1.689393|      0.3100573| 0.1514779|  -11.15273|      0| 26|FALSE      |
```

fit: Abundance ~ 0 + Condition + (1|Mixture)
```
|Contrast |    log2FC| percentControl|        SE| Tstatistic| Pvalue| DF|isSingular |
|:--------|---------:|--------------:|---------:|----------:|------:|--:|:----------|
|F7       | -1.689393|      0.3100573| 0.1514779|  -11.15273|      0| 26|FALSE      |
```

fit: Abundance ~ 1 + Condition + (1|Mixture)
```
|Contrast |    log2FC| percentControl|        SE| Tstatistic| Pvalue| DF|isSingular |
|:--------|---------:|--------------:|---------:|----------:|------:|--:|:----------|
|F7       | -1.689393|      0.3100573| 0.1514779|  -11.15273|      0| 26|FALSE      |
```

### F10

fit: Abundance ~ Condition + (1|Mixture)
```
|Contrast |    log2FC| percentControl|        SE| Tstatistic| Pvalue| DF|isSingular |
|:--------|---------:|--------------:|---------:|----------:|------:|--:|:----------|
|F10      | -1.834294|      0.2804287| 0.1514779|  -12.10931|      0| 26|FALSE      |
```

fit: Abundance ~ 0 + Condition + (1|Mixture)
```
|Contrast |    log2FC| percentControl|        SE| Tstatistic| Pvalue| DF|isSingular |
|:--------|---------:|--------------:|---------:|----------:|------:|--:|:----------|
|F10      | -1.834294|      0.2804287| 0.1514779|  -12.10931|      0| 26|FALSE      |
```

fit:Abundance ~ 1 + Condition + (1|Mixture)
```
|Contrast |    log2FC| percentControl|        SE| Tstatistic| Pvalue| DF|isSingular |
|:--------|---------:|--------------:|---------:|----------:|------:|--:|:----------|
|F10      | -1.834294|      0.2804287| 0.1514779|  -12.10931|      0| 26|FALSE      |
```


## Assessing 'Mutant-Control' comparison.

fit: Abundance ~ Condition + (1|Mixture)
```
|Contrast       |    log2FC| percentControl|        SE| Tstatistic| Pvalue| DF|isSingular |
|:--------------|---------:|--------------:|---------:|----------:|------:|--:|:----------|
|Mutant-Control | -1.516956|      0.3494225| 0.0572533|  -26.49552|      0| 26|FALSE      |
```

fit: Abundance ~ 0 + Condition + (1|Mixture)
```
|Contrast       |    log2FC| percentControl|        SE| Tstatistic| Pvalue| DF|isSingular |
|:--------------|---------:|--------------:|---------:|----------:|------:|--:|:----------|
|Mutant-Control | -1.516956|      0.3494225| 0.0572533|  -26.49552|      0| 26|FALSE      |
```

fit: Abundance ~ 1 + Condition + (1|Mixture)
```
|Contrast       |    log2FC| percentControl|        SE| Tstatistic| Pvalue| DF|isSingular |
|:--------------|---------:|--------------:|---------:|----------:|------:|--:|:----------|
|Mutant-Control | -1.516956|      0.3494225| 0.0572533|  -26.49552|      0| 26|FALSE      |
```


## Fitting module-level mixed-models.

fit: Abundance ~ Condition + Protein + (1|Mixture)
```
Linear mixed model fit by REML. t-tests use Satterthwaite's method [
lmerModLmerTest]
Formula: fx
   Data: msstats_prot %>% filter(Protein %in% washc_prots)

REML criterion at convergence: 93.4

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-3.2580 -0.5918  0.0858  0.6241  2.7300 

Random effects:
 Groups   Name        Variance Std.Dev.
 Mixture  (Intercept) 0.003304 0.05748 
 Residual             0.071832 0.26801 
Number of obs: 210, groups:  Mixture, 3

Fixed effects:
                     Estimate Std. Error        df t value Pr(>|t|)    
(Intercept)           7.49356    0.08520  43.19599  87.957  < 2e-16 ***
ConditionControl.F4  -0.87286    0.09787 190.00000  -8.919 3.92e-16 ***
ConditionControl.F5  -0.58882    0.09787 190.00000  -6.017 8.98e-09 ***
ConditionControl.F6  -0.29095    0.09787 190.00000  -2.973  0.00333 ** 
ConditionControl.F7  -0.26114    0.09787 190.00000  -2.668  0.00828 ** 
ConditionControl.F8  -0.42892    0.09787 190.00000  -4.383 1.94e-05 ***
ConditionControl.F9  -0.61773    0.09787 190.00000  -6.312 1.90e-09 ***
ConditionMutant.F10  -1.70151    0.09787 190.00000 -17.386  < 2e-16 ***
ConditionMutant.F4   -2.02720    0.09787 190.00000 -20.714  < 2e-16 ***
ConditionMutant.F5   -1.82345    0.09787 190.00000 -18.632  < 2e-16 ***
ConditionMutant.F6   -1.71322    0.09787 190.00000 -17.506  < 2e-16 ***
ConditionMutant.F7   -1.67446    0.09787 190.00000 -17.110  < 2e-16 ***
ConditionMutant.F8   -1.82874    0.09787 190.00000 -18.686  < 2e-16 ***
ConditionMutant.F9   -1.85850    0.09787 190.00000 -18.990  < 2e-16 ***
ProteinQ6PGL7         0.35219    0.05849 190.00000   6.022 8.75e-09 ***
ProteinQ8C2E7         0.14150    0.05849 190.00000   2.419  0.01649 *  
ProteinQ8VDD8         0.08571    0.05849 190.00000   1.465  0.14446    
ProteinQ9CR27         0.73278    0.05849 190.00000  12.529  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

fit: Abundance ~ 0 + Genotype + BioFraction + Protein + (1|Mixture)
```
Linear mixed model fit by REML. t-tests use Satterthwaite's method [
lmerModLmerTest]
Formula: fx0
   Data: msstats_prot %>% filter(Protein %in% washc_prots)

REML criterion at convergence: 98.5

Scaled residuals: 
     Min       1Q   Median       3Q      Max 
-2.76559 -0.65135  0.04919  0.68597  2.55047 

Random effects:
 Groups   Name        Variance Std.Dev.
 Mixture  (Intercept) 0.003228 0.05681 
 Residual             0.077166 0.27779 
Number of obs: 210, groups:  Mixture, 3

Fixed effects:
                 Estimate Std. Error        df t value Pr(>|t|)    
GenotypeMutant    5.36020    0.07406  26.74555  72.373  < 2e-16 ***
GenotypeControl   6.72686    0.07406  26.74555  90.826  < 2e-16 ***
BioFractionF5     0.24389    0.07172 196.00000   3.400 0.000815 ***
BioFractionF6     0.44794    0.07172 196.00000   6.245 2.57e-09 ***
BioFractionF7     0.48223    0.07172 196.00000   6.723 1.89e-10 ***
BioFractionF8     0.32120    0.07172 196.00000   4.478 1.28e-05 ***
BioFractionF9     0.21191    0.07172 196.00000   2.954 0.003515 ** 
BioFractionF10    0.59927    0.07172 196.00000   8.355 1.18e-14 ***
ProteinQ6PGL7     0.35219    0.06062 196.00000   5.810 2.49e-08 ***
ProteinQ8C2E7     0.14150    0.06062 196.00000   2.334 0.020593 *  
ProteinQ8VDD8     0.08571    0.06062 196.00000   1.414 0.158992    
ProteinQ9CR27     0.73278    0.06062 196.00000  12.088  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

fit: Abundance ~ 1 + Genotype + BioFraction + Protein + (1|Mixture)
```
Linear mixed model fit by REML. t-tests use Satterthwaite's method [
lmerModLmerTest]
Formula: fx1
   Data: msstats_prot %>% filter(Protein %in% washc_prots)

REML criterion at convergence: 98.5

Scaled residuals: 
     Min       1Q   Median       3Q      Max 
-2.76559 -0.65135  0.04919  0.68597  2.55047 

Random effects:
 Groups   Name        Variance Std.Dev.
 Mixture  (Intercept) 0.003228 0.05681 
 Residual             0.077166 0.27779 
Number of obs: 210, groups:  Mixture, 3

Fixed effects:
                 Estimate Std. Error        df t value Pr(>|t|)    
(Intercept)       5.36020    0.07406  26.74555  72.373  < 2e-16 ***
GenotypeControl   1.36666    0.03834 196.00000  35.647  < 2e-16 ***
BioFractionF5     0.24389    0.07172 196.00000   3.400 0.000815 ***
BioFractionF6     0.44794    0.07172 196.00000   6.245 2.57e-09 ***
BioFractionF7     0.48223    0.07172 196.00000   6.723 1.89e-10 ***
BioFractionF8     0.32120    0.07172 196.00000   4.478 1.28e-05 ***
BioFractionF9     0.21191    0.07172 196.00000   2.954 0.003515 ** 
BioFractionF10    0.59927    0.07172 196.00000   8.355 1.18e-14 ***
ProteinQ6PGL7     0.35219    0.06062 196.00000   5.810 2.49e-08 ***
ProteinQ8C2E7     0.14150    0.06062 196.00000   2.334 0.020593 *  
ProteinQ8VDD8     0.08571    0.06062 196.00000   1.414 0.158992    
ProteinQ9CR27     0.73278    0.06062 196.00000  12.088  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```


fit: Abundance ~ Condition + Protein + (1|Mixture)
```
|Contrast       |    log2FC| percentControl|        SE| Tstatistic| Pvalue|  DF|isSingular |
|:--------------|---------:|--------------:|---------:|----------:|------:|---:|:----------|
|Mutant-Control | -1.366663|      0.3877871| 0.0369896|  -36.94728|      0| 190|FALSE      |
```

fit: Abundance ~ 0 + Genotype + BioFraction + Protein + (1|Mixture)
```
|Contrast       |    log2FC| percentControl|        SE| Tstatistic| Pvalue|  DF|isSingular |
|:--------------|---------:|--------------:|---------:|----------:|------:|---:|:----------|
|Mutant-Control | -1.366663|      0.3877871| 0.0383383|  -35.64748|      0| 196|FALSE      |
```

fit: Abundance ~ 1 + Genotype + BioFraction + Protein + (1|Mixture)
```
|Contrast       |    log2FC| percentControl|        SE| Tstatistic| Pvalue|  DF|isSingular |
|:--------------|---------:|--------------:|---------:|----------:|------:|---:|:----------|
|Mutant-Control | -1.366663|      0.3877871| 0.0383383|  -35.64748|      0| 196|FALSE      |
```

## MSstatsTMT results:
fit reduced model with a single Mixture and multiple Runs:
Abundance ~ 1 + (1 | Run) + Condition

```
|Protein |Label                  |    log2FC|        SE|       DF| pvalue| adj.pvalue|issue |
|:-------|:----------------------|---------:|---------:|--------:|------:|----------:|:-----|
|Q3UMB9  |Mutant.F4-Control.F4   | -1.306658| 0.1514779| 26.00002|      0|          0|NA    |
|Q3UMB9  |Mutant.F5-Control.F5   | -1.378141| 0.1514779| 26.00002|      0|          0|NA    |
|Q3UMB9  |Mutant.F6-Control.F6   | -1.599893| 0.1514779| 26.00002|      0|          0|NA    |
|Q3UMB9  |Mutant.F7-Control.F7   | -1.689393| 0.1514779| 26.00002|      0|          0|NA    |
|Q3UMB9  |Mutant.F8-Control.F8   | -1.636860| 0.1514779| 26.00002|      0|          0|NA    |
|Q3UMB9  |Mutant.F9-Control.F9   | -1.173450| 0.1514779| 26.00002|      0|          0|NA    |
|Q3UMB9  |Mutant.F10-Control.F10 | -1.834294| 0.1514779| 26.00002|      0|          0|NA    |
|Q3UMB9  |Mutant-Control         | -1.516956| 0.0572533| 26.00002|      0|          0|NA    |
```
