## Model0 -- Intrafraction Comparisons

```
Linear mixed model fit by REML. t-tests use Satterthwaite's method 
[lmerModLmerTest]
Formula: Abundance ~ 0 + Genotype.BioFraction + (1|Mixture)
   Data: msstats_prot %>% filter(Protein == swip)

REML criterion at convergence: 3.7

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-1.4998 -0.6099 -0.1471  0.7489  1.4237 

Random effects:
 Groups   Name        Variance Std.Dev.
 Mixture  (Intercept) 0.009667 0.09832 
 Residual             0.034474 0.18567 
Number of obs: 42, groups:  Mixture, 3

Fixed effects:
                     Estimate Std. Error      df t value Pr(>|t|)    
ConditionControl.F10   7.6192     0.1213 17.2481   62.81   <2e-16 ***
ConditionControl.F4    6.7117     0.1213 17.2481   55.33   <2e-16 ***
ConditionControl.F5    6.9462     0.1213 17.2481   57.26   <2e-16 ***
ConditionControl.F6    7.2407     0.1213 17.2481   59.69   <2e-16 ***
ConditionControl.F7    7.3216     0.1213 17.2481   60.36   <2e-16 ***
ConditionControl.F8    7.1298     0.1213 17.2481   58.78   <2e-16 ***
ConditionControl.F9    6.9549     0.1213 17.2481   57.34   <2e-16 ***
ConditionMutant.F10    5.7850     0.1213 17.2481   47.69   <2e-16 ***
ConditionMutant.F4     5.4051     0.1213 17.2481   44.56   <2e-16 ***
ConditionMutant.F5     5.5681     0.1213 17.2481   45.90   <2e-16 ***
ConditionMutant.F6     5.6415     0.1213 17.2481   46.51   <2e-16 ***
ConditionMutant.F7     5.6331     0.1213 17.2481   46.44   <2e-16 ***
ConditionMutant.F8     5.4930     0.1213 17.2481   45.28   <2e-16 ***
ConditionMutant.F9     5.7813     0.1213 17.2481   47.66   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Goodness of fit

|       R2m|       R2c|
|---------:|---------:|
| 0.9351423| 0.9493463|

_Marginal and conditional R-squared for mixed effect model (Nakagawa et al., 2014)._


### Statistical Results

|protein |contrast             |    log2FC| Pvalue| Tstatistic|        SE|       DF|
|:-------|:--------------------|---------:|------:|----------:|---------:|--------:|
|Q3UMB9  |Mutant.F7-Control.F7 | -1.688576|      0|  -11.13831| 0.1516008| 25.99987|

----

## Model1 -- Mutant-Control Comparison

```
Linear mixed model fit by REML. t-tests use Satterthwaite's method 
[lmerModLmerTest]
Formula: Abundance ~ 0 + Genotype + BioFraction + (1|Subject)
   Data: msstats_prot %>% filter(Protein == swip)

REML criterion at convergence: -7.2

Scaled residuals: 
     Min       1Q   Median       3Q      Max 
-2.12532 -0.55464 -0.09399  0.58495  1.98324 

Random effects:
 Groups   Name        Variance Std.Dev.
 Subject  (Intercept) 0.03412  0.1847  
 Residual             0.02293  0.1514  
Number of obs: 42, groups:  Subject, 6

Fixed effects:
                Estimate Std. Error       df t value Pr(>|t|)    
GenotypeControl  6.81676    0.12546  6.32017  54.332 1.11e-09 ***
GenotypeMutant   5.30003    0.12546  6.32017  42.243 5.44e-09 ***
BioFractionF5    0.19875    0.08742 30.00000   2.273 0.030328 *  
BioFractionF6    0.38272    0.08742 30.00000   4.378 0.000134 ***
BioFractionF7    0.41895    0.08742 30.00000   4.792 4.19e-05 ***
BioFractionF8    0.25304    0.08742 30.00000   2.894 0.007018 ** 
BioFractionF9    0.30969    0.08742 30.00000   3.542 0.001319 ** 
BioFractionF10   0.64373    0.08742 30.00000   7.363 3.34e-08 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Goodness of fit

|       R2m|       R2c|
|---------:|---------:|
| 0.9162647| 0.9663459|

### Statistical Results

|protein |contrast       |    log2FC|    Pvalue| Tstatistic|        SE|       DF|
|:-------|:--------------|---------:|---------:|----------:|---------:|--------:|
|Q3UMB9  |Mutant-Control | -1.516727| 0.0006566|  -9.605806| 0.1578969| 3.999999|
