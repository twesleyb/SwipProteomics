>Loading SwipProteomics
`lmer: ~Abundance+(1 | BioFraction:Genotype) + Genotype`

|Protein | P.value|   Log2.FC| t.statistic|  variance|
|:-------|-------:|---------:|-----------:|---------:|
|Q3UMB9  |       0| -8.649613|   -44.00426| 0.0386371|
> is log2FC overestimated?

```
Linear mixed model fit by REML. t-tests use Satterthwaite's method 
[lmerModLmerTest]

Formula: Abundance+(1 | BioFraction:Genotype) + Genotype
   Data: subdat

REML criterion at convergence: 10.4

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-1.8695 -0.3797  0.1429  0.6688  1.4651 

Random effects:
 Groups               Name        Variance Std.Dev.
 BioFraction:Genotype (Intercept) 0.03938  0.1984  
 Residual                         0.04413  0.2101  
Number of obs: 42, groups:  BioFraction:Genotype, 14

Fixed effects:
               Estimate Std. Error       df t value Pr(>|t|)    
(Intercept)     7.13269    0.08791 12.00000   81.14  < 2e-16 ***
GenotypeMutant -1.51693    0.12432 12.00000  -12.20 4.01e-08 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            (Intr)
GenotypMtnt -0.707
```
