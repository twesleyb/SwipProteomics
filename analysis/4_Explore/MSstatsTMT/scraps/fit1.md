# SwipProteomics

* fit: ~Abundance0 + (1 | BioFraction) + Genotype

```
Linear mixed model fit by REML. 
t-tests use Satterthwaite's method [lmerModLmerTest]

Formula: fx1
   Data: subdat

REML criterion at convergence: 9.9

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-1.8636 -0.7523  0.1214  0.6825  1.9166 

Random effects:
 Groups      Name        Variance Std.Dev.
 BioFraction (Intercept) 0.03169  0.1780  
 Residual                0.05105  0.2259  
Number of obs: 42, groups:  BioFraction, 7

Fixed effects:
                Estimate Std. Error      df t value Pr(>|t|)    
GenotypeControl  7.13269    0.08342 8.73972   85.51 4.42e-14 ***
GenotypeMutant   5.61576    0.08342 8.73972   67.32 3.57e-13 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            GntypC
GenotypMtnt 0.651 
```
