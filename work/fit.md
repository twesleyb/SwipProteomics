# Protein-level model

`log2(rel_Intensity) ~ 0 + Condition + (1 | Mixture)`

```
Linear mixed model fit by REML. t-tests use Satterthwaite's method [
lmerModLmerTest]
Formula: fx0
   Data: swip_tmt %>% subset(Protein == swip)

REML criterion at convergence: -42.7

Scaled residuals: 
     Min       1Q   Median       3Q      Max 
-2.20225 -0.44617 -0.06154  0.47135  1.61037 

Random effects:
 Groups   Name        Variance Std.Dev.
 Mixture  (Intercept) 0.000000 0.0000  
 Residual             0.007362 0.0858  
Number of obs: 42, groups:  Mixture, 3

Fixed effects:
                     Estimate Std. Error       df t value Pr(>|t|)    
ConditionControl.F10 -0.07476    0.04954 28.00000  -1.509    0.142    
ConditionMutant.F10  -1.80336    0.04954 28.00000 -36.404  < 2e-16 ***
ConditionControl.F4  -1.00095    0.04954 28.00000 -20.206  < 2e-16 ***
ConditionMutant.F4   -2.17772    0.04954 28.00000 -43.961  < 2e-16 ***
ConditionControl.F5  -0.56207    0.04954 28.00000 -11.346 5.52e-12 ***
ConditionMutant.F5   -1.90642    0.04954 28.00000 -38.484  < 2e-16 ***
ConditionControl.F6  -0.23737    0.04954 28.00000  -4.792 4.90e-05 ***
ConditionMutant.F6   -1.72952    0.04954 28.00000 -34.913  < 2e-16 ***
ConditionControl.F7  -0.23383    0.04954 28.00000  -4.720 5.96e-05 ***
ConditionMutant.F7   -1.77800    0.04954 28.00000 -35.892  < 2e-16 ***
ConditionControl.F8  -0.40120    0.04954 28.00000  -8.099 8.10e-09 ***
ConditionMutant.F8   -1.75399    0.04954 28.00000 -35.407  < 2e-16 ***
ConditionControl.F9  -0.57757    0.04954 28.00000 -11.659 2.93e-12 ***
ConditionMutant.F9   -1.75181    0.04954 28.00000 -35.363  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation matrix not shown by default, as p = 14 > 12.
Use print(x, correlation=TRUE)  or
    vcov(x)        if you need it

optimizer (nloptwrap) convergence code: 0 (OK)
boundary (singular) fit: see ?isSingular
```



|Contrast       |    log2FC| percentControl|        SE| Tstatistic| Pvalue| DF|       S2|isSingular |
|:--------------|---------:|--------------:|---------:|----------:|------:|--:|--------:|:----------|
|Mutant-Control | -1.401866|      0.3784393| 0.0264791|  -52.94235|      0| 28| 0.007362|TRUE       |

## Module-level Model

`log2(rel_Intensity) ~ 0 + Condition + (1 | Protein)`


```
Linear mixed model fit by REML. t-tests use Satterthwaite's method [
lmerModLmerTest]
Formula: fx1
   Data: swip_tmt %>% subset(Protein %in% washc_prots)

REML criterion at convergence: 59.6

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-3.5264 -0.4525 -0.0063  0.5300  2.3597 

Random effects:
 Groups   Name        Variance Std.Dev.
 Protein  (Intercept) 0.03751  0.1937  
 Residual             0.06457  0.2541  
Number of obs: 168, groups:  Protein, 4

Fixed effects:
                     Estimate Std. Error      df t value Pr(>|t|)    
ConditionControl.F10  -0.1323     0.1215  6.8223  -1.089 0.312963    
ConditionMutant.F10   -1.8102     0.1215  6.8223 -14.901 1.86e-06 ***
ConditionControl.F4   -0.9709     0.1215  6.8223  -7.992 0.000105 ***
ConditionMutant.F4    -2.1107     0.1215  6.8223 -17.375 6.69e-07 ***
ConditionControl.F5   -0.5311     0.1215  6.8223  -4.372 0.003473 ** 
ConditionMutant.F5    -1.8258     0.1215  6.8223 -15.029 1.76e-06 ***
ConditionControl.F6   -0.2199     0.1215  6.8223  -1.810 0.114291    
ConditionMutant.F6    -1.6551     0.1215  6.8223 -13.624 3.37e-06 ***
ConditionControl.F7   -0.1727     0.1215  6.8223  -1.422 0.199148    
ConditionMutant.F7    -1.6723     0.1215  6.8223 -13.766 3.15e-06 ***
ConditionControl.F8   -0.3751     0.1215  6.8223  -3.088 0.018206 *  
ConditionMutant.F8    -1.7523     0.1215  6.8223 -14.424 2.31e-06 ***
ConditionControl.F9   -0.6418     0.1215  6.8223  -5.283 0.001242 ** 
ConditionMutant.F9    -1.8748     0.1215  6.8223 -15.433 1.47e-06 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation matrix not shown by default, as p = 14 > 12.
Use print(x, correlation=TRUE)  or
    vcov(x)        if you need it
```


|Contrast       |    log2FC| percentControl|        SE| Tstatistic| Pvalue|  DF|        S2|isSingular |
|:--------------|---------:|--------------:|---------:|----------:|------:|---:|---------:|:----------|
|Mutant-Control | -1.379633|      0.3843165| 0.0392109|  -35.18497|      0| 151| 0.0645747|FALSE      |
