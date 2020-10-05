> Loading SwipProteomics

`lmer: ~Abundance1 + (1 | Run) + Condition`

|Comparison             |Protein | Log2 Fold Change| P-value| P-adjust|        SE|       DF|
|:----------------------|:-------|----------------:|-------:|--------:|---------:|--------:|
|Mutant.F10-Control.F10 |Q3UMB9  |        -1.835485|       0|        0| 0.1513823| 26.00002|
|Mutant.F4-Control.F4   |Q3UMB9  |        -1.307498|       0|        0| 0.1513823| 26.00002|
|Mutant.F5-Control.F5   |Q3UMB9  |        -1.377369|       0|        0| 0.1513823| 26.00002|
|Mutant.F6-Control.F6   |Q3UMB9  |        -1.598311|       0|        0| 0.1513823| 26.00002|
|Mutant.F7-Control.F7   |Q3UMB9  |        -1.688940|       0|        0| 0.1513823| 26.00002|
|Mutant.F8-Control.F8   |Q3UMB9  |        -1.637331|       0|        0| 0.1513823| 26.00002|
|Mutant.F9-Control.F9   |Q3UMB9  |        -1.173553|       0|        0| 0.1513823| 26.00002|

```
Linear mixed model fit by REML. t-tests use Satterthwaite's method 
[lmerModLmerTest]
Formula: Abundance1 + (1 | Run) + Condition
   Data: msstats_prot %>% filter(Protein == prot)

REML criterion at convergence: 3.7

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-1.5019 -0.6087 -0.1490  0.7500  1.4234 

Random effects:
 Groups   Name        Variance Std.Dev.
 Run      (Intercept) 0.009755 0.09877 
 Residual             0.034375 0.18540 
Number of obs: 42, groups:  Run, 3

Fixed effects:
                    Estimate Std. Error      df t value Pr(>|t|)    
(Intercept)           7.6219     0.1213 17.1227  62.843  < 2e-16 ***
ConditionControl.F4  -0.9095     0.1514 26.0000  -6.008 2.41e-06 ***
ConditionControl.F5  -0.6757     0.1514 26.0000  -4.464 0.000138 ***
ConditionControl.F6  -0.3812     0.1514 26.0000  -2.518 0.018302 *  
ConditionControl.F7  -0.3003     0.1514 26.0000  -1.983 0.057978 .  
ConditionControl.F8  -0.4917     0.1514 26.0000  -3.248 0.003198 ** 
ConditionControl.F9  -0.6658     0.1514 26.0000  -4.398 0.000165 ***
ConditionMutant.F10  -1.8355     0.1514 26.0000 -12.125 3.33e-12 ***
ConditionMutant.F4   -2.2170     0.1514 26.0000 -14.645 4.51e-14 ***
ConditionMutant.F5   -2.0531     0.1514 26.0000 -13.562 2.66e-13 ***
ConditionMutant.F6   -1.9795     0.1514 26.0000 -13.076 6.11e-13 ***
ConditionMutant.F7   -1.9892     0.1514 26.0000 -13.140 5.46e-13 ***
ConditionMutant.F8   -2.1290     0.1514 26.0000 -14.064 1.15e-13 ***
ConditionMutant.F9   -1.8394     0.1514 26.0000 -12.150 3.18e-12 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation matrix not shown by default, as p = 14 > 12.
```
