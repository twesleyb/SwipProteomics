Example of a protein-wise lmer model:

* NOTE: Random effect of Run: (1|Run)
* `Run` = c(M1_1, M2_1, M3_1) ~ Experiment/Batch

```
Linear mixed model fit by REML. t-tests use Satterthwaite's method 
[lmerModLmerTest]

Formula: Abundance ~ 1 + (1 | Run) + Condition
   Data: data

REML criterion at convergence: 25.5

Scaled residuals: 
     Min       1Q   Median       3Q      Max 
-1.52171 -0.44833  0.02691  0.54551  1.95544 

Random effects:
 Groups   Name        Variance Std.Dev.
 Run      (Intercept) 0.003653 0.06044 
 Residual             0.081285 0.28511 
Number of obs: 42, groups:  Run, 3

Fixed effects:
                    Estimate Std. Error       df t value Pr(>|t|)    
(Intercept)          6.36282    0.16826 27.34251  37.815  < 2e-16 ***
ConditionControl.F4 -1.60807    0.23279 26.00000  -6.908 2.47e-07 ***
ConditionControl.F5 -1.51464    0.23279 26.00000  -6.507 6.75e-07 ***
ConditionControl.F6 -1.31989    0.23279 26.00000  -5.670 5.79e-06 ***
ConditionControl.F7 -1.11641    0.23279 26.00000  -4.796 5.76e-05 ***
ConditionControl.F8 -0.58089    0.23279 26.00000  -2.495 0.019259 *  
ConditionControl.F9 -0.01111    0.23279 26.00000  -0.048 0.962312    
ConditionMutant.F10 -0.27203    0.23279 26.00000  -1.169 0.253181    
ConditionMutant.F4  -1.19421    0.23279 26.00000  -5.130 2.39e-05 ***
ConditionMutant.F5  -1.44978    0.23279 26.00000  -6.228 1.37e-06 ***
ConditionMutant.F6  -1.15739    0.23279 26.00000  -4.972 3.62e-05 ***
ConditionMutant.F7  -0.95663    0.23279 26.00000  -4.109 0.000351 ***
ConditionMutant.F8  -0.52631    0.23279 26.00000  -2.261 0.032363 *  
ConditionMutant.F9   0.03264    0.23279 26.00000   0.140 0.889571    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
