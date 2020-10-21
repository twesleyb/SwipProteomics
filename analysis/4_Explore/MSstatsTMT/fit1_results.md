Loading SwipProteomics

```
Linear mixed model fit by REML. t-tests use Satterthwaite's method 
[lmerModLmerTest]

Formula: fx1
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

Correlation of Fixed Effects:
            GntypC GntypM BFrcF5 BFrcF6 BFrcF7 BFrcF8 BFrcF9
GenotypMtnt  0.208                                          
BioFractnF5 -0.348 -0.348                                   
BioFractnF6 -0.348 -0.348  0.500                            
BioFractnF7 -0.348 -0.348  0.500  0.500                     
BioFractnF8 -0.348 -0.348  0.500  0.500  0.500              
BioFractnF9 -0.348 -0.348  0.500  0.500  0.500  0.500       
BioFrctnF10 -0.348 -0.348  0.500  0.500  0.500  0.500  0.500

```

|protein |contrast                       |    log2FC| percentControl| Pvalue| Tstatistic|        SE|     DF|isSingular |
|:-------|:------------------------------|---------:|--------------:|------:|----------:|---------:|------:|:----------|
|Q3UMB9  |GenotypeMutant-GenotypeControl | -1.516727|      0.3494779|      0|  -9.605806| 0.1578969| 160.44|FALSE      |


|protein |symbol |contrast                       |     log2FC| percentControl| Pvalue| Padjust| Tstatistic|        SE|        DF|
|:-------|:------|:------------------------------|----------:|--------------:|------:|-------:|----------:|---------:|---------:|
|Q8C2E7  |Washc5 |GenotypeMutant-GenotypeControl | -1.3010547|      0.4058294|      0|       0| -21.351209| 0.0609359| 1077.2423|
|Q6PGL7  |Washc2 |GenotypeMutant-GenotypeControl | -1.0675369|      0.4771329|      0|       0| -12.042081| 0.0886505|  508.9758|
|P41234  |Abca2  |GenotypeMutant-GenotypeControl | -0.4095020|      0.7528832|      0|       0|  -9.302825| 0.0440191| 2064.3142|
|Q3UMB9  |Washc4 |GenotypeMutant-GenotypeControl | -1.5167270|      0.3494779|      0|       0|  -9.605806| 0.1578969|  160.4400|
|P35282  |Rab21  |GenotypeMutant-GenotypeControl | -0.5448780|      0.6854494|      0|       0|  -8.600514| 0.0633541|  996.5758|
|Q8K4Z0  |Lgi2   |GenotypeMutant-GenotypeControl | -0.4477726|      0.7331739|      0|       0|  -8.314603| 0.0538538| 1379.2022|

Total number of significant proteins: 98
