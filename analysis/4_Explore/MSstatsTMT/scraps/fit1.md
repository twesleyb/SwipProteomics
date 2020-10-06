```
Linear mixed model fit by REML ['lmerModLmerTest']
Formula: Abundance ~ 0 + (1|BioFraction) + Genotype
   Data: subdat
REML criterion at convergence: 9.9247
Random effects:
 Groups      Name        Std.Dev.
 BioFraction (Intercept) 0.1780  
 Residual                0.2259  
Number of obs: 42, groups:  BioFraction, 7
Fixed Effects:
GenotypeControl   GenotypeMutant  
          7.133            5.616  
```



Contrast:

| GenotypeControl| GenotypeMutant|
|---------------:|--------------:|
|              -1|              1|

Results:

|Protein | P.value|   Log2.FC| t.statistic|  variance|
|:-------|-------:|---------:|-----------:|---------:|
|Q3UMB9  |       0| -1.516927|   -21.75455| 0.0048622|
