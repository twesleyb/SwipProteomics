# Module Level analysis

The function (fx) to be fit:
```
fx <- formula("Abundance ~ 0 + Genotype + BioFraction + (1|Protein)")
```

Fit the model for a given module:
```
fm <- lmerTest::lmer(fx,msstats_filt %>% filter(Module=="M1")

summary(fm,ddf="Satterthwaite")

```

#### Results
```
Linear mixed model fit by REML. t-tests use Satterthwaite's method 
[lmerModLmerTest]
Formula: fx
   Data: subdat

REML criterion at convergence: 13327.6

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-4.2553 -0.5179  0.0413  0.6088  3.6346 

Random effects:
 Groups   Name        Variance Std.Dev.
 Protein  (Intercept) 0.4510   0.6715  
 Residual             0.1877   0.4333  
Number of obs: 10416, groups:  Protein, 248

Fixed effects:
                  Estimate Std. Error         df t value Pr(>|t|)    
GenotypeMutant   8.523e+00  4.430e-02  2.821e+02  192.38   <2e-16 ***
GenotypeControl  8.559e+00  4.430e-02  2.821e+02  193.21   <2e-16 ***
BioFractionF5   -8.246e-01  1.588e-02  1.016e+04  -51.91   <2e-16 ***
BioFractionF6   -1.755e+00  1.588e-02  1.016e+04 -110.48   <2e-16 ***
BioFractionF7   -2.261e+00  1.588e-02  1.016e+04 -142.33   <2e-16 ***
BioFractionF8   -2.545e+00  1.588e-02  1.016e+04 -160.19   <2e-16 ***
BioFractionF9   -2.832e+00  1.588e-02  1.016e+04 -178.28   <2e-16 ***
BioFractionF10  -2.668e+00  1.588e-02  1.016e+04 -167.98   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            GntypM GntypC BFrcF5 BFrcF6 BFrcF7 BFrcF8 BFrcF9
GentypCntrl  0.982                                          
BioFractnF5 -0.179 -0.179                                   
BioFractnF6 -0.179 -0.179  0.500                            
BioFractnF7 -0.179 -0.179  0.500  0.500                     
BioFractnF8 -0.179 -0.179  0.500  0.500  0.500              
BioFractnF9 -0.179 -0.179  0.500  0.500  0.500  0.500       
BioFrctnF10 -0.179 -0.179  0.500  0.500  0.500  0.500  0.500
```


|       R2m|       R2c|
|---------:|---------:|
| 0.6022799| 0.8831005|

data(partition)
data(msstats_prot)
data(fx0)
data(fx1)




