## Pairwise comparisons between conditions:
Model-based testing for differentially abundant protins between pairs of
conditions was carried out through a contrast of condition means.
Beta 1 .. Bc the parameters associated with Condition (i.e. the levels of
Condition) and l1 .. 1c a vector of coefficients where lc = 0
Contrast vector: l = (1,-1,0,0) 
null = lcBc = 0

estimate the parameters with REML [see 7 8 9] to obtain Beta hat,
the contrast, and cooresponding t-statistic [see 10]

```
t = l^T * Beta^hat / sqrt(l * s^2 * V^hat * l^T)
```

* `s^2` is the estimate of theta2 [eq S1]
* `V^hat` is unscaled variance-covariance matrix of B^hat
* `s^2 * V^hat` is the (scaled) variance-covariance matrix of B^hat
* `sqrt(l * s^2 * V^hat * l^T)` is the standard error of the contrast
> I need to learn how to write latex eq.

estimates of s2 and Vhat are obtained by REML the matrix V is a function of
estimates sigma^2hat_M for random effect of mixtures sigma^2hat_T for random
effect of technical replicates
sigma^2_S for the random effect of subject

standard error of the contrast takes into account both technical variancee
sigma^2_M, sigma^2_T, sigma^2 and biological variance sigma^2_S

the degrees of freedom for the t-statstic are derifed by Satterhwaite approx
(see 11)

The valculation of VAR(s^2*l*V^hat*l^T) see 10 and 12

when the number of biological replicates in each condition is small,
the Empirical Bayes moderation as implemented by the R package limma.
Estimate of the error variance sigma^2 follows a scaled chi-squared distribution
with v degrees of freem
s^2|sigma^2  ~ sigma^2/v CHI^2_d0

dof is estimated by lmerTEST [see 12]

parameters d0 and s0^2 are estimated from the distribution of observed sigma^2
of all proteins using an EB approach (limma)

posterior variance estimate is incorporated into residual variance of each
protein

moderate t statistic = t^~ = l^T * beta^hat / sqrt(l * s^2~ * V^hat * l^T)

with df + d0 dof

[8] Patterson, H. D. and Thompson, R. (1971). Recovery of inter-block
information when block sizes are unequal.
Biometrika, 58, 545–554.

REF REML:

[9] Harville, D. A. (1977). Maximum likelihood approaches to variance component
estimation and to related problems.
Journal of the American Statistical Association, 72, 320–338.

[10] Wang, T. and Merkle, E. C. (2018). merderiv: Derivative computations for
linear mixed effects models with
application to robust standard errors. Journal of Statistical Software, Code
Snippets, 87, 1–16.

[11] Satterthwaite, F. E. (1946). An approximate distribution of estimates of
variance components. Biometrics bulletin,
2, 110–114.
