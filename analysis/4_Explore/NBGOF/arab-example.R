#!/usr/bin/env Rscript

# adapted from arab.Rd documentation

# imports
renv::load()
library(NBGOF)

# load the data
data(arab) # 26222 by 6 matrix

## simulation set-up:
seed = 539
sim = 999
conf.env = 0.95
nc = detectCores() - 1
set.seed(seed)

## subset read counts to exclude all zero rows within each treatment:
good.arab = arab[(rowSums(arab[ ,1:3]) != 0 & rowSums(arab[ ,4:6]) != 0), ]
m = 500    # number of genes retained

# consider a single group with 3 replicates
samp.idx = sample(1:dim(good.arab)[1], m)
arab.sub = good.arab[samp.idx,1:3]
lib.sizes = colSums(arab.sub)
y = arab.sub

## model matrix for arab.sub:
x = as.matrix(rep(1,3)) # a single group

## GOF tests for different dispersion models:
fnbp.arab = nb.gof.m(counts=y, x=x, sim=sim, model="NBP", seed=1, ncores=nc)

summary(fnbp.arab)

# ERROR:
# Error in if (method == "MAPL") { : argument is of length zero
# Calls: nb.gof.m -> model.nbp.m -> estimate.dispersion
# Execution halted
