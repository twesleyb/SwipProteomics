#!/usr/bin/env Rscript

## MSstatsTMT calls to lmerTest::lmer in MSstatsTMT::linearModel.functions.R

# Mixture: a TMT experiment -- the multiplex mixture of TMT-labeled samples
# Mixture:TechRepMixture: technical replicates of a mixture
# Subject:Group:Mixture: subject, a biological replicate, e.g. WT1:WT:Exp1

# Fixed effects:
# Group: the contrast of interest? 

# response: (protein) Abundance
# mixed-effects:
#   * Mixture
#   * Mixture:TechRepMixture
#   * Subject:Group:Mixture
#
#   * Run - a TMT run, e.g. Exp1.126
#   * Subject:Group - e.g. WT1.Control
# 
#   * Mixture
#   * Mixture:TechRepMixture

# model 1
lmer(Abundance ~ 1 + (1|Mixture) + (1|Mixture:TechRepMixture) + Group + (1|Subject:Group:Mixture), data = data)

# model 2
lmer(Abundance ~ 1 + (1|Run) + Group + (1|Subject:Group), data)

# model 3
lmer(Abundance ~ 1 + (1|Mixture) + (1|Mixture:TechRepMixture) + Group, data = data)

# model 4: random effect of TMT run? e.g. Exp1.126 v Exp1.126N
lmer(Abundance ~ 1 + (1|Run) + Group, data = data)

# model 5: simple comparison
lm(Abundance ~ 1 + Group, data = data)
