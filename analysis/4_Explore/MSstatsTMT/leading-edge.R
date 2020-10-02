 #!/usr/bin/env Rscript 

# 'fitted.models' contains: 
# NOTE: Outline for lmer case only!
# [1] protein name = UniProt Accession
# [*] model = list():
#     [2] fm = fitted lm or lmer object
#     [3] fixEffs --| lme4::fixef(fm) -----------|
#     [4] sigma ----| stats::sigma(fm) ----------| All from .rhoInit
#     [5] thopt ----| lme4::getME(fm, "theta") --|
#     [6] A -- rho$A = calcApvar(rho)
# [7] s2 = sigma^2
# [8] s2_df = degrees of freedom
# [9] coeff = coefficients = same as fm$coeff


# load renv
root <- "~/projects/SwipProteomics"
renv::load(root)


# load dev functions
load_fun <- function() {
	# ./dev is assumed to be in cwd
	devfun <- list.files("./dev",pattern="*.R",full.names=TRUE)
	invisible(sapply(devfun,source))
}

load_fun()

# check results against MSstatsTMT output:
load(file.path(root,"rdata","fitted.models.rda"))
# == fitted.models # len(fitted.models) == 5

# these steps reproduce the objects contained within the list 
# returned by .linear.model.fitting(data) with alot less munge
# fitted.models <- .linear.model.fitting(data)

# load msstats preprocessed protein data
load(file.path(root,"rdata","msstats_prot.rda"))
data_prot <- msstats_prot


#df <- data_prot %>% filter(Protein == prot)
#lmerTest::lmer(Abundance ~ Condition + (1 | BioReplicate), data=df)
# Fix levels

# [1] Protein Accession
proteins <- unique(data_prot$Protein)
prot <- sample(proteins,1)
message("fitting protein: ", prot)
# lmer: Abundance ~ Condition + (1 | BioReplicate) 

# [2] Linear model (lmer)
fm = fitted.models$model[[1]]$model

# given fm, we can calculate...

# [3] Fixed effects (e.g. Condition)
fixed_effects = lme4::fixef(fm)
stopifnot(all(fixed_effects == fitted.models$model[[1]]$fixEffs))

# [4] Sigma - residual standard deviation
sigma_ <- stats::sigma(fm)
stopifnot(sigma_ == fitted.models$model[[1]]$sigma)

# [5] thopt -  extract theta the random-effects parameter estimates
# parameterized as the relative Cholesky factors of each random effect
thopt = lme4::getME(fm, "theta") 

# [6] A - .calcApvar(rho) 
calcApvar <- function(fm,thopt,sigma_) {
  # alternative: to .calcApvar
  # requires .devfunTheta 
  # requires .myhess
  #dd <- .devfunTheta(fm) # generate a deviance function = devfun 
  dd = lme4::devfun2(fm)
  h = .myhess(dd, c(thopt, sigma_)) # hessian given devfun and params
  ch = try(chol(h)) # cholesky
  A = 2 * chol2inv(ch)
  # check
  eigval <- eigen(h, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigval) < sqrt(.Machine$double.eps)) {
	  warning("Asymptotic covariance matrix A is not positive!")
  }
  return(A)
}

# FIXME: not working!
A = .calcApvar(fm,thopt,sigma_)
#stopifnot(A == fitted.models$model[[1]]$A)
# functions calling update model are not working with error:
# Error in list2env(data) : first argument must be a named list

# [7] s2 = sigma^2
av = anova(fm)
#show_tests(av)
s2 = av$"Mean Sq"/av$"F value"
stopifnot(s2 == fitted.models$s2)

# [8] s2_df = degrees of freedom
s2_df == av$DenDF
stopifnot(s2_df == fitted.models$s2_df)

# [9] coeff = coefficients = same as fm$coeff
coeff = lme4::fixef(fm)
stopifnot(all(coeff == fitted.models$coeff[[1]]))
