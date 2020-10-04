 #!/usr/bin/env Rscript 

# description: Working through MSstatsTMT protein-level models and 
# statistical comparisons

# input:
# * preprocessed protein-level data from PDtoMSstatsTMTFormat()
contrast = "pairwise"
moderated = TRUE 
padj_method = "BH"

## key statistics:
# [*] protein name - UniProt Accession
# [*] fm - fitted lm or lmer object
# [*] sigma
# [*] theta (thopt)
# [*] A (Apvar) - asymptotic variance-covariance matrix
# [*] s2 - sigma^2
# [*] s2_df - degrees of freedom
# [*] coeff - coefficients 

## prepare environment --------------------------------------------------------

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root)

# imports
suppressPackageStartupMessages({
	library(dplyr) # all other calls should be in the form pac::fun
})


## functions ------------------------------------------------------------------

# function to source MSstatsTMT's internal functions in dev
# consider moving MSstatsTMT to source or something
load_fun <- function(funcdir="~/projects/SwipProteomics/src/MSstatsTMT") {
	# these are core MSstats functions utilized by groupComparisons()
	fun <- list.files(funcdir, pattern="*.R$", full.names=TRUE)
	n <- length(fun)
	if (n==0) { 
		warning("No R files in 'funcdir'.") 
	} else {
		invisible(sapply(fun,source))
	}
} #EOF


## load the data ---------------------------------------------------------------

# load MSstatsTMT's guts
load_fun()

# load msstats preprocessed protein data
myfile <- file.path(root,"rdata","msstats_prot.rda")
load(myfile)
# msstats_prot
data_prot <- msstats_prot


## begin protein-level modeling -------------------------------------------------

# remove rows with NA intensities
if (any(is.na(data))) { 
	data <- data[!is.na(data$Abundance), ]
	warning("Missing values were removed from input 'data'.")
}

# for intra-fraction comparisons MSstatsTMT fits the model:
#    >>>    lmerTest::lmer(Abundance ~ 1 + (1|Run) + Condition, data)


#function <- fits model, calculates rho,
#df.prior and s2.prior are 0
# should call subsequent moderation function to estimate df.post and s2.post ==
# requires all fits

fx <- formula("Abundance ~ 1 + (1|Run) + Condition")

#protein = sample(proteins,1)
protein = "Q3UMB9"

# fixme: how to catch warnings?
fm <- lmerTest::lmer(fx, data = msstats_prot %>% filter(Protein == protein))

# compute fit model statistics
rho <- getRho(fm) 

#names(rho)
# "coeff" "sigma" "thopt" "df"    "s2"    "model" "A"

# generate all potential pairwise contrasts
pairwise_contrasts <- .makeContrast(levels(data$Condition))

# the comparisons we are interested in are:
# "Mutant.F4-Control.F4"
# "Mutant.F5-Control.F5"
# ...

# pairwise contrasts sorted the levels aphabetically, so the contrasts above are
# now of the form:
# "Control.F4-Mutant.F4"
# "Control.F5-Mutant.F5"

# make inverse
new_rows <- sapply(strsplit(rownames(pairwise_contrasts),"-"), function(x) {
		   paste(x[2],x[1],sep="-") })
new_contrasts <- -1*pairwise_contrasts
rownames(new_contrasts) <- new_rows

comp <- paste(paste("Mutant",paste0("F",seq(4,10)),sep="."),
      paste("Control",paste0("F",seq(4,10)),sep="."),
      sep="-")

#head(comp)
#[1] "Mutant.F4-Control.F4" "Mutant.F5-Control.F5" "Mutant.F6-Control.F6"
#[4] "Mutant.F7-Control.F7" "Mutant.F8-Control.F8" "Mutant.F9-Control.F9"

# subset contrast_matrix, keeping the comparisons we are interested in
contrast_matrix <- new_contrasts[comp,]

n <- dim(contrast_matrix)[1] # == 7, the total number of Biological Fractions

stopifnot(n==7)

# each row of the contrast matrix is what MSstatsTMT calls a
# contrast.matrix.single
# coerce these to the form MSstatsTMT wants using a loop
subdat <- data %>% filter(Protein == protein)

contrast_list <- lapply(seq(nrow(contrast_matrix)), function(x) {
				cm <- .make.contrast.single(fm, 
							    contrast_matrix[x,],
							    subdat)
	  return(cm)
      })

# set names 
names(contrast_list) <- rownames(contrast_matrix)

# lets add contrasts_list to rho
rho$contrasts <- contrast_list

#head(contrast_list[[1]],4)
#(Intercept) ConditionControl.F4 ConditionControl.F5 ConditionControl.F6
#          0                  -1                   0                   0
# NOTE: ConditionMutant.F4 == 1

# for each of these contrasts, compute the fold change using
# the coefficients from the fit model
# basically, this is a fancy way of doing coeff(+1) - coeff(-1)
rho$FC <- lapply(contrast_list, function(cm) {
			      FC <- (cm %*% rho$coeff)[, 1]
			      return(FC)
})

#head(fc_list,1)
#$`Mutant.F4-Control.F4`
#[1] -1.306407

## The above can be done for every protein ------------------------------------

# inputs:
data <- msstats_prot
proteins <- unique(data$Protein)
# * contrasts_list
# NOTE: contrasts_list is used to generate the comparisons tested for each 
# protein. Its the format MSstatsTMT expects
# the model function to fit:
fx <- formula("Abundance ~ 1 + (1|Run) + Condition")

# there should be no missing values
# remove rows with NA intensities
if (any(is.na(data))) { 
	data <- data[!is.na(data$Abundance), ]
	warning("Missing values were removed from input 'data'.")
}

# empty list for output of loop
protein_fits <- list()

for (protein in proteins) {

# fit the model to the proteins data
subdat <- data %>% filter(Protein == protein)
fm <- lmerTest::lmer(fx, data = msstats_prot %>% filter(Protein == protein))
# compute model statistics
rho <- getRho(fm) 


# A function to compute unmoderated posterior s2.post
cp <- function(s2.prior=0,df.prior=0,s2,s2_df) {
	s2.post <- (s2.prior*df.prior + s2 * s2_df) / (df.prior + s2_df)
	return(s2.post)
}

# the statistics for all contrasts are stored as a list 'stats' in rho

rho$stats <- list()
for (comp in names(contrast_list)){
	# the formatted contrast
	message("Tested contrast: ",comp)
	cm <- contrast_list[[comp]]
	# we store the proteins statistics in a list
	stats_list <- list()
	stats_list$Comparison <- comp
	stats_list$Protein <- protein
	# compute stuff (from MSstatsTMT)
	s2.post <- cp(s2=rho$s2,s2_df=rho$s2_df)
	vss <- .vcovLThetaL(fm)
	varcor <- vss(t(cm),c(rho$thopt, rho$sigma))
	vcov <- varcor$unscaled.varcor * rho$s2
	se2 <- as.matrix(t(cm) %*% as.matrix(vcov) %*% cm)
	# calculate variance
        vcov.post <- varcor$unscaled.varcor * s2.post
        variance <- as.matrix(t(cm) %*% as.matrix(vcov.post) %*% cm) 
	# calculate degrees of freedom
	g <- .mygrad(function(x) vss(t(cm), x)$varcor, c(rho$thopt, rho$sigma))
	denom <- t(g) %*% rho$A %*% g
	df.post <- 2 * (se2)^2 / denom + df.prior
	# compute fold change and the t-statistic
	FC <- (cm %*% rho$coeff)[, 1]
	t <- FC / sqrt(variance)
	# compute the p-value
	p <- 2 * pt(-abs(t), df = df.post)
	# put it all together:
	stats_list$"Log2 Fold Change" <- FC
	stats_list$"P-value" <- p
	stats_list$"SE" <- as.numeric(sqrt(variance))
	stats_list$"DF" = as.numeric(df.post)
	# add to rho
	rho$stats[[comp]] <- stats_list
 } #EOL for each contrast

#length(rho$stats)
#[1] 7

# take a look at the stats for WASH for all intrafraction comparisons:
bind_rows(rho$stats) %>% knitr::kable()

# a function to do the above:
testContrasts <- function(rho,contrast_list) {
	rho$stats <- list()
	for (comp in names(contrast_list)){
	  # the formatted contrast
	  cm <- contrast_list[[comp]]
	  # we store the proteins statistics in a list
	  stats_list <- list()
	  stats_list$Comparison <- comp
	  stats_list$Protein <- protein
	  # compute stuff (from MSstatsTMT)
	  s2.post <- cp(s2=rho$s2,s2_df=rho$s2_df)
	  vss <- .vcovLThetaL(fm)
	  varcor <- vss(t(cm),c(rho$thopt, rho$sigma))
	  vcov <- varcor$unscaled.varcor * rho$s2
	  se2 <- as.matrix(t(cm) %*% as.matrix(vcov) %*% cm)
	  # calculate variance
          vcov.post <- varcor$unscaled.varcor * s2.post
          variance <- as.matrix(t(cm) %*% as.matrix(vcov.post) %*% cm) 
	  # calculate degrees of freedom
	  g <- .mygrad(function(x) vss(t(cm), x)$varcor, c(rho$thopt, rho$sigma))
	  denom <- t(g) %*% rho$A %*% g
	  df.post <- 2 * (se2)^2 / denom + df.prior
	  # compute fold change and the t-statistic
	  FC <- (cm %*% rho$coeff)[, 1]
	  t <- FC / sqrt(variance)
	  # compute the p-value
	  p <- 2 * pt(-abs(t), df = df.post)
	  # put it all together:
	  stats_list$"Log2 Fold Change" <- FC
	  stats_list$"P-value" <- p
	  stats_list$"SE" <- as.numeric(sqrt(variance))
	  stats_list$"DF" = as.numeric(df.post)
	  # add to rho
	  rho$stats[[comp]] <- stats_list
	} #EOL for each contrast
	# return rho updated with protein stats
	return(rho)
}

rho <- testContrasts(rho,contrast_list)

# now we have:
# do all of the the unmoderated work for a single protein:
# just need contrast_list the data and a function
# FIXME: how to capture boundary fit message?
# FIXME: need to parallelize for speed


silence <- function(x,...){
sink(tempfile())
on.exit(sink())
invisible(force(x))
}

protein_fits <- list()
for (protein in proteins) {
	subdat <- data %>% filter(Protein == protein)

	fm <- tryCatch(
		       lmerTest::lmer(fx, data = subdat)
		       warning = function(w) {
		       }

	rho <- getRho(fm) 
	rho <- testContrasts(rho,contrast_list)
	protein_fits[[protein]] <- rho
}



## moderation -----------------------------------------------------------------

## NOTE: for moderated t-statistic, all models must be fit.


# calc posterior s2
# function: generate contrsats: given data generate contrasts to be tested for
# all proteins

getContrasts <- function(data, fit_proteins){
	# do once
	pairwise_contrasts <- .makeContrast(levels(data$Condition))
	# make inverse
	#new_rows <- sapply(strsplit(rownames(pairwise_contrasts),"-"), function(x) {
	#		   paste(x[2],x[1],sep="-") })
	#new_contrasts <- -1*pairwise_contrasts
	#rownames(new_contrasts) <- new_rows
	#all_contrasts <- rbind(pairwise_contrasts,new_contrasts)
	all_contrasts <- pairwise_contrasts
	# generate contrasts to be tested
	# test all possible pairwise comparisons:
	# NOTE: order of levels determines sign of fold change: default is
	# alphabetically sorted

	contrast_list <- list()
	for (contrast in rownames(pairwise_contrasts)){
		contrast.matrix.single <- pairwise_contrasts[contrast,]
		# should not depend upon fm! names of coeff for each protein is
		# the same bc each protein is fit with the same model!
	  	contrast_list[[contrast]] <- .make.contrast.single(fm, contrast.matrix.single, data)
	} else if (contrast %in% rownames(all_contrasts)) {
		contrast_list[[contrast]] <- .make.contrast.single(fm,all_contrasts[contrast,],data)
	} else {
		stop("Input 'contrast' not valid.")
	} #EIS

# given all fits and contrasts, for every fit.contrast do:
# fold change
# variance
# t
# p

# now we can compute fold change
cm <- contrast_list[[1]]
fold_change <- lapply(contrast_list, function(cm) {
			      FC <- (cm %*% rho$coeff)[, 1]
			      return(FC)
})


#if (inherits(fit$model, "lm")) {
#  variance <- diag(t(cm) %*% summary(fit$model)$cov.unscaled %*% cm) * s2.post
#  df.post <- s2_df + df.prior
#} else {

# for lmer:
vss <- .vcovLThetaL(fm)
varcor <- vss(t(cm), c(rho$thopt, rho$sigma)) ## for the theta and sigma parameters
vcov <- varcor$unscaled.varcor * rho$s2
se2 <- as.matrix(t(cm) %*% as.matrix(vcov) %*% cm)
# what should be added to rho?
#varcor

## calculate df
g <- .mygrad(function(x) vss(t(cm), x)$varcor, c(rho$thopt, rho$sigma))

rho$s2.prior = 0
rho$df.prior = 0

rho$s2.post <- (rho$s2.prior * rho$df.prior + rho$s2 * s2_df) / (df.prior + s2_df)

denom <- try(t(g) %*% rho$A %*% g, silent = TRUE)

if (inherits(denom, "try-error")) {
	df.post <- s2_df + df.prior
} else {
	df.post <- 2 * (se2)^2 / denom + df.prior
}

  ## calculate the t statistic
  t <- FC / sqrt(variance)

  ## calculate p-value
  p <- 2 * pt(-abs(t), df = df.post)



## NOTES about fx - formula
# Where 'Abundance' is the response, 'Run' is a mixed-effect, and 'Condition'
# is a fixed-effect. The model is fit by a call to lmerTest::lmer().
# lmerTest is a package that wraps 'lme4'.

# 'Run' cooresponds to a single multi-plex TMT experiment.
# This is the affect of experimental Batch.

# 'Condition' for intrafraction comparisons is Treatment::BioFraction 
# (e.g. Control.F7 and Mutant.F9)


# create constrast matrix for pairwise comparisons between conditions
## all of this is done to 
#conditions <- as.character(unique(data$Condition))
#all_contrasts <- .makeContrast(conditions)

# fixme: we are not really interested in all comparisons!
#ncomp <- nrow(all_contrasts)


# fit a linear model for each protein
# this loops to fit models for all proteins given the experimental design
# automatically determined from the formatting of the MSstatsTMT input
#fitted.models <- .linear.model.fitting(data)
# fixme: what happens if lm case is passed to lmer?
# otherwise its as simlple as:



# fitted.models is a list that contains the fm and things calculated from fm
# for every protein


## .linear.model.fitting ------------------------------------------------------

#.linear.model.fitting <- function(data) {

#proteins <- unique(as.character(data$Protein))
#num.protein <- length(proteins)

# objects to store the output from loop:
#s2.all <- NULL # sigma^2
#s2_df.all <- NULL # degree freedom of sigma^2
#pro.all <- NULL # testable proteins
#coeff.all <- list() # coefficients
#linear.models <- list() 

#sub_data <- data %>% dplyr::filter(Protein == prot) ## data for swip

## Record the annotation information
#sub_annot <- unique(sub_data[, c(
#"Run", "Channel", "BioReplicate",
#"Condition", "Mixture", "TechRepMixture"
#)])

## check the experimental design
#sub_singleSubject <- .checkSingleSubject(sub_annot)
#sub_TechReplicate <- .checkTechReplicate(sub_annot)
#sub_bioMixture <- .checkMulBioMixture(sub_annot)
#sub_singleRun <- .checkSingleRun(sub_annot)

# apply appropriate model:
# lmer(A ~ 1 + (1|R) + C)
# FIXME: why not just pass a formula?
# and fit the model?
#fit <- fit_reduced_model_mulrun(sub_data) 
# equivalent to :
#fm = formula(Abundance ~ 1 + (1|Run) + Condition)
#fx = lmerTest::lmer(fm,sub_data)

## estimate variance and df from linear models
# pass

# rho <- function() {

# populate list rho with (fixEffs (coeff), sigma, and thopt(theta))
# these things can be done with simple calls to lme4

########################
# FIXME: unmask what rho does
#rho <- .rhoInit(rho, fit, TRUE)  # can all be easily calculated from fm

## asymptotic variance-covariance matrix for theta and sigma
#rho$A <- .calcApvar(rho)  # tricky because .updateModel
# FIXME: replace with a function that takes fm
# can replace all with single function instead of two
# given fit, get: A, sigma, coeff, theta

# add to step above:
av <- anova(fit)
#av <- anova(fx)
#coeff <- lme4::fixef(rho$model) # done in step above
s2_df <- av$DenDF
s2 <- av$"Mean Sq" / av$"F value"

##-----------------------------------------------------------------------------
## perform empirical bayes moderation
  #if (moderated) { ## moderated t statistic
    ## NOTE: this step takes into acount all fits with nobs > 1!
    ## Estimate the prior variance and degree freedom
    # vector of s2 for all fits with > 1 obs
    # vector of dof (npeptides) for all fits with > 1 obs
    if (moderated) {
    eb_input_s2 <- fitted.models$s2[fitted.models$s2_df != 0]
    eb_input_df <- fitted.models$s2_df[fitted.models$s2_df != 0]
    eb_fit <- limma::squeezeVar(eb_input_s2, eb_input_df)
    if (is.infinite(eb_fit$df.prior)) {
	    df.prior <- 0
	    s2.prior <- 0
    }
	    df.prior <- eb_fit$df.prior
	    s2.prior <- eb_fit$var.prior
    }
  
## else ordinary t statistic
s2.prior <- 0
df.prior <- 0

# now go back to protein level stuff: assessing contrasts between groups

## get the data for protein i
sub_data <- data %>% dplyr::filter(Protein == prot) 

# create constrast matrix for pairwise comparisons between conditions
## all of this is done to 
conditions <- as.character(unique(data$Condition))

all_contrasts <- .makeContrast(conditions)

# fixme: we are not really interested in all comparisons!
ncomp <- nrow(all_contrasts)

## record the contrast matrix for each protein
sub.contrast.matrix <- contrast.matrix

# can we subset here?
all_contrasts

## get the linear model for proteins[i]
#fit
#s2 <- s2.all[proteins[i]]
#s2_df <- s2_df.all[proteins[i]]
#coeff 

# calculate posterior:
s2.post <- (s2.prior * df.prior + s2 * s2_df) / (df.prior + s2_df)

## example: one contrast
j = 1
count = 0

# i.e. for a given contrast or row of the contrast.matrix

# groups with positive coefficients
positive.groups <- colnames(sub.contrast.matrix)[sub.contrast.matrix[j, ] > 0]

# groups with negative coefficients
negative.groups <- colnames(sub.contrast.matrix)[sub.contrast.matrix[j, ] < 0]

# make sure at least one group from each side of the contrast exist
contrast.matrix.single <- as.vector(sub.contrast.matrix[j, ])
names(contrast.matrix.single) <- colnames(sub.contrast.matrix)

cm <- .make.contrast.single(fit, contrast.matrix.single, sub_data)

# munge to reverse
cm[cm == -1] <- 0
cm[names(cm) == "ConditionControl.F4"] <- -1

# negative index all else 0
#
# why do we need fit here?
# used to get names of coeff

## return contrast matrix

# logFC
FC <- (cm %*% coeff)[, 1]

# alot simpler:
a = coeff["ConditionControl.F10"] # 1
b = coeff["ConditionControl.F4"] # -1
# FC = a - b

## variance and df
# for lm case:
#variance <- diag(t(cm) %*% summary(fit$model)$cov.unscaled %*% cm) * s2.post
#df.post <- s2_df + df.prior

## for the mixed effects case:
thopt = rho$thopt
sigma = rho$sigma

# do something with theta and sigma parameters
fit = rho$model
vss <- .vcovLThetaL(fit) # returns a function, includes call to .updateModel
varcor <- vss(t(cm), c(thopt, sigma)) 

# NOTE: this is not equivalent to:
#vss <- update(fit,devFunOnly=TRUE) # get deviance function
#varcor <- vss(theta)

vcov <- varcor$unscaled.varcor * s2

se2 <- as.matrix(t(cm) %*% as.matrix(vcov) %*% cm)

## calculate variance
vcov.post <- varcor$unscaled.varcor * s2.post
variance <- as.matrix(t(cm) %*% as.matrix(vcov.post) %*% cm)

## calculate df
# now this is crazy
g <- .mygrad(function(x) vss(t(cm), x)$varcor, c(thopt, sigma))

# HERE's A!
denom <- try(t(g) %*% A %*% g, silent = TRUE) 

# when does error arrise?
if (inherits(denom, "try-error")) {
        df.post <- s2_df + df.prior 
} else {
	df.post <- 2 * (se2)^2 / denom + df.prior
}

## calculate the t statistic
t <- FC / sqrt(variance) # variance is from posterior distribution

# calculate p-value
p <- 2 * pt(-abs(t), df = df.post)

## key statistics:
#log2FC = FC
#SE = sqrt(variance)
#DF = df.post
# if (s2_df == 0) "SingleMeasurePerCondition"

# p 
# [1,] 2.400825e-06 

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
  res <- as.data.frame(res[seq_len(count), ])
  res$Protein <- as.factor(res$Protein)
  res$log2FC <- as.numeric(as.character(res$log2FC))
  res$pvalue <- as.numeric(as.character(res$pvalue))
  res$adjusted.pvalue <- NA
  comps <- unique(res$Comparison)

  ## Adjust multiple tests for each comparison
  for (i in seq_along(comps)) {
    res[res$Comparison == comps[i], "adjusted.pvalue"] <- p.adjust(res[res$Comparison == comps[i], "pvalue"], adj.method)
  }

  res <- res[, c(
    "Protein",
    "Comparison",
    "log2FC",
    "SE",
    "DF",
    "pvalue",
    "adjusted.pvalue",
    "issue"
  )]
  return(res)
}















# fit the model:
fx <- formula(Abundance ~ 1 + (1|Run) + Condition)
fm <- lmerTest::lmer(fx, data) 

# calc fixed effects
fixed_effects <- lme4::fixef(fm)

# calc Sigma
sigma_ <- stats::sigma(fm)

# calc theta
theta_ <- lme4::getME(fm, "theta")

# calc s2  ~ sigma_^2
av = anova(fm)
s2 = av$"Mean Sq"/av$"F value" 

# calc s2_df
s2_df = av$DenDF

# calc coeff
coeff = lme4::fixef(fm) # == fixed_effects

# all we are missing is A!
A = .calcApvar(rho)

dd = .devfunTheta(rho$model) # this is the problem child
#is.function(dd) == TRUE

.devfunTheta(fm)

# load previously fit models
#load(file.path(root,"rdata","fitted.models.rda"))
#old_fm = fitted.models$model[[1]]$model
#.devfunTheta(old_fm)

# we already have thopt (theta) and sigma!
h = .myhess(dd,c(rho$thopt,sigma=rho$sigma))

# now do this:

# priors
if (moderated) { 
	# need all fits

    eb_input_s2 <- fitted.models$s2[fitted.models$s2_df != 0]

    eb_input_df <- fitted.models$s2_df[fitted.models$s2_df != 0]
    # perform empirical bayes moderation
    eb_fit <- limma::squeezeVar(eb_input_s2, eb_input_df)
    if (is.infinite(eb_fit$df.prior)) {
      df.prior <- 0
      s2.prior <- 0
    } else {
      df.prior <- eb_fit$df.prior
      s2.prior <- eb_fit$var.prior
    }
  # ordinary t-statistic
  } else { 
    s2.prior <- 0
    df.prior <- 0
  }


# if model is fittable, then do this
s2.prior = df.prior = 0 # non-moderated t-statistic
s2.post <- (s2.prior * df.prior + s2 * s2_df) / (df.prior + s2_df)

          ## logFC 
          # cm = numeric
          # (Intercept) ConditionMutant
          # 0           NA

cm = setNames(c(0,NA),nm=c("(Intercept)", "ConditionControl.F4"))

coeff = lmer::fixef(fm)
          FC <- (cm %*% coeff)[, 1]



















#####################################################
#data = data
moderated
contrast.matrix
adj.method

#' @import statmod
#' @importFrom limma squeezeVar
#' @importFrom dplyr %>% group_by filter
#' @importFrom stats aggregate anova coef lm median medpolish 
#' @import from stats model.matrix p.adjust pt t.test xtabs
#' @keywords internal

groups <- as.character(unique(data$Condition)) # groups
if (length(groups) < 2) {
  stop("Please check the 'Condition' column in 'annotation' data.frame",
       "There must be at least two conditions!")
}

## contrast matrix can be matrix or character vector
if (is.matrix(contrast.matrix)) {
  # comparison come from contrast.matrix
  if (!all(colnames(contrast.matrix) %in% groups)) {
    stop("Please check input 'contrast.matrix'. ",
  "Column names of 'contrast.matrix' must match Conditions!")
    }
  } else {
    # create constrast matrix for pairwise comparison
    contrast.matrix <- .makeContrast(groups)
  }

  # The number of comparisons
  ncomp <- nrow(contrast.matrix)

  ## fit the linear model for each protein
  fitted.models <- .linear.model.fitting(data)

  ## estimate the prior variance and degree freedom
  # moderated t-statistic 
  if (moderated) { 
    eb_input_s2 <- fitted.models$s2[fitted.models$s2_df != 0]
    eb_input_df <- fitted.models$s2_df[fitted.models$s2_df != 0]
    # perform empirical bayes moderation
    eb_fit <- limma::squeezeVar(eb_input_s2, eb_input_df)
    if (is.infinite(eb_fit$df.prior)) {
      df.prior <- 0
      s2.prior <- 0
    } else {
      df.prior <- eb_fit$df.prior
      s2.prior <- eb_fit$var.prior
    }
  # ordinary t-statistic
  } else { 
    s2.prior <- 0
    df.prior <- 0
  }

  # extract the linear model fitting results
  proteins <- fitted.models$protein # proteins
  s2.all <- fitted.models$s2 # group variance
  s2_df.all <- fitted.models$s2_df # degree freedom of s2
  lms <- fitted.models$model # linear models
  coeff.all <- fitted.models$coeff # coefficients

  num.protein <- length(proteins)

  ## store the inference results
  res <- as.data.frame(matrix(rep(NA, 7 * num.protein * ncomp), ncol = 7)) 
  colnames(res) <- c("Protein", "Comparison", "log2FC", 
		     "pvalue", "SE", "DF", "issue")
  data$Condition <- as.factor(data$Condition) # make sure group is factor
  data$Run <- as.factor(data$Run)

  # check the number of MS runs in the data
  nrun <- length(unique(data$Run)) 
  count <- 0
  for (i in seq_along(proteins)) {

    ## get the data for protein i
    sub_data <- data %>% dplyr::filter(Protein == proteins[i]) 

    ## record the contrast matrix for each protein
    sub.contrast.matrix <- contrast.matrix
    sub_groups <- as.character(unique(sub_data[, c("Condition")]))
    sub_groups <- sort(sub_groups) # sort the groups based on alphabetic order

    ## get the linear model for proteins[i]
    fit <- lms[[proteins[i]]]
    s2 <- s2.all[proteins[i]]
    s2_df <- s2_df.all[proteins[i]]
    coeff <- coeff.all[[proteins[i]]]

    if (!is.character(fit)) { ## check the model is fittable

      s2.post <- (s2.prior * df.prior + s2 * s2_df) / (df.prior + s2_df)

      ## Compare one specific contrast
      # perform testing for required contrasts
      for (j in seq_len(nrow(sub.contrast.matrix))) {
        count <- count + 1
        res[count, "Protein"] <- proteins[i] ## protein names
        res[count, "Comparison"] <- row.names(sub.contrast.matrix)[j] ## comparison

        # groups with positive coefficients
        positive.groups <- colnames(sub.contrast.matrix)[sub.contrast.matrix[j, ] > 0]
        # groups with negative coefficients
        negative.groups <- colnames(sub.contrast.matrix)[sub.contrast.matrix[j, ] < 0]
        # make sure at least one group from each side of the contrast exist
        if (any(positive.groups %in% sub_groups) &
          any(negative.groups %in% sub_groups)) {
          contrast.matrix.single <- as.vector(sub.contrast.matrix[j, ])
          names(contrast.matrix.single) <- colnames(sub.contrast.matrix)

          cm <- .make.contrast.single(fit$model, contrast.matrix.single, sub_data)

          ## logFC
          FC <- (cm %*% coeff)[, 1]

          ## variance and df
          if (inherits(fit$model, "lm")) {
            variance <- diag(t(cm) %*% summary(fit$model)$cov.unscaled %*% cm) * s2.post
            df.post <- s2_df + df.prior
          } else {
            vss <- .vcovLThetaL(fit$model)
            varcor <- vss(t(cm), c(fit$thopt, fit$sigma)) ## for the theta and sigma parameters
            vcov <- varcor$unscaled.varcor * s2
            se2 <- as.matrix(t(cm) %*% as.matrix(vcov) %*% cm)

            ## calculate variance
            vcov.post <- varcor$unscaled.varcor * s2.post
            variance <- as.matrix(t(cm) %*% as.matrix(vcov.post) %*% cm)

            ## calculate df
            g <- .mygrad(function(x) vss(t(cm), x)$varcor, c(fit$thopt,fit$sigma))
            denom <- try(t(g) %*% fit$A %*% g, silent = TRUE)
            if (inherits(denom, "try-error")) {
              df.post <- s2_df + df.prior
            } else {
              df.post <- 2 * (se2)^2 / denom + df.prior
            }
          }

          ## calculate the t statistic
          t <- FC / sqrt(variance)

          ## calculate p-value
          p <- 2 * pt(-abs(t), df = df.post)
          res[count, "pvalue"] <- p

          ## save testing results
          res[count, "log2FC"] <- FC
          res[count, "SE"] <- sqrt(variance)
          res[count, "DF"] <- df.post

          if (s2_df == 0) {
            res[count, "issue"] <- "SingleMeasurePerCondition"
          } else {
            res[count, "issue"] <- NA
          }
        } else {
          # at least one condition is missing
          out <- .issue.checking(
            data = sub_data,
            contrast.matrix = sub.contrast.matrix[j, ]
          )

          res[count, "log2FC"] <- out$logFC
          res[count, "pvalue"] <- NA
          res[count, "SE"] <- NA
          res[count, "DF"] <- NA
          res[count, "issue"] <- out$issue
        }
      } # for constrast matrix
    } else {
      # very few measurements so that the model is unfittable
      for (j in 1:nrow(sub.contrast.matrix)) {
        count <- count + 1
        res[count, "Protein"] <- proteins[i] ## protein names
        res[count, "Comparison"] <- row.names(sub.contrast.matrix)[j] ## comparison

        out <- .issue.checking(
          data = sub_data,
          contrast.matrix = sub.contrast.matrix[j, ]
        )

        res[count, "log2FC"] <- out$logFC
        res[count, "pvalue"] <- NA
        res[count, "SE"] <- NA
        res[count, "DF"] <- NA
        res[count, "issue"] <- out$issue
      } # end loop for comparison
    } # if the linear model is fittable
  } # for each protein

  res <- as.data.frame(res[seq_len(count), ])
  res$Protein <- as.factor(res$Protein)
  res$log2FC <- as.numeric(as.character(res$log2FC))
  res$pvalue <- as.numeric(as.character(res$pvalue))
  res$adjusted.pvalue <- NA
  comps <- unique(res$Comparison)

  ## Adjust multiple tests for each comparison
  for (i in seq_along(comps)) {
    res[res$Comparison == comps[i], "adjusted.pvalue"] <- p.adjust(res[res$Comparison == comps[i], "pvalue"], adj.method)
  }

  res <- res[, c(
    "Protein",
    "Comparison",
    "log2FC",
    "SE",
    "DF",
    "pvalue",
    "adjusted.pvalue",
    "issue"
  )]
  return(res)
} #EOF .proposed.model

###############################################################################
## Resume .proposed.model
  ### check column name in order to use groupComparisonPlot from MSstats
  colnames(result)[colnames(result) == "Comparison"] <- "Label"
  colnames(result)[colnames(result) == "adjusted.pvalue"] <- "adj.pvalue"

  return(result)
}





























###############################################################################
# THIS IS THE MODEL FOR ADJUSTED COMPARISONS
# lmer: Abundance ~ Condition + (1 | BioReplicate) 
proteins <- unique(data_prot$Protein)
prot <- sample(proteins,1)
df <- data_prot %>% filter(Protein == prot)

# munge BioReplicate should be Condition.Fraction
df$BioReplicate <- gsub("\\.R[1,2,3]","",as.character(df$BioReplicate))

# munge Condition
df$Condition <- gsub("\\.F[0-9]{1,2}","",as.character(df$BioReplicate))

# munge BioFraction == subcellular fraction
df$BioFraction <- sapply(strsplit(df$BioReplicate,"\\."),"[",2)

fm0 = lmerTest::lmer(Abundance ~ Condition + (1|BioFraction), data=df) # how to understand if model makes sense?
fm1 = lmerTest::lmer(Abundance ~ Condition + (1|BioReplicate:BioFraction), data=df) # how to understand if model makes sense?
anova(fm0,fm1) # p=1, use fm0?


# Fix levels

# [1] Protein Accession
message("fitting protein: ", prot)
# lmer: Abundance ~ Condition + (1 | BioReplicate) 

# [2] Linear model (lmer)
# extract from list
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

	# not working.... i think you need to be in the same env that the model
	# was generated in for some reason...

  # alternative: to .calcApvar
  # requires .devfunTheta 
  # requires .myhess
  dd <- .devfunTheta(fm) # generate a deviance function = devfun 
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
A = calcApvar(fm,thopt,sigma_) # list2env(data) first arg must be a named list

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
