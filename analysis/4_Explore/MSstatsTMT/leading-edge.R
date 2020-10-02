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

# load the MSstats preprocessed protein-level data
load(file.path(root,"rdata","msstats_prot.rda"))

# function to source functions in dev
load_fun <- function() {
	# NOTE: ./dev is assumed to be in cwd
	# these are core MSstats functions utilized by groupComparisons()
	fun <- list.files("./dev",pattern="*.R$",full.names=TRUE)
	invisible(sapply(fun,source))
}

load_fun()

# Workthrough 

# check results against MSstatsTMT output:
#load(file.path(root,"rdata","fitted.models.rda"))
# == fitted.models # len(fitted.models) == 5

# these steps reproduce the objects contained within the list 
# returned by .linear.model.fitting(data) with alot less munge
# fitted.models <- .linear.model.fitting(data)

# load msstats preprocessed protein data
load(file.path(root,"rdata","msstats_prot.rda"))
data_prot <- msstats_prot
#proteins <- unique(data_prot$Protein)
#prot <- sample(proteins,1)
prot <- "Q3UMB9"

data = data_prot[data_prot$Protein == prot,]

contrast.matrix = "pairwise"
moderated = FALSE
adj.method = "BH"
remove_norm_channel = TRUE
remove_empty_channel = TRUE

#groupComparisonTMT <- function(data,
#                               contrast.matrix = "pairwise",
#                               moderated = FALSE,
#                               adj.method = "BH",
#                               remove_norm_channel = TRUE,
#                               remove_empty_channel = TRUE) {

## check input data
required.info <- c(
"Protein", "BioReplicate", "Abundance", "Run", "Channel",
"Condition", "TechRepMixture", "Mixture"
)

if (!all(required.info %in% colnames(data))) {
  missedAnnotation <- which(!(required.info %in% colnames(data)))
  missedAnnotation.comb <- paste(required.info[missedAnnotation], collapse = ", ")
  if (length(missedAnnotation) == 1) {
    stop(paste(
      "Please check the required input. ** columns :",
      missedAnnotation.comb, "is missed."
    ))
  } else {
    stop(paste(
      "Please check the required input. ** columns :",
      missedAnnotation.comb, ", are missed."
    ))
  }
}

## remove 'Empty' column : It should not used for further analysis
if (remove_empty_channel & is.element("Empty", unique(data$Condition))) {
  data <- data[data$Condition != "Empty", ]
  data$Condition <- factor(data$Condition)
}

## remove 'Norm' column : It should not used for further analysis
if (remove_norm_channel & is.element("Norm", unique(data$Condition))) {
  data <- data[data$Condition != "Norm", ]
  data$Condition <- factor(data$Condition)
}

## remove the rows with NA intensities
data <- data[!is.na(data$Abundance), ]

## Inference
load_fun()

library(dplyr)
# fit: A ~ 1 + (1|Run) + C 
#result <- .proposed.model(data, moderated, contrast.matrix, adj.method)

moderated = TRUE
contrast.matrix = "pairwise"
adj.method = "BH"

  # create constrast matrix for pairwise comparisons between conditions
  conditions <- as.character(unique(data$Condition))
  contrast.matrix <- .makeContrast(conditions)
  ncomp <- nrow(contrast.matrix)

  ## fit the linear model for each protein
  fitted.models <- .linear.model.fitting(data)

  # fitted.models is a list that contains the fm and things calculated from fm
  # we are struggling to reproduce A asymptoptic var-covar matrix
  #.linear.model.fitting <- function(data) {

  data$Protein <- as.character(data$Protein) 
  proteins <- as.character(unique(data$Protein))
  num.protein <- length(proteins)
  linear.models <- list() # linear models

  s2.all <- NULL # sigma^2
  s2_df.all <- NULL # degree freedom of sigma^2
  pro.all <- NULL # testable proteins
  coeff.all <- list() # coefficients

  sub_data <- data %>% dplyr::filter(Protein == prot) ## data for swip

    ## Record the annotation information
    sub_annot <- unique(sub_data[, c(
      "Run", "Channel", "BioReplicate",
      "Condition", "Mixture", "TechRepMixture"
    )])

    ## check the experimental design
    sub_singleSubject <- .checkSingleSubject(sub_annot)
    sub_TechReplicate <- .checkTechReplicate(sub_annot)
    sub_bioMixture <- .checkMulBioMixture(sub_annot)
    sub_singleRun <- .checkSingleRun(sub_annot)

    # apply appropriate model:
    # lmer(A ~ 1 + (1|R) + C)
    fit <- fit_reduced_model_mulrun(sub_data) 
    # equivalent to :
    #fm = formula(Abundance ~ 1 + (1|Run) + Condition)
    #fx = lmerTest::lmer(fm,sub_data)

    ## estimate variance and df from linear models
    # pass

    ## Estimate the coeff from lmerTest model
    rho <- list() 
    rho <- .rhoInit(rho, fit, TRUE)  # can all be easily calculated from fm

    ## asymptotic variance-covariance matrix for theta and sigma
    rho$A <- .calcApvar(rho)  # tricky because .updateModel

   #.calcApvar <- function(rho=list(fm,thopt,sigma)) {
    fit = rho$model
    thopt = rho$thopt
    sigma = rho$sigma

    dd <- .devfunTheta(fit)

    h <- .myhess(dd, c(thopt, sigma))
    ch <- try(chol(h), silent = TRUE)
    
    A <- 2 * chol2inv(ch)

    eigval <- eigen(h, symmetric = TRUE, only.values = TRUE)$values

  if (min(eigval) < sqrt(.Machine$double.eps)) { # tolerance
    warning("Asymptotic covariance matrix A is not positive!")
  }

    av <- anova(rho$model)
    #av <- anova(fx)
    coeff <- lme4::fixef(rho$model)
    s2_df <- av$DenDF
    s2 <- av$"Mean Sq" / av$"F value"

# picking up in ______________ .proposed.model?

  ## perform empirical bayes moderation
  #if (moderated) { ## moderated t statistic
    ## Estimate the prior variance and degree freedom
    #eb_input_s2 <- fitted.models$s2[fitted.models$s2_df != 0]
    #eb_input_df <- fitted.models$s2_df[fitted.models$s2_df != 0]
    #eb_fit <- limma::squeezeVar(eb_input_s2, eb_input_df)
    #if (is.infinite(eb_fit$df.prior)) {
    #  df.prior <- 0
    #  s2.prior <- 0
    #} else {
    #  df.prior <- eb_fit$df.prior
    #  s2.prior <- eb_fit$var.prior
    #}
  
    ## ordinary t statistic
    s2.prior <- 0
    df.prior <- 0

    ## get the data for protein i
    sub_data <- data %>% dplyr::filter(Protein == prot) 

    ## record the contrast matrix for each protein
    sub.contrast.matrix <- contrast.matrix

    ## get the linear model for proteins[i]
    #fit
    #s2 <- s2.all[proteins[i]]
    #s2_df <- s2_df.all[proteins[i]]
    #coeff 

      s2.post <- (s2.prior * df.prior + s2 * s2_df) / (df.prior + s2_df)

      ## Compare one specific contrast
    j = 1
    count = 0

# groups with positive coefficients
positive.groups <- colnames(sub.contrast.matrix)[sub.contrast.matrix[j, ] > 0]

# groups with negative coefficients
negative.groups <- colnames(sub.contrast.matrix)[sub.contrast.matrix[j, ] < 0]

# make sure at least one group from each side of the contrast exist
contrast.matrix.single <- as.vector(sub.contrast.matrix[j, ])
names(contrast.matrix.single) <- colnames(sub.contrast.matrix)

cm <- .make.contrast.single(fit, contrast.matrix.single, sub_data)

# logFC
FC <- (cm %*% coeff)[, 1]

a = coeff["ConditionControl.F10"] # 1
b = coeff["ConditionControl.F4"] # -1
# FC = a - b

## variance and df
# for lm case:
#variance <- diag(t(cm) %*% summary(fit$model)$cov.unscaled %*% cm) * s2.post
#df.post <- s2_df + df.prior

## for the theta and sigma parameters
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
