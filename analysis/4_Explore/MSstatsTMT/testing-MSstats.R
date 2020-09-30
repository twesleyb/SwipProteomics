#!/usr/bin/env Rscript

# title: SwipProteomics
# description: analysis of Swip TMT spatial proteomics data with MSstatsTMT
# author: Tyler W Bradshaw <twesleyb10@gmail.com>

## Input ----------------------------------------------------------------------
# data in root/rdata:
# * data_prot.rda

## Prepare the working environment --------------------------------------------

root <- "~/projects/SwipProteomics"
datadir <- file.path(root,"data")
rdatdir <- file.path(root,"rdata") # temporary/large files not tracked by git


## Prepare the R environment ---------------------------------------------------

# load renv
renv::load(root,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
	library(MSstatsTMT) # twesleyb/MSstatsTMT
})

# load functions in root/R
suppressPackageStartupMessages({ devtools::load_all() })


## Perform protein-level statistical tests --------------------------------------
# Tests for significant changes in protein abundance across conditions based on
# a family of linear mixed-effects models in TMT experiment.

data(msstats_prot); data_prot <- msstats_prot

# test for all the possible pairs of conditions	
#all_results <- groupComparisonTMT(data_prot)	

# groupComparisonTMT
data = data_prot
contrast.matrix = 'pairwise'
moderated = FALSE
adj.method = 'BH'

# 1 ---------------------------------------------------------------------------

groupComparisonTMT <- function(data,
                               contrast.matrix = "pairwise",
                               moderated = FALSE,
                               adj.method = "BH") {

  # check, the  input data should contain the following columns:
  required_cols <- c("Protein", "BioReplicate", "Abundance", "Run", "Channel",
		     "Condition", "TechRepMixture", "Mixture")
  if (!all(required_cols %in% colnames(data))) {
  	stop("Check input 'data'. The following columns are required:", 
  	     paste(required_cols,collapse=", "))
  }

  # remove 'Empty' column : It should not used for further analysis
  if ("Empty" %in% levels(data$Condition)) {
  	stop("Column 'Condition' cannot contain the level 'Empty'.")
  }

  # remove 'Norm' column : It should not used for further analysis
  if ("Norm" %in% levels(data$Condition)) {
  	stop("Column 'Condition' cannot contain the level 'Norm'.")
  }

  # remove the rows with NA intensities
  if (any(is.na(data$Abundance))) {
	  warning("Input 'data' cannot contain missing values (NA). ",
		  "Missing values will be removed.")
	  data <- data[!is.na(data$Abundance), ]
  }

  # protein-level inference see [2]
  result <- .proposed.model(data, moderated, contrast.matrix, adj.method)

  # check column name in order to use groupComparisonPlot from MSstats
  colnames(result)[colnames(result) == "Comparison"] <- "Label"
  colnames(result)[colnames(result) == "adjusted.pvalue"] <- "adj.pvalue"

  # return results
  return(result)
}

# 2 ---------------------------------------------------------------------------

.proposed.model <- function(data,
                            moderated = TRUE,
                            contrast.matrix = "pairwise",
                            adj.method = "BH") {

  groups <- as.character(unique(data$Condition))
  if (length(groups) < 2) {
    stop("Input 'data' has fewer than two 'Conditions'. ",
	 "Please check annotation file input.")
  }

  # check input contrast.matrix
  if (is.matrix(contrast.matrix)) {
    if (!all(colnames(contrast.matrix) %in% groups)) {
      stop("Please check input 'contrast.matrix'. ",
	   "Column names should match levels in 'data$Condition'.")
    }
  } else if (is.character(contrast.matrix)) {
    # create constrast matrix for pairwise comparisons see [3]
    contrast.matrix <- .makeContrast(groups)
  } else {
	  stop("Problem with input 'contrast.matrix'. ",
	       "Expected matrix or character vector.")
  }

  ncomp <- nrow(contrast.matrix)

  # fit the linear model for each protein, see [4]
  fitted.models <- .linear.model.fitting(data)

  ## perform empirical bayes moderation
  if (moderated) { ## moderated t statistic
    ## Estimate the prior variance and degree freedom

    eb_input_s2 <- fitted.models$s2[fitted.models$s2_df != 0]
    eb_input_df <- fitted.models$s2_df[fitted.models$s2_df != 0]

    eb_fit <- limma::squeezeVar(eb_input_s2, eb_input_df)

    if (is.infinite(eb_fit$df.prior)) {
      df.prior <- 0
      s2.prior <- 0
    } else {
      df.prior <- eb_fit$df.prior
      s2.prior <- eb_fit$var.prior
    }
  } else { ## ordinary t statistic
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
  res <- as.data.frame(matrix(rep(NA, 7 * num.protein * ncomp), ncol = 7)) ## store the inference results
  colnames(res) <- c("Protein", "Comparison", "log2FC", "pvalue", "SE", "DF", "issue")
  data$Group <- as.factor(data$Group) # make sure group is factor
  data$Run <- as.factor(data$Run)
  nrun <- length(unique(data$Run)) # check the number of MS runs in the data
  count <- 0
  for (i in seq_along(proteins)) {
    message(paste("Testing for Protein :", proteins[i], "(", i, " of ", num.protein, ")"))

    ## get the data for protein i
    sub_data <- data %>% dplyr::filter(Protein == proteins[i]) ## data for protein i
    ## record the contrast matrix for each protein
    sub.contrast.matrix <- contrast.matrix

    sub_groups <- as.character(unique(sub_data[, c("Group")]))
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
            g <- .mygrad(function(x) vss(t(cm), x)$varcor, c(fit$thopt, fit$sigma))
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
}

#------------------------------------------------------------------------------


# 3 ---------------------------------------------------------------------------

.makeContrast <- function(groups) {
  ncomp <- length(groups) * (length(groups) - 1) / 2 # Number of comparison
  contrast.matrix <- matrix(rep(0, length(groups) * ncomp),
    ncol = length(groups)
  )
  colnames(contrast.matrix) <- groups
  count <- 0
  contrast.matrix.rownames <- NULL
  for (j in seq_len(length(groups) - 1)) {
    for (k in (j + 1):length(groups)) {
      count <- count + 1
      # save row name
      contrast.matrix.rownames <- c(
        contrast.matrix.rownames,
        paste(groups[j], groups[k], sep = "-")
      )
      # set constrast value
      contrast.matrix[count, groups[j]] <- 1
      contrast.matrix[count, groups[k]] <- -1
    }
  }
  rownames(contrast.matrix) <- contrast.matrix.rownames
  return(contrast.matrix)
}

#------------------------------------------------------------------------------


# 4 ---------------------------------------------------------------------------
## FIT MIXED EFFECT LINEAR MODELS APPROPRIATE FOR THE EXPERIMENTAL DESIGN

.linear.model.fitting <- function(data) {

  # check input
  # FIXME: move to singular input checking function!
  data$Protein <- as.character(data$Protein) 
  proteins <- as.character(unique(data$Protein))
  num.protein <- length(proteins)

  # outputs from LOOP:
  s2.all <- NULL          # sigma^2
  pro.all <- NULL         # testable proteins
  s2_df.all <- NULL       # degree freedom of sigma^2
  coeff.all <- list()     # coefficients
  linear.models <- list() # linear models

  # LOOP to do inference for each protein individually
  # FIXME: parallelize
  for (prot in proteins) {

    sub_data <- data %>% dplyr::filter(Protein == prot) ## data for protein_i

    # collect the protein's metadata
    sub_annot <- unique(sub_data[, c(
      "Run", "Channel", "BioReplicate",
      "Condition", "Mixture", "TechRepMixture"
    )])

    # check the experimental design, see [5]
    sub_singleSubject <- .checkSingleSubject(sub_annot)
    sub_TechReplicate <- .checkTechReplicate(sub_annot)
    sub_BioMixture <- .checkMultiBioMixture(sub_annot)
    sub_singleRun <- .checkSingleRun(sub_annot)

    ## NASTY - do not enter
    ## got condition f one way anova??
    if (sub_singleSubject) { # no biological variation within each condition and mixture
      if (sub_TechReplicate & sub_BioMixture) { # multiple mixtures and technical replicates
        # fit the full model with mixture and techrep effects for spiked-in data
	message("\ncondition: a")
        fit <- fit_full_model_spikedin(sub_data) # see [6.1]
        if (is.null(fit)) { # if the full model is not applicable,
          # then fit the reduced model with only run effect
	  message("\ncondition: b")
          fit <- fit_reduced_model_mulrun(sub_data) # see [6.2]
        }
        if (is.null(fit)) { # if the second model is not applicable,
          # then fit one-way anova model
	  message("\ncondition: c")
          fit <- fit_reduced_model_onerun(sub_data) # see [6.3]
        }
      } else { # sub_singleSubject == TRUE 
        if (sub_TechReplicate | sub_BioMixture) { # multiple TechReplicates or BioMixtures
          # fit the reduced model with only run effect
	  message("\ncondition: d")
          fit <- fit_reduced_model_mulrun(sub_data) # see [6.2]
          if (is.null(fit)) { # if the second model is not applicable,
            # then fit one-way anova model
	    message("\ncondition: e")
            fit <- fit_reduced_model_onerun(sub_data) # see [6.3] # not null
          }
        } else { # single run case
          # fit one-way anova model
	  message("\ncondition: f")
          fit <- fit_reduced_model_onerun(sub_data) # see [6.3]
        }
      }
    } else { # biological variation exists within each condition and mixture
      if (sub_bioMixture) { # multiple biological mixtures
        if (sub_TechReplicate) { # multiple technical replicate MS runs
          # fit the full model with mixture, techrep, subject effects
	  message("\ncondition: g")
          fit <- fit_full_model(sub_data) # see [6.4]
          if (is.null(fit)) { # full model is not applicable
            # fit the reduced model with run and subject effects
	    message("\ncondition: h")
            fit <- fit_reduced_model_techrep(sub_data) # see [6.5]
          }
          if (is.null(fit)) { # if the second model is not applicable
            # then fit one-way anova model
	    message("\ncondition: i") 
            fit <- fit_reduced_model_onerun(sub_data) # see [6.3]
          }
        } else { # single technical replicate MS run
          # fit the reduced model with only run effect
	  message("\ncondition: j") 
          fit <- fit_reduced_model_mulrun(sub_data) # see [6.2]
          if (is.null(fit)) { # second model is not applicable
            # fit one-way anova model
	    message("\ncondition: k") 
            fit <- fit_reduced_model_onerun(sub_data) # see [6.3]
          }
        }
      } else { # single biological mixture
        if (sub_TechReplicate) { # multiple technical replicate MS runs
          # fit the reduced model with run and subject effects
	  message("\ncondition: l") 
          fit <- fit_reduced_model_techrep(sub_data) # see [6.5]
          if (is.null(fit)) { # second model is not applicable
            # fit one-way anova model
	    message("\ncondition: m") 
            fit <- fit_reduced_model_onerun(sub_data) # see [6.3]
          }
        } else { # single run
          # fit one-way anova model
	  message("\ncondition: n") 
          fit <- fit_reduced_model_onerun(sub_data) # see [6.3]
        } # single technical replicate MS run
      } # single biological mixture
    } # biological variation

    ## estimate variance and df from linear models
    if (!is.null(fit)) { # the model is fittable
      if (inherits(fit, "lm")) { # single run case
        ## Estimate the coeff from fixed model
        av <- anova(fit)
        coeff <- coef(fit)

        s2_df <- av["Residuals", "Df"]

        if (s2_df == 0) {
          s2 <- 0
        } else {
          # use error variance for testing
          s2 <- av["Residuals", "Mean Sq"]
        }

        linear.models[[proteins[i]]] <- list(model = fit)
      } else {
        ## Estimate the coeff from lmerTest model
        rho <- list() ## environment containing info about model
        rho <- .rhoInit(rho, fit, TRUE) ## save lmer outcome in rho envir variable
        rho$A <- .calcApvar(rho) ## asymptotic variance-covariance matrix for theta and sigma

        av <- anova(rho$model)
        coeff <- lme4::fixef(rho$model)
        s2_df <- av$DenDF
        s2 <- av$"Mean Sq" / av$"F value"

        linear.models[[proteins[i]]] <- rho
      }

      pro.all <- c(pro.all, proteins[i])
      s2.all <- c(s2.all, s2)
      s2_df.all <- c(s2_df.all, s2_df)
      coeff.all[[proteins[i]]] <- coeff
    } else { # the model is not fittble
      # message(proteins[i], " is untestable due to no enough measurements.")
      linear.models[[proteins[i]]] <- "unfittable"
      pro.all <- c(pro.all, proteins[i])
      s2.all <- c(s2.all, NA)
      s2_df.all <- c(s2_df.all, NA)
      coeff.all[[proteins[i]]] <- NA
    }
  } # for each protein
  names(s2.all) <- proteins
  names(s2_df.all) <- proteins

  return(list(
    protein = pro.all,
    model = linear.models,
    s2 = s2.all,
    s2_df = s2_df.all,
    coeff = coeff.all
  ))
}

#------------------------------------------------------------------------------


# 5 ----------------------------------------------------------------------------
# functions to check the experimental design

## check single subject within each condition in each mixture
.checkSingleSubject <- function(annotation) {
  temp <- unique(annotation[, c("Mixture", "Condition", "BioReplicate")])
  temp$Condition <- factor(temp$Condition)
  temp$Mixture <- factor(temp$Mixture)
  temp1 <- xtabs(~ Mixture + Condition, data = temp)
  singleSubject <- all(temp1 <= "1")
  return(singleSubject)
}

## check .checkTechReplicate
.checkTechReplicate <- function(annotation) {
  temp <- unique(annotation[, c("Mixture", "Run")])
  temp$Mixture <- factor(temp$Mixture)
  temp1 <- xtabs(~Mixture, data = temp)
  TechReplicate <- all(temp1 != "1")
  return(TechReplicate)
}

## check whether there are multiple biological mixtures
.checkMultiBioMixture <- function(annotation) {
  temp <- unique(annotation[, "Mixture"])
  temp <- as.vector(as.matrix(temp))
  return(length(temp) > 1)
}

## check whether there is only single run
.checkSingleRun <- function(annotation) {
  temp <- unique(annotation[, "Run"])
  temp <- as.vector(as.matrix(temp))
  return(length(temp) == 1)
}

# 6 ----------------------------------------------------------------------------
# mixed effect linear models (lmerTest::lmer, lme4)

## 6.1
fit_full_model_spikedin <- function(data) {
  #' @importFrom lmerTest lmer
  #' fit the whole plot and subplot model if the data has no biological variation,
  #' multiple mixtures with multiple technical replicate runs
  fit <- suppressMessages(try(lmerTest::lmer(Abundance ~ 1 + (1 | Mixture) + (1 | Mixture:TechRepMixture)
    + Condition, data = data), TRUE))
  if (!inherits(fit, "try-error")) {
    return(fit)
  } else { # the parameters are not estimable, return null
    return(NULL)
  }
}


## 6.2 
fit_reduced_model_mulrun <- function(data) {
  #' @importFrom lmerTest lmer
  #' fit the reduced with only run effect
  #' fit the whole plot and subplot model if the data has no biological variation,
  #' multiple mixtures or multiple technical replicate runs
  #' or if the data has multiple mixtures but single technical replicate MS run
  fit <- suppressMessages(try(lmerTest::lmer(Abundance ~ 1 + (1 | Run) + Condition, data = data), TRUE))
  if (!inherits(fit, "try-error")) {
    return(fit)
  } else { # the parameters are not estimable, return null
    return(NULL)
  }
}


## 6.3
fit_reduced_model_onerun <- function(data) {
  #' @importFrom lmerTest lmer
  #' fit one-way anova model
  #' fit the whole plot and subplot model if the data has single run
  fit <- suppressMessages(try(lm(Abundance ~ 1 + Condition, data = data), TRUE))
  if (!inherits(fit, "try-error")) {
    return(fit)
  } else { # the parameters are not estimable, return null
    return(NULL)
  }
}


## 6.4
fit_full_model <- function(data) {
  #' @importFrom lmerTest lmer
  #' fit the full model with mixture, techrep and subject effects
  #' fit the whole plot and subplot model if the data has
  #' multiple mixtures, multiple technical replicate runs per mixture and biological variation
  fit <- suppressMessages(try(lmerTest::lmer(Abundance ~ 1 + (1 | Mixture) + (1 | Mixture:TechRepMixture) + # whole plot
    Condition + # subplot
    (1 | BioReplicate:Condition:Mixture), data = data), TRUE))
  if (!inherits(fit, "try-error")) {
    return(fit)
  } else { # if the parameters are not estimable, return null
    return(NULL)
  }
}


## 6.5
fit_reduced_model_techrep <- function(data) {
  #' @importFrom lmerTest lmer
  #' @keywords internal
  #' fit the reduced model with run and subject effects
  #' fit the whole plot and subplot model if the data has
  #' single mixture with multiple technical replicate runs
  fit <- suppressMessages(try(lmerTest::lmer(Abundance ~ 1 + (1 | Run) + # whole plot
    Condition + # subplot
    (1 | BioReplicate:Condition), data = data), TRUE))
  if (!inherits(fit, "try-error")) {
    return(fit)
  } else { # the parameters are not estimable, return null
    return(NULL)
  }
}


#------------------------------------------------------------------------------
