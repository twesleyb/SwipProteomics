#!/usr/bin/env Rscript

# title: SwipProteomics
# description: working through MSstats groupComparisons function
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

# subset the data for testing
rand_prots <- sample(data_prot$Protein,10)

# test for all the possible pairs of conditions	
#all_results <- groupComparisonTMT(data_prot)	

# args for groupComparisonTMT
data = data_prot %>% filter(Protein %in% rand_prots)
contrast.matrix = 'pairwise'
moderated = FALSE
adj.method = 'BH'

# 1 ---------------------------------------------------------------------------
# init descent

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
	  n_out <- length(data$Protein[is.na(data$Abundance)])
	  data <- data[!is.na(data$Abundance), ]
	  warning("Input 'data' cannot contain missing values (NA).\n",
	    paste(n_out,"proteins with missing values were removed from 'data'."))
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
# going deeper

# mostly figured out... 
# main things: fit the protein wise-models

.proposed.model <- function(data,
                            moderated = FALSE,
                            contrast.matrix = "pairwise",
                            adj.method = "BH") {
  # check that there are some Conditions to compare (>2)
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
  # fit the linear model for each protein, see [4 +++]
  fitted.models <- .linear.model.fitting(data)
  #save(fitted.models,file="fitted.models.rda",version=2)
  # names(fitted.models)
  # [1] "protein" "model"   "s2"      "s2_df"   "coeff"
  ## if moderated, then perform empirical bayes moderation
  if (moderated) {
    ## moderated t-statistic
    # Estimate the prior variance and degree freedom
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
  } else { 
    ## ordinary t-statistic
    # NOTE: this is the default
    s2.prior <- 0
    df.prior <- 0
  } # EIS
  ## extract the linear model fitting results from list
  proteins <- fitted.models$protein # proteins
  s2.all <- fitted.models$s2 # group variance
  s2_df.all <- fitted.models$s2_df # degree freedom of s2
  lms <- fitted.models$model # linear models
  coeff.all <- fitted.models$coeff # coefficients
  ## init empty df to store the inference results
  num.protein <- length(proteins)
  res <- as.data.frame(matrix(rep(NA, 7 * num.protein * ncomp), ncol = 7)) 
  colnames(res) <- c("Protein", "Comparison", "log2FC", 
		     "pvalue", "SE", "DF", "issue")
  data$Condition <- as.factor(data$Condition) # make sure group is factor
  data$Run <- as.factor(data$Run)
  nrun <- length(unique(data$Run)) # check the number of MS runs in the data
  # LOOP to perform statistical comparisons using the protein-wise models:
  # ARG for testing
  #i = 1
  count <- 0
  for (i in seq_along(proteins)) {
    ## get the data for protein i
    prot = proteins[i]
    sub_data <- data %>% dplyr::filter(Protein == prot) 
    ## record the contrast matrix for each protein
    sub.contrast.matrix <- contrast.matrix
    sub_groups <- as.character(unique(sub_data[, c("Condition")]))
    # sort the groups based on alphabetic order
    sub_groups <- sort(sub_groups) 
    ## get protein's linear model
    fit <- lms[[prot]]
    s2 <- s2.all[prot]
    s2_df <- s2_df.all[prot]
    coeff <- coeff.all[[prot]]
    ## if the model is fittable
    if (!is.character(fit)) { 
      s2.post <- (s2.prior * df.prior + s2 * s2_df) / (df.prior + s2_df)
      # LOOP to perform statistical testing for every contrast
      for (j in seq(nrow(sub.contrast.matrix))) {
        count <- count + 1
        res[count, "Protein"] <- proteins[i]
        res[count, "Comparison"] <- row.names(sub.contrast.matrix)[j] 
        # groups with positive coefficients
	idy <- sub.contrast.matrix[j, ] > 0
        positive.groups <- colnames(sub.contrast.matrix)[idy]
        # groups with negative coefficients
	idy <- sub.contrast.matrix[j, ] < 0
        negative.groups <- colnames(sub.contrast.matrix)[idy]
        # make sure at least one group from each side of the contrast exist
        if (any(positive.groups %in% sub_groups) &
          any(negative.groups %in% sub_groups)) {
          contrast.matrix.single <- as.vector(sub.contrast.matrix[j, ])
          names(contrast.matrix.single) <- colnames(sub.contrast.matrix)
	  # class(contrast.matrix.single) = "numeric"
	  # Condition A = 1; Condition B = -1; other = 0
	  # generate contrast matrix
          cm <- .make.contrast.single(fit$model, # see [12]
				      contrast.matrix.single, sub_data)
          # calculate logFC from coefficients
          FC <- (cm %*% coeff)[, 1]
          ## variance and df
          if (inherits(fit$model, "lm")) {
	    # lm
            variance <- diag(t(cm) %*% 
			     summary(fit$model)$cov.unscaled %*% cm) * s2.post
            df.post <- s2_df + df.prior
          } else {
	    # lmer
	    # returns a function that takes Lc and thpars
            vss <- .vcovLThetaL(fit$model) # see [13]
	    ## for the theta and sigma parameters:
            varcor <- vss(t(cm), c(fit$thopt, fit$sigma)) 
            vcov <- varcor$unscaled.varcor * s2
            se2 <- as.matrix(t(cm) %*% as.matrix(vcov) %*% cm)
            ## calculate variance
            vcov.post <- varcor$unscaled.varcor * s2.post
            variance <- as.matrix(t(cm) %*% as.matrix(vcov.post) %*% cm)
            ## calculate df see [14]
            g <- .mygrad(function(x) vss(t(cm),x)$varcor,c(fit$thopt,fit$sigma))
            denom <- try(t(g) %*% fit$A %*% g, silent = TRUE)
            if (inherits(denom, "try-error")) {
              df.post <- s2_df + df.prior
            } else {
              df.post <- 2 * (se2)^2 / denom + df.prior
            }
          } # EIS
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
          # continued from:
	  # if (any(positive.groups %in% sub_groups) &
          #  any(negative.groups %in% sub_groups)) {
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
        } #EIS
      } # ENDS LOOP
    # continued from:
    # if (!is.character(fit)) { 
    } else {
    # very few measurements so that the model is unfittable
      for (k in 1:nrow(sub.contrast.matrix)) {
        count <- count + 1
        res[count, "Protein"] <- proteins[i]
        res[count, "Comparison"] <- row.names(sub.contrast.matrix)[j] 
        out <- .issue.checking( # see [15]
          data = sub_data,
          contrast.matrix = sub.contrast.matrix[j, ]
        )
        res[count, "log2FC"] <- out$logFC
        res[count, "pvalue"] <- NA
        res[count, "SE"] <- NA
        res[count, "DF"] <- NA
        res[count, "issue"] <- out$issue
      } # EOL end loop for comparison
    } # EIS the linear model is fittable
  } # EOL for each protein
  # init empty data frame for results
  res <- as.data.frame(res[seq_len(count), ])
  res$Protein <- as.factor(res$Protein)
  res$log2FC <- as.numeric(as.character(res$log2FC))
  res$pvalue <- as.numeric(as.character(res$pvalue))
  res$adjusted.pvalue <- NA
  comps <- unique(res$Comparison)
  ## Adjust multiple tests for each comparison
  for (i in seq_along(comps)) {
	  q <- p.adjust(res[res$Comparison == comps[i], "pvalue"], adj.method)
    res[res$Comparison == comps[i], "adjusted.pvalue"] <- q
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


# 3 ---------------------------------------------------------------------------
# short dive

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

## 4 ---------------------------------------------------------------------------
# fit protein-wise mixed linear models appropriate for the experimental design

.linear.model.fitting <- function(data) {
  # check input
  data$Protein <- as.character(data$Protein) 
  proteins <- as.character(unique(data$Protein))
  num.protein <- length(proteins)
  # empty objects for output from LOOP:
  s2.all <- NULL          # sigma^2
  pro.all <- NULL         # testable proteins
  s2_df.all <- NULL       # degree freedom of sigma^2
  coeff.all <- list()     # coefficients
  linear.models <- list() # linear models
  # LOOP to do inference for each protein individually see [5-7 +8 +9-11]
  for (prot in proteins) {
    # collect the proteins data and metadata
    sub_data <- data %>% dplyr::filter(Protein == prot)
    sub_annot <- unique(sub_data[, c(
      "Run", "Channel", "BioReplicate",
      "Condition", "Mixture", "TechRepMixture")])
    # check the experimental design, see [5]
    sub_singleSubject <- .checkSingleSubject(sub_annot)
    sub_TechReplicate <- .checkTechReplicate(sub_annot)
    sub_BioMixture <- .checkMultiBioMixture(sub_annot)
    sub_singleRun <- .checkSingleRun(sub_annot)
    ## LOGIC TREE TO AUTOMATICALLY FIT APPROPRIATE LM
    if (sub_singleSubject) { # no bio variation within each condition and mix
      if (sub_TechReplicate & sub_BioMixture) { # multiple mix and tech rep
        # fit the full model with mixture and techrep effects for spiked-in data
	#message("\n(a) Abundance ~ 1 + (1|Mixture) + ",
	#	"(1|Mixture:TechRepMixture) + Group")
        fit <- fit_full_model_spikedin(sub_data) # see [6.1]
        if (is.null(fit)) { # if the full model is not applicable,
          # then fit the reduced model with only run effect
	  #message("\n(b) Abundance ~ 1 + (1|Run) + Group")
          fit <- fit_reduced_model_mulrun(sub_data) # see [6.2]
        }
        if (is.null(fit)) { # if the second model is not applicable,
          # then fit one-way anova model
	  #message("\n(c) Abundance ~ 1 + Group")
          fit <- fit_reduced_model_onerun(sub_data) # see [6.3]
        }
      } else { # sub_singleSubject == TRUE 
        if (sub_TechReplicate | sub_BioMixture) { # multiple TechRep or BioMix
          # fit the reduced model with only run effect
	  #message("\n(d) Abundance ~ 1 + (1|Mixture) + ",
	  #	  "(1|Mixture:TechRepMix) + Group + ",
	  # 	  "(1|Subject:Group:Mixture)")
          fit <- fit_reduced_model_mulrun(sub_data) # see [6.2]
          if (is.null(fit)) { # if the second model is not applicable,
            # then fit one-way anova model
	    #message("\n(e) Abundance ~ 1 + Group")
            fit <- fit_reduced_model_onerun(sub_data) # see [6.3]
          }
        } else { # single run case
          # fit one-way anova model
	  #message("\n(f) Abundance ~ 1 + Group")
          fit <- fit_reduced_model_onerun(sub_data) # see [6.3]
        }
      }
    } else { # biological variation exists within each condition and mixture
      if (sub_bioMixture) { # multiple biological mixtures
        if (sub_TechReplicate) { # multiple technical replicate MS runs
          # fit the full model with mixture, techrep, subject effects
	#  message("\n(g) Abundance ~ 1 + (1|Mixture) + ",
	# 	  "(1|Mixture:TechRepMix) + Group + ",
	# 	  "(1|Subject:Group:Mixture)")
          fit <- fit_full_model(sub_data) # see [6.4]
          if (is.null(fit)) { # full model is not applicable
            # fit the reduced model with run and subject effects
	#    message("\n(h) Abundance ~ 1 + (1|Run) + Group + ",
	# 	    "(1|Subject:Group")
            fit <- fit_reduced_model_techrep(sub_data) # see [6.5]
          }
          if (is.null(fit)) { # if the second model is not applicable
            # then fit one-way anova model
	    #message("\n(i) Abundance ~ 1 + Group")
            fit <- fit_reduced_model_onerun(sub_data) # see [6.3]
          }
        } else { # single technical replicate MS run
          # fit the reduced model with only run effect
	#  message("\n(j) Abundance ~ 1 + (1|Mixture) + ",
	#	  "(1|Mixture:TechRepMix) + Group + ",
	#	  "(1|Subject:Group:Mixture)")
          fit <- fit_reduced_model_mulrun(sub_data) # see [6.2]
          if (is.null(fit)) { # second model is not applicable
            # fit one-way anova model
	    #message("\n(k) Abundance ~ 1 + Group")
            fit <- fit_reduced_model_onerun(sub_data) # see [6.3]
          }
        }
      } else { # single biological mixture
        if (sub_TechReplicate) { # multiple technical replicate MS runs
          # fit the reduced model with run and subject effects
	  #message("\n(L) Abundance ~ 1 + (1|Run) + Group + ",
	#	  "(1|Subject:Group")
          fit <- fit_reduced_model_techrep(sub_data) # see [6.5]
          if (is.null(fit)) { # second model is not applicable
            # fit one-way anova model
	#    message("\n(m) Abundance ~ 1 + Group")
            fit <- fit_reduced_model_onerun(sub_data) # see [6.3]
          }
        } else { # single run
          # fit one-way anova model
	#  message("\n(n) Abundance ~ 1 + Group")
          fit <- fit_reduced_model_onerun(sub_data) # see [6.3]
        } # single technical replicate MS run
      } # single biological mixture
    } # biological variation
    # you now have a fit.
    #writeLines(capture.output(summary(fit)),"fit.txt")
    # NOTE: see [9-11] before IFS
    if (is.null(fit)) {
	    # the model is not fittable
	    linear.models[[prot]] <- "unfittable"
	    pro.all <- c(pro.all, proteins[i])
	    s2.all <- c(s2.all, NA)
	    s2_df.all <- c(s2_df.all, NA)
	    coeff.all[[prot]] <- NA
    } else {
      # extract some basic parameters from the fit
      # if single run case, 
      if (inherits(fit, "lm")) {
        # then estimate variance and df from linear model
        av <- anova(fit)
        coeff <- coef(fit)
        s2_df <- av["Residuals", "Df"]
        if (s2_df == 0) {
          s2 <- 0
        } else {
          # use error variance for testing
          s2 <- av["Residuals", "Mean Sq"]
        }
        linear.models[[prot]] <- list(model = fit)
      # if not lm, then lmer
      } else {
        rho <- list() 
        rho <- .rhoInit(rho, fit, FALSE) # see [7]; was TRUE 
	# names(rho) - key info about the model (fit)
	# [1] "model"   "fixEffs" "sigma"   "thopt"
	# FIXME: .rhoInit should default to FALSE for expediancy
	# add asymptotic var-covar matrix for theta and sigma to list rho
        rho$A <- .calcApvar(rho) # see [9-11]
	# rho$A
	#       [,1]        [,2]
	# [1,]  0.119417491 -0.0010351106
	# [2,] -0.001035111  0.0001256123 
	# anova
        av <- anova(rho$model)
        coeff <- lme4::fixef(rho$model) # didnt we do this somewhere else?
        s2_df <- av$DenDF # DF?
        s2 <- av$"Mean Sq" / av$"F value" # mean^2 / F = Sigma^2
	# store protein model in list
	#message("stored model")
        linear.models[[prot]] <- rho
      } # EIS 
    } #EIS
    # loop munge
    pro.all <- c(pro.all, prot)
    s2.all <- c(s2.all, s2)
    s2_df.all <- c(s2_df.all, s2_df)
    coeff.all[[prot]] <- coeff
  } # EOL for each protein
  # wrap-up
  names(s2.all) <- proteins
  names(s2_df.all) <- proteins
  results_list <- list(protein = pro.all,
			 model = linear.models,
			 s2 = s2.all,
			 s2_df = s2_df.all,
			 coeff = coeff.all)
  return(results_list)
} # EOF


## 5 ----------------------------------------------------------------------------
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
# functions to mixed effect linear models (calls lmerTest::lmer)

## 6.1
fit_full_model_spikedin <- function(data) {
  #' @importFrom lmerTest lmer
  #' fit the whole plot and subplot model if the data has no biological variation
  #' multiple mixtures with multiple technical replicate runs
  fit <- suppressMessages(try(lmerTest::lmer(Abundance ~ 1 + (1 | Mixture) + 
					     (1 | Mixture:TechRepMixture) + 
					     Condition, data = data), TRUE))
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
  #' fit whole plot and subplot model if the data has no biological variation
  #' multiple mixtures or multiple technical replicate runs
  #' or if the data has multiple mixtures but single technical replicate MS run
  fit <- suppressMessages(try(lmerTest::lmer(Abundance ~ 1 + (1 | Run) + 
					     Condition, data = data), TRUE))
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
  #' multiple mixtures, multiple technical replicate runs per mixture, 
  #' and biological variation
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
  fit <- suppressMessages(try(lmerTest::lmer(Abundance ~ 1 + (1 | Run) + # plot
    Condition + # subplot
    (1 | BioReplicate:Condition), data = data), TRUE))
  if (!inherits(fit, "try-error")) {
    return(fit)
  } else { # the parameters are not estimable, return null
    return(NULL)
  }
}


# 7 ---------------------------------------------------------------------------

.rhoInit <- function(rho, model, change.contr = FALSE, mf.final = NULL) {
  #' @importFrom lme4 fixef getME
  #' @importFrom stats sigma
  # Create rho vector containing info about mixed model
  # model = fit
  if (change.contr) { # updating model, this looks expensive!
    rho$model <- .updateModel(model, mf.final = mf.final, 
			      change.contr = change.contr) # see [8]
  } else {
    # model = fit
    rho$model <- model
  }
  # done for both -- add to empty list rho
  # NOTE: THIS IS THE IMPORTANT STUFF!
  rho$fixEffs <- lme4::fixef(rho$model)
  rho$sigma <- stats::sigma(rho$model)
  rho$thopt <- lme4::getME(rho$model,"theta")
  # NOTE: THIS IS THE IMPORTANT STUFF! 
  # For each protein, given its lm fit, 
  # * extract the fixed effect estimates (lme4::fixef)
  # * extract the stdev of residuals (stats::sigma)
  # * extract generalized components from fit (lme4::getME)
  return(rho)
}


# 8 ---------------------------------------------------------------------------
# too deep
# FIXME: consider removing
# NO, can't scrap this. Utilized by some of the complex envs used in later fun

.updateModel <- function(object, mf.final = NULL, ..., change.contr = FALSE) {
  #' @importFrom stats formula getCall terms update.formula
  # useful things:
  # * call(fit) | fit = model yields the lmer formula
  # * model.matrix(fit) yields the contrast matrix
  if (is.null(call <- getCall(object))) { # extract the function call
  # e.g. lmerTest::lmer(formula = Abundance ~ 1 + Condition
    stop("object should contain a 'call' component")
  }
  extras <- match.call(expand.dots = FALSE)$...
  if (!is.null(mf.final)) {
    call$formula <- update.formula(formula(object), mf.final)
  }
  if (any(grepl("sample", call))) {
    call <- as.list(call)[-which(names(as.list(call)) %in% c("data", "subset"))]
    call[["data"]] <- quote(model.frame(object))
  }
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
    }
  }
  if (change.contr) {
    mm <- model.matrix(object)
    contr <- attr(mm, "contrasts")
    ## change contrasts for F tests calculations
    ## list of contrasts for factors
    if (change.contr && length(which(unlist(contr) != "contr.SAS")) > 0) {
      names.facs <- names(contr)
      l.lmerTest.private.contrast <- as.list(rep("contr.SAS", 
						 length(names.facs)))
      names(l.lmerTest.private.contrast) <- names(contr)
      call[["contrasts"]] <- l.lmerTest.private.contrast
    }
    else if (!is.null(contr)) {
      call[["contrasts"]] <- contr[names(contr) %in% attr(terms(call$formula), 
							  "term.labels")]
    }
  }
  call <- as.call(call)
  ff <- environment(formula(object))
  pf <- parent.frame() ## save parent frame in case we need it
  sf <- sys.frames()[[1]]
  ff2 <- environment(object)
  tryCatch(eval(call, envir = ff),
    error = function(e) {
      tryCatch(eval(call, envir = sf),
        error = function(e) {
          tryCatch(eval(call, envir = pf),
            error = function(e) {
              eval(call, envir = ff2)
            }
          )
        }
      )
    }
  )
}

# 9 ---------------------------------------------------------------------------
# utilizes .devFunTheta and .myhess doing some math I don't understand
# seems like .devFunTheta is this thing: `ans` = function(thpars) { }
# which does some math with the lm = fit = model = envff

# maths: 

# dev <- envff$pp$ldL2() + 
#          ( envff$resp$wrss() + 
#           envff$pp$sqrL(1) ) / sigsq + n * log(2*pi*sigsq)

# p <- ncol(envff$pp$RX())
# dev <- dev + 2 * determinant(envff$pp$RX())$modulus - p * log(2 * pi * sigsq)

# envff$resp$wrss = fun return("the weighted residual sum of squares")
# envff$pp$sqrL(1) = 0.006859864 # some sort of constant

# NOTE: depends upon 10 and 11

.calcApvar <- function(rho) {
  # calc asymptotic variance covariance matrix of variance parameters 
  # based on theta parameters and sigma
  #' @keywords internal
  # returns devFunTheta a class with attributes:
  # srcref, thopt, and class
  # attr(dd,"thopt") = 0.5316409
  # attr(dd,"srcref") = a class including: function(thpars) { ff(thpars...
  dd <- .devfunTheta(rho$model) # see [10]
  # calculate 2x2 hessian matrix
  h <- .myhess(dd, c(rho$thopt, sigma = rho$sigma)) # see [11]
  #       [,1]        [,2]
  # [1,]  18.03627   148.6282
  # [2,] 148.62825 17146.7792
  # calculate Choleski factorization of h (decomposition)
  ch <- try(base::chol(h), silent = TRUE)
  #       [,1]        [,2]
  # [1,] 4.246914  34.99677
  # [2,] 0.000000 126.18243
  if (inherits(ch, "try-error")) {
    return(rho)
  }
  # calculate the inverse of the Choleski matrix
  A <- 2 * base::chol2inv(ch)
  #       [,1]        [,2]
  # [1,]  0.119417491 -0.0010351106
  # [2,] -0.001035111  0.0001256123
  # calculate eigenvalue of hessian matrix h (spectral decomposition)
  eigval <- eigen(h, symmetric = TRUE, only.values = TRUE)$values
  # [1] 17148.06878    16.74671 
  # tol ~ sqrt(.Machine$double.eps)
  isposA <- !(min(eigval) < sqrt(.Machine$double.eps))
  if (!isposA) {
    warning("Asymptotic covariance matrix A is not positive!")
  }
  # return inverse of cholesky matrix of hessian d
  # NOTE: important bits:
  # looks like .myhess is some sort of wrapper around hessian, using devFunTheta
  #  >>>  dd <- .devfunTheta(rho$model) 
  #  >>>  h <- .myhess(dd, c(rho$thopt, sigma = rho$sigma))
  #  >>>  ch <- try(base::chol(h), silent = TRUE)
  #  >>>  A <- 2 * base::chol2inv(ch)
  return(A)
}

# 10 ---------------------------------------------------------------------------
# this is confusing... another call to updateModule, but with sneaky param
# devFunOnly, ... do some math ...

.devfunTheta <- function(fm) {
  # fm = rho$model = fit
  # devfun function as a function of optimal parameters
  #' @importFrom lme4 isGLMM isLMM getME
  #' @importFrom methods is
  #' @keywords internal
  stopifnot(is(fm, "merMod"))
  np <- length(fm@pp$theta) # e.g. 1
  nf <- length(lme4::fixef(fm)) # e.g. 14
  is_glm <- lme4::isGLMM(fm) # is the model a glm?
  # if the model is not a glm...
  if (!is_glm) {
    np <- np + 1L # e.g. np + 1 = 2
  }
  n <- nrow(fm@pp$V)
  # annother call to update model...
  # includes sneaky arg devFunOnly
  ff <- .updateModel(fm, devFunOnly = TRUE) 
  ## ff is now a function!
  reml <- lme4::getME(fm, "is_REML") # e.g. TRUE
  envff <- environment(ff)
  # if a linear model,
  is_lm <- lme4::isLMM(fm) # TRUE
  if (is_lm) {
    # then declare a function ans
    ans <- function(thpars) {
      stopifnot(is.numeric(thpars), length(thpars) == np)
      ff(thpars[-np])
      sigsq <- thpars[np]^2
      dev <- envff$pp$ldL2() + (envff$resp$wrss() + envff$pp$sqrL(1)) / sigsq + n * log(2 * pi * sigsq)
      if (reml) {
        p <- ncol(envff$pp$RX())
        dev <- dev + 2 * determinant(envff$pp$RX())$modulus - p * log(2 * pi * sigsq)
      }
      return(dev)
    } # EOF
  } # EIS
  attr(ans, "thopt") <- fm@pp$theta # update thopt parameter with
  class(ans) <- ".devfunTheta" # update a class attribute
  return(ans)
}


# 11 ---------------------------------------------------------------------------

.myhess <- function(fun, x, fx = NULL, delta = 1e-4) {
  # usage: h <- .myhess(dd, c(rho$thopt, sigma = rho$sigma))
  # fun = dd -- the function, c(thopt,sigma), delta)
  #fun = dd
  #x = c(rho$thopt,rho$sigma)
  #fx = NULL
  #delta = 1e-4
  # calculate hessian matrix
  #' @keywords internal
  nx <- length(x) 
  #fx <- if (!is.null(fx)) fx else fun(x, ...)
  fx <- if (!is.null(fx)) fx else fun(x) # function(thpars)
  H <- array(NA, dim = c(nx, nx)) # empty
  # loop to fill H
  for (j in 1:nx) {
    ## Diagonal elements:
    xadd <- xsub <- x
    xadd[j] <- x[j] + delta
    xsub[j] <- x[j] - delta
    H[j, j] <- (fun(xadd) - 2 * fx +
      fun(xsub)) / delta^2
    ## Upper triangular (off diagonal) elements:
    for (i in 1:nx) {
      if (i >= j) break
      xaa <- xas <- xsa <- xss <- x
      xaa[c(i, j)] <- x[c(i, j)] + c(delta, delta)
      xas[c(i, j)] <- x[c(i, j)] + c(delta, -delta)
      xsa[c(i, j)] <- x[c(i, j)] + c(-delta, delta)
      xss[c(i, j)] <- x[c(i, j)] - c(delta, delta)
      H[i, j] <- (fun(xaa) - fun(xas) -
        fun(xsa) + fun(xss)) /
        (4 * delta^2)
    }
  }
  ## Fill in lower triangle:
  H[lower.tri(H)] <- t(H)[lower.tri(H)]
  return(H)
}

# 12 --------------------------------------------------------------------------
#cm <- .make.contrast.single(fit$model, contrast.matrix.single, sub_data)

#fit = fit$model
#contrast = contrast.matrix.single
#sub_data = sub_data

.make.contrast.single <- function(fit, contrast, sub_data) {
  ## make constrast
  #' @importFrom stats coef
  #' @importFrom lme4 fixef
  #' @keywords internal
  ## when there are some groups which are all missing
  sub_groups <- as.character(levels(sub_data[, c("Condition")]))
  # groups with positive coefficients
  positive.groups <- names(contrast)[contrast > 0]
  # groups with negative coefficients
  negative.groups <- names(contrast)[contrast < 0]
  # if some groups not exist in the protein data
  check <- all(positive.groups %in% sub_groups) & 
	  all(negative.groups %in% sub_groups)
  # do something if some groups are missing
  if (!check) {
    contrast.single <- contrast[sub_groups]
    ## tune the coefficients of positive groups so that their summation is 1
    temp <- contrast.single[contrast.single > 0]
    temp <- temp * (1 / sum(temp, na.rm = TRUE))
    contrast.single[contrast.single > 0] <- temp
    ## tune the coefficients of positive groups so that their summation is 1
    temp2 <- contrast.single[contrast.single < 0]
    temp2 <- temp2 * abs(1 / sum(temp2, na.rm = TRUE))
    contrast.single[contrast.single < 0] <- temp2
    ## set the coefficients of non-existing groups to zero
    contrast[] <- 0
    contrast[sub_groups] <- contrast.single
  }
  ##
  if (inherits(fit, "lm")) {
    coef_name <- names(stats::coef(fit))
  } else {
    coef_name <- names(lme4::fixef(fit))
  }
  ## intercept
  temp <- coef_name[grep("Intercept", coef_name)]
  intercept_c <- rep(0, length(temp))
  names(intercept_c) <- temp
  ##
  if (length(temp) == 0) {
    intercept_c <- NULL
  }
  ## group
  temp <- coef_name[grep("Condition", coef_name)]
  tempcontrast <- contrast[sub_groups]
  ##
  group_c <- tempcontrast[gsub("Condition", "", temp)]
  names(group_c) <- temp
  ##
  if (length(temp) == 0) {
    group_c <- NULL
  }
  ## combine all
  newcontrast <- c(intercept_c, group_c)
  if (inherits(fit, "lm")) {
    contrast1 <- newcontrast[!is.na(stats::coef(fit))]
  } else {
    # lmer
    contrast1 <- newcontrast[!is.na(lme4::fixef(fit))]
  }
  return(contrast1)
}


# 13 ---------------------------------------------------------------------------

.vcovLThetaL <- function(fm) {
  ## returns Lc %*% vcov as a function of theta parameters %*% t(Lc)
  #' @importFrom lme4 isGLMM isLMM fixef
  #' @importFrom methods is
  #' @keywords internal
  stopifnot(is(fm, "merMod"))
  np <- length(fm@pp$theta)
  nf <- length(lme4::fixef(fm))
  if (!lme4::isGLMM(fm)) { # if not a GLM
    np <- np + 1L
  }
  # .updateModel returns a function with some stuff in it?
  ff2 <- .updateModel(fm, devFunOnly = TRUE)
  # convert to an environment
  envff2 <- environment(ff2)
  #ls(envff2)
  #[1] "lmer_Deviance" "lower"         "pp"            "resp" 
  #is.environment(envff2) == TRUE
  if (lme4::isLMM(fm)) {
    ans <- function(Lc, thpars) {
      stopifnot(is.numeric(thpars), length(thpars) == np)
      sigma2 <- thpars[np]^2
      ff2(thpars[-np])
      # what is RXi? pp is of class merPredD from lme4- some sort of generator
      vcov_unscaled <- tcrossprod(envff2$pp$RXi())
      vcov_out <- sigma2 * vcov_unscaled
      return(list(
        varcor = as.matrix(Lc %*% as.matrix(vcov_out) %*% t(Lc)),
        unscaled.varcor = vcov_unscaled,
        sigma2 = sigma2
      ))
    }
  }
  class(ans) <- ".vcovLThetaL"
  return(ans)
}


# 14 ---------------------------------------------------------------------------

.mygrad <- function(fun, x, delta = 1e-4,
                    method = c("central", "forward", "backward"), ...) {
  ## calculate gradient
  #' @keywords internal
  method <- match.arg(method)
  nx <- length(x)
  if (method %in% c("central", "forward")) {
    Xadd <- matrix(rep(x, nx), nrow = nx, byrow = TRUE) + diag(delta, nx)
    fadd <- apply(Xadd, 1, fun, ...)
  }
  if (method %in% c("central", "backward")) {
    Xsub <- matrix(rep(x, nx), nrow = nx, byrow = TRUE) - diag(delta, nx)
    fsub <- apply(Xsub, 1, fun, ...) ## eval.parent perhaps?
  }
  res <- switch(method,
    "forward" = (fadd - fun(x, ...)) / delta,
    "backward" = (fun(x, ...) - fsub) / delta,
    "central" = (fadd - fsub) / (2 * delta)
  )
  res
}


# 15 --------------------------------------------------------------------------

.issue.checking <- function(data,
                            contrast.matrix) {
  ## check the reason for results with NA
  #' @keywords internal
  #' check the possible reason for untestable comparison
  ## choose each comparison
  contrast.matrix.sub <- contrast.matrix
  # groups in the sub data
  sub_groups <- as.character(unique(data$Group))
  # groups with positive coefficients
  positive.groups <- names(contrast.matrix.sub)[contrast.matrix.sub > 0]
  # groups with negative coefficients
  negative.groups <- names(contrast.matrix.sub)[contrast.matrix.sub < 0]
  if (is.null(positive.groups) | is.null(negative.groups)) {
    stop("Please check the contrast.matrix. 
         Each row must have both positive and negative values,
         and their sum must be 1!")
  }
  if (any(positive.groups %in% sub_groups) &
    any(negative.groups %in% sub_groups)) {
    logFC <- NA
    issue <- "unfittableModel"
  } else {
    # more than one condition
    if (all(!positive.groups %in% sub_groups) &
      any(negative.groups %in% sub_groups)) {
      logFC <- (-Inf)
      issue <- "oneConditionMissing"
    } else {
      if (any(positive.groups %in% sub_groups) &
        all(!negative.groups %in% sub_groups)) {
        logFC <- Inf
        issue <- "oneConditionMissing"
      } else {
        logFC <- NA
        issue <- "completeMissing"
      }
    }
  }
  # return
  return(list(logFC = logFC, issue = issue))
}
