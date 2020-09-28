#!/usr/bin/env Rscript

# title: SwipProteomics
# description: analysis of Swip TMT spatial proteomics data with MSstatsTMT
# author: Tyler W Bradshaw <twesleyb10@gmail.com>


## Input ----------------------------------------------------------------------

# input data in root/data:
input_dir = "data/PSM.zip"

# PSM.zip contains:
input_psm = "PSM/5359_PSM_Report.xlsx"
input_samples = "PSM/5359_Sample_Report.xlsx"

## Prepare the working environment --------------------------------------------

root <- "~/projects/SwipProteomics"
datadir <- file.path(root,"data")
rdatdir <- file.path(root,"rdata") # temporary/large files not tracked by git

# Prepare the R environment ---------------------------------------------------

# load renv
renv::load(root,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
	library(MSstatsTMT)
})

# load functions in root/R
suppressPackageStartupMessages({ devtools::load_all() })


## unzip the directory containing the raw data ---------------------------------

unzip(file.path(root,input_dir)) # unzip root/PSM.zip into current dir


# load sample data in root/rdata --------------------------------------
# this excel spreadsheet was from Greg, exported from PD(?)

# pass colnames to read_excel
col_names <- c("Sample","Mixture","MS.Channel","drop",
	       "Channel","Proteomics ID","ConditionFraction","Experiment")
samples <- readxl::read_excel(input_samples,col_names=col_names)

unlink(input_samples)


## read PSM data from excel ------------------------------------------------

# read raw data, NOTE: this takes several minutes
raw_pd <- readxl::read_excel(input_psm,progress=FALSE)

unlink(input_psm) 

unlink(tools::file_path_sans_ext(basename(input_dir))) # rmdir ./PSM

## re-format sample annotations for MSstats -----------------------------------
# Extract relevant annotations for MSstatsTMT
# Format colums for MSstatsTMT

f1 <- function(x) { # x = samples$ConditionFraction
	# extract sample 'Condition'
	paste0("F",as.numeric(sapply(strsplit(x,"Control|Mutant|SPQC"),"[",2)))
}

f2 <- function(x) { # x = samples$ConditionFraction 
	# extract sample 'Fraction'
	sapply(strsplit(x,"[0-9]{1,2}"),"[",1)
}

# Munge ConditionFraction Column to Fraction and Condition cols
samples$Fraction <- f1(samples$ConditionFraction)
samples$Condition <- f2(samples$ConditionFraction)
samples$Condition <- as.character(interaction(samples$Condition,samples$Fraction)) # added
samples$ConditionFraction <- NULL

# Set the formatting of the channel for normalization for MSstats
samples$Condition[grepl("SPQC",samples$Condition)] <- "Norm" 

# clean-up Mixture column 
samples$Mixture <- gsub("F","M",samples$Mixture)

# Remove un-needed col
samples$drop <- NULL

# FIXME: how should BioReplicate be defined?
samples$BioReplicate <- paste(samples$Condition,paste0("R",as.numeric(as.factor(samples$Experiment))),sep=".")


## re-format PSM data for MSstatsTMT ---------------------------------------------

# make columns look like MSstats by replacing special characters with '.'
chars <- c(" ","\\[","\\]","\\:","\\(","\\)","\\/","\\+","\\#","\\-")
new_names <- gsub(paste(chars,collapse="|"),".",colnames(raw_pd))
colnames(raw_pd) <- new_names

# add X if starts column starts with ".."
colnames(raw_pd) <- gsub("^\\.\\.","X..",colnames(raw_pd))


## map Spectrum.Files to MS.Channels ------------------------------------------

# munge extra

# collect all Spectrum.Files grouped by Experiment
# split 'Spectrum.File' at first "_", ID##### is experimental identifier.
all_files <- raw_pd$Spectrum.File
exp_files <- lapply(split(all_files, sapply(strsplit(all_files,"_"),"[",1)),
		    unique)
files_dt <- data.table("Experiment ID" = rep(names(exp_files),
					     times=sapply(exp_files,length)),
	               "Run" = unlist(exp_files))

# collect all MS.Channels, grouped by Experiment
all_channels <- samples$MS.Channel
exp_channels <- split(all_channels, sapply(strsplit(all_channels,"_"),"[",1))
exp_dt <- data.table("Experiment ID" = rep(names(exp_channels),
					times=sapply(exp_channels,length)),
	             "MS.Channel" = unlist(exp_channels))

# use exp_channels to create numeric ID for MS Fraction
exp_fraction_list <- lapply(exp_channels, function(x) {
	       setNames(as.numeric(as.factor(x)),x)
	   })
x <- unlist(exp_fraction_list)
named_fractions <- setNames(x,sapply(strsplit(names(x),"\\."),"[",2))


# build annotation file -------------------------------------------------------

# annotation: data.frame which contains 
# Run - match Spectrum.File
# Fraction - TMT mixture may be fractionated into multiple fractions
# TechRepMixture - all 1, no repeats of a mixture
# Mixture - concatenation of TMT labeled samples - an MS experiment
# Channel - the TMT channels
# BioReplicate - individual subjects (mice)
# Condition - WT MUT norm (SPQC)

# create annotation_dt from Spectrum.Files and MS.Runs
annotation_dt <- left_join(files_dt,exp_dt,by="Experiment ID")
idx <- match(annotation_dt$"MS.Channel", samples$MS.Channel)
annotation_dt$Fraction <- named_fractions[annotation_dt$MS.Channel]
annotation_dt$TechRepMixture <- rep(1,length(idx))
annotation_dt$Mixture <- samples$Mixture[idx]
annotation_dt$Condition <- samples$Condition[idx]
#annotation_dt$Subject <- samples$BioReplicate[idx]
annotation_dt$Channel <- samples$Channel[idx]

# Add BioFraction column 
# FIXME: how to pass additional covariates to MSstats?
annotation_dt$BioFraction <- samples$Fraction[idx]

# how to handle repeated measures design? 
annotation_dt$BioReplicate <- samples$BioReplicate[idx]
#as.character(interaction(samples$BioReplicate[idx],
#						       samples$Fraction[idx]))

# whoops, fix norm#.F# -- should be just Norm - the SPQC sample which is
# analyzed in Technical duplicate for every experiment?
annotation_dt$BioReplicate[grepl("Norm",annotation_dt$BioReplicate)] <- "Norm"

# Remove un-needed cols
annotation_dt$"Experiment ID" <- NULL
annotation_dt$"MS.Channel" <- NULL

# check the annotation file
# this basic design is repeated 3x (M = 3):
knitr::kable(annotation_dt[c(1:16),])

# save to file
myfile <- file.path(datadir,"annotation_dt.rda")
save(annotation_dt,file=myfile,version=2)
message(paste("saved",myfile))


## coerce data to MSstats format ---------------------------------

# NOTE: this takes a considerable amount of time
data_pd <- PDtoMSstatsTMTFormat(raw_pd, annotation_dt)

# ** Shared PSMs (assigned in multiple proteins) are removed.
# ** 705 features have 1 or 2 intensities across runs and are removed.
# ** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.
# ** For peptides overlapped between fractions of M1_1, use the fraction with maximal average abundance.
# ** For peptides overlapped between fractions of M2_1, use the fraction with maximal average abundance.
# ** For peptides overlapped between fractions of M3_1, use the fraction with maximal average abundance.
# ** Fractions belonging to same mixture have been combined.

#knitr::kable(head(data_pd))

# save to file
myfile <- file.path(rdatdir,"data_pd.rda")
save(data_pd,file=myfile,version=2)
message(paste("saved",myfile))


# Protein summarize and normalization -----------------------------------------
# use MSstats for protein summarization	
# Sample Summary does not look correct, should be 

# NOTE: this takes a considerable amount of time
# FIXME: Joining, by = ("Run", "Channel") # unexpected output
# FIXME: remove extremely long message about normalization between runs
# Normalization between MS runs for Protein : Q9DA08 ( 8525  of  8556 )
data_prot <- proteinSummarization(data_pd,
				  method="msstats",	
                                  global_norm=TRUE,	
                                  reference_norm=TRUE,	
                                  remove_norm_channel = TRUE)

# ** Protein-level summarization done by MSstats.

# save to file
myfile <- file.path(rdatdir,"data_prot.rda")
save(data_prot,file=myfile,version=2)
message(paste("saved",myfile))
	

# Protein-level statistical testing -------------------------------------------
# Tests for significant changes in protein abundance across conditions based on
# a family of linear mixed-effects models in TMT experiment. Experimental
# design of case-control study (patients are not repeatedly measured) is
# automatically determined based on proper statistical model.	

# test for all the possible pairs of conditions	
all_results <- groupComparisonTMT(data_prot)	

#------------------------------------------------------------------------------
# get boundary fit error from  lme... something is not correct with our design
# or MSstats's interpretation of the design.
# ^ fixed BioReplicate Condition.BioFraction.Replicate (e.g. Control.F4.R1)
# not all proteins get good fits.

# groupComparisonTMT

library(dplyr)

load_args <- function(){
load(file.path(rdatdir,"data_prot.rda"))
data = data_prot
contrast.matrix = 'pairwise'
moderated = FALSE
adj.method = 'BH'
remove_norm_channel = TRUE
remove_empty_channel = TRUE
}

groupComparisonTMT <- function(data,
                               contrast.matrix = "pairwise",
                               moderated = FALSE,
                               adj.method = "BH") {

  # check: input data should contain the following columns:
  required_cols <- c("Protein", "BioReplicate", "Abundance", "Run", "Channel",
		     "Condition", "TechRepMixture", "Mixture")

  # stop if any required columns are missing
  if (!all(required_cols %in% colnames(data))) {
    is_missing <- !(required.info %in% colnames(data))
      stop("Please check input 'data'.",
	   "The following column(s) are required:",
	   paste(required_cols[is_missing], collapse = ", "))
  }

  # stop if "Empty" in data$Condition
  if ("Empty" %in% unique(data$Condition)) {
	  stop("Column 'Condition' cannot contain 'Empty'.")
  }

  # stop if "Norm" in data$Condition
  if ("Norm" %in% data$Condition) {
	  stop("Column 'Condition' cannot contain 'Norm'.")
  }

  # check for rows with NA intensities
  if (any(is.na(data$Abundance))) {
	  warning("Missing values not tolerated.",
		  " Missing values (NA) will be removed.")
	  data <- data[!is.na(data$Abundance), ]
  }

  ## change some column names as used in group comparison function
  ## Ting: need to change later for time course design
  # colnames(data)[colnames(data) == "BioReplicate"] <- "Subject"
  # colnames(data)

  # data$runchannel <- paste(as.numeric(as.factor(data$Run)),
  # as.numeric(as.factor(data$Channel)),
  # sep = "_")
  # colnames(data)[colnames(data) == 'runchannel'] <- 'Subject'
  #colnames(data)[colnames(data) == "Condition"] <- "Group"

  #########################################################################
  ## Inference
  result <- .proposed.model(data, moderated, contrast.matrix, adj.method)
  #########################################################################

  ### check column name in order to use groupComparisonPlot from MSstats
  colnames(result)[colnames(result) == "Comparison"] <- "Label"
  colnames(result)[colnames(result) == "adjusted.pvalue"] <- "adj.pvalue"

  return(result)
}


.makeContrast <- function(conditions) {
  # generate contrast matrix given vector of conditions
  ncomp <- length(conditions) * (length(conditions) - 1) / 2
  contrast.matrix <- matrix(rep(0, length(conditions) * ncomp),
    ncol = length(conditions)
  )
  colnames(contrast.matrix) <- conditions
  count <- 0
  contrast.matrix.rownames <- NULL
  for (j in seq_len(length(conditions) - 1)) {
    for (k in (j + 1):length(conditions)) {
      count <- count + 1
      # save row name
      contrast.matrix.rownames <- c(
        contrast.matrix.rownames,
        paste(conditions[j], conditions[k], sep = "-")
      )
      # set constrast value
      contrast.matrix[count, conditions[j]] <- 1
      contrast.matrix[count, conditions[k]] <- -1
    }
  }
  rownames(contrast.matrix) <- contrast.matrix.rownames
  return(contrast.matrix)
}

###############################################################################
## functions to check the experimental design

.checkSingleSubject <- function(annotation) {
  # check single subject within each Mixture.Condition
  df <- unique(annotation[, c("Mixture", "Condition", "BioReplicate")])
  df$Group <- factor(df$Condition)
  df$Mixture <- factor(df$Mixture)
  xtab <- xtabs(~ Mixture + Group, data = df)
  singleSubject <- all(xtab <= "1")
  return(singleSubject)
}

.checkTechReplicate <- function(annotation) {
  # check .checkTechReplicate
  df <- unique(annotation[, c("Mixture", "Run")])
  df$Mixture <- factor(df$Mixture)
  xtab <- xtabs(~Mixture, data = df)
  isTechRep <- all(xtab != "1")
  return(isTechRep)
}

.checkMultiBioMixture <- function(annotation) {
  # check whether there are multiple biological mixtures
  mixtures <- unique(as.character(annotation[, "Mixture"]))
  return(length(mixtures) > 1)
}

.checkSingleRun <- function(annotation) {
  # check whether there is only single run
  runs <- unique(as.character(annotation[, "Run"]))
  return(length(runs) == 1)
}

###############################################################################


###############################################################################
## lmer functions

fit_full_model <- function(data) {
  ## 1. fit the full model with mixture, techrep and subject effects
  #' @importFrom lmerTest lmer
  #' @keywords internal
  #' fit the whole plot and subplot model if the data has
  #' multiple mixtures, multiple technical replicate runs per mixture and 
  #' biological variation
  fit <- suppressMessages({
	  try(lmerTest::lmer(Abundance ~ 1 + (1 | Mixture) + 
			     (1 | Mixture:TechRepMixture) + # whole plot
			     Group + # subplot
		             (1 | Subject:Group:Mixture), data = data), 
	      silent = TRUE) })
  if (!inherits(fit, "try-error")) {
    return(fit)
  } else { 
  # the parameters are not estimable, return null
    return(NULL)
  }
}

fit_reduced_model_techrep <- function(data) {
  ## 2. fit the reduced model with run and subject effects
  #' @importFrom lmerTest lmer
  #' @keywords internal
  #' fit the whole plot and subplot model if the data has
  #' single mixture with multiple technical replicate runs
  fit <- suppressMessages({
	  try(lmerTest::lmer(Abundance ~ 1 + (1 | Run) + # whole plot
			     Group + # subplot
			     (1 | Subject:Group), data = data), silent=TRUE) })
  if (!inherits(fit, "try-error")) {
    return(fit)
  } else { 
    # else the parameters are not estimable, return null
    return(NULL)
  }
}

fit_full_model_spikedin <- function(data) {
  ## 3. fit the reduced model with mixture and techrep effects
  #' @importFrom lmerTest lmer
  #' @keywords internal
  #' fit the whole plot and subplot model if the data has no biological variation,
  #' multiple mixtures with multiple technical replicate runs
  fit <- suppressMessages({
	  try(lmerTest::lmer(Abundance ~ 1 + (1 | Mixture) + 
			     (1 | Mixture:TechRepMixture) + 
			     Group, data = data), silent = TRUE) })
  if (!inherits(fit, "try-error")) {
    return(fit)
  } else { 
  # the parameters are not estimable, return null
    return(NULL)
  }
}

fit_reduced_model_multirun <- function(data) {
  ## 4. fit the reduced with only run effect
  #' @importFrom lmerTest lmer
  #' @keywords internal
  #' fit the whole plot and subplot model if the data has no biological variation,
  #' multiple mixtures or multiple technical replicate runs
  #' or if the data has multiple mixtures but single technical replicate MS run
  fit <- suppressMessages({
	  try(lmerTest::lmer(Abundance ~ 1 + (1 | Run) + Group, data = data), 
	      silent = TRUE) })
  if (!inherits(fit, "try-error")) {
    return(fit)
  } else {
  # the parameters are not estimable, return null
  return(NULL)
  }
}

fit_reduced_model_onerun <- function(data) {
  # fit one-way anova model
  # fit the whole plot and subplot model if the data has single run
  #' @importFrom lmerTest lmer
  fit <- suppressMessages(try(lm(Abundance ~ 1 + Group, data = data), TRUE))
  if (!inherits(fit, "try-error")) {
    return(fit)
  } else {
    # the parameters are not estimable, return null
    return(NULL)
  }
}

###############################################################################
.linear.model.fitting <- function(data) {

  require(dplyr)

  data$Protein <- as.character(data$Protein) 
  proteins <- as.character(unique(data$Protein))
  num.protein <- length(proteins)

  linear.models <- list() # linear models
  s2.all <- NULL # sigma^2
  s2_df.all <- NULL # degree freedom of sigma^2
  pro.all <- NULL # testable proteins
  coeff.all <- list() # coefficients

  # loop to do inference for each protein individually
  for (i in seq_along(proteins)) {

    # data for protein i
    sub_data <- data %>% dplyr::filter(Protein == proteins[i]) 

    # Record the annotation information
    sub_annot <- unique(sub_data[, c(
      "Run", "Channel", "BioReplicate",
      "Condition", "Mixture", "TechRepMixture"
    )])

    # check the experimental design
    sub_singleSubject <- .checkSingleSubject(sub_annot)
    sub_TechReplicate <- .checkTechReplicate(sub_annot)
    sub_BioMixture <- .checkMultiBioMixture(sub_annot)
    sub_singleRun <- .checkSingleRun(sub_annot)

    # All NULL?
 fit <- fit_full_model_spikedin(sub_data)

 fit <- fit_reduced_model_multirun(sub_data)

 fit <- fit_reduced_model_onerun(sub_data)

 fit <- fit_full_model(sub_data)

 fit <- fit_reduced_model_techrep(sub_data)

    # estimate variance and df from linear models
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


###############################################################################
data = data
if (!exists("moderated")) { moderated = TRUE }
if (!exists("contrast.matrix")) { contrast.matrix = "pairwise" }
if (!exists("adj.method")) { adj.method = "BH" }

.proposed.model <- function(data,
                            moderated = TRUE,
                            contrast.matrix = "pairwise",
                            adj.method = "BH") {

  # stop if there are fewer than 2 conditions	
  conditions <- as.character(unique(data$Condition)) # all conditions
  if (length(conditions) < 2) {
    stop("Please check the Condition column in annotation file.",
	 "There must be at least two conditions!")
  }

  ## contrast matrix can be matrix or character vector.

  if (is.matrix(contrast.matrix)) {
    # if matrix, then colnames must be in conditions 
    if (!all(colnames(contrast.matrix) %in% conditions)) {
      stop("Please check input 'contrast.matrix'.",
	   "Column names must be in unique(data$Condition)!")
    }
  } else if (is.character(contrast.matrix)) {
    # create constrast matrix for pairwise comparison
    contrast.matrix <- .makeContrast(conditions)
  }

  # number of comparisons
  ncomp <- nrow(contrast.matrix)

  # fit the linear model for each protein
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
