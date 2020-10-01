#!/usr/bin/env Rscript

fit_full_model <- function(data) {
#############################################
## fit the full model with mixture, techrep and subject effects
#############################################
#' @importFrom lmerTest lmer
#' @keywords internal
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

fit_reduced_model_techrep <- function(data) {
#############################################
## fit the reduced model with run and subject effects
#############################################
#' @importFrom lmerTest lmer
#' @keywords internal
#' fit the whole plot and subplot model if the data has
#' single mixture with multiple technical replicate runs
  fit <- suppressMessages(try(lmerTest::lmer(Abundance ~ 1 + (1 | Run) + # whole plot
    Condition + # subplot
    (1 | BioReplicate:Condition), data = data), TRUE))

  if (!inherits(fit, "try-error")) {
    return(fit)
  } else { # if the parameters are not estimable, return null
    return(NULL)
  }
}

fit_full_model_spikedin <- function(data) {
#############################################
## fit the reduced model with mixture and techrep effects
#############################################
#' @importFrom lmerTest lmer
#' @keywords internal
#' fit the whole plot and subplot model if the data has no biological variation,
#' multiple mixtures with multiple technical replicate runs
  fit <- suppressMessages(try(lmerTest::lmer(Abundance ~ 1 + (1 | Mixture) + (1 | Mixture:TechRepMixture)
    + Condition, data = data), TRUE))

  if (!inherits(fit, "try-error")) {
    return(fit)
  } else { # if the parameters are not estimable, return null
    return(NULL)
  }
}

fit_reduced_model_mulrun <- function(data) {
#############################################
## fit the reduced with only run effect
#############################################
#' @importFrom lmerTest lmer
#' @keywords internal
#' fit the whole plot and subplot model if the data has no biological variation,
#' multiple mixtures or multiple technical replicate runs
#' or if the data has multiple mixtures but single technical replicate MS run
  fit <- suppressMessages(try(lmerTest::lmer(Abundance ~ 1 + (1 | Run) + Condition, data = data), TRUE))

  if (!inherits(fit, "try-error")) {
    return(fit)
  } else { # if the parameters are not estimable, return null
    return(NULL)
  }
}

fit_reduced_model_onerun <- function(data) {
#############################################
## fit one-way anova model
#############################################
#' @importFrom lmerTest lmer
#' @keywords internal
#' fit the whole plot and subplot model if the data has single run
  fit <- suppressMessages(try(lm(Abundance ~ 1 + Condition, data = data), TRUE))

  if (!inherits(fit, "try-error")) {
    return(fit)
  } else { # if the parameters are not estimable, return null
    return(NULL)
  }
}
