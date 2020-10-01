#############################################
## fit the proper linear model for each protein
#############################################
#' @importFrom lme4 fixef
#' @import lmerTest
#' @importFrom stats vcov
#' @importFrom dplyr filter
#' @keywords internal
#' fit the proper linear model for each protein
.linear.model.fitting <- function(data) {
  Abundance <- Condition <- Protein <- NULL

  data$Protein <- as.character(data$Protein) ## make sure protein names are character
  proteins <- as.character(unique(data$Protein)) ## proteins
  num.protein <- length(proteins)
  linear.models <- list() # linear models
  s2.all <- NULL # sigma^2
  s2_df.all <- NULL # degree freedom of sigma^2
  pro.all <- NULL # testable proteins
  coeff.all <- list() # coefficients
  ## do inference for each protein individually
  for (i in seq_along(proteins)) {
    #message(paste("Model fitting for Protein :", proteins[i], "(", i, " of ", num.protein, ")"))
    sub_data <- data %>% dplyr::filter(Protein == proteins[i]) ## data for protein i
    # sub_groups <- as.character(unique(sub_data$Condition))
    # if(length(sub_groups) == 1){
    #   stop("Only one condition!")
    # }
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

    if (sub_singleSubject) { # no biological variation within each condition and mixture
      if (sub_TechReplicate & sub_bioMixture) { # multiple mixtures and technical replicates
        # fit the full model with mixture and techrep effects for spiked-in data
        fit <- fit_full_model_spikedin(sub_data)

        if (is.null(fit)) { # full model is not applicable
          # fit the reduced model with only run effect
          fit <- fit_reduced_model_mulrun(sub_data)
        }

        if (is.null(fit)) { # the second model is not applicable
          # fit one-way anova model
          fit <- fit_reduced_model_onerun(sub_data)
        }
      } else {
        if (sub_TechReplicate | sub_bioMixture) { # multiple mixtures or multiple technical replicates
          # fit the reduced model with only run effect
          fit <- fit_reduced_model_mulrun(sub_data)

          if (is.null(fit)) { # the second model is not applicable
            # fit one-way anova model
            fit <- fit_reduced_model_onerun(sub_data)
          }
        } else { # single run case
          # fit one-way anova model
          fit <- fit_reduced_model_onerun(sub_data)
        }
      }
    } else { # biological variation exists within each condition and mixture
      if (sub_bioMixture) { # multiple biological mixtures
        if (sub_TechReplicate) { # multiple technical replicate MS runs
          # fit the full model with mixture, techrep, subject effects
          fit <- fit_full_model(sub_data)

          if (is.null(fit)) { # full model is not applicable
            # fit the reduced model with run and subject effects
            fit <- fit_reduced_model_techrep(sub_data)
          }

          if (is.null(fit)) { # second model is not applicable
            # fit one-way anova model
            fit <- fit_reduced_model_onerun(sub_data)
          }
        } else { # single technical replicate MS run
          # fit the reduced model with only run effect
          fit <- fit_reduced_model_mulrun(sub_data)

          if (is.null(fit)) { # second model is not applicable
            # fit one-way anova model
            fit <- fit_reduced_model_onerun(sub_data)
          }
        }
      } else { # single biological mixture
        if (sub_TechReplicate) { # multiple technical replicate MS runs
          # fit the reduced model with run and subject effects
          fit <- fit_reduced_model_techrep(sub_data)

          if (is.null(fit)) { # second model is not applicable
            # fit one-way anova model
            fit <- fit_reduced_model_onerun(sub_data)
          }
        } else { # single run
          # fit one-way anova model
          fit <- fit_reduced_model_onerun(sub_data)
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
