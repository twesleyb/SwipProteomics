#!/usr/bin/env Rscript

#' linear.model.fitting(data)
#'
#' fit the proper linear model for each protein
#'
#' @importFrom lme4 fixef
#' @import lmerTest
#' @importFrom stats vcov
#' @importFrom dplyr filter
#' @keywords internal

.linear.model.fitting <- function(data) {

  proteins <- as.character(unique(as.character(data$Protein)))
  num.protein <- length(proteins)

  # empty objects for output of LOOP:
  s2.all <- NULL          # sigma^2
  pro.all <- NULL         # testable proteins
  s2_df.all <- NULL       # degree freedom of sigma^2
  coeff.all <- list()     # coefficients
  linear.models <- list() # lmers

  ## LOOP to do inference for each protein individually:
  for (i in seq_along(proteins)) {
    sub_data <- data %>% dplyr::filter(Protein == proteins[i]) 

    # Record the annotation information
    sub_annot <- unique(sub_data[, c(
      "Run", "Channel", "BioReplicate",
      "Condition", "Mixture", "TechRepMixture"
    )])

    # check the experimental design
    sub_singleSubject <- .checkSingleSubject(sub_annot)
    sub_TechReplicate <- .checkTechReplicate(sub_annot)
    sub_bioMixture <- .checkMulBioMixture(sub_annot)
    sub_singleRun <- .checkSingleRun(sub_annot)

    ## parse the experimental design:
    if (sub_singleSubject) { # if no Bio variation w/in ea Condition and Mix

      if (sub_TechReplicate & sub_bioMixture) { # if multi Mix and TechRep

	# fit: lmer(A ~ 1 + (1|M) + (1|M|T) + C)
        fit <- fit_full_model_spikedin(sub_data) # [1]

        if (is.null(fit)) { # if full model is not applicable
          # then fit:  lmer(A ~ 1 + (1|R) + C)
          fit <- fit_reduced_model_mulrun(sub_data) # [2]
        }

        if (is.null(fit)) { # if the second model is not applicable
          # then fit: lm(A ~ 1 + C)
          fit <- fit_reduced_model_onerun(sub_data) # [3]
        }

      } else { 

        if (sub_TechReplicate | sub_bioMixture) { # if multi Mix | TechRep

          # fit [2] 
          fit <- fit_reduced_model_mulrun(sub_data)

          if (is.null(fit)) { # if the second model is not applicable
            # fit [3]
            fit <- fit_reduced_model_onerun(sub_data)
          }

        } else { # single run case

          # fit [3]
          fit <- fit_reduced_model_onerun(sub_data)
        }
      }

    } else { # biological variation exists within each Condition and Mixture

      if (sub_bioMixture) { # if multiple biological mixtures

        if (sub_TechReplicate) { # if multiple technical replicate MS runs

          # fit: lmer(A ~ 1 + (1|M) + (1|M:TechRepMix) + C + (1|BioRep:C:M))
          fit <- fit_full_model(sub_data) # [4]

          if (is.null(fit)) { # if the full model is not applicable
            # then fit: lmer(A ~ 1 + (1|Run) + C + (1|BioRep:C)
            fit <- fit_reduced_model_techrep(sub_data) # [5]
          }

          if (is.null(fit)) { # if model 5 is not applicable
            # fit [3]
            fit <- fit_reduced_model_onerun(sub_data)
          }

        } else { # single technical replicate MS run
	  # fit [2]
          fit <- fit_reduced_model_mulrun(sub_data)

          if (is.null(fit)) { # if model 2 is not applicable
            # then fit [3]
            fit <- fit_reduced_model_onerun(sub_data)
          }
        }

      } else { # single biological mixture

        if (sub_TechReplicate) { # multiple technical replicate MS runs
          # fit [5]
          fit <- fit_reduced_model_techrep(sub_data) 

          if (is.null(fit)) { # if model 5 is not applicable
	    # fit [3]
            fit <- fit_reduced_model_onerun(sub_data)
          }

        } else { # single run
          # fit [3]
          fit <- fit_reduced_model_onerun(sub_data)

        } # single technical replicate MS run
      } # single biological mixture
    } # biological variation

    ## estimate variance and df from linear models

    if (!is.null(fit)) { # the model is fittable

      if (inherits(fit, "lm")) { # if single run case [3]

        # estimate the coeff from fixed-effects model
        av <- anova(fit)
        coeff <- stats::coef(fit)
        s2_df <- av["Residuals", "Df"]
        if (s2_df == 0) {
          s2 <- 0
        } else {
          # use error variance for testing
          s2 <- av["Residuals", "Mean Sq"]
        }

	# store the fitted model in linear.models list
        linear.models[[proteins[i]]] <- list(model = fit)

      } else { # mixed effects models [1, 2, 4, 5]

        # estimate the coeff from mixed-effects model 
        rho <- list() 
        rho <- .rhoInit(rho, fit, TRUE) # save lmer outcome in rho envir var
        rho$A <- .calcApvar(rho) # asymptotic var-covar matrix for theta & sigma

        av <- anova(rho$model) # utilize lmerTest added capability
        coeff <- lme4::fixef(rho$model)
        s2_df <- av$DenDF
        s2 <- av$"Mean Sq" / av$"F value"

	# store the fitted model in linear.models list
        linear.models[[proteins[i]]] <- rho
      }

      pro.all <- c(pro.all, proteins[i])
      s2.all <- c(s2.all, s2)
      s2_df.all <- c(s2_df.all, s2_df)
      coeff.all[[proteins[i]]] <- coeff

    } else { # the model is not fittble
      linear.models[[proteins[i]]] <- "unfittable"
      pro.all <- c(pro.all, proteins[i])
      s2.all <- c(s2.all, NA)
      s2_df.all <- c(s2_df.all, NA)
      coeff.all[[proteins[i]]] <- NA
    }
  } # EOL for each protein

  names(s2.all) <- proteins
  names(s2_df.all) <- proteins

  return(list(
    protein = pro.all,
    model = linear.models,
    s2 = s2.all,
    s2_df = s2_df.all,
    coeff = coeff.all
  ))
} #EOF
