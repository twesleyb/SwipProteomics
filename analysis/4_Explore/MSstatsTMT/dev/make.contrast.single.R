.make.contrast.single <- function(fit, contrast, sub_data) {
#############################################
## make constrast
#############################################
# 	MSstats
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
  if (!(all(positive.groups %in% sub_groups) &
    all(negative.groups %in% sub_groups))) {
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

  if (inherits(fit, "lm")) {
    coef_name <- names(stats::coef(fit))
  } else {
    coef_name <- names(lme4::fixef(fit))
  }

  ## intercept
  temp <- coef_name[grep("Intercept", coef_name)]
  intercept_c <- rep(0, length(temp))
  names(intercept_c) <- temp
  if (length(temp) == 0) {
    intercept_c <- NULL
  }

  ## group
  temp <- coef_name[grep("Condition", coef_name)]
  tempcontrast <- contrast[sub_groups]
  group_c <- tempcontrast[gsub("Condition", "", temp)]
  names(group_c) <- temp
  if (length(temp) == 0) {
    group_c <- NULL
  }

  ## combine all
  newcontrast <- c(intercept_c, group_c)
  if (inherits(fit, "lm")) {
    contrast1 <- newcontrast[!is.na(stats::coef(fit))]
  } else {
    contrast1 <- newcontrast[!is.na(lme4::fixef(fit))]
  }

  return(contrast1)
}
