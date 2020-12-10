#' lmerFit

#' @importFrom lme4 lmerControl
#' @importFrom lmerTest lmer

lmerFit <- function(args_list,control=TRUE) {

  # fit lmerTest lmer with some lmerControls
  if (control) {
    args_list[["control"]] <- lme4::lmerControl(check.conv.singular = "ignore",
   					        calc.derivs = FALSE, 
					        check.rankX = "stop.deficient")
  }

# fit the model
  fm <- do.call(lmerTest::lmer, args_list)

  return(fm)
} #EOF
