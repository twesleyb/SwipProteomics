#!/usr/bin/env Rscript

.checkSingleSubject <- function(annotation) {
###############################################################
## check single subject within each condition in each mixture
###############################################################
#' @keywords internal
  temp <- unique(annotation[, c("Mixture", "Condition", "BioReplicate")])
  temp$Condition <- factor(temp$Condition)
  temp$Mixture <- factor(temp$Mixture)
  temp1 <- xtabs(~ Mixture + Condition, data = temp)
  singleSubject <- all(temp1 <= "1")

  return(singleSubject)
}

.checkTechReplicate <- function(annotation) {
#############################################
## check .checkTechReplicate
#############################################
#' @keywords internal
  temp <- unique(annotation[, c("Mixture", "Run")])
  temp$Mixture <- factor(temp$Mixture)
  temp1 <- xtabs(~Mixture, data = temp)
  TechReplicate <- all(temp1 != "1")

  return(TechReplicate)
}

.checkMulBioMixture <- function(annotation) {
#############################################
## check whether there are multiple biological mixtures
#############################################
#' @keywords internal
  temp <- unique(annotation[, "Mixture"])
  temp <- as.vector(as.matrix(temp))

  return(length(temp) > 1)
}

.checkSingleRun <- function(annotation) {
#############################################
## check whether there is only single run
#############################################
#' @keywords internal
  temp <- unique(annotation[, "Run"])
  temp <- as.vector(as.matrix(temp))

  return(length(temp) == 1)
}
