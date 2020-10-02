#!/usr/bin/env Rscript

#' check.design
#' @keywords internal

.checkSingleSubject <- function(annotation) {
  # check if there is a single Biological subject
  temp <- unique(annotation[, c("Mixture", "Condition", "BioReplicate")])
  temp$Condition <- factor(temp$Condition)
  temp$Mixture <- factor(temp$Mixture)
  temp1 <- xtabs(~ Mixture + Condition, data = temp)
  singleSubject <- all(temp1 <= "1")
  return(singleSubject)
}


.checkTechReplicate <- function(annotation) {
  # check if there is Technical Replication of Mixture
  temp <- unique(annotation[, c("Mixture", "Run")])
  temp$Mixture <- factor(temp$Mixture)
  temp1 <- xtabs(~Mixture, data = temp)
  TechReplicate <- all(temp1 != "1")
  return(TechReplicate)
}


.checkMulBioMixture <- function(annotation) {
  # check whether there are multiple biological Mixtures
  temp <- unique(annotation[, "Mixture"])
  temp <- as.vector(as.matrix(temp))
  return(length(temp) > 1)
}


.checkSingleRun <- function(annotation) {
  # check whether there is only single run
  temp <- unique(annotation[, "Run"])
  temp <- as.vector(as.matrix(temp))
  return(length(temp) == 1)
}
