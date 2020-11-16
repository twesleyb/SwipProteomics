#!/usr/bin/env Rscript

#' @export rmOutlierPSM

rmOutlierPSM <- function(mixture, psm_df, nbins=5, nSD=4, nComplete=2) {

  # this function identifies QC PSM-level outliers
  # NOTE: PSM with incomplete (n != nComplete) features are removed
  # ^we cannot assess reproducibility of a PSM if it was only quantified once

  # collect psm with incomplete (missing) data
  is_incomplete <- psm_df %>% 
	  dplyr::filter(Mixture == mixture & Condition == "Norm") %>%
	  group_by(PSM) %>% mutate(nObs=length(PSM)) %>% ungroup() %>% 
	  filter(nObs != nComplete) %>%
	  select(PSM) %>% unlist() %>% unique()

  # subset, keep psm with complete features
  filt_df <- psm_df %>% 
	  dplyr::filter(Mixture == mixture & Condition == "Norm") %>% 
	  group_by(PSM) %>% mutate(nObs=length(PSM)) %>% ungroup()  %>%
	  # drop PSM with missing features
	  filter(nObs == nComplete) %>% 
	  # calculate mean of log2 Intensity for binning psm
	  group_by(Mixture,PSM) %>% mutate(meanQC=mean(log2(Intensity))) %>% 
	  ungroup()

  # group ratio data into intensity bins
  breaks <- stats::quantile(filt_df$meanQC, seq(0, 1, length.out=nbins+1),
  		   names=FALSE,na.rm=TRUE)
  filt_df <- filt_df %>% 
  	mutate(bin = cut(meanQC, breaks, labels=FALSE, include.lowest=TRUE))

  # compute log ratio (difference) of QC measurments
  filt_df <- filt_df %>% group_by(Mixture,PSM) %>% 
  	mutate(ratioQC = diff(log2(Intensity))) %>% ungroup()

  # summarize mean and SD of each bin
  filt_df <- filt_df %>% group_by(bin) %>% 
  	mutate(meanRatio = mean(ratioQC,na.rm=TRUE),
  	       binSD = sd(ratioQC), .groups="drop") %>% ungroup()

  # flag outliers
  filt_df <- filt_df %>% mutate(isLow = ratioQC < meanRatio - nSD*binSD)
  filt_df <- filt_df %>% mutate(isHigh = ratioQC > meanRatio + nSD*binSD)
  filt_df <- filt_df %>% mutate(isOutlier = isLow | isHigh)

  # summary
  summary_df <- filt_df %>% group_by(bin) %>%
  	summarize(meanIntensity = mean(Intensity),
  		  meanRatio = mean(ratioQC,na.rm=TRUE),
  		  binSD = sd(ratioQC), 
  		  nOut = sum(isOutlier), .groups="drop")

  # status
  if (length(is_incomplete) > 0) {
    warning(length(is_incomplete),
	    " PSM with incomplete features in mixture ", mixture,
	    " were removed.", call.=FALSE)
  }

  out_df <- filt_df %>% select(Mixture,PSM, bin, meanQC, ratioQC, binSD, isOutlier)
  return(out_df %>% filter(isOutlier))
} # EOF
