#' lmerTestProtein
#' @import dplyr lmerTest data.table
#' @export lmerTestProtein


lmerTestProtein <- function(protein, fx, msstats_prot, contrasts) {

  suppressPackageStartupMessages({
    library(lme4)
    library(dplyr)
    library(lmerTest)
    library(data.table)
  })

  getIndex <- function(namen,dm=lme4::fixef(fm)) {
    # a helper function to find column index 
    idy <- which(grepl(namen,names(dm)))
  }

  # check input args
  stopifnot(inherits(fx,"formula"))
  stopifnot(inherits(protein,"character"))
  stopifnot(inherits(msstats_prot,"data.frame"))

  # input contrasts should be astats_df matrix, numeric vector, or list of such 
  if (inherits(contrasts,"matrix")) {
	  # rows of matrix are contrasts
	  # rownames of matrix are contrast names (comparisons)
	  contrast_list <- unlist(apply(contrasts,1,list),recursive=FALSE)
	  stopifnot(all(sapply(contrast_list,sum)==0))
  } else if (inherits(contrasts,"numeric")) {
	  # the sum of contrast vector should be 0
	  # FIXME: why is sum sometimes wrong?!?!
	  stopifnot(sum(contrasts[contrasts<0], contrasts[contrasts>0]) == 0)

	  contrast_list <- list(contrasts)
  } else if (inherits(contrasts,"list")) {
          contrast_list <- contrasts
	  # comparison names will be generated from pos and neg coefficients
  } else {
	  stop("problem parsing input 'contrasts'.")
  }

  # subset the data
  subdat <- msstats_prot %>% as.data.table() %>% dplyr::filter(Protein == protein)
  if (any(is.na(subdat))) {
	# the data cannot contain missing values
  	warning("The data cannot contain missing values.")
        return(NULL)
  }

  # fit the model
  fm <- lmerTest::lmer(fx, data=subdat)

  # evaluate statistical comparisons for all contrasts
  stats_list <- list()

  for (comparison in contrast_list) {

    # the contrast to be tested, ensure it matches names(fixef(fm))
    contrast <- comparison[names(sort(sapply(names(comparison),getIndex)))]
    ## FIXME: we do this work twice!

    pos_group <- paste0(names(contrast[contrast>0]),collapse="+")
    neg_group <- paste(names(contrast[contrast<0]),collapse="+")
    comparison <- paste(pos_group,neg_group,sep="-")

    # assess contrast
    test_results <- lmerTestContrast(fm, contrast)

    # collect results
    test_results$Protein <- protein 
    stats_list[[comparison]] <- test_results
  } # EOL for every comparison

  # compile results
  stats_df <- do.call(rbind, stats_list)
  rownames(stats_df) <- NULL

  # sort cols
  stats_df <- stats_df[,c("Protein","Contrast","log2FC","percentControl",
	      "Tstatistic","Pvalue","SE", "DF", "isSingular")]

  return(stats_df)
} #EOF
