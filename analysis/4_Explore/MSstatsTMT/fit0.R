#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: fit lmer model for comparisons between Control and Mutant mice
# adjusted for differences in subcellular fraction.

## prepare the env
root <- "~/projects/SwipProteomics"
renv::load(root); devtools::load_all(root)

## inputs:
data(msstats_prot)
data(gene_map)
data(swip)

## other imports
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
})


## Functions ------------------------------------------------------------------
# Modified from internal MSstatsTMT functions.

# * lmerTestProtein - model fitting and statistical testing for a given protein

lmerTestProtein <- function(protein, fx, msstats_prot, 
			    contrast_matrix, s2_prior=0, df_prior=0) {
  # subset the data
  subdat <- msstats_prot %>% filter(Protein == protein)
  if (any(is.na(subdat))) {
	# the data cannot contain missing values
  	warning("The data cannot contain missing values.")
        return(NULL)
  }
  # fit the model
  fm <- lmerTest::lmer(fx, data=subdat)
  # compute Satterthwaite degrees of freedom and other key statistics
  model_summary <- summary(fm, ddf = "Satterthwaite")
  s2_df <- as.numeric(model_summary$coefficients[,"df"][1]) 
  coeff <- model_summary$coeff[,"Estimate"] # == lme4::fixef(fm)
  sigma <- model_summary$sigma # == stats::sigma(fm)
  theta <- model_summary$optinfo$val # aka thopt == lme4::getME(fm, "theta")
  vcov <- model_summary$vcov # variance-covariance matrix
  se2 <- as.numeric(contrast_matrix %*% vcov %*% contrast_matrix) # == variance
  s2_prior = sigma^2
  # calculate posterior s2
  s2_post <- (s2_prior * df_prior + sigma^2 * s2_df) / (df_prior + s2_df)
  # calcuate symtoptic var-covar matrix
  A <- fm@vcov_varpar
  # FIXME: there might not always be two?
  g <- c(contrast_matrix %*% fm@Jac_list[[1]] %*% contrast_matrix,
	 contrast_matrix %*% fm@Jac_list[[2]] %*% contrast_matrix)
  denom <- as.numeric(t(g) %*% A %*% g)
  # NOTE: which is correct?
  df_post <- (2 * se2) / (denom + df_prior) # or se2^2 ???
  #df_post <- (2 * se2^2) / (denom + df_prior)
  # compute fold change and the t-statistic
  FC <- (contrast_matrix %*% coeff)[, 1]
  t <- FC / sqrt(se2) 
  # compute the p-value given t-statistic and df.post
  p <- 2 * pt(-abs(t), df = df_post) 
  # compile results
  rho <- list()
  rho$protein <- protein
  rho$model <- fx
  comparison <- paste(names(contrast_matrix)[contrast_matrix == +1], 
		      names(contrast_matrix)[contrast_matrix == -1],sep="-")
  rho$stats <- data.frame(protein=protein,contrast=comparison,
  		 log2FC=FC, percentControl=2^FC, Pvalue=p,
  		 Tstatistic=t, SE=sqrt(se2), DF=df_post, 
		 isSingular=lme4::isSingular(fm))
  return(rho)
} #EOF


## check Swip's fit -----------------------------------------------------------

# demonstrate fit:
fx0 <- formula("Abundance ~ 0 + Condition + (1|Mixture)")
fm0 <- lmerTest::lmer(fx0, msstats_prot %>% filter(Protein == swip))

# ddf options: c("Satterthwaite", "Kenward-Roger", "lme4"))
model_summary <- summary(fm0, ddf = "Satterthwaite")
model_summary

# FIXME: need to replace with more reproducible code
knitr::kable(r.squaredGLMM.merMod(fm0))

# build a contrast matrix:
cm0 <- lme4::fixef(fm0)
cm0[] <- 0
cm0["ConditionControl.F7"] <- -1
cm0["ConditionMutant.F7"] <- +1 

# test a comparison defined by contrast_matrix
model0 <- lmerTestProtein(swip, fx0, msstats_prot, cm0)
model0$stats %>% knitr::kable()


## loop to fit all proteins ----------------------------------------------------

n_cores <- parallel::detectCores() - 1
BiocParallel::register(BiocParallel::SnowParam(n_cores))

prots = unique(as.character(msstats_prot$Protein))

results_list <- foreach(protein = prots) %dopar% {
	suppressMessages({
	  try(lmerTestProtein(protein, fx0, msstats_prot, cm0),silent=TRUE)
	})
} # EOL


## collect results ------------------------------------------------------------

idx <- unlist(sapply(results_list,class)) != "try-error"
filt_list <- results_list[which(idx)]
results_df <- bind_rows(sapply(filt_list,"[[","stats"))

# drop singular
results_df <- results_df %>% filter(!isSingular)
results_df$isSingular <- NULL

## annotate with gene symbols
idx <- match(results_df$protein,gene_map$uniprot)
results_df <- tibble::add_column(results_df,
  				 symbol=gene_map$symbol[idx],
  				 .after="protein")

## adjust pvals 
results_df <- tibble::add_column(results_df, 
			 Padjust=p.adjust(results_df$Pvalue,"BH"),
			 .after="Pvalue")

## sort
results_df <- results_df %>% arrange(Pvalue)

# examine top results
results_df %>% head() %>% knitr::kable()

# status
message("Total number of significant proteins: ",
	sum(results_df$Padjust < 0.05))


