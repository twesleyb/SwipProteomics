#!/usr/bin/env Rscript

# title:
# author:
# description: inspired by edgeR::gof function

# prep the env
root = "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

# create dir for output
figsdir <- file.path(root,"figs","TMT")
if (!dir.exists(figsdir)) { dir.create(figsdir); message("mkdir ",figsdir) }

# set plotting theme and font
fontdir <- file.path(root,"fonts")
ggtheme(); set_font("Arial", font_path=fontdir)

# load the data
data(swip)
data(swip_tmt)
data(msstats_prot)

# imports
library(dplyr)
library(ggplot2)
library(doParallel)


# fit the model to swip
df <- swip_tmt %>% filter(Accession == swip)
df$Abundance <- log2(df$Intensity)

# BECAUSE WE HAVE ALREADY DEALT WITH MIXTURE IT SHOULD NOT BE IN THE MODEL
fx <- formula("Abundance ~ 0 + Genotype:Fraction + (1|DOB)")
fm = lmerTest::lmer(fx,df)
contrast<-lme4::fixef(fm)
contrast[] <- 0
contrast[grepl("WT",names(contrast))] <- -1/sum(grepl("WT",names(contrast)))
contrast[grepl("MUT",names(contrast))] <- +1/sum(grepl("MUT",names(contrast)))
lmerTestContrast(fm,contrast) %>% mutate(Contrast='Mutant-Control') %>% unique() %>% knitr::kable()



# functions -------------------------------------------------------------------

lmerGOF <- function(fm) {
  # calc Nakagawa R2 and get deviance and residual df from a fit model
  # we will use these stats to evaluate the overall goodness of fit
  r2 <- setNames(as.numeric(r.squaredGLMM.merMod(fm)),nm=c("R2.fixef","R2.total"))
  dev <- setNames(as.numeric(stats::deviance(fm,REML=FALSE)),"deviance")
  ddf  <- setNames(as.numeric(stats::df.residual(fm)),"residual.DF")
  return(c(r2,dev,ddf))
} #EOL

  
## main -----------------------------------------------------------------------

## fit all models

# the model to be fit:
fx <- formula("Abundance ~ 0 + Condition + (1|Mixture)")

# register parallel backend
n_cores <- parallel::detectCores() -1
doParallel::registerDoParallel(cores=n_cores)

# loop through all proteins, fit model and calc gof stats
proteins <- unique(as.character(msstats_prot$Protein))
fit_list <- foreach(prot = proteins) %dopar% {
  input <- list(fx,msstats_prot %>% filter(Protein == prot))
  suppressMessages({
	  try(do.call(lmerTest::lmer,input),silent=TRUE)
  })
}
names(fit_list) <- proteins

# remove errors 
idx <- sapply(fit_list,inherits,"try-error")
filt_list <- fit_list[!idx]

# drop singular
idx <- sapply(filt_list,lme4::isSingular)
filt_list <- filt_list[!idx]


## calculate goodness-of-fit

## loop through fits and calculate goodness of fit stats
gof_stats <- foreach(i = seq(filt_list)) %dopar% {
	lmerGOF(filt_list[[i]])
}

# collect results
gof_results <- data.table::as.data.table(do.call(rbind,gof_stats))
gof_results$Protein <- names(filt_list)

# calc p-values and identify outliers
gof_results$Pvalue <- pchisq(gof_results$deviance,
			     df = gof_results$residual.DF, 
			     lower.tail = FALSE, log.p = FALSE)
gof_results$isOutlier <- p.adjust(gof_results$Pvalue,method="holm") < 0.05

message("N Outlier Protein: ", sum(gof_results$isOutlier))

gof_results %>% filter(isOutlier) %>% knitr::kable()


## save data
myfile <- file.path(root,"data","gof_results.rda")
save(gof_results,file=myfile,version=2)

## generate a plot
print(nrow(gof_results))

# inspired by edgeR gof function
myfile <- file.path(figsdir,"lmer-gof.pdf")
pdf(myfile)
x <- gof_results$deviance
n <- length(x)
col <- rep_len("black", n)
col[gof_results$isOutlier] <- "blue"
z <- (x-mean(x))/sd(x)
pch <- rep_len(1,n)
pch[gof_results$isOutlier] <- 16
qqnorm(z,col=col,pch=pch)
abline(0,1)
dev.off()
