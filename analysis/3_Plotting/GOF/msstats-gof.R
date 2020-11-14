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
data(msstats_prot)

# imports
library(dplyr)
library(data.table)
library(ggplot2)
library(data.table)
library(doParallel)


# fit the model to swip
fx <- formula("Abundance ~ 1 + Condition + (1|Mixture)")
fit <- lmerTest::lmer(fx, msstats_prot %>% subset(Protein == swip))

(fit)

# create a contrast
contrast<-lme4::fixef(fit)
contrast[] <- 0
contrast[grepl("Control",names(contrast))] <- -1/sum(grepl("Mutant",names(contrast)))
contrast[grepl("Mutant",names(contrast))] <- +1/sum(grepl("Mutant",names(contrast)))


lmerTestContrast(fm,contrast) %>% mutate(Contrast='Mutant-Control') %>% unique() %>% knitr::kable()



# functions -------------------------------------------------------------------

devStats <- function(fit) {
	# calculate residual deviance and residual degrees of freedom
  d <- setNames(as.numeric(stats::deviance(fit,REML=FALSE)),"deviance")
  df  <- setNames(as.numeric(stats::df.residual(fit)),"residual.DF")
  return(c(d,df))
}
  
## main -----------------------------------------------------------------------

## fit all models

# register parallel backend
n_cores <- parallel::detectCores() -1
doParallel::registerDoParallel(cores=n_cores)

# loop through all proteins, fit model and calc gof stats
proteins <- unique(as.character(msstats_prot$Protein))
fit_list <- foreach(prot = proteins) %dopar% {
  input <- list(fx, msstats_prot %>% filter(Protein == prot))
  suppressMessages({
	  try(do.call(lmerTest::lmer,input),silent=TRUE)
  })
}
names(fit_list) <- proteins

# remove errors 
idx <- sapply(fit_list,inherits,"try-error")
filt_list <- fit_list[!idx]

# check singular
idx <- sapply(filt_list,lme4::isSingular)
sum(idx)

## calculate goodness-of-fit
## loop through fits and calculate goodness of fit stats
gof_df <- as.data.table(t(sapply(fit_list,devStats)))

# collect results
gof_df$Protein <- names(filt_list)

# calc p-values and identify outliers
gof_df$Pvalue <- pchisq(gof_df$deviance,
			     df = gof_df$residual.DF, 
			     lower.tail = FALSE, log.p = FALSE)
gof_df$isOutlier <- p.adjust(gof_df$Pvalue,method="holm") < 0.05

message("N Outlier Protein: ", sum(gof_df$isOutlier))

gof_results %>% filter(isOutlier) %>% knitr::kable()


## save data
myfile <- file.path(root,"data","gof_results.rda")
save(gof_results,file=myfile,version=2)

## generate a plot
print(nrow(gof_df))

# inspired by edgeR gof function
myfile <- file.path(figsdir,"lmer-gof.pdf")
pdf(myfile)

x <- gof_df$deviance
n <- length(x)
col <- rep_len("black", n)
col[gof_df$isOutlier] <- "blue"
z <- (x-mean(x))/sd(x)
pch <- rep_len(1,n)
pch[gof_df$isOutlier] <- 16
qqnorm(z,col=col,pch=pch)
abline(0,1)

dev.off()
