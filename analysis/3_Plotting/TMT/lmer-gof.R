#!/usr/bin/env Rscript

# title:
# author:
# description: inspired by edgeR::gof function

## prep the env
root = "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)


## create dir for output
figsdir <- file.path(root,"figs","TMT")
if (!dir.exists(figsdir)) { dir.create(figsdir); message("mkdir ",figsdir) }


## set plotting theme and font
fontdir <- file.path(root,"fonts")
ggtheme(); set_font("Arial", font_path=fontdir)


## data
data(swip)
data(msstats_prot)


## imports
library(dplyr)
library(ggplot2)
library(doParallel)


## functions 

lmerGOF <- function(fm) {
  # calculate Nakagawa R2 and get deviance and residual df from a fit model
  # we will use these stats to evaluate the overall goodness of fit
  r2 <- setNames(as.numeric(r.squaredGLMM.merMod(fm)),nm=c("R2.fixef","R2.total"))
  dev <- setNames(as.numeric(stats::deviance(fm,REML=FALSE)),"deviance")
  ddf  <- setNames(as.numeric(stats::df.residual(fm)),"residual.DF")
  return(c(r2,dev,ddf))
} #EOL


## main

# the model to be fit:
fx <- formula("Abundance ~ 0 + Condition + (1|Mixture)")

# register parallel backend
n_cores <- parallel::detectCores() -1
doParallel::registerDoParallel(cores=n_cores)

# loop through all proteins, fit model, calc gof stats
proteins <- unique(as.character(msstats_prot$Protein))
fit_list <- foreach(prot = proteins) %dopar% {
  input <- list(fx,msstats_prot %>% filter(Protein == prot))
  suppressMessages({
	  try(do.call(lmerTest::lmer,input),silent=TRUE)
  })
}
names(fit_list) <- proteins

# drop errors 
idx <- sapply(fit_list,inherits,"try-error")
filt_list <- fit_list[!idx]

# drop singular
idx <- sapply(filt_list,lme4::isSingular)
filt_list <- filt_list[!idx]


## loop through fits and calculate goodness of fit stats
gof_stats <- foreach(i = seq(filt_list)) %dopar% {
	lmerGOF(filt_list[[i]])
}


## collect results
gof_results <- as.data.table(do.call(rbind,gof_stats))
gof_results$Protein <- names(filt_list)

# calcualte p-values and identify outliers
gof_results$Pvalue <- pchisq(gof_results$deviance, 
			     df = gof_results$residual.DF, 
			     lower.tail = FALSE, log.p = FALSE)
gof_results$outlier <- p.adjust(gof_results$Pvalue,method="holm") < 0.05

message("Outliers: ", sum(gof_results$outlier))
message("Percent: ", round(100*sum(gof_results$outlier)/nrow(gof_results),3))

## save data
myfile <- file.path(root,"data","gof_results.rda")
save(gof_results,file=myfile,version=2)


## generate a plot
plot <- ggplot(gof_results, aes(sample = deviance))
plot <- plot + stat_qq() + stat_qq_line(color="darkred")
plot <- plot + ggtitle("qq-plot of residual deviance")
plot <- plot + ylab("Sample Quantiles")
plot <- plot + xlab("Theoretical Quantiles")
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(axis.line.y=element_line())
plot <- plot + theme(axis.line.x=element_line())


# FIXME: add color for outliers!
myfile <- file.path(figsdir,"lmer-GOF.pdf")
ggsave(myfile,plot)
