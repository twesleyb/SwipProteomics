#!/usr/bin/env Rscript

## Illustrate the Module-level analysis

root <- "~/projects/SwipProteomics"
renv::load(root)

<<fit the model, eval=TRUE, echo=TRUE>>=

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

## load SwipProteomics data
# devtools::install_github("twesleyb/SwipProteomics")
library(SwipProteomics)

data(swip)
data(gene_map)
data(partition)
data(msstats_prot)

## formula to be fit:
fx0 <- formula("Abundance ~ 0 + Genotype:BioFraction + (1|Mixture) + (1|Protein)")

wash_prots <- mapID("Washc")

#wash_module <- partition[swip] # M17
#prots <- names(partition[partition==wash_module])

# fit the model
fm <- lmerTest::lmer(fx0, msstats_prot %>% filter(Protein %in% wash_prots))

# create contrast

# evaluate contrast
lmerTestContrast(fm,contrast) %>% knitr::kable()

@

#########################################



library(doParallel)

lmerGOF <- function(fm) {
  dev <- setNames(as.numeric(stats::deviance(fm,REML=FALSE)),"deviance")
  ddf <- setNames(as.numeric(stats::df.residual(fm)),"residual.DF")
  return(c(dev,ddf))
}

myfile <- file.path(root,"rdata","ne_surprise_partition.csv")
part_df <- fread(myfile,drop=1)
partition <- unlist(part_df) + 1
partition[partition %in% as.numeric(names(which(table(partition)<5)))] <- 0

modules <- split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))
proteins = unlist(modules)

n = parallel::detectCores() - 1
doParallel::registerDoParallel(n)

stats_list <- foreach(module = names(modules)) %dopar% {
	subdat <- msstats_prot %>% filter(Protein %in% modules[[module]])
	lmer_args <- list(fx0,subdat)
	fm <- do.call(lmerTest::lmer,lmer_args)
	return(lmerGOF(fm))
} #EOL
names(stats_list) <- names(modules)

df <- data.table(do.call(rbind,stats_list))

# yikes
x <- df$deviance
n <- length(x)
col <- rep_len("black",n)
z <- (x-mean(x))/sd(x)
pch <- rep_len(1,n)
qqnorm(z,col=col,pch=pch)
abline(0,1)


dm <- msstats_prot %>% reshape2::dcast(Protein ~ Mixture + Channel + Condition,
				       value.var = "Abundance") %>% 
                       as.data.table() %>% as.matrix(rownames="Protein")
plot_data <- vsn::meanSdPlot(dm)
plot <- plot_data$gg

library(ggplot2)

plot <- plot + ggtitle("vsn::meanSdPlot")
plot <- plot + theme(title = element_text(colour="red"))
