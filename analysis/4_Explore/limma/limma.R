## Define function: tmt_limma


# Load the library

library(limma)


# create a basic design matrix

group <- as.factor(c(rep("QC",3),rep("WT",4),rep("KO",4)))

group <- factor(group, levels(group)[c(3,1,2)]) # set the factor order

design <- model.matrix(~ 0 + group)

colnames(design) <- c("WT", "KO", "QC")


# the contrast is where the KO vs WT sub-comparison is selected

contrast <- makeContrasts(KO-WT, levels = design)

contrast


# do the linear model fitting

data_limma <- sub

fit <- lmFit(data_limma, design)



# get the fit for the contrast of interest

fit2 <- contrasts.fit(fit, contrast)



# do the empirical Bayes moderation of the test statistic to asses differential expression.

fit2 <- eBayes(fit2, trend = TRUE)



# grab the information in topTable so we can get the data to plot candidates

# the coef parameter has to do with the contrast of interest

# specify no sorting of results and a number that is longer than the data table

tt_limma <- topTable(fit2, coef = 1, sort.by = "none", number = Inf)



# let's see what columns we have

head(tt_limma)

summary(decideTests(fit2))
