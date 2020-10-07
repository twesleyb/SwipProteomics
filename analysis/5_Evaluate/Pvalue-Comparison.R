#!/usr/bin/env Rscript

# title: SwipProteomics
# description: 
# author: Tyler W Bradshaw <twesleyb10@gmail.com>
# os: windows linux subsystem (WSL)

## Input ----------------------------------------------------------------------
# * all_results.rda - MSstats intrafraction results

# specify the projects root directory:
root <- "~/projects/SwipProteomics"

## Output ---------------------------------------------------------------------

## FUNCTIONS -----------------------------------------------------------------

mkdir <- function(...) {
	# create a new directory
	newdir <- file.path(...)
	if (!dir.exists(newdir)) { 
		dir.create(newdir)
		message(paste("created",newdir))
	}
}


fix_colname <- function(df,old_colname,new_colname) {
	# change a column's name in a data.frame
	colnames(df)[which(colnames(df) == old_colname)] <- new_colname
	return(df)
}

## Prepare the working environment --------------------------------------------

# project directories
datadir <- file.path(root,"data"); mkdir(datadir)
rdatdir <- file.path(root,"rdata"); mkdir(rdatdir)
downdir <- file.path(root,"downloads"); mkdir(downdir)


# Prepare the R environment ---------------------------------------------------

# load renv
renv::load(root, quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	suppressWarnings({ library(getPPIs) }) # FIXME: annoying warnings!
	library(data.table)
	library(MSstatsTMT)
})

# load functions in root/R
suppressPackageStartupMessages({ devtools::load_all() })

# load the data in root/data
data(swip_tmt) # the preprocessed data (contains edgeR stats)
data(swip_dep) # the DEP stats

# load the data in root/rdata
myfile <- file.path(root,"rdata","all_results.rda")
load(myfile)

# set global plotting settings
ggtheme(); set_font("Arial", font_path = fontdir)

## compare pvalues to DEP and EdgeR pipelines ----------------------------------

# 1. combine EdgeR and DEP stats
# collect edgeR stats for intra-fraction comparisons
tmpdf <-  swip_tmt %>% 
	group_by(Fraction) %>% 
	dplyr::select(Fraction,Accession,PValue) %>% 
	unique()
# combine edgeR and DEP stats
stats_df <- swip_dep %>% 
	dplyr::select(Fraction,Accession,p.val) %>% 
	left_join(tmpdf,by=c("Fraction","Accession"))
# clean-up colnames
stats_df <- fix_colname(stats_df,"PValue","EdgeR")
stats_df <- fix_colname(stats_df,"p.val","DEP")

# 2. add MSstatsTMT stats
# annotate with fraction
fraction <- regmatches(all_results$Label,
		       regexpr("\\F[0-9]{1,2}",all_results$Label))
all_results$Fraction <- fraction
# get relevant cols
df <- all_results %>% dplyr::select(Fraction,Protein,pvalue)
# fix colnames
colnames(df) <- c("Fraction","Accession","MSstatsTMT")

# combine with other stats
all_stats <- left_join(df,stats_df,by=c("Accession","Fraction"))

# drop rows in which pvals from all methods is NA
idx <- apply(all_stats, 1, function(x) all(is.na(x[c(3,4,5)])))
filt_stats <- all_stats[!idx,]

# the spearman rank correlation is very high
rho1 <- cor(x=filt_stats$MSstatsTMT, y=filt_stats$DEP, 
	    method="spearman", use='pairwise.complete.obs')
rho2 <- cor(x=filt_stats$MSstatsTMT, y=filt_stats$EdgeR, 
	    method="spearman", use='pairwise.complete.obs')
rho3 <- cor(x=filt_stats$DEP, y=filt_stats$EdgeR, 
	    method="spearman", use='pairwise.complete.obs')

df <- data.table('cor(MSstatsTMT, DEP)'=rho1,
		 'cor(MSstatsTMT, EdgeR)'=rho2,
		 'cor(DEP, EdgeR)' = rho3)
knitr::kable(df)


## ----------------------------------------------------------------------------

# generate a plot examining coorelation between edgeR and DEP pvalues
plot <- ggplot(data=stats_df,aes(x=EdgeR,y=DEP))
plot <- plot + geom_point()
plot <- plot + xlab("P-Value (edgeR)")
plot <- plot + ylab("P-Value (DEP)")
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(panel.border=element_rect(colour="black",fill=NA,size=1))

# save as pdf
myfile <- file.path(figsdir,"PValue_correlation_scatterplot.pdf")
ggsave(myfile,plot,height=5,width=5)
message(paste("\nSaved",myfile))

## ----------------------------------------------------------------------------

# tidy the data
df <- reshape2::melt(stats_df,id=c("Fraction","Accession"),
		     value.name = "PValue",
		     variable.name="Method")

# loop generate p-value histograms for every intra-fraction comparison
plots <- list()
for (fraction in unique(df$Fraction)) {
	plot <- ggplot(df %>% filter(Fraction == fraction),
		       aes(x=PValue,color=Method))
	plot <- plot + geom_histogram(bins=100)
	plot <- plot + theme(panel.background = element_blank())
	plot <- plot + theme(panel.border=element_rect(colour="black",
						       fill=NA,size=1))
	plot <- plot + ggtitle(paste("Fraction:",fraction))
	plots[[fraction]] <- plot
}

# save as a single pdf
myfile <- file.path(figsdir,"PValue_Histograms.pdf")
ggsavePDF(plots,myfile)
message(paste("\nSaved",myfile))

if (file.exists("Rplots.pdf")) { unlink("Rplots.pdf") }

# DONE!

# We have:

# EdgeR lvl 0
# DEP lvl 0
# MSstats lvl 0

# EdgeR lvl 1
# MSstats lvl 1

# EdgeR modules
# MSstats modules
