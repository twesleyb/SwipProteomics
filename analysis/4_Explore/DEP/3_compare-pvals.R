#!/usr/bin/env Rscript

# title: SwipProteomics
# description: 
# author: twab <twesleyb10@gmail.com>
# os: windows linux subsystem (WSL)

## INPUT ----------------------------------------------------------------------
# specify project's root directory
ROOT <- "~/projects/SwipProteomics"

## OUTPUT --------------------------------------------------------------------
# * plot comparing edgeR and DEP p-values

## FUNCTIONS -----------------------------------------------------------------

fix_colname <- function(df,old_colname,new_colname) {
	# change a column's name in a data.frame
	colnames(df)[which(colnames(df) == old_colname)] <- new_colname
	return(df)
}


mkdir <- function(...) {
	# create a new directory
	newdir <- file.path(...)
	if (!dir.exists(newdir)) { 
		dir.create(newdir)
		message(paste("created",newdir))
	}
}


## Prepare the working environment ---------------------------------------------

# directory for output tables:
tabsdir <- file.path(ROOT,"tables"); mkdir(tabsdir)

# directory for output data:
datadir <- file.path(ROOT,"data"); mkdir(datadir)

# directory for temporary data:
rdatdir <- file.path(ROOT,"rdata"); mkdir(rdatdir)

# directory for figures:
figsdir <- file.path(ROOT,"figs","DEP"); mkdir(figsdir)


## Prepare the R environment -------------------------------------------------

# load renv
renv::load(ROOT,quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(DEP) # twesleyb/DEP 
	library(dplyr) # for manipulating data
	library(data.table) # for working with data.tables
	library(ggplot2) # for plotting
})

# load functions in root/R
devtools::load_all(ROOT)

# set global plotting settings
ggtheme(); set_font("Arial", font_path = fontdir)

# load the data in root/data
data(swip_tmt) # the preprocessed data (contains edgeR stats)
data(swip_dep) # the DEP stats

## ----------------------------------------------------------------------------

# collect edgeR stats for intra-fraction comparisons
tmpdf <-  swip_tmt %>% 
	group_by(Fraction) %>% 
	select(Fraction,Accession,PValue) %>% 
	unique()
# combine edgeR and DEP (limma) stats
stats_df <- swip_dep %>% 
	select(Fraction,Accession,p.val) %>% 
	left_join(tmpdf,by=c("Fraction","Accession"))

stats_df <- fix_colname(stats_df,"PValue","EdgeR")
stats_df <- fix_colname(stats_df,"p.val","DEP")

# the spearman rank correlation is very high
rho <- cor(x=stats_df$EdgeR, y=stats_df$DEP, method="spearman")
knitr::kable(as.data.table(setNames(list(rho),nm="Spearman R")))

# FIXME: what about foldchange?

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
