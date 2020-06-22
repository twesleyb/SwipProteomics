#!/usr/bin/env Rscript

#' ---
#' title: Swip Proteomics
#' description:
#' authors: Tyler W Bradshaw
#' ---

## INPUTS:

## OPTIONS:

## OUTPUTS:

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. 

# Load renv -- use renv::load NOT activate!
rootdir <- getrd()
renv::load(rootdir,quiet=TRUE) # NOTE: getrd is a f(x) in .Rprofile.

# Load required packages and functions.
suppressPackageStartupMessages({
	library(dplyr) # For manipulating data.
	library(ggplot2) # For making plots.
	library(data.table) # For working with tables.
})

# Load project specific functions and data.
suppressWarnings({ devtools::load_all() })

# Project directories:
fontdir <- file.path(rootdir, "fonts") # Arial font for plots.
figsdir <- file.path(rootdir, "figs","TMT") # Output figures.

# Set global plotting settings.
ggtheme()
set_font("Arial", font_path = fontdir)

#---------------------------------------------------------------------
## Load the data.
#---------------------------------------------------------------------

#----------------------------------------------------------------------
## Plot commonly dysregulated prots, adjusted for fraction differences.
#----------------------------------------------------------------------

# Calculate protein abundance, adjusted for fraction differences.
dm <- tmt_protein

logCPM <- edgeR::cpm(dge, log=TRUE)

# Remove effect of fraction.
dm <- limma::removeBatchEffect(logCPM,batch=dge$samples$fraction,
				   design=model.matrix(~treatment,
						       data=dge$samples))

# Collect results.
dt <- as.data.table(dm,keep.rownames="Accession") %>% 
	reshape2::melt(id.var="Accession",
		       variable.name="Sample",
		       value.name="Intensity")
dt$Sample <- as.character(dt$Sample)

# Combine with additional meta data.
idy <- which(colnames(tmt_protein)=="Intensity")
temp_prot <- as.data.table(tmt_protein)[,-..idy] %>% filter(Treatment != "SPQC")
adjusted_prot <- left_join(temp_prot,dt, by=c("Accession","Sample"))

# Wash proteins + control.
prots <- c("Washc4","Washc1","Washc2","Washc5","Tubb4a")
names(prots) <- gene_map$uniprot[match(prots,gene_map$symbol)]

# Combine adjusted data and stats, subset to keep proteins of interest.
df <- adjusted_prot %>% filter(Accession %in% names(prots))

# Colors for WT and mutant groups.
colors = c(Control="#47b2a4",Mutant="#B86FAD")

# Labels will simply be WT and Mutant.
xlabels <- rep(c("WT","Mutant"),times=length(prots))

# Order of the factors.
factor_order <- paste(rep(names(prots),each=2), c("Control","Mutant"),sep=".")
df$Accession.Treatment <- as.character(interaction(df$Accession,df$Treatment))
df$Accession.Treatment <- factor(df$Accession.Treatment,levels=factor_order)

# Generate a plot.
plot <- ggplot(df, aes(x=Accession.Treatment, y=Intensity,
		       fill=Treatment)) + 
	geom_boxplot() + geom_point(aes(fill=Treatment,shape=Treatment)) +
	theme(axis.text.x=element_text(angle=45))
plot <- plot + scale_fill_manual(name="Genotype",values=colors)
plot <- plot + theme(legend.position="none")
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(panel.border=element_rect(colour="black",fill="NA",size=1))
plot <- plot + theme(axis.title.x = element_blank())
plot <- plot + scale_x_discrete(labels=xlabels)
plot <- plot + ylab("log2(Adjusted Intensity)")

# Add some lines to break up the data.
plot <- plot + geom_vline(xintercept=seq(2.5,length(proteins)*2,by=2),
			  linetype="dotted",size=0.5)

# Perform t-tests.
# FIXME: Substitute with glm pvalues.
prot_list <- df %>% group_by(Accession) %>% group_split()
data_ttests <- lapply(prot_list,function(subdf) {
			      x <- subdf$Intensity
			      y <- subdf$Treatment
			      result <- t.test(x~y,paired=FALSE,
					       alternative="greater")
			      ttest_dt <- as.data.table(t(unlist(result)))
			      return(ttest_dt)
			  })
names(data_ttests) <- names(prots)
ttest_dt <- bind_rows(data_ttests,.id="Accession")
ttest_dt$p.adjust <- p.adjust(ttest_dt$p.value,method="bonferroni")

# Annotate the plot with stats.
stats <- df %>% filter(Treatment=="Control") %>% 
	group_by(Accession.Treatment) %>% 
	dplyr::summarize(ypos = 1.02*max(Intensity),.groups="drop")
stats$Accession.Treatment <- as.character(stats$Accession.Treatment)
stats$Accession <- sapply(strsplit(stats$Accession.Treatment,"\\."),"[",1)
stats <- stats %>% dplyr::select(Accession.Treatment,Accession,ypos)
stats$xpos <- seq(1.5,by=2,length.out=5)
stats$symbol <- "ns"
stats <- left_join(stats,ttest_dt,by="Accession")
stats$symbol[stats$p.adjust < 0.05] <- "*"
stats$symbol[stats$p.adjust < 0.005] <- "**"
stats$symbol[stats$p.adjust < 0.0005] <- "***"

# Add significance stars.
plot <- plot + 
	annotate("text",x=stats$xpos,y=stats$ypos,label=stats$symbol,size=4)

# Annotate with protein names.
symbols <- prots 
build <- ggplot_build(plot)
ymax <- build$layout$panel_params[[1]][["y.range"]][2]
plot <- plot + annotate("text",x=seq(1.5,length(prots)*2,by=2),
			y=ymax,label=symbols,size=5)

# Save.
if (save_plots) {
	myfile <- file.path(figsdir,
			    "Select_Proteins_Adjusted_Abundance.pdf")
	ggsave(myfile,plot, height = fig_height, width = fig_width)
}
