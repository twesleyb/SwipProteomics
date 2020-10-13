#!/usr/bin/env Rscript

# author: twab
# title: SwipProteomics
# description: examining lmer models

# references:

# * Hlavac, Marek (2018). stargazer: Well-Formatted Regression and Summary 
#   Statistics Tables. R package version 5.2.2. 
#   https://CRAN.R-project.org/package=stargazer 

# * MuMin package

# * Nakagawa, S., Schielzeth, H. (2013) A general and simple method
#   for obtaining R<U+00B2> from Generalized Linear Mixed-effects Models. Methods
#   in Ecology and Evolution 4: 133<U+2013>

## inputs 
root = "~/projects/SwipProteomics"


## functions ------------------------------------------------------------------

subprot <- function(protein,required=c("msstats_prot")) {
	# subset the protein-level data for a given protein
	# FIXME: error message about env is not informative
	#stopifnot(all(sapply(required,exists)))
	require(dplyr, quietly=TRUE)
	return(msstats_prot %>% dplyr::filter(Protein == protein))
} #EOF

submod <- function(module,required=c("msstats_prot")) {
	# subset the module-level data for a given module
	# FIXME: error message about env is not informative
	require(dplyr, quietly=TRUE)
	#stopifnot(all(sapply(required,exists)))
	stopifnot("Module" %in% colnames(msstats_prot))
	return(msstats_prot %>% dplyr::filter(Module == module))
} #EOF


## Prepare environment --------------------------------------------------------

# load renv
if (dir.exists(file.path(root,"renv"))) { renv::load(root) }

# load project
devtools::load_all(quiet=TRUE)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

# dir for output
figsdir <- file.path(root,"figs","Modules", "Groups")
if (!dir.exists(figsdir)) { dir.create(figsdir); message("mkdir ",figsdir) }

# set font and ggplot theme
fontdir <- file.path(root, "fonts")
ggtheme(); set_font("Arial",font_path=fontdir)

# load the data
data(swip)
data(gene_map)
data(partition)
data(sig_modules)
data(msstats_prot)
data(module_colors)
data(msstats_results)


## plot protein summary --------------------------------------------------------

plot_protein_summary <- function(protein) {
  # a function to generate a proteins summary plot
  # requires SwipProteomics for col2hex and data
  require(dplyr,quietly=TRUE)
  require(ggplot2,quietly=TRUE)
  # defaults
  wt_color = "#47b2a4"
  mut_color <- col2hex(c("R"=148,"G"=33,"B"=146))
  # prepare the data
  msstats_df <- left_join(msstats_prot,msstats_results,by=c("Protein","BioFraction"))
  msstats_df$Module <- paste0("M",partition[msstats_df$Protein])
  gene <- gene_map$symbol[match(protein,gene_map$uniprot)]
  df <- subset(msstats_df,msstats_df$Protein == protein) %>% 
	select(Mixture,Channel,Condition,Protein,Abundance,SE,Module) %>% 
	unique() %>% 
	group_by(Condition,Protein) %>% 
	summarize(mean_Abundance = mean(Abundance),
		  SD = sd(Abundance),
		  SE = unique(SE),.groups="drop")
  df <- df %>% mutate(CV = SD/mean_Abundance)
  df <- df %>% mutate(norm_Abundance = mean_Abundance/max(mean_Abundance))
  df$BioFraction <- sapply(strsplit(as.character(df$Condition),"\\."),"[",2)
  df$Genotype <- sapply(strsplit(as.character(df$Condition),"\\."),"[",1)
  df$BioFraction <- factor(df$BioFraction,
			   levels=c("F4","F5","F6","F7","F8","F9","F10"))
  # Generate the plot.
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction)
  plot <- plot + aes(y = norm_Abundance)
  plot <- plot + aes(group = Genotype, colour = Genotype, 
		     shape=Genotype, fill=Genotype,shade=Genotype)
  plot <- plot + aes(ymin=norm_Abundance - CV)
  plot <- plot + aes(ymax=norm_Abundance + CV)
  plot <- plot + geom_line()
  plot <- plot + geom_ribbon(alpha=0.1, linetype="blank")
  plot <- plot + geom_point(size=2)
  plot <- plot + ggtitle(paste(gene,protein,sep=" | "))
  plot <- plot + ylab("Normalized Abundance")
  plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))
  plot <- plot + theme(axis.text.x = element_text(color="black", size=11, 
						  angle = 0, hjust = 1, 
						  family = "Arial"))
  plot <- plot + theme(axis.text.y = element_text(color="black",size=11, 
						  angle = 0, hjust = 1, 
						  family = "Arial"))
  plot <- plot + theme(panel.background = element_blank())
  plot <- plot + theme(axis.line.x=element_line())
  plot <- plot + theme(axis.line.y=element_line())
  plot <- plot + scale_colour_manual(values=c(wt_color,mut_color))
  plot <- plot + scale_fill_manual(values=c(wt_color,mut_color))
  plot <- plot + theme(legend.position = "none")
  return(plot)
} #EOF


## main ------------------------------------------------------------------------

# Loop to do work.
modules <- split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))

for (module in names(modules)) {
  # Get the modules color
  wt_color <- "#47b2a4" # teal
  mut_color <- module_colors[module]
  # 1. Generate all the protein plots for a given module using 
  #    plot_protein_summary
  plots <- list()
  message("\nWorking on: ", module)
  proteins <- modules[[module]]
  pbar <- txtProgressBar(max=length(proteins),style=3)
  for (protein in proteins) {
    # FIXME: colors do not seem correct
    plot <- plot_protein_summary(protein)
    plot_label <- paste("Module:", partition[protein])
    yrange <- plot$data %>% dplyr::filter(Protein == protein) %>% 
  	    select(norm_Abundance) %>% range()
    ypos <- yrange[1] - 0.1* diff(yrange)
    plot <- plot + annotate(geom="label",x=7, y=ypos, label=plot_label)
    plots[[protein]] <- plot
    setTxtProgressBar(pbar,value=match(protein,proteins))
  } # EOL for proteins
  close(pbar)
  # 2. combine the data containing the normalized data  for all proteins
  # (each protein scaled to maximum within a plot) into a single df
  df <- bind_rows(sapply(plots,"[", "data"))
  fx <- "norm_Abundance ~ 0 + (1|BioFraction) + (1|Protein) + Genotype"
  fm <- lmerTest::lmer(fx,df)
  # goodness of fit
  #rho = r.squaredGLMM.merMod(fm) # R2marginal (fixef) and R2conditional (total)
  # 3. get fitted values for all proteins, and then
  # collect medians at each time point --> this will be plot as a summary
  # of the module
  df$fit_Abundance <- fitted(fm)
  df <- df %>% group_by(Condition,BioFraction) %>% 
	  mutate(best_fit=median(fit_Abundance))
  # 4. Generate combined plot with 'best fit' line
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction)
  plot <- plot + aes(y = norm_Abundance)
  plot <- plot + aes(group = interaction(Protein,Genotype))
  plot <- plot + aes(colour = Genotype)
  plot <- plot + aes(fill=Genotype, shade=Genotype)
  plot <- plot + aes(ymin=norm_Abundance - CV)
  plot <- plot + aes(ymax=norm_Abundance + CV)
  plot <- plot + geom_ribbon(alpha=0.09, linetype="blank")
  plot <- plot + geom_line(aes(y=best_fit),size=1.0,linetype="solid")
  plot <- plot + ggtitle(paste0(module," (n = ",length(proteins),")"))
  plot <- plot + ylab("Normalized Protein Abundance")
  plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))
  plot <- plot + theme(axis.text.x = element_text(color="black"))
  plot <- plot + theme(axis.text.x = element_text(size=11))
  plot <- plot + theme(axis.text.x = element_text(angle=0))
  plot <- plot + theme(axis.text.x = element_text(hjust=0))
  plot <- plot + theme(axis.text.x = element_text(family="Arial"))
  plot <- plot + theme(axis.text.y = element_text(color="black"))
  plot <- plot + theme(axis.text.y = element_text(size=11))
  plot <- plot + theme(axis.text.y = element_text(angle=0))
  plot <- plot + theme(axis.text.y = element_text(hjust=0))
  plot <- plot + theme(axis.text.y = element_text(family="Arial"))
  plot <- plot + theme(panel.background = element_blank())
  plot <- plot + theme(axis.line.x=element_line())
  plot <- plot + theme(axis.line.y=element_line())
  plot <- plot + scale_colour_manual(values=c(wt_color,mut_color))
  plot <- plot + scale_fill_manual(values=c(wt_color,mut_color))
  plot <- plot + theme(legend.position = "none")
  plot <- plot + scale_x_discrete(expand = c(0.05,0))
  # save
  myfile = file.path(figsdir, paste0(module,"_all_Protein_Summary.pdf"))
  ggsavePDF(plot,file=myfile)
  message("\nSaved: ", module)
} # EOL for modules
