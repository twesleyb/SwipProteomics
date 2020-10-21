#!/usr/bin/env Rscript

# title: 
# description: plot protein abundance
# authors: Tyler W Bradshaw

## Inputs --------------------------------------------------------------------

# options:

# Input data in root/data/
root = "~/projects/SwipProteomics"


## Prepare environment --------------------------------------------------------

renv::load(root,quiet=TRUE)
suppressWarnings({ devtools::load_all() })

# load the data
data(swip)
data(gene_map)
data(msstats_prot)
data(module_colors)
data(leidenalg_partition)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

# project dirs
fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs", "Modules","Groups")
if (! dir.exists(figsdir)) {
	dir.create(figsdir,recursive=TRUE)
}

# set plotting theme and font
ggtheme(); set_font("Arial",font_path=fontdir)

# wash complex proteins
washc_prots = gene_map$uniprot[grep("Washc*",gene_map$symbol)]


## plot protein summary --------------------------------------------------------

data(fit1_results)
data(fit0_results)

plot_protein_summary <- function(protein,fit0_results=NULL) {
  # a function to generate a proteins summary plot
  # requires SwipProteomics for col2hex and data
  require(dplyr,quietly=TRUE)
  require(ggplot2,quietly=TRUE)
  # defaults
  wt_color = "#47b2a4"
  mut_color <- col2hex(c("R"=148,"G"=33,"B"=146))
  # prepare the data
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
#partition <- setNames(rep(1,length(washc_prots)),nm=washc_prots)
#modules <- split(names(partition),partition)
#names(modules) <- paste0("M",names(modules))

modules <- split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))

# drop modules that are less than min_size
min_size = 5
too_small <- names(which(sapply(modules,length) < min_size))
modules <- modules[names(modules) %notin% too_small]


for (module in names(modules)) {

  # Get the modules color
  wt_color = "#47b2a4"
  #mut_color <- col2hex("purple")
  mut_color <- module_colors[module]

  # 1. Generate all the plots for a given module using plot_protein_summary
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

  # 2. combine the data containing the normalized data 
  # (each protein scaled to maximum within a plot) into a single df
  # and generate combined plot
  df <- bind_rows(sapply(plots,"[", "data"))

  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction)
  plot <- plot + aes(y = norm_Abundance)
  plot <- plot + aes(group = interaction(Protein,Genotype), colour = Genotype, 
		     shape=Genotype, fill=Genotype,shade=Genotype)
  plot <- plot + aes(ymin=norm_Abundance - CV)
  plot <- plot + aes(ymax=norm_Abundance + CV)
  plot <- plot + geom_line()
  plot <- plot + geom_ribbon(alpha=0.1, linetype="blank")
  plot <- plot + geom_point(size=2)
  plot <- plot + ggtitle(paste0(module," (n = ",length(proteins),")"))
  plot <- plot + ylab("Normalized Protein Abundance")
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

  # save
  myfile = file.path(figsdir, paste0(module,"_all_Proteins.pdf"))
  ggsavePDF(plot,file=myfile)
  message("\nSaved: ", module)
} # EOL for modules


###############################################################################

# our task is to summarize a module. should this be done for each fraction or
# each mixture?

# proteins = washc_prots
subdat <- msstats_prot %>% filter(Protein %in% proteins)

data_list <- subdat %>% group_by(BioFraction) %>% group_split()
# length == 7

df = data_list[[1]]

#fx = lm(Abundance ~ BioFraction.Genotype + Run)
myfile <- file.path(root,"rdata","msstats_psm.rda")
load(myfile) # msstats_psm
data(leidenalg_partition) # partition

# collect the data
msstats_psm <- msstats_psm %>% mutate(Module = partition[ProteinName])
msstats_psm <- msstats_psm %>% filter(!is.na(Module))
msstats_psm <- msstats_psm %>% filter(!is.na(Intensity))
# now, for each run...
#msstats_psm <- msstats_psm %>% mutate(Feature = paste(Mixture,Channel,Module,ProteinName,PeptideSequence, collapse="_"))


#################################################
## another attempt...

# how to summarize module abundance...
df = plot$data
df = df %>% group_by(Genotype,BioFraction) %>% mutate(Fit = median(norm_Abundance))

fm <- lm(norm_Abundance ~ 0 + Condition + Protein, data = plot$data)
cf <- as.data.table(coef(fm),keep.rownames="Condition")
colnames(cf) <- c("Condition", "Fit")
cf <- cf %>% mutate(Condition = gsub("Condition","", Condition))
cf <- cf %>% filter(Condition %in% plot$data$Condition)

# combine fit data and plot data
df = left_join(plot$data,cf,by="Condition")

  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction)
  plot <- plot + aes(y = norm_Abundance)
  plot <- plot + aes(group = interaction(Protein,Genotype), colour = Genotype, 
		     shape=Genotype, fill=Genotype,shade=Genotype)
  plot <- plot + aes(ymin=norm_Abundance - CV)
  plot <- plot + aes(ymax=norm_Abundance + CV)
  plot <- plot + geom_line()
  plot <- plot + geom_ribbon(alpha=0.1, linetype="blank")
  plot <- plot + geom_point(size=2)
  plot <- plot + geom_line(aes(x=BioFraction,y=Fit),size=2,colour="red")


  plot <- plot + ggtitle(paste0(module," (n = ",length(proteins),")"))
  plot <- plot + ylab("Normalized Protein Abundance")
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

