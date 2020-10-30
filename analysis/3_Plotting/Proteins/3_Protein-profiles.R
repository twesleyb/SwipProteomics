#!/usr/bin/env Rscript

# title: 
# description: plot protein abundance
# authors: Tyler W Bradshaw

## Inputs --------------------------------------------------------------------

# options:
#mut purple: col2hex(c("R"=148,"G"=33,"B"=146))

# Input data in root/data/
root = "~/projects/SwipProteomics"

## Prepare environment --------------------------------------------------------

renv::load(root,quiet=TRUE)
suppressWarnings({ devtools::load_all() })

# load the data
data(swip)
data(gene_map)
data(partition)
data(msstats_prot)
data(module_colors)
data(msstats_results)

# all protein modules
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))

# generate named vector of protein module colors
color_list <- sapply(names(module_colors),function(x) {
		       setNames(rep(module_colors[x],length(modules[[x]])),
				nm=modules[[x]]) })
prot_colors <- unlist(color_list)
names(prot_colors) <- sapply(strsplit(names(prot_colors),"\\."),"[",2)

# other imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

# project dirs
fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs","Proteins")

if (! dir.exists(figsdir)) {
	dir.create(figsdir,recursive=TRUE)
}

ggtheme(); set_font("Arial",font_path=fontdir)


## prepare the data -----------------------------------------------------------

# combine msstats_results (stats) and msstats_prot (normalized protein data)
msstats_results <- msstats_results %>% filter(Contrast != "Mutant-Control")
biofraction <- sapply(strsplit(msstats_results$Contrast, "\\."),"[",3)
msstats_results$BioFraction <- biofraction
shared_cols <- intersect(colnames(msstats_results),colnames(msstats_prot))
prot_df <- left_join(msstats_prot,msstats_results,by=shared_cols)

# annotate with module membership
prot_df <- prot_df %>% filter(Protein %in% names(partition)) %>%
	mutate(Module = paste0("M",partition[Protein])) 

# this is a nice tidy table of the final data!
#swip_tmt <- prot_df
#myfile <- file.path(root,"rdata","swip_tmt.rda")
#save(swip_tmt,file=myfile,version=2)


## Function ------------------------------------------------------------------

# a function to generate protein profile plot:
plot_profile <- function(prot_df, protein,
			 gene = gene_map$symbol[match(protein,gene_map$uniprot)],
			 wt_color = "#47b2a4",
			 mut_color = prot_colors[protein]) {
  # prepare stats for labeling
  stats_df <- prot_df %>% filter(Protein == protein) %>% 
  	select(Protein,BioFraction,FDR) %>% unique()
  stats_df$stars <- ""
  stats_df$stars[stats_df$FDR < 0.05] <- "*"
  stats_df$stars[stats_df$FDR < 0.005] <- "**"
  stats_df$stars[stats_df$FDR < 0.0005] <- "***"
  # prepare the data 
  df <- prot_df %>% 
  # subset -- keep data from single protein
  filter(Protein == protein) %>% 
  # group by Genotype.BioFraction.Protein
  group_by(Condition,Protein) %>%
  # calculate the average of the three mixtures
  summarize(mean_Abundance = mean(Abundance),
	    SD = sd(Abundance),
	    SE = unique(SE),
	    N = length(Abundance),
	    .groups="drop") %>%
  # calculate coefficient of variation
  mutate(CV = SD/mean_Abundance) %>%
  # scale profile
  mutate(norm_Abundance = mean_Abundance/max(mean_Abundance))
  # munge to annotate with Genotype and BioFraction
  condition <- as.character(df$Condition)
  df$Genotype <- factor(sapply(strsplit(condition,"\\."),"[",1),
		        levels=c("Control","Mutant"))
  df$BioFraction <- factor(sapply(strsplit(condition,"\\."),"[",2),
		        levels=c("F4","F5","F6","F7","F8","F9","F10"))
  # Generate the plot
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction)
  plot <- plot + aes(y = mean_Abundance)
  plot <- plot + aes(group = Genotype)
  plot <- plot + aes(colour = Genotype)
  plot <- plot + aes(shape = Genotype)
  plot <- plot + aes(fill = Genotype)
  plot <- plot + aes(shade = Genotype)
  #plot <- plot + aes(ymin=norm_Abundance - CV)
  #plot <- plot + aes(ymax=norm_Abundance + CV)
  plot <- plot + aes(ymin=mean_Abundance - mean_Abundance * CV)
  plot <- plot + aes(ymax=mean_Abundance + mean_Abundance * CV)
  plot <- plot + geom_line()
  plot <- plot + geom_ribbon(alpha=0.1, linetype="blank")
  plot <- plot + geom_point(size=2)
  plot <- plot + ggtitle(paste(gene,protein,sep=" | "))
  plot <- plot + ylab("Scaled Protein Abundance")
  plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))
  plot <- plot + theme(axis.text.x = element_text(color="black", size=11))
  plot <- plot + theme(axis.text.x = element_text(angle = 0, hjust = 1)) 
  plot <- plot + theme(axis.text.x = element_text(family = "Arial"))
  plot <- plot + theme(axis.text.y = element_text(color="black", size=11))
  plot <- plot + theme(axis.text.y = element_text(angle = 0, hjust = 1)) 
  plot <- plot + theme(axis.text.y = element_text(family = "Arial"))
  plot <- plot + theme(panel.background = element_blank())
  plot <- plot + theme(axis.line.x=element_line())
  plot <- plot + theme(axis.line.y=element_line())
  plot <- plot + scale_colour_manual(values=c(wt_color,prot_colors[protein]))
  plot <- plot + scale_fill_manual(values=c(wt_color,prot_colors[protein]))
  plot <- plot + theme(legend.position = "none")
  # annotate with module membership and significance stars
  yrange <- setNames(range(df$mean_Abundance),nm=c("min","max"))
  ypos1 <- yrange["min"]  + 0.05 * diff(yrange) # Module#
  ypos2 <- yrange["max"] + 0.05 * diff(yrange) # stars
  module_label <- paste0("M",partition[protein])
  plot <- plot + annotate(geom="label",x=7.25,y=ypos1,label=module_label)
  plot <- plot + annotate("text",x=stats_df$BioFraction,y=ypos2,
		label=stats_df$stars,size=7)
  return(plot)
} #EOF


## get swip's fit

# save plots for a subset of proteins indivually
washc_prots <- gene_map$uniprot[grepl("Washc*",gene_map$symbol)]

for (protein in washc_prots) {
  gene = gene_map$symbol[match(protein,gene_map$uniprot)]
  wt_color = "#47b2a4"

  plot <- plot_profile(prot_df,protein)
  fx0 <- formula("Abundance ~ 0 + Condition + (1|Mixture)")
  fm0 <- lmerTest::lmer(fx0,msstats_prot %>% filter(Protein == protein))
  fx1 <- formula("Abundance ~ 0 + Genotype + BioFraction + (1|Mixture)")
  fm1 <- lmerTest::lmer(fx1,msstats_prot %>% filter(Protein == protein))
  wt_yint <- lme4::fixef(fm1)["GenotypeControl"]
  mut_yint <- lme4::fixef(fm1)["GenotypeMutant"]
  plot <- plot + geom_hline(yintercept = wt_yint,linetype="dashed",color=wt_color)
  plot <- plot + geom_hline(yintercept = mut_yint,linetype="dashed",color=prot_colors[protein])

  (summary(fm1,ddf="Satterthwait"))
  myfile <- file.path(figsdir,paste(protein,gene,"profile.pdf",sep="_"))
  ggsave(myfile,plot)
} # EOL


quit()
stop()

## generate plots -------------------------------------------------------------

# protein names sorted by module membership
sorted_prots <- as.character(unlist(split(names(partition),partition)[-1]))

# loop for all proteins
plots <- list()
pbar <- txtProgressBar(max=length(sorted_prots),style=3)
for (protein in sorted_prots) {
        plots[[protein]] <- plot_profile(prot_df,protein)
        setTxtProgressBar(pbar,value=match(protein,sorted_prots))
} # EOL for proteins
close(pbar)

# save
myfile = file.path(figsdir,"Protein_profiles.pdf")
ggsavePDF(plots, file=myfile)
