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
data(protein_gof)
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
sub_results <- msstats_results %>% filter(Contrast != "Mutant-Control")
biofraction <- sapply(strsplit(sub_results$Contrast, "\\."),"[",3)
sub_results$BioFraction <- biofraction
shared_cols <- intersect(colnames(sub_results),colnames(msstats_prot))
prot_df <- left_join(msstats_prot,sub_results,by=shared_cols)

# annotate with module membership
prot_df <- prot_df %>% filter(Protein %in% names(partition)) %>%
	mutate(Module = paste0("M",partition[Protein])) 

# proteins with overall sig change
sigprots <- msstats_results %>% 
	filter(Contrast == 'Mutant-Control' & FDR<0.05) %>% 
	select(Protein) %>% unlist() %>% as.character() %>% unique()


## Function ------------------------------------------------------------------

# a function to generate protein profile plot:
plot_profile <- function(prot_df, protein, sigprots, protein_gof,
			 gene = gene_map$symbol[match(protein,gene_map$uniprot)],
			 wt_color = "#47b2a4",
			 mut_color = prot_colors[protein]) {
	r2 <- protein_gof %>% filter(Protein == protein) %>% 
		select(R2.total) %>% as.numeric()
	title_anno <- paste0("(R2 = ",round(r2,3),")")
	title_colors <- c("darkred"=TRUE,"black"=FALSE)
	title_color <- names(which(title_colors==protein %in% sigprots))
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
  summarize(mean_Abundance = mean(norm_Abundance), # using batch corrected norm_Abundance!
	    SD = sd(norm_Abundance),
	    SE = unique(SE),
	    N = length(Abundance),
	    .groups="drop") %>%
  # calculate coefficient of variation
  mutate(CV = SD/mean_Abundance) %>%
  # scale profile
  mutate(scale_Abundance = mean_Abundance/max(mean_Abundance))
  # munge to annotate with Genotype and BioFraction
  condition <- as.character(df$Condition)
  df$Genotype <- factor(sapply(strsplit(condition,"\\."),"[",1),
		        levels=c("Control","Mutant"))
  df$BioFraction <- factor(sapply(strsplit(condition,"\\."),"[",2),
		        levels=c("F4","F5","F6","F7","F8","F9","F10"))
  # Generate the plot
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction)
  plot <- plot + aes(y = scale_Abundance)
  plot <- plot + aes(group = Genotype)
  plot <- plot + aes(colour = Genotype)
  plot <- plot + aes(shape = Genotype)
  plot <- plot + aes(fill = Genotype)
  plot <- plot + aes(shade = Genotype)
  plot <- plot + aes(ymin=scale_Abundance - CV)
  plot <- plot + aes(ymax=scale_Abundance + CV)
  plot <- plot + geom_line()
  plot <- plot + geom_ribbon(alpha=0.1, linetype="blank")
  plot <- plot + geom_point(size=2)
  plot <- plot + ggtitle(paste(gene," | ",protein,title_anno))
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
  plot <- plot + theme(plot.title = element_text(color=title_color))
  plot <- plot + theme(axis.line.y=element_line())
  plot <- plot + scale_colour_manual(values=c(wt_color,prot_colors[protein]))
  plot <- plot + scale_fill_manual(values=c(wt_color,prot_colors[protein]))
  plot <- plot + theme(legend.position = "none")
  # annotate with module membership and significance stars
  yrange <- setNames(range(df$scale_Abundance),nm=c("min","max"))
  ypos1 <- yrange["min"]  + 0.05 * diff(yrange) # Module#
  ypos2 <- yrange["max"] + 0.05 * diff(yrange) # stars
  module_label <- paste0("M",partition[protein])
  plot <- plot + annotate(geom="label",x=7.25,y=ypos1,label=module_label)
  plot <- plot + annotate("text",x=stats_df$BioFraction,y=ypos2,
		label=stats_df$stars,size=7)
  return(plot)
} #EOF


#plot_profile(prot_df,swip,sigprots,protein_gof)

## generate plots -------------------------------------------------------------

# protein names sorted by module membership
sorted_prots <- as.character(unlist(split(names(partition),partition)[-1]))

# loop for all proteins
plots <- list()
pbar <- txtProgressBar(max=length(sorted_prots),style=3)
for (protein in sorted_prots) {
        plots[[protein]] <- plot_profile(prot_df,protein,sigprots,protein_gof)
        setTxtProgressBar(pbar,value=match(protein,sorted_prots))
} # EOL for proteins
close(pbar)

# save
myfile = file.path(figsdir,"Protein_profiles.pdf")
ggsavePDF(plots, file=myfile)
