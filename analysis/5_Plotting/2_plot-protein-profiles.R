#!/usr/bin/env Rscript

# title: SwipProteomics
# description: plot protein abundance
# author: Tyler W Bradshaw

# plot median Abundance

## ---- Inputs

# input data in root/data/
root = "~/projects/SwipProteomics"

input_part = "ne_surprise_partition"
input_colors = "ne_surprise_colors"


## ---- Prepare environment 

renv::load(root,quiet=TRUE)

# library(SwipProteomics)
devtools::load_all(root, quiet=TRUE)

# load the data
data(swip)
data(gene_map)
data(protein_gof)
data(msstats_prot)
data(msstats_results)

data(list=input_part) # partition
data(list=input_colors) # module colors

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
	dir.create(figsdir, recursive=TRUE)
}

# set plotting theme
ggtheme()
set_font("Arial",font_path=fontdir)


## ---- prepare the data

# all protein modules
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))

# generate named vector of protein module colors
#color_list <- sapply(names(module_colors),function(x) {
#		       setNames(rep(module_colors[x],length(modules[[x]])),
#				nm=modules[[x]]) })
#prot_colors <- unlist(color_list)
#names(prot_colors) <- sapply(strsplit(names(prot_colors),"\\."),"[",2)

# drop 'Mutant-Control' contrasts from results and then
# combine msstats_results with msstats_prot
filt_results <- msstats_results %>% filter(Contrast != "Mutant-Control") %>%
	mutate(BioFraction = sapply(strsplit(Contrast, "\\."),"[",3))
shared_cols <- intersect(colnames(msstats_prot),colnames(filt_results))
prot_df <- left_join(msstats_prot,filt_results,by=shared_cols)

# annotate with module membership
prot_df <- prot_df %>% filter(Protein %in% names(partition)) %>%
	mutate(Module = paste0("M",partition[Protein])) 


## ---- Function 

# a function to generate protein profile plot:

protein = sample(unique(msstats_prot$Protein),1)

plotProfile <- function(protein, gene_map, msstats_prot, msstats_results, 
			 protein_gof) {

  wt_color <- "#47b2a4"
  mut_color <- "#b671af"

  gene <- gene_map$symbol[match(protein,gene_map$uniprot)]

  # proteins with overall sig change
  sig_prots <- msstats_results %>% 
	filter(Contrast == 'Mutant-Control' & FDR<0.05) %>%  ungroup() %>%
	select(Protein) %>% unlist() %>% as.character() %>% unique()

  # plot title and color
  r2 <- protein_gof %>% filter(Protein == protein) %>% 
	  select(R2.fixef) %>% as.numeric()
  title_anno <- paste0("(R2_fixef = ",round(r2,3),")")

  title_colors <- c("darkred"=TRUE,"black"=FALSE)
  title_color <- names(which(title_colors==protein %in% sig_prots))

  # prepare stats for labeling
  stats_df <- prot_df %>% filter(Protein == protein) %>% 
  	select(Protein,BioFraction,FDR) %>% unique()
  stats_df$stars <- ""
  stats_df$stars[stats_df$FDR < 0.05] <- "*"
  stats_df$stars[stats_df$FDR < 0.005] <- "**"
  stats_df$stars[stats_df$FDR < 0.0005] <- "***"

  # prepare the data 
  df <- prot_df %>% 
    filter(Protein == protein) %>% 
    group_by(Condition, Protein) %>%
    summarize(med_Abundance = median(Abundance),
	      ymin = min(Abundance),
	      ymax = max(Abundance),
	      .groups="drop") 

    # munge to annotate with Genotype and BioFraction
    condition <- as.character(df$Condition)
    df$Genotype <- factor(sapply(strsplit(condition,"\\."),"[",1),
		        levels=c("Control","Mutant"))
    df$BioFraction <- factor(sapply(strsplit(condition,"\\."),"[",2),
		        levels=c("F4","F5","F6","F7","F8","F9","F10"))

    # generate the plot
    plot <- ggplot(df)
    plot <- plot + aes(x = BioFraction)
    plot <- plot + aes(y = med_Abundance)
    plot <- plot + aes(group = Genotype)
    plot <- plot + aes(colour = Genotype)
    plot <- plot + aes(shape = Genotype)
    plot <- plot + aes(fill = Genotype)
    plot <- plot + aes(shade = Genotype)
    plot <- plot + aes(ymin=ymin)
    plot <- plot + aes(ymax=ymax)
    plot <- plot + geom_line()
    plot <- plot + geom_ribbon(alpha=0.1, linetype="blank")
    plot <- plot + geom_point(size=2)
    plot <- plot + ggtitle(paste(gene," | ",protein)) # title_anno?
    plot <- plot + ylab("log2(Protein Intensity)")
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
    plot <- plot + scale_colour_manual(values=c(wt_color,mut_color))
    plot <- plot + scale_fill_manual(values=c(wt_color,mut_color))
    plot <- plot + theme(legend.position = "none")

    # annotate with module membership and significance stars
    yrange <- setNames(range(df$med_Abundance),nm=c("min","max"))
    ypos1 <- yrange["min"]  + 0.05 * diff(yrange) # Module#
    ypos2 <- yrange["max"] + 0.05 * diff(yrange) # stars
    module_label <- paste0("M",partition[protein])
    plot <- plot + annotate(geom="label",x=7.25,y=ypos1,label=module_label)
    plot <- plot + annotate("text",x=stats_df$BioFraction,y=ypos2,
		  label=stats_df$stars,size=7)

  return(plot)
} #EOF


## ---- generate plots 

# protein names sorted by module membership
sorted_prots <- as.character(unlist(split(names(partition),partition)[-1]))

message("\nGenerating plots for ",
	formatC(length(sorted_prots),big.mark=",")," proteins.")

# loop for all proteins
plots <- list()
pbar <- txtProgressBar(max=length(sorted_prots),style=3)
for (protein in sorted_prots) {
	plots[[protein]] <- plotProfile(protein, gene_map, msstats_prot, 
				         msstats_results, protein_gof)
        setTxtProgressBar(pbar,value=match(protein,sorted_prots))
} # EOL for proteins
close(pbar)


## ---- save results 

message("\nSaving plots as a single pdf, this will take several minutes.")
myfile = file.path(figsdir,"protein_profiles.pdf")
ggsavePDF(plots, file=myfile)
