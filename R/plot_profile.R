#' plot_profile

plot_profile <- function(protein, wt_color = "#47b2a4",
				 mut_color = col2hex(c("R"=148,"G"=33,"B"=146)),
				 scale_abundance = FALSE) {
  # assumes msstats_prot and msstats_results are available in the current env
  # * gene_map
  # * msstats_prot
  # * msstats_results
  # require(dplyr)
  # require(ggplot2)
  # get protein's gene symbol
  gene <- gene_map$symbol[match(protein,gene_map$uniprot)]
  # annotate msstats_results with BioFraction
  biofraction <- sapply(strsplit(msstats_results$Label,"\\."),"[",3)
  prot_df <- msstats_prot
  results_df <- msstats_results
  results_df$BioFraction <- factor(biofraction,levels=c("F4","F5","F6","F7","F8","F9","F10"))
  # prepare the data 
  df <- prot_df %>% left_join(results_df,by=c("Protein","BioFraction")) %>%
	  # subset -- only keep clusterd proteins
	  filter(Protein %in% names(partition)) %>%  
	  # annotate with module membership
	  mutate(Module = paste0("M",partition[Protein])) %>%
	  # subset -- keep data from single protein
	  filter(Protein == protein) %>% 
	  # collect data
	  select(Protein, Mixture, Channel, Condition, BioFraction, Abundance, SE, Module) %>%
          # group by Genotype.BioFraction.Protein
	  group_by(Condition,Protein) %>%
	  # calculate the average of the three mixtures
	  summarize(mean_Abundance = mean(Abundance),
		    SD = sd(Abundance),
		    SE = unique(SE),.groups="drop") %>% 
	  mutate(CV = SD/mean_Abundance) %>%
	  # scale profile
	  mutate(norm_Abundance = mean_Abundance/max(mean_Abundance)) %>%
	  # munge
	  mutate(BioFraction = sapply(strsplit(as.character(Condition),"\\."),"[",2)) %>%
	  mutate(Genotype = sapply(strsplit(as.character(Condition),"\\."),"[",1)) %>%
	  mutate(BioFraction = factor(BioFraction, levels=c("F4","F5","F6","F7","F8","F9","F10")))
  # Generate the plot
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction)
  # switch - scale = TRUE/FALSE
  if (scale_abundance) {
	  plot <- plot + aes(y = norm_Abundance)
  } else {
	  plot <- plot + aes(y = mean_Abundance)
  }
  plot <- plot + aes(group = Genotype)
  plot <- plot + aes(colour = Genotype)
  plot <- plot + aes(shape = Genotype)
  plot <- plot + aes(fill = Genotype)
  plot <- plot + aes(shade = Genotype)
  # switch - scale = TRUE/FALSE
  if (scale_abundance) {
    plot <- plot + aes(ymin=norm_Abundance - CV)
    plot <- plot + aes(ymax=norm_Abundance + CV)
  } else {
    plot <- plot + aes(ymin=mean_Abundance - mean_Abundance * CV)
    plot <- plot + aes(ymax=mean_Abundance + mean_Abundance * CV)
  }
  plot <- plot + geom_line()
  plot <- plot + geom_ribbon(alpha=0.1, linetype="blank")
  plot <- plot + geom_point(size=2)
  plot <- plot + ggtitle(paste(gene,protein,sep=" | "))
  plot <- plot + ylab("Protein Abundance")
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
  plot <- plot + scale_colour_manual(values=c(wt_color,mut_color))
  plot <- plot + scale_fill_manual(values=c(wt_color,mut_color))
  plot <- plot + theme(legend.position = "none")
  return(plot)
} #EOF
