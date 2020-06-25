#FIXME: Incomplete

# Generate missing value density plots.
tp <- sl_peptide

df <- tp %>% filter(Treatment != "QC") %>% filter(Experiment == "Exp1") %>%
	group_by(Accession,Sequence) %>%
	dplyr::summarize(Mean_Intensity = mean(Intensity,na.rm=TRUE), 
			      Contains_Missing = any(is.na(Intensity)),
			      .groups = "drop") %>%
	filter(!is.na(Mean_Intensity)) # Remove peptides all(missing values).

plot <- ggplot(df)
plot <- plot + aes(x=log2(Mean_Intensity),
		   fill=Contains_Missing,
		   colour=Contains_Missing)
plot <- plot + geom_density(alpha=0.1,size=1) + ggtitle("Experiment 1")

# NOTE: From this plot, the distributions overlap. Missing values are
# missing at random or even missing completely at random.

if (save_plots) {
	myfile <- file.path(figsdir,"Peptide_MV_Density.pdf")
	ggsave(myfile,plot,height=fig_height,width=fig_width)
}

