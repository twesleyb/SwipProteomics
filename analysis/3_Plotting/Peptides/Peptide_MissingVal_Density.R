
#!/usr/bin/env Rscript

#' ---
#' title: Swip Proteomics Plotting
#' description: 
#' authors: Tyler W Bradshaw
#' ---

## OPTIONS:

## OUTPUT:

#---------------------------------------------------------------------
## Misc functions
#---------------------------------------------------------------------

# Get the repository's root directory.
getrd <- function(here=getwd(), dpat= ".git") {
	in_root <- function(h=here, dir=dpat) { 
		check <- any(grepl(dir,list.dirs(h,recursive=FALSE))) 
		return(check)
	}
	# Loop to find root.
	while (!in_root(here)) { 
		here <- dirname(here) 
	}
	root <- here
	return(root)
}

# Load renv.
root <- getrd()
renv::load(root, quiet=TRUE)

# Generate missing value density plots.
# FIXME: tidy_peptide is not saved, so we dont have the necessary data.

quit()
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

