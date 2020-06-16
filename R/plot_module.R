# Function to generate module plots.
#' plot_module
#' @export 
plot_module <- function(module_stats,module) {

	suppressPackageStartupMessages({
		library(dplyr)
		library(data.table)
		library(ggplot2)
	})

	# Define colors for WT and Swip mutants.
	colors <- c("#000000","#303030","#5E5E5E", # WT Blacks
		    "#942192","#B847B4","#DC6AD7") # Swip Purples

	# Subset the data.
	df <- module_stats %>% filter(Module == module)

	# Insure levels of Fraction are in correct order.
	df$Fraction <- factor(df$Fraction, 
			      levels=c("F4","F5","F6","F7","F8","F9","F10"))

	# Extract some summary info for annotating the plot.
	pval <- formatC(unique(df$KW.padj),format="e",digits=2)
	pve <- formatC(unique(df$PVE),digits=2)
	nProts <- unique(df$nProteins)

	# Collect statistics for annotating the plot with significance stars.
	stats <- df %>% group_by(Fraction) %>% 
		summarize(ypos = 1.2*max(ME),DT.pval = unique(DT.pval))
	stats$symbol <- ""
	stats$symbol[stats$DT.pval<0.1] <- "."
	stats$symbol[stats$DT.pval<0.05] <- "*"
	stats$symbol[stats$DT.pval<0.005] <- "**"
	stats$symbol[stats$DT.pval<0.0005] <- "***"

	# Generate the plot.
	plot <- ggplot(df, aes(x=Fraction,y=ME,
		       group=interaction(Experiment,Treatment),
		       colour=interaction(Experiment,Treatment)))

	# Add geoms.
	plot <- plot + geom_point() + geom_path()

	# Customize labels.
	mytitle <- paste0("Module",module,"\n(P.Adj = ",pval,"; PVE = ",pve,")")
	plot <- plot + ggtitle(mytitle) + ylab("Module Summary Abundance")

	# Add significance stars annotation.
	plot <- plot + 
		annotate("text",x=stats$Fraction,
			 y=max(stats$ypos),label=stats$symbol,size=11)

	# Customize colors and modify legend title and labels.
	mylabs <- paste(c(rep('Control',3),rep('Mutant',3)),c(1,2,3))
	plot <- plot + scale_colour_manual(name="Replicate",
					   values=colors,
					   labels=mylabs) 

	return(plot)
}
