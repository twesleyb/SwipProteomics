plot_protein <- function(data,protein,id.col="Accession") {

	# Generate protein plot.
	# Subset the data.
	suppressPackageStartupMessages({
		library(dplyr)
		library(ggplot2)
		library(data.table)
	})

	# Colors for the plot.
	# Generated shade of black online: https://mycolor.space/
	graphpad_purple <- c("R"=148,"G"=33,"B"=146)
	colors <- c("#000000","#303030","#5E5E5E", # WT Blacks
	   	    "#942192","#B847B4","#DC6AD7") # Swip Purples
	
	# Subset the data.
	df <- subset(data,data[[id.col]] == protein)
	gene <- unique(df$Gene)

	# Insure Fraction is a factor, and levels are in correct order.
	df$Fraction <- factor(df$Fraction,
			      levels=c("F4","F5","F6","F7","F8","F9","F10"))
	df$"Cfg Force (xg)" <- factor(df$"Cfg Force (xg)")
	levels(df$"Cfg Force (xg)") <- c("5,000","9,000","12,000","15,000",
					 "30,000", "79,000","120,000")

	# Collect FDR stats.
	stats <- df %>% group_by(Treatment,Fraction) %>% 
		summarize(Intensity = max(Intensity), FDR = unique(FDR))
	stats$ypos <- 1.02 * log2(max(stats$Intensity))
	stats <- stats %>% filter(Treatment == "Control")
	stats$Cfg.Force <- levels(df$"Cfg Force (xg)")
	stats$symbol <- ""
	stats$symbol[stats$FDR<0.1] <- "."
	stats$symbol[stats$FDR<0.05] <- "*"
	stats$symbol[stats$FDR<0.005] <- "**"
	stats$symbol[stats$FDR<0.0005] <- "***"
	#stats$ypos <- 1.02*log2(stats$Intensity)

	# Generate the plot.
	plot <- ggplot(df, aes(x = Fraction, y = log2(Intensity),
			       group = interaction(Experiment,Treatment),
			       colour = interaction(Experiment,Treatment))) + 
                geom_point(aes(shape=Treatment,
			       fill=Treatment),size=2) + 
		geom_line() + 
		ggtitle(protein)
	# Annotate with significance stars.
	plot <- plot + annotate("text",x=stats$Fraction,
				y=max(stats$ypos),label=stats$symbol,size=7)
	# Add Custom colors and modify legend title and labels.
	mylabs <- paste(c(rep('Control',3),rep('Mutant',3)),c(1,2,3))
	plot <- plot + scale_colour_manual(name="Replicate",
				           values=colors,
					   labels=mylabs) 
	plot <- plot + scale_x_discrete(labels=stats$Cfg.Force)
	plot <- plot + xlab("Force (xg)")
	plot <- plot + ggtitle(paste(gene,protein,sep=" | "))

	# Edit y axis.
	plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))

	# Make x and y axes consistent.
	plot <- plot + theme(axis.text.x = element_text(color="black",size=11,
							angle = 0, hjust = 1, 
							family = "Arial"))
	plot <- plot + theme(axis.text.y = element_text(color="black",size=11,
							angle = 0, hjust = 1, 
							family = "Arial"))

	# Remove background.
	plot <- plot + theme(panel.background = element_blank())

	# Add x and y axes.
	plot <- plot + theme(axis.line.x=element_line())
	plot <- plot + theme(axis.line.y=element_line())

	# Remove legend.
	plot <- plot + theme(legend.position = "none")

	return(plot)
}
