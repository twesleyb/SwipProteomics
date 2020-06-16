#!/usr/bin/env Rscript

ggplotPCAprot <- function(data,...){
	# Generate PCA plot of proteins.
	# param Data is a normalized data matrix.
	# Perform PCA.
	pca <- prcomp(data,...)
	pca_summary <- as.data.frame(t(summary(pca)$importance))
	idx <- order(pca_summary[["Proportion of Variance"]],decreasing=TRUE)
	pca_summary <- pca_summary[idx,]
	top2_pc <- head(pca_summary[["Proportion of Variance"]],2)
	names(top2_pc) <- head(rownames(pca_summary),2)
	# Plot axis labels:
	x_label <- paste0(names(top2_pc)[1],
			  " (PVE: ",round(100*top2_pc[1],2)," %)")
	y_label <- paste0(names(top2_pc)[2],
			  " (PVE: ",round(100*top2_pc[2],2)," %)")
	# Collect data for plotting.
	df <- as.data.frame(pca$x[,names(top2_pc)])
	colnames(df) <- c("x","y")
	# Generate the plot.
	plot <- ggplot(df, aes(x,y)) + geom_point()
	plot <- plot + xlab(x_label)
	plot <- plot + ylab(y_label)
	plot <- plot + theme(axis.title.x = element_text(color = "black")) 
	plot <- plot + theme(axis.title.x = element_text(size = 11))
       	plot <- plot + theme(axis.title.x = element_text(face = "bold"))
	plot <- plot + theme(axis.title.y = element_text(color = "black")) 
	plot <- plot + theme(axis.title.y = element_text(size = 11))
       	plot <- plot + theme(axis.title.y = element_text(face = "bold"))
	plot <- plot + theme(panel.background = element_blank())
	plot <- plot + theme(panel.border = element_rect(fill=NA))
	plot <- plot + scale_x_continuous(expand = c(0,0))
	plot <- plot + scale_y_continuous(expand = c(0,0))
	# Return the plot
	return(plot)
}
