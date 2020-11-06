#!/usr/bin/env Rscript

root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

# load the raw psm data
myfile <- file.path(root,"rdata","msstats_psm.rda")
load(myfile)

library(grid)
library(dplyr)
library(gtable)
library(ggplot2)
library(data.table)

## functions -------------------------------------------------------------------

plot_corQC <- function(msstats_prot, mix, nbins=5, nSD=4) {
	# subset the data for QC PSM from a single mixture
	subdf <- msstats_psm %>%
		filter(Mixture == mix, Condition == "Norm") %>%
		group_by(Mixture,PSM) %>% 
		mutate(nObs = sum(!is.na(Intensity))) %>% 
		filter(nObs==2) %>% 
		mutate(log.Intensity = log2(Intensity)) 
	# calculate log ratio and cast into a wide data.frame
	df <- subdf %>% 
		reshape2::dcast(PSM ~ Channel + Condition, 
				value.var="log.Intensity") %>%
	        as.data.table() %>% as.matrix(rownames="PSM") %>% as.data.frame()
	colnames(df) <- c("QC1","QC2")
	# calculate ratio and mean of QC replicates
	ratioQC <- apply(df,1,diff)
	meanQC <- rowMeans(df)
	df$ratioQC <- ratioQC
	df$meanQC <- meanQC
	# bin by mean log.Intensity
	breaks <- stats::quantile(df$meanQC, 
				  seq(0,1, length.out=nbins+1), 
				  names=FALSE)
	df <- df %>% 
		mutate(bin = cut(meanQC, breaks, 
				 lables=FALSE, include.lowest=TRUE))
	df$bin <- as.factor(as.numeric(df$bin))
	# calculate bin stats
	df <- df %>% group_by(bin) %>% 
		summarize(meanRatio=mean(ratioQC),
			  SD=sd(ratioQC),.groups="drop") %>% 
		left_join(df,by="bin")
	# annotate outliers
	isHigh <- df$ratioQC > df$meanRatio + nSD*df$SD 
	isLow <- df$ratioQC < df$meanRatio - nSD*df$SD
	df$isOutlier <- isHigh | isLow
	# linear fit
	fit <- lm(QC1 ~ QC2, df)
	rho <- cor.test(~ QC1 + QC2, df, method = "pearson", conf.level = 0.95)
	m <- paste("Slope =", round(coef(fit)[2], 4))
	r2 <- paste("R2 =", round(rho$estimate, 4))
	## trendlines for outlier thresholds
	#l1 <- geom_abline(intercept = coef(fit)[1] + 1, 
	#	       slope = coef(fit)[2],
	#	       color = "black", linetype = "dashed")
	#l2 <- geom_abline(intercept = coef(fit)[1] - 1, 
	#	       slope = coef(fit)[2],
	#	       color = "black", linetype = "dashed")
	# colors
	palette <- "Blues"
	colors <- head(rev(RColorBrewer::brewer.pal(nbins+1,palette)),nbins)
	colors <- setNames(colors,levels(df$bin))
	df$color <- as.character(colors[df$bin])
	df$color[df$isOutlier] <- "red"
	# generate the plot
	plot <- ggplot(df)
	plot <- plot + aes(x = QC1, y = QC2, group=isOutlier)
	plot <- plot + geom_point(colour=df$color)
	plot <- plot + geom_abline(intercept = coef(fit)[1], 
				   slope = coef(fit)[2],
				   color = "black", linetype = "dashed")
	plot <- plot + scale_fill_manual(values=df$color)
	plot <- plot + theme(panel.background=element_blank())
	plot <- plot + theme(axis.line.x=element_line())
	plot <- plot + theme(axis.line.y=element_line())
	plot <- plot + theme(axis.text.x = element_text(size=11,family="Arial"))
	plot <- plot + theme(axis.text.y = element_text(size=11,family="Arial"))
	# create annotation layer.
	mytable <- rbind(r2, m)
	build <- ggplot_build(plot)
	yrange <- unlist(build$layout$panel_params[[1]][8])
	xrange <- range(plot$data$QC1)
	xmin <- min(xrange)
	xmax <- max(xrange)
	xdelta <- xmax - xmin
	ymin <- min(yrange)
	ymax <- max(yrange)
	ydelta <- ymax - ymin
	tt <- ttheme_default(base_size = 11, 
			     core = list(bg_params=list(fill="white")))
	tab <- tableGrob(mytable, rows = NULL, theme = tt)
	g <- gtable_add_grob(tab, grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
			     t = 1, b = nrow(tab), l = 1, r = ncol(tab))
	## add annotation
	plot <- plot + annotation_custom(g, xmin = xmin - 0.65 * xdelta, xmax,
				         ymin = ymin + 0.8 * ydelta, ymax)
	plot <- plot + ggtitle(mix)
	return(plot)
}

## main ------------------------------------------------------------------------

p1 <- plot_corQC(msstats_prot,"M1")
p2 <- plot_corQC(msstats_prot,"M2")
p3 <- plot_corQC(msstats_prot,"M3")

myfile <- file.path(root,"figs","Outliers")
if (!dir.exists(basename(myfile))) { dir.create(basename(myfile)) }

pdf(myfile,onefile=TRUE)
p1
p2
p3
invisible(dev.off())
