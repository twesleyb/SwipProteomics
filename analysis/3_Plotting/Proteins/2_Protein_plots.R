#!/usr/bin/env Rscript

# title: 
# description: plot protein abundance
# authors: Tyler W Bradshaw

## Inputs --------------------------------------------------------------------

# options:
wt_color = "#47b2a4"
mut_color <- col2hex(c("R"=148,"G"=33,"B"=146))
colors <- c("#000000","#303030","#5E5E5E", # WT Blacks
	    "#942192","#B847B4","#DC6AD7") # Swip Purples

# Input data in root/data/
# * msstats_prot


## Output ---------------------------------------------------------------------
# * a single pdf with plots of all proteins


## Functions -------------------------------------------------------------------

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


plot_protein <- function(prot_df, gene_map, protein,legend=FALSE) {
  # function to generate a proteins plot
  # Subset the data.
  gene <- gene_map$symbol[match(protein,gene_map$uniprot)]
  df <- subset(prot_df,prot_df$Protein == protein)
  # Insure Fraction is a factor, and levels are in correct order.
  df$BioFraction <- factor(df$BioFraction,
  			      levels=c("F4","F5","F6","F7","F8","F9","F10"))
  # Collect FDR stats.
  stats <- df %>% group_by(Genotype,BioFraction) %>% 
  		summarize(`Max Abundance` = max(Abundance),
  			  FDR = unique(adj.pvalue),.groups="drop")
  stats$ypos <- 1.02 * max(stats$`Max Abundance`)
  stats <- stats %>% filter(Genotype == "Control")
  stats$symbol <- ""
  stats$symbol[stats$FDR<0.1] <- "."
  stats$symbol[stats$FDR<0.05] <- "*"
  stats$symbol[stats$FDR<0.005] <- "**"
  stats$symbol[stats$FDR<0.0005] <- "***"
  # Generate the plot.
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction, y = Abundance, 
  		   group = interaction(Mixture,Genotype),
  		   colour = interaction(Mixture,Genotype))
  plot <- plot + geom_point(aes(shape=Genotype, fill=Genotype),size=2)
  plot <- plot + geom_line()
  plot <- plot + ggtitle(protein)
  # Annotate with significance stars.
  if (any(stats$FDR<0.1)) {
    plot <- plot + annotate("text",x=stats$BioFraction,y=max(stats$ypos),label=stats$symbol,size=7)
  }
  # Add Custom colors and modify legend title and labels.
  mylabs <- paste(c(rep('Control',3),rep('Mutant',3)),c(1,2,3))
  plot <- plot + scale_colour_manual(name="Replicate", values=colors, labels=mylabs) 
  #plot <- plot + scale_x_discrete(labels=stats$Cfg.Force)
  #plot <- plot + xlab("Force (xg)")
  plot <- plot + ggtitle(paste(gene,protein,sep=" | "))
  # Edit y axis.
  plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))
  # Make x and y axes consistent.
  plot <- plot + theme(axis.text.x = element_text(color="black",size=11, angle = 0, hjust = 1, family = "Arial"))
  plot <- plot + theme(axis.text.y = element_text(color="black",size=11, angle = 0, hjust = 1, family = "Arial"))
  # Remove background.
  plot <- plot + theme(panel.background = element_blank())
  # Add x and y axes.
  plot <- plot + theme(axis.line.x=element_line())
  plot <- plot + theme(axis.line.y=element_line())
  # Remove legend.
  if (!legend) {
	  plot <- plot + theme(legend.position = "none")
  }
  return(plot)
} #EOF


norm_to_max <- function(df) {
	norm <- function(x) { x*(1/max(x)) }
	df$Normalized.Abundance <- norm(df$Abundance)
	return(df)
}

## Set-up the workspace -------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

# Load additional functions in root/R/
suppressWarnings({ devtools::load_all() })

# Project directories:
fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs", "Proteins")

# If necessary, create figsdir.
if (! dir.exists(figsdir)) {
	dir.create(figsdir,recursive=TRUE)
}

# Set theme for the plots:
ggtheme(); set_font("Arial",font_path=fontdir)

# load the data
data(swip)
data(gene_map)
data(partition)
data(sig_modules)
data(msstats_prot)
data(module_colors)
data(msstats_results)

# combine protein data and statistical results
prot_df <- left_join(msstats_prot,msstats_results,by=c("Protein","BioFraction"))

# annotate with module membership
prot_df$Module <- paste0("M",partition[prot_df$Protein])


## Generate the plots ---------------------------------------------------------

# sort proteins by module membership, drop M0
sorted_prots <- unlist(split(names(partition),partition)[-1])

# Loop to generate plots for all_proteins.
message("\nGenerating plots for ", 
	formatC(length(sorted_prots),big.mark=","), " proteins.")
plots <- list()
pbar <- txtProgressBar(max=length(sorted_prots),style=3)
for (protein in sorted_prots) {
	# generate a proteins plot
	plot <- plot_protein(prot_df,gene_map,protein)
	# annotate with module assignment
	plot_label <- paste("Module:", partition[protein])
	yrange <- plot$data %>% dplyr::filter(Protein == protein) %>% 
		select(Abundance) %>% range()
	ypos <- yrange[1] - 0.1* diff(yrange)
	plot <- plot + annotate(geom="label",x=7, y=ypos, label=plot_label)
	setTxtProgressBar(pbar,value=match(protein,sorted_prots))
}
close(pbar)

plots[[swip]]

rand_prot <- names(sample(partition[partition==partition[swip]],1))
plots[[rand_prot]]

# Generate a legend ------------------------------------------------------------

# Generate a plot with a legend.
plot <- plot_protein(prot_df,gene_map, swip, legend = TRUE)
plot_legend <- cowplot::get_legend(plot) 

# save
myfile <- file.path(figsdir,"Protein_plots_legend.pdf")
ggsave(plot_legend,file=myfile,width=4.5,height=4.5)





########################################################################
## New script
#!/usr/bin/env Rscript

# title: 
# description: plot protein abundance
# authors: Tyler W Bradshaw

## Inputs --------------------------------------------------------------------

# options:
wt_color = "#47b2a4"
mut_color <- col2hex(c("R"=148,"G"=33,"B"=146))

# Input data in root/data/

## Output ---------------------------------------------------------------------


## Functions -------------------------------------------------------------------

## Prepare environment --------------------------------------------------------

root <- getrd()
renv::load(root,quiet=TRUE)

suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

suppressWarnings({ devtools::load_all() })

fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs", "Proteins")

if (! dir.exists(figsdir)) {
	dir.create(figsdir,recursive=TRUE)
}

ggtheme(); set_font("Arial",font_path=fontdir)

# load the data
data(swip)
data(gene_map)
data(partition)
data(sig_modules)
data(msstats_prot)
data(module_colors)
data(msstats_results)

# plot protein summary
# function to generate a proteins plot

plot_protein_summary(swip,wt_color,mut_color)

plot_protein_summary <- function(protein,wt_color,mut_color) {
	require(dplyr,quietly=TRUE)
	require(ggplot2,quietly=TRUE)
  msstats_df <- left_join(msstats_prot,msstats_results,by=c("Protein","BioFraction"))
  msstats_df$Module <- paste0("M",partition[msstats_df$Protein])
  gene <- gene_map$symbol[match(protein,gene_map$uniprot)]
  df <- subset(msstats_df,msstats_df$Protein == protein) %>% 
	select(Mixture,Channel,Condition,Protein,Abundance,SE,Module) %>% unique() %>% 
	group_by(Condition,Protein) %>% summarize(mean_Abundance = mean(Abundance),
						  SD = sd(Abundance),
						  SE = unique(SE),.groups="drop")
  df <- df %>% mutate(CV = SD/mean_Abundance)
  df <- df %>% mutate(norm_Abundance = mean_Abundance/max(mean_Abundance))
  df$BioFraction <- sapply(strsplit(as.character(df$Condition),"\\."),"[",2)
  df$Genotype <- sapply(strsplit(as.character(df$Condition),"\\."),"[",1)
  df$BioFraction <- factor(df$BioFraction,
			   levels=c("F4","F5","F6","F7","F8","F9","F10"))
  # Generate the plot.
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction)
  plot <- plot + aes(y = norm_Abundance)
  plot <- plot + aes(group = Genotype, colour = Genotype, shape=Genotype, fill=Genotype,shade=Genotype)
  plot <- plot + aes(ymin=norm_Abundance - CV)
  plot <- plot + aes(ymax=norm_Abundance + CV)
  plot <- plot + geom_line()
  plot <- plot + geom_ribbon(alpha=0.1, linetype="blank")
  plot <- plot + geom_point(size=2)
  plot <- plot + ggtitle(paste(gene,protein,sep=" | "))
  plot <- plot + ylab("Normalized Abundance")
  plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))
  plot <- plot + theme(axis.text.x = element_text(color="black",size=11, angle = 0, hjust = 1, family = "Arial"))
  plot <- plot + theme(axis.text.y = element_text(color="black",size=11, angle = 0, hjust = 1, family = "Arial"))
  plot <- plot + theme(panel.background = element_blank())
  plot <- plot + theme(axis.line.x=element_line())
  plot <- plot + theme(axis.line.y=element_line())
  plot <- plot + scale_colour_manual(values=c(wt_color,mut_color))
  plot <- plot + scale_fill_manual(values=c(wt_color,mut_color))
  plot <- plot + theme(legend.position = "none")
  return(plot)
} #EOF


quit()


## Combine plots for all proteins from a module together ----------------------

# Loop to do work.
message("\nAligning module prortein plots.")
grouped_plots <- list()
modules <- unique(partition[partition!=0])
sorted_modules <- modules[order(modules)] # Sort
sorted_modules <- paste0("M",sorted_modules)

pbar <- txtProgressBar(max=length(sorted_modules),style=3)
for (module in sorted_modules){
	# Get the protein data for a module.
	subdf <- prot_df %>% filter(Module == module)
	prots <- unique(as.character(subdf$Protein))
	prot_plots <- plots[prots]
	protein_list <- lapply(prot_plots,function(x) x$data)
	# Remove any null data arising from pesky prot.
	idx <- sapply(protein_list,is.null)
	if (any(idx)) { 
		message(paste("Warning: NULL data is removed."))
		protein_list <- protein_list[-which(idx)] 
	}
	# Get modules colors.
	module_color <- module_colors[module]
	# Normalize/scale.
	norm_prot <- lapply(protein_list,norm_to_max) # protein-wise normalization to max
	# Combine data for all proteins together.
	norm_df <- bind_rows(norm_prot)
	# To simplify plot, calculate protein-wise mean.
	norm_df <- norm_df %>% group_by(Protein,BioFraction,Genotype) %>%
		summarize(Mean.Normalized.Abundance=mean(Normalized.Abundance),
			  .groups="drop")
	# Insure Fraction and Cfg force are factors.
	# Sort factor levels in a logical order.
	norm_df$BioFraction <- factor(norm_df$BioFraction,
		      levels=c("F4","F5","F6","F7","F8","F9","F10"))
	# Fit with lm, add fitted values to df.
	#### stop:
	#fit <- lmerTest::lmer(Mean.Normalized.Abundance ~ (1|BioFraction) + Genotype, data = prot_df)
	#prot_df$Fitted.Intensity <- fit$fitted.values
	# Generate plot.
	plot <- ggplot(norm_df)
	plot <- plot + aes(x = BioFraction)
	plot <- plot + aes(y = Mean.Normalized.Abundance)
	plot <- plot + aes(group = interaction(Genotype,Protein))
	plot <- plot + aes(color = Genotype)
	plot <- plot + geom_line(alpha=0.31)
	#plot <- plot + geom_line(aes(y=Fitted.Intensity),size=1.5,alpha=0.5)
	plot <- plot + geom_point(aes(shape=Genotype, fill=Genotype),size=1.5)
	plot <- plot + scale_colour_manual(name="Replicate",
					   values = c(wt_color,module_color))
	plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))
	plot <- plot + theme(axis.text.x = element_text(color="black",size=11,
							angle = 0, hjust = 1, 
							family = "Arial"))
	plot <- plot + theme(axis.text.y = element_text(color="black",size=11,
							angle = 0, hjust = 1, 
							family = "Arial"))
	plot <- plot + theme(panel.background = element_blank())
	plot <- plot + theme(axis.line.x=element_line())
	plot <- plot + theme(axis.line.y=element_line())
	plot <- plot + theme(legend.position = "none")
	plot <- plot + ggtitle(paste("Module:",module))
	# Add module annotations.
	yrange <- range(plot$data$Mean.Normalized.Abundance)
	ymax <- yrange[1] + 0.10 * diff(yrange)
        # Add plot to list.
	grouped_plots[[module]] <- plot
	setTxtProgressBar(pbar,value=match(module,sorted_modules))
} # EOL
close(pbar)


## Save the data --------------------------------------------------------------

# Save all proteins.
message("\nSaving all plots, this will take several minutes.")
myfile <- file.path(figsdir,"Protein_plots.pdf")
ggsavePDF(plots, myfile)

# Save.
message("\nSaving all modules, this will take several minutes.")
myfile <- file.path(figsdir,"Module_Protein_plots.pdf")
ggsavePDF(grouped_plots, myfile)
