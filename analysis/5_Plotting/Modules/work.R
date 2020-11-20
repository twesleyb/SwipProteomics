#!/usr/bin/env Rscript

# title: 
# description: plot protein abundance
# authors: Tyler W Bradshaw

## Inputs --------------------------------------------------------------------

# options:

# Input data in root/data/
root = "~/projects/SwipProteomics"


## Prepare environment --------------------------------------------------------

renv::load(root,quiet=TRUE)
suppressWarnings({ devtools::load_all() })

# load the data
data(swip)
data(gene_map)
data(ne_surprise_surprise_partition)
data(module_gof)
#data(msstats_prot)
data(swip_tmt)
data(module_colors)
data(msstats_results)


# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
	library(doParallel)
})

# project dirs
fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs", "Modules")
if (! dir.exists(figsdir)) {
	dir.create(figsdir,recursive=TRUE)
}

# set plotting theme and font
ggtheme(); set_font("Arial",font_path=fontdir)


## Function ------------------------------------------------------------------

#  # subset the data
#df <- reshape2::melt(norm_prot[prots,])
#colnames(df) <- c("Protein","Sample","Abundance")
#df <- df %>% 
#	mutate(Mixture = sapply(strsplit(as.character(Sample),"\\_"),"[",1),
#	       Genotype = sapply(strsplit(as.character(Sample),"\\_"),"[",2),
#	       BioFraction = sapply(strsplit(as.character(Sample),"\\_"),"[",3))
#df$Sample <- NULL

modules= split(names(partition),partition)
names(modules) <- paste0("M",names(modules))

partition[swip]

module_prots <- modules[["M33"]]

plot_profile <- function(module, msstats_prot, partition,
			 module_colors, module_gof, wt_color = "#47b2a4") {

	if (module %notin% module_gof$Module) {
	       warning(module," is not in 'module_gof'.")
	       return(NULL)
	}

# protein-level fits
# NOTE: we scale each protein to its maximum
fm_list <- list()
args_list <- list()
#args_list[["formula"]] <- scale_Abundance ~ 0 + Genotype:BioFraction + (1|Mixture)
args_list[["formula"]] <- scale_Abundance ~ 0 + Genotype:BioFraction 

for (prot in module_prots) {
	#args_list[["data"]] <- msstats_prot %>% subset(Protein == prot) %>%
	#	mutate(scale_Abundance=Abundance/max(Abundance))
	args_list[["data"]] <- swip_tmt %>% subset(Protein == prot) %>%
		mutate(scale_Abundance=log2(Intensity)/max(log2(Intensity)))
	fm_list[[prot]] <- do.call(lm,args_list)
	#fm_list[[prot]] <- lmerFit(args_list)
}

# protein-level plots
# plot fit proteins
plot_list <- list()
for (prot in names(fm_list)){
  df <- data.table("coeff" = names(coef(fm_list[[prot]])),
  		   "beta" = coef(fm_list[[prot]]))
  df <- df %>% mutate(coeff = as.character(coeff)) %>%
        mutate(Genotype=gsub("Genotype","",
			     sapply(strsplit(coeff,":"),"[",1)),
               BioFraction=gsub("BioFraction","",
				sapply(strsplit(coeff,":"),"[",2)))
  levels(df$BioFraction) <- c("F4","F5","F6","F7","F8","F9","F10")
  levels(df$Genotype) <- c("Control","Mutant")
  plot <- ggplot(df)
  plot <- plot + aes(x=BioFraction,y=beta,group=Genotype)
  plot <- plot + aes(color=Genotype)
  plot <- plot + geom_path()
  plot <- plot + ggtitle(prot)
  plot_list[[prot]] <- plot
}

# module-level fit
# fit model to scaled protein abundance
fx1 <- scale_Abundance ~ 0 + Genotype:BioFraction + (1|Protein)
args_list <- list()
args_list[["formula"]] <- fx1
args_list[["data"]] <- swip_tmt %>% subset(Protein %in% module_prots) %>%
	group_by(Protein) %>% mutate(scale_Abundance=log2(Intensity)/max(log2(Intensity)))
fm1 <- lmerFit(args_list)
fm_list[["fit"]] <- fm1

# collect coefficients
df <- reshape2::melt(do.call(rbind,lapply(fm_list,fixef)))
colnames(df) <- c("Protein","coeff","beta")
df <- df %>% mutate(coeff = as.character(coeff)) %>%
      mutate(Genotype=gsub("Genotype","",sapply(strsplit(coeff,":"),"[",1)),
             BioFraction=gsub("BioFraction","",sapply(strsplit(coeff,":"),"[",2)))
levels(df$BioFraction) <- c("F4","F5","F6","F7","F8","F9","F10")
levels(df$Genotype) <- c("Control","Mutant")

df$group <- as.character(df$Protein) == "fit"
plot <- ggplot(df)
plot <- plot + aes(x=BioFraction,y=beta,group=interaction(Protein,Genotype))
plot <- plot + aes(color=Genotype)
plot <- plot + geom_path()
plot
