#!/usr/bin/env Rscript

# title: SwipProteomics
# description:
# author: twab

## ---- Inputs

# input data in root/data/
root = "~/projects/SwipProteomics"

input_part = "ne_surprise_surprise_partition"
input_colors = "ne_surprise_surprise_colors"


## ---- Prepare the R environment

renv::load(root, quiet=TRUE)
devtools::load_all(root,quiet=TRUE)

# load the data
data(gene_map)
data(msstats_prot)
data(list=input_part)
data(list=input_colors)

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
ggtheme()
set_font("Arial",font_path=fontdir)


## ---- load geneLists

# FIXME: rm notin and write_excel from geneLists!
library(geneLists)

data(corum) 


## ---- convert corum pathways to uniprot IDs

# create df mapping partition uniprot to entrez
df <- data.table(uniprot = names(partition))
idx <- match(df$uniprot,gene_map$uniprot)
df <- df %>% mutate(entrez = gene_map$entrez[idx])

# loop through corum complexes, map entrez to uniprot, drop NA
corum_prots <- lapply(corum, function(x) {
			      idx <- match(x, df$entrez)
			      uniprot <- df$uniprot[idx[!is.na(idx)]]
			      return(uniprot)
})


# insure proteins are all in the same module
prot_list <- list()
for (i in c(1:length(corum_prots))) {
  part <- partition[corum_prots[[i]]]
  if (all(table(part) == 1)) {
	  warning("complex is split among multiple modules.")
	  prot_list[[i]] <- NA
  } else {
    idy <- as.numeric(names(which(table(part)==max(table(part)))))
    prots <- names(part[part %in% idy])
    prot_list[[i]] <- prots
  }
}
names(prot_list) <- names(corum_prots)
idx <- unlist(sapply(prot_list,is.na))
prot_list <- prot_list[!idx]

# drop small complexes
sizes <- sapply(prot_list,length)
idx <- sizes < 3
filt_list <- prot_list[!idx]


results_list <- list()
pbar <- txtProgressBar(max=length(filt_list),style=3)

for (path in names(filt_list)){
  prots <- filt_list[[path]]
  lmer_args <- list()
  fx <- scale_Abundance ~ 0 + Condition + (1|Mixture) + (1|Protein)
  lmer_args[["data"]] <- msstats_prot %>% subset(Protein %in% prots) %>% 
	  group_by(Protein) %>% mutate(scale_Abundance = Abundance/max(Abundance))
  lmer_args[["formula"]] <- fx
  lmer_args[["control"]] <- lme4::lmerControl(check.conv.singular="ignore")
  fm <- do.call(lmerTest::lmer,lmer_args)
  LT <- getContrast(fm, "Mutant","Control")
  res <- lmerTestContrast(fm, LT) %>% mutate(Contrast="Mutant-Control") %>% 
	  unique()
  res$nProts <- length(prots)
  results_list[[path]] <- res
  setTxtProgressBar(pbar,value=match(path,names(filt_list)))
}
close(pbar)

# collect results
results_df <- bind_rows(results_list,.id="Pathway") %>% 
	mutate(Padjust = p.adjust(Pvalue, method="bonferroni")) %>% arrange(Pvalue)


## save to file
myfile <- file.path(root,"rdata","complex_results.csv")
fwrite(results_df, myfile)


## ---- Function


plot_complex <- function(path, prots, msstats_prots) {
  # color for Control condition
  wt_color = "#47b2a4"
  mut_color = "#b671af"
  # Subset
  subdat <- msstats_prot %>% subset(Protein %in% prots)
  # number of proteins in module
  nprots <- length(unique(subdat$Protein))
  # set factor order (levels)
  subdat$Genotype <- factor(subdat$Genotype,levels= c("Control","Mutant"))
  subdat$BioFraction <- factor(subdat$BioFraction,
			 levels=c("F4","F5","F6","F7","F8","F9","F10"))
  # scale to max, take median of three replicates
  df <- subdat %>% group_by(Protein) %>%
	  mutate(scale_Abundance = Abundance/max(Abundance)) %>%
	  group_by(Protein, Genotype, BioFraction) %>% 
	  summarize(med_Abundance = median(scale_Abundance), 
	          SD = sd(scale_Abundance),
	          N = length(scale_Abundance),
	          .groups="drop")
  # calculate coefficient of variation (CV == unitless error) and scale to max
  df <- df %>% mutate(CV = SD/med_Abundance)
  # get module fitted data by fitting linear model to scaled Abundance
  fm <- lmerTest::lmer(med_Abundance ~ 0 + Genotype:BioFraction + (1|Protein), df)
  # collect coefficients
  fit_df <- data.table("coef" = names(lme4::fixef(fm)),
		       "fit_y" = lme4::fixef(fm)) %>%
    mutate(Genotype = gsub("Genotype","",sapply(strsplit(coef,"\\:"),"[",1))) %>%
    mutate(BioFraction=gsub("BioFraction","",sapply(strsplit(coef,"\\:"),"[",2)))
  # combine module data and fitted values
  df <- left_join(df, fit_df,by=c("Genotype","BioFraction"))
  # again, insure factor order is correct
  df$Genotype <- factor(df$Genotype,levels=c("Control","Mutant"))
  df$BioFraction <- factor(df$BioFraction,
			   levels=c("F4","F5","F6","F7","F8","F9","F10"))
  # get marginal r2 for annot plot title
  r2 <- r.squaredGLMM.merMod(fm)[,"R2m"]
  r2_anno <- paste("(",paste(paste(c("R2.Fixef = "),
			     round(r2,3)),collapse=" | "),")")
  # Generate the plot
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction)
  plot <- plot + aes(y = med_Abundance)
  plot <- plot + aes(group = interaction(Genotype,Protein))
  plot <- plot + aes(colour = Genotype)
  plot <- plot + aes(shape = Genotype)
  plot <- plot + aes(fill = Genotype)
  plot <- plot + aes(shade = Genotype)
  plot <- plot + aes(ymin=med_Abundance - CV)
  plot <- plot + aes(ymax=med_Abundance + CV)
  plot <- plot + geom_line(alpha=0.25)
  plot <- plot + theme(legend.position = "none")
  plot <- plot + ggtitle(paste0(path, " (n = ",nprots,")\n", r2_anno))
  plot <- plot + ylab("Scaled Abundance")
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
  # add fitted lines
  plot <- plot + geom_line(aes(y=fit_y, group=interaction("fit",Genotype)),
			   linetype="dashed",alpha=1,size=0.75)
  # set colors
  plot <- plot + scale_colour_manual(values=c(wt_color,mut_color))
  return(plot)
} #EOF


## ---- generate plots 

plot_list <- list()
for (path in names(filt_list)){
	prots <- filt_list[[path]]
	plot <- plot_complex(path, prots, msstats_prot)
	plot_list[[path]] <- plot
} #EOL


## ---- save as pdf

myfile <- file.path(root,"figs","Modules","complex_profiles.pdf")
ggsavePDF(plot_list, myfile)

warnings()
