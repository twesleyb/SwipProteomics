#!/usr/bin/env Rscript

# title: SwipProteomics
# description: examine some protein complexes
# author: twab

## ---- Inputs

# input data in root/data/
root = "~/projects/SwipProteomics"
input_part = "ne_surprise_partition"
input_colors = "ne_surprise_colors"


## ---- Prepare the R environment

renv::load(root, quiet=TRUE)
devtools::load_all(root, quiet=TRUE)

# load the data
data(gene_map)
data(washc_prots)
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
set_font("Arial", font_path=fontdir)


## ---- load geneLists

library(geneLists) # twesleyb/geneLists

# load corum protein complexes (mapped to mouse entrez)
data(corum) 


## ---- convert corum pathways to uniprot IDs

# create df mapping partition uniprot to entrez
df <- data.table(uniprot = names(partition))
idx <- match(df$uniprot,gene_map$uniprot)
df <- df %>% mutate(entrez = gene_map$entrez[idx])

# loop through corum complexes, map entrez to uniprot, drop NA
# NOTE: this takes a couple seconds... 
corum_prots <- lapply(corum, function(x) {
			      idx <- match(x, df$entrez)
			      uniprot <- df$uniprot[idx[!is.na(idx)]]
			      stopifnot(all(uniprot %in% msstats_prot$Protein))
			      return(uniprot)
})


# insure proteins are all in the same module
prot_list <- list()
multi_list <- list()
for (i in c(1:length(corum_prots))) {
  part <- partition[corum_prots[[i]]]
  if (all(table(part) == 1)) {
	  #warning("complex is split among multiple modules.")
	  prot_list[[i]] <- NA
  } else {
    idy <- as.numeric(names(which(table(part)==max(table(part)))))
    prots <- names(part[part %in% idy])
    if (length(idy) > 1){
    multi_list[[names(corum_prots)[i]]] <- part
      } else {
        prot_list[[i]] <- setNames(rep(idy,length(prots)),nm=prots)
    }
  }
}
names(prot_list) <- names(corum_prots)
idx <- sapply(prot_list,function(x) any(is.na(x)))
prot_list <- prot_list[!idx]


# Warnings are bc complexes may be split between modules. We collect the
# proteins that are all assigned to the same module

# drop small complexes
sizes <- sapply(prot_list,length)
idx <- sizes < 3
filt_list <- prot_list[!idx]


## ---- add RAVE complex just for kicks

rave <- mapID(c("Rogdi","Dmxl2","Wdr7"))
filt_list[["RAVE"]] <- partition[rave]

filt_list[["WASHC"]] <- partition[washc_prots]

cbind(k = length(filt_list), med_size = median(sapply(filt_list,length))) %>%
	knitr::kable()


## ---- loop to do work

# loop to do module-level lmer statistical analysis
results_list <- list()
pbar <- txtProgressBar(max=length(filt_list),style=3)
for (path in names(filt_list)){
  # build input args list for lmerTest
  lmer_args <- list()
  prots <- names(filt_list[[path]])
  # formula for statistical comparisons
  fx <- log2(scale_Intensity) ~ 1 + Condition + (1|Protein)
  lmer_args[["formula"]] <- fx
  lmer_args[["data"]] <- msstats_prot %>% subset(Protein %in% prots) %>% 
	  group_by(Protein) %>% mutate(Intensity = 2^Abundance) %>%
	  mutate(scale_Intensity=Intensity/sum(Intensity))
  lmer_args[["control"]] <- lme4::lmerControl(check.conv.singular="ignore")
  fm <- do.call(lmerTest::lmer,lmer_args)
  LT <- getContrast(fm, "Mutant","Control")
  res <- lmerTestContrast(fm, LT) %>% 
	  mutate(Contrast="Mutant-Control") %>% unique()
  r2 <- r.squaredGLMM.merMod(fm)
  res$R2.fixef <- r2[,"R2m"]
  res$R2.total <- r2[,"R2c"]
  res$nProts <- length(prots)
  res$Module <- paste0("M",unique(filt_list[[path]]))
  results_list[[path]] <- res
  setTxtProgressBar(pbar,value=match(path,names(filt_list)))
}
close(pbar)


# check the WASH Complex result
results_list[["WASHC"]] %>% knitr::kable()


# collect results
results_df <- bind_rows(results_list,.id="Pathway") %>% 
	mutate(Padjust = p.adjust(Pvalue, method="bonferroni")) %>% arrange(Pvalue)

# summary sig
results_df %>% 
	summarize(nSig= sum(Padjust < 0.05), Total = length(unique(Pathway))) %>% 
	knitr::kable()


## ---- Function


plotComplex <- function(path, prots, msstats_prot) {
  #############################################################################
  #path = sample(names(filt_list),1)
  #prots = filt_list[[path]]
  #plotComplex(path,prots,msstats_prot)
  #############################################################################
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
  # prepare the data for plotting
  df <- subdat %>% 
	  mutate(Intensity = 2^Abundance) %>% 
	  group_by(Protein, Genotype, BioFraction) %>%
	  summarize(med_Intensity = median(Intensity), .groups="drop") %>%
	  group_by(Protein) %>% 
	  mutate(scale_Intensity = scale01(log2(med_Intensity/sum(med_Intensity))))
  # get module fitted data by fitting linear model to scaled Abundance
  fx <- scale_Intensity ~ 0 + Genotype:BioFraction + (1|Protein)
  lmer_control <- lme4::lmerControl(check.conv.singular="ignore")
  fm <- lmerTest::lmer(fx, df, control = lmer_control)
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
  # Generate the plot
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction)
  plot <- plot + aes(y = scale_Intensity)
  plot <- plot + aes(group = interaction(Genotype,Protein))
  plot <- plot + aes(colour = Genotype)
  plot <- plot + aes(shape = Genotype)
  plot <- plot + aes(fill = Genotype)
  plot <- plot + aes(shade = Genotype)
  plot <- plot + geom_line(alpha=0.25)
  plot <- plot + theme(legend.position = "none")
  plot <- plot + ggtitle(paste0(path, " (n = ",nprots,")"))
  plot <- plot + ylab("Scaled Protein Intensity")
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
  plot <- plot + geom_line(aes(y=fit_y, group=interaction("fit",Genotype)),
			   linetype="dashed",alpha=1,size=0.75)
  plot <- plot + scale_colour_manual(values=c(wt_color,mut_color))
  return(plot)
} #EOF


## ---- generate plots 

plot_list <- list()
pbar <- txtProgressBar(max=length(filt_list),style=3)
for (path in names(filt_list)){
	prots <- names(filt_list[[path]])
	plot <- plotComplex(path, prots, msstats_prot)
	plot_list[[path]] <- plot
	setTxtProgressBar(pbar, value=match(path,names(filt_list)))
} #EOL
close(pbar)


## ---- save results to file

# save stats
myfile <- file.path(root,"tables","S7_SWIP-TMT_Complex_Results.xlsx")
write_excel(list("complexes"=results_df), myfile)


# save plots
myfile <- file.path(root,"figs","Modules","complex_profiles.pdf")
ggsavePDF(plot_list, myfile)


## ---- examine multi


plotComplex2 <- function(path,prots, msstats_prot){
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
    # prepare the data for plotting
    df <- subdat %>% 
  	  mutate(Intensity = 2^Abundance) %>% 
  	  group_by(Protein, Genotype, BioFraction) %>%
  	  summarize(med_Intensity = median(Intensity), .groups="drop") %>%
  	  group_by(Protein) %>% 
  	  mutate(scale_Intensity = scale01(log2(med_Intensity/sum(med_Intensity))))
    # get module fitted data by fitting linear model to scaled Abundance
    fx <- scale_Intensity ~ 0 + Genotype:BioFraction + (1|Protein)
    lmer_control <- lme4::lmerControl(check.conv.singular="ignore")
    fm <- lmerTest::lmer(fx, df, control = lmer_control)
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
    # Generate the plot
    plot <- ggplot(df)
    plot <- plot + aes(x = BioFraction)
    plot <- plot + aes(y = scale_Intensity)
    plot <- plot + aes(group = interaction(Genotype,Protein))
    plot <- plot + aes(colour = Genotype)
    plot <- plot + aes(shape = Genotype)
    plot <- plot + aes(fill = Genotype)
    plot <- plot + aes(shade = Genotype)
    plot <- plot + geom_line(alpha=0.25)
    plot <- plot + theme(legend.position = "none")
    plot <- plot + ggtitle(paste0(path, " (n = ",nprots,")"))
    plot <- plot + ylab("Scaled Protein Intensity")
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
    plot <- plot + geom_line(aes(y=fit_y, group=interaction("fit",Genotype)),
  			   linetype="dashed",alpha=1,size=0.75)
    plot <- plot + scale_colour_manual(values=c(wt_color,mut_color))
    return(plot)
} #EOF

mylist <- list()
for (path in names(multi_list)){
	prots = names(multi_list[[path]])
	mylist[[path]] <- plotComplex2(path,prots,msstats_prot)
}

# save
myfile <- file.path(root,"figs","Modules","split-complexes.pdf")
ggsavePDF(mylist,myfile)
