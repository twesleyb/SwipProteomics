#!/usr/bin/env Rscript

#' ---
#' title: Swip Proteomics Plotting
#' description: Plot an overview of the network.
#' authors: Tyler W Bradshaw
#' ---

## ---- input 
root <- "~/projects/SwipProteomics"


## ---- Prepare the workspace

# Load renv
renv::load(root, quiet=TRUE)

# Imports
suppressPackageStartupMessages({
	library(gridExtra)
	library(data.table)
	library(ggplot2)
	library(grid)
})

devtools::load_all()

# Set plotting theme.
set_font("Arial",font_path=file.path(root,"fonts"))

# Load the data
data(msstats_prot)
data(ne_surprise_partition)

# Summarize key module statistics
n_nodes <- formatC(length(unique(msstats_prot$Protein)),big.mark=",")
n_modules <- length(unique(partition))-1 # M0 is not a module
p_clustered <- round(100*sum(partition!=0)/length(partition),2)
median_size <- median(sapply(split(partition,partition)[-1],length))

df <- data.table("Nodes" = n_nodes,
		 "Modules" = n_modules,
		 "Percent Clustered" = p_clustered,
		 "Median Module Size" = median_size)

# Create table.
gtab <- tableGrob(df, rows=NULL, theme=ttheme_default())

# Save.
myfile <- file.path(root,"figs","Networks","Network_Summary.pdf")
ggsaveTable(gtab,myfile)




data(swip)
data(msstats_results)
df = msstats_results %>% filter(Protein == swip) %>% filter(Contrast == 'Mutant-Control')

df %>% tibble::add_column(percentControl = 2^df$log2FC, .after="log2FC")

df = df %>% dplyr::mutate(Pvalue = formatC(Pvalue,digits=3))
df = df %>% dplyr::mutate(SE = formatC(SE,digits=3))
df = df %>% dplyr::mutate(DF = formatC(DF,digits=3))
df = df %>% dplyr::mutate(FDR = formatC(FDR,digits=3))
df = df %>% dplyr::mutate(log2FC = formatC(log2FC,digits=3))
df = df %>% dplyr::select(-Padjust, -Entrez)
gtab <- tableGrob(df, rows=NULL, theme=ttheme_default())
plot(gtab)


m = "M310"
df = module_results %>% filter(Module == m)
colnames(df)[colnames(df)=="nProts"] <- "n"
df = df %>% dplyr::mutate(Pvalue = formatC(Pvalue,digits=3))
df = df %>% dplyr::mutate(SE = formatC(SE,digits=3))
df = df %>% dplyr::mutate(percentControl = formatC(percentControl,digits=3))
df = df %>% dplyr::mutate(Padjust = formatC(Padjust,digits=3))
df = df %>% dplyr::mutate(DF = formatC(DF,big.mark=","))
df = df %>% dplyr::mutate(log2FC = formatC(log2FC,digits=3))
df = df %>% dplyr::mutate(S2 = formatC(S2,digits=3))
df = df %>% dplyr::mutate(Tstatistic = formatC(Tstatistic,digits=3))
df = df %>% select(Contrast, log2FC, Tstatistic, DF, Padjust)
gtab <- tableGrob(df, rows=NULL, theme=ttheme_default())
ggsaveTable(gtab,"tab.pdf")

plot(gtab)


data(ne_surprise2_partition)
data(gene_map)

prots=lopit_prots[["PROTEASOME"]]
part = partition[prots]
p = part[part==310]
mapID(names(p[!is.na(p)]),"uniprot","symbol")


modules = split(names(partition),partition)
names(modules) <- paste0("M", names(modules))
prots = modules[[m]]
length(prots)

lmer_args[["formula"]] <- Abundance ~ 1 + Condition + (1|Protein)
lmer_args[["data"]] <- msstats_prot %>% subset(Protein %in% prots)


fm = do.call(lmerTest::lmer, lmer_args)
LT = getContrast(fm, "Mutant", "Control")

lmerTestContrast(fm, LT) %>% mutate(Contrast = "Mutant-Control") %>% 
	unique() %>% knitr::kable()


###############################################################################
## THESE STATISTICS MATCH WHAT WE PLOT

data(washc_prots)

## Input
m = "M310"
#prots <- washc_prots
prots <- modules[[m]]
data.table(nProts = length(prots)) %>% knitr::kable()

lmer_args <- list()
lmer_args[["formula"]] <- log2(rel_Intensity) ~ 0 + Condition + (1|Protein)
lmer_args[["data"]] <- msstats_prot %>% subset(Protein %in% prots) %>% 
	mutate(Intensity = 2^Abundance) %>%
	group_by(Protein) %>%
	mutate(rel_Intensity = Intensity / sum (Intensity) )
fm <- do.call(lmerTest::lmer, lmer_args)
#summary(fm, ddf="Satterthwaite")
LT <- getContrast(fm, "Mutant","Control")
res <- lmerTestContrast(fm, LT) %>% cleanRes() 
res %>% knitr::kable()

cleanRes <- function(data) {
	data %>% mutate(Contrast = 'Mutant-Control') %>% 
		unique()
}

qqnorm(residuals(fm))
qqline(residuals(fm))

r.squaredGLMM.merMod(fm)

## i fear the the variance estimate for protein may be shrunk

## we can fit subtly differnt models to the same data to evaluate different
## things... that doesnt mean each fit/model is appropriate to answer all
## questions

fx <- scale01(log2(rel_Intensity)) ~ (1|Genotype) + (1|BioFraction) + (1|Protein) + (1|Mixture)
fx <- log2(rel_Intensity) ~ (1|Genotype) + (1|BioFraction) + (1|Protein) + (1|Mixture)
lmer_args[["formula"]] <- fx
lmer_args[["data"]] <- msstats_prot %>% subset(Protein %in% prots) %>% 
	mutate(Intensity = 2^Abundance) %>%
	group_by(Protein) %>%
	mutate(rel_Intensity = Intensity / sum (Intensity) )
fm = do.call(lme4::lmer, lmer_args)
vp <- getVariance(fm)
as.data.table(t(vp/sum(vp))) %>% knitr::kable()


# We scale the evalues to be in the range 0 - 1 for convienence when plotting,
# this may give a visual effect of making things look better than they really
# are
