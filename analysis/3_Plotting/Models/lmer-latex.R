#!/usr/bin/env Rscript

# author: twab
# title: SwipProteomics
# description: examining lmer models

# references:

# * Hlavac, Marek (2018). stargazer: Well-Formatted Regression and Summary 
#   Statistics Tables. R package version 5.2.2. 
#   https://CRAN.R-project.org/package=stargazer 

# * MuMin package

# * Nakagawa, S., Schielzeth, H. (2013) A general and simple method
#   for obtaining R<U+00B2> from Generalized Linear Mixed-effects Models. Methods
#   in Ecology and Evolution 4: 133<U+2013>

## inputs 
root = "~/projects/SwipProteomics"


## functions ------------------------------------------------------------------

subprot <- function(protein,required=c("msstats_prot")) {
	# subset the protein-level data for a given protein
	# FIXME: error message about env is not informative
	#stopifnot(all(sapply(required,exists)))
	require(dplyr, quietly=TRUE)
	return(msstats_prot %>% dplyr::filter(Protein == protein))
} #EOF

submod <- function(module,required=c("msstats_prot")) {
	# subset the module-level data for a given module
	# FIXME: error message about env is not informative
	require(dplyr, quietly=TRUE)
	#stopifnot(all(sapply(required,exists)))
	stopifnot("Module" %in% colnames(msstats_prot))
	return(msstats_prot %>% dplyr::filter(Module == module))
} #EOF

# NOTE: this doesn't work for some reason... see note below.
#generate_report <- function(protein,gene_map,partition,msstats_prot,figsdir) {
#  # generate a latex report for a given protein
#  # latex tables showing R2 of marginal and conditional types
#  # * marginal represents variance explained by fixed effects, e.g. Genotype
#  # * conditional represents the variance explained by the entire model 
#  #   including both (fixef and mixef).
#  # imports
#  # requires SwipProteomics -- needs R2 functions 
#  require(dplyr, quietly=TRUE)
#  require(stargazer, quietly=TRUE)
#  # annotate msstats with module membership
#  msstats_prot$Module <- paste0("M",partition[msstats_prot$Protein])
#  # get protein's module and gene name
#  module <- paste0("M",partition[protein])
#  gene <- gene_map$symbol[match(protein,gene_map$uniprot)]
#  # fit models
#  fm0 <- lmerTest::lmer("Abundance ~ (1|Mixture) + Condition",subprot(protein)) # | Condition = Genotype.BioFraction
#  fm1 <- lmerTest::lmer("Abundance ~ 0 + (1|BioFraction) + Genotype", subprot(protein))
#  fm2 <- lmerTest::lmer("Abundance ~ 0 + (1|BioFraction) + (1|Protein) + Genotype", submod(module))
#  # if lmerTest was used, then coerce class to something stargazer can work with
#  class(fm0) <- "lmerMod"; class(fm1) <- "lmerMod"; class(fm2) <- "lmerMod"
#  # NOTE: ^fixes error about $ and S4 class (https://stackoverflow.com/questions/31319030)
#  # init latex document
#  output_file <- file.path(figsdir,paste(protein,gene,"report.tex",sep="_"))
#  if (file.exists(output_file)) { unlink(output_file); warning("rm ",output_file) }
#  latex_file <- file(output_file, open="a")
#  writeLines("\\documentclass[11pt]{report}",latex_file)
#  writeLines("\\begin{document}",latex_file)
#  # title
#  writeLines(paste0("\\title{",gene,"|",protein,"}"),latex_file)
#  writeLines("\\author{Tyler W. A. Bradshaw}",latex_file)
#  writeLines("\\maketitle",latex_file)
#  # latex tables of model summaries
#  writeLines(stargazer(fm0,title=as.character(attr(fm0,"call"))[2]),latex_file) # fm0
#  writeLines(stargazer(fm1,title=as.character(attr(fm1,"call"))[2]),latex_file) # fm1
#  writeLines(stargazer(fm2,title=as.character(attr(fm2,"call"))[2]),latex_file) # fm2
#  writeLines("marginal (fixed effects) and conditional (total variance) R2:")
#  writeLines(knitr::kable(r.squaredGLMM.merMod(fm0),format="latex"),latex_file)
#  writeLines(knitr::kable(r.squaredGLMM.merMod(fm1),format="latex"),latex_file)
#  writeLines(knitr::kable(r.squaredGLMM.merMod(fm2),format="latex"),latex_file)
#  # end
#  writeLines("\\end{document}",latex_file)
#  close(latex_file)
#} #EOF

## main ----------------------------------------------------------------------------

# load renv
if (dir.exists(file.path(root,"renv"))) { renv::load(root) }

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(stargazer)
})

# load project
devtools::load_all(quiet=TRUE)

# dir for output
figsdir <- file.path(root,"figs","Models")
if (!dir.exists(figsdir)) { dir.create(figsdir); message("mkdir ",figsdir) }

# load the required data
data(swip)
data(gene_map)
data(partition) #FIXME: where is partition?
data(msstats_prot)

washc = gene_map$uniprot[grep("Washc*",gene_map$symbol)]
partition = setNames(rep(1,length(washc)),nm=washc)
proteins = names(partition)

# NOTE: function doesnt work, environment error, use of NULL env is defunct
# generate_report(swip,gene_map,partition,msstats_prot,figsdir)
#proteins = unique(as.character(msstats_prot$Protein))

pbar <- txtProgressBar(max=length(proteins),style=3)
for (protein in proteins) {
  # generate a latex report for a given protein
  # latex tables showing R2 of marginal and conditional types
  # * marginal represents variance explained by fixed effects, e.g. Genotype
  # * conditional represents the variance explained by the entire model 
  #   including both (fixef and mixef).
  # imports
  # requires SwipProteomics -- needs R2 functions 
  #require(dplyr, quietly=TRUE)
  #require(stargazer, quietly=TRUE)
  # annotate msstats with module membership
  msstats_prot$Module <- paste0("M",partition[msstats_prot$Protein])
  # get protein's module and gene name
  module <- paste0("M",partition[protein])
  gene <- gene_map$symbol[match(protein,gene_map$uniprot)]
  # fit models
  ## FIXME: if any error => pass
  fm0 <- lmerTest::lmer("Abundance ~ (1|Mixture) + Condition",subprot(protein)) 
  #^ | Condition = Genotype.BioFraction
  fm1 <- lmerTest::lmer("Abundance ~ 0 + (1|BioFraction) + Genotype", 
			subprot(protein))
  fm2 <- lmerTest::lmer("Abundance ~0+ (1|BioFraction) + (1|Protein) + Genotype",
			submod(module))
  ## FIXME: breaks if error
  # if lmerTest was used, then coerce class to something stargazer can work with
  class(fm0) <- "lmerMod"; class(fm1) <- "lmerMod"; class(fm2) <- "lmerMod"
  # NOTE: ^fixes error about $ and S4 class 
  # From: (https://stackoverflow.com/questions/31319030)
  # init latex document
  output_file <- file.path(figsdir,paste(protein,gene,"report.tex",sep="_"))
  if (file.exists(output_file)) { unlink(output_file); warning("rm ",output_file) }
  latex_file <- file(output_file, open="a")
  writeLines("\\documentclass[11pt]{report}",latex_file)
  writeLines("\\begin{document}",latex_file)
  # title
  writeLines(paste0("\\title{",gene,"|",protein,"}"),latex_file)
  writeLines("\\author{Tyler W. A. Bradshaw}",latex_file)
  writeLines("\\maketitle",latex_file)
  # latex tables of model summaries
  writeLines(stargazer(fm0,title=as.character(attr(fm0,"call"))[2]),latex_file) 
  # fm0
  writeLines(stargazer(fm1,title=as.character(attr(fm1,"call"))[2]),latex_file) 
  # fm1
  writeLines(stargazer(fm2,title=as.character(attr(fm2,"call"))[2]),latex_file) 
  # fm2
  writeLines("marginal (fixed effects) and conditional (total variance) R2:",
	     latex_file)
  writeLines(knitr::kable(r.squaredGLMM.merMod(fm0),format="latex"),latex_file)
  writeLines(knitr::kable(r.squaredGLMM.merMod(fm1),format="latex"),latex_file)
  writeLines(knitr::kable(r.squaredGLMM.merMod(fm2),format="latex"),latex_file)
  # end
  writeLines("\\end{document}",latex_file)
  close(latex_file)
  setTxtProgressBar(pbar,match(protein,proteins))
} # EOL
