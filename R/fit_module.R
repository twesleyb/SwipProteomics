#' fit_module
#'
fit_module <- function(prots, msstats_prot, gene_map) {

  # subset data for all prots in a module
  subdat <- msstats_prot %>% subset(Protein %in% prots)

  # number of proteins in module
  nprots <- length(unique(subdat$Protein))

  # set factor order (levels)
  # always do this explicitly before plotting to save yourself trouble
  subdat$Genotype <- factor(subdat$Genotype,levels= c("Control","Mutant"))
  subdat$BioFraction <- factor(subdat$BioFraction,
			 levels=c("F4","F5","F6","F7","F8","F9","F10"))

  # for each protein: scale Abundance to max and take median of three* replicates
  #* NOTE: there may be less than three BioReplicates if a protein was identified in
  # less than 3 of the TMT mixtures. 
  df <- subdat %>% group_by(Protein) %>%
	  mutate(scale_Abundance = Abundance/max(Abundance)) %>%
	  group_by(Protein, Genotype, BioFraction) %>% 
	  summarize(med_Abundance = median(scale_Abundance), 
	          SD = sd(scale_Abundance),
	          N = length(scale_Abundance),
	          .groups="drop")

  # calculate coefficient of variation (CV == unitless error)
  df <- df %>% mutate(CV = SD/med_Abundance)

  # get module fitted data by fitting linear model to scaled Abundance
  fx <- "med_Abundance ~ 0 + Genotype:BioFraction + (1|Protein)"
  fm <- lmerTest::lmer(formula=fx, data=df)

  # assess overall contrast
  LT <- getContrast(fm,"Mutant","Control")
  result <- lmerTestContrast(fm,LT) %>% 
	  mutate(Contrast='Mutant-Control') %>% unique() 

  # annotate results with protein info
  idx <- match(prots, gene_map$uniprot)
  gene_prot <- paste(paste(gene_map$symbol[idx],prots,sep="|"),collapse="; ")
  result <- result %>% mutate("Proteins" = gene_prot) %>% mutate('nProts'=nprots) 

  return(result)
} #EOF
