#!/usr/bin/env Rscript

#' ---
#' title: Swip TMT Proteomics
#' description: analysis of modules for enrichment of WASH proteome
#' authors: Tyler W Bradshaw
#' ---

## Optional parameters:
BF_alpha = 0.10 # Significance threshold.

#---------------------------------------------------------------------
## Misc function - getrd().
#---------------------------------------------------------------------

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

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global imports.
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
	library(geneLists) # for gene lists (pathways)
})

# Load functions in root/R and data in root/data.
suppressWarnings({ devtools::load_all() })

# Directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Load the gene lists from geneLists.
#data(list="iPSD") # iPSD (Arhgef9, Gphn, InSyn1) proteome from Uezu et al.
#data(list="spence2019") # Nascent synapse proteome from Spence et al.
#data(list="gao2020") # Interneuron proteome from Gao et al. -- can we use this?
#data(list="dube2020") # Presynapse proteome from Dube et al. -- can we use this?
data(list="mshallmark") # Mouse hallmark genes from BROAD.
data(list="hlgd") # Human (+mouse) lysosome genes.
#data(list=c("sfariAnimal","sfariGene")) # SFARI ASD database
data(list="synGO") # SynGO Neuron paper from Koopmanns et al.
#data(list="mitocarta2") # Mitochondrial gene/protein database.
#data(list=c("iossifov2012ASD","iossifov2014ASD")) # ASD genes from Iossifov+.
#data(list="wilkinson2017gapgef") # Syngap1 IP 
#data(list="wang2017Epilepsy") # Epilepsy genes
#data(list="weingarten2014AZ") # Active zone proteome
#data(list="boyken2013presynapse") # Presynaptic genes from Boyken et al.
#data(list="lee2017shank3proteome") # Shank3 proteome from Lee et al.,
#data(list="han2013shank3proteome") # Shank3 proteome from Lee et al.,
#data(list="takamori2006SV") # Takamori SV proteome.
#data(list="synsysnet") # SynSysNet synaptic protein database.
data(list="lopitDCpredictions") # Predicted subcellular locations from Geledaki.
#data(list="msgene") # MS genes
#data(list="alsgene") # ALS genes
#data(list="alsodgene") # ALS online database
#data(list="alzgene") # Alzheimer genes
data(list="disgeneLSD") # LSD genes from DisGeneNet.
#data(list="pdgene") # PD genes
#data(list="lelieveld2016ID") # ID genes from levivel et al.,
#data(list="rauch2012ID") # ID genes
#data(list="deligt2012ID") # ID genes
data(list="corum") # CORUM protein complexes.
#data(list="peng2004psd") # Peng 2004 PSD proteome.
#data(list="chen2014AMPARs") # AMPAR proteome.
#data(list="biesemann2014sortedsynaptosomes") # FACS (e/in) sorted synaptosomes.
#data(list="dosemeci2007PSD95proteome") # DLG4 proteome.
data(list="itzhak2017") # cultured neuron spatial proteomics

# Clean up peng data.
#peng2004psd$All <- NULL
#names(peng2004psd) <- paste("Peng et al., 2004:",names(peng2004psd))

# Load the data from root/data.
data(gene_map) # gene mapping data
data(partition) # graph partition
data(tmt_protein) # the proteomics data
data(wash_interactome) # WASH1 BioID from this study, Courtland et al., 2020.
data(sig_modules) # modules with sig DA.

#--------------------------------------------------------------------
## Do work.
#--------------------------------------------------------------------

# Collect list of modules, map Uniprot accession to Entrez..
modules <- split(names(partition),partition)[-1]
module_entrez <- lapply(modules,function(x) { 
				gene_map$entrez[match(x,gene_map$uniprot)] })
all_entrez <- unlist(module_entrez,use.names=FALSE)

# Collect WASH BioID genes.
wash_prots <- unique(wash_interactome$Accession)
wash_genes <- na.omit(gene_map$entrez[match(wash_prots,gene_map$uniprot)])

# Collect Iossifov et al., ASD genes (2 studies).
#iossifov_genes <- c(unique(iossifov2012ASD$iossifov2012ASD,iossifov2014ASD$ASD))

# Clean-up takamori lists.
#takamori2006SV[["All"]] <- NULL
#names(takamori2006SV) <- paste("Takamoir et al., Presynapse:",
#			       names(takamori2006SV))

# Clean up names of synsysnet lists.
#names(synsysnet) <- paste("SynSysNet:",names(synsysnet))

# Clean up lopit dc predictions.
names(lopitDCpredictions) <- paste("LopitDC:",names(lopitDCpredictions))

# Collect list of entrez ids for pathways of interest.
gene_lists <- list(
		   #"ARHGEF9-BioID"=iPSD$Arhgef9,
		   #"GPHNN-BioID" = iPSD$Gphn,
		   #"INSYN1-BioID" = iPSD$InSyn1,
		   #"SYP-BioID" = dube2020$"syp-presynapse",
		   "WASH1-BioID" = unique(wash_genes),
		   "HLGD: Lysosome"=hlgd$LGD,
		   #"SFARI" = c(unique(sfariAnimal$ASD,sfariGene$ASD)),
		   #"MitoCarta2: Mitochondria" = mitocarta2$mitocarta2,
		   #"Wilkinson: Syngap1-IP" = wilkinson2017gapgef$Syngap1,
		   #"Wilkinson: Agap2-IP" =   wilkinson2017gapgef$Agap2,
	       	   #"Wilkinson: Kalrn-IP" =   wilkinson2017gapgef$Kalrn,
		   #"Weingarten Presynaptic AZ" = weingarten2014AZ$weingarten2014AZ,
		   #"Lee et al., Shank3" = unique(unlist(lee2017shank3proteome[c(1,2)])),
		   #"Han et al., Shank3" = unique(han2013shank3proteome$IP),
		   "Lysosome Storage Disorder" = disgeneLSD$"lysosomal storage diseases"
		   #"Multiple Sclerosis" = msgene$MS,
		   #"Alzheimers" = alzgene$ALZ,
		   #"ALS" = unique(c(alsgene$ALS,alsodgene$ALS)),
		   #"Parkisons" = pdgene$PD,
		   #"Chen et al., 2014: AMPAR" = chen2014AMPARs$"AMPAR Proteome",
		   #"Biesmann Sorted Synaptosomes" = biesemann2014sortedsynaptosomes$"FASS Enriched (>2)",
		   #"Dosemeci: DLG4-Proteome" = dosemeci2007PSD95proteome$"Affinity Purified PSD95-Complex"
		   )

# Clarify corum names:
names(corum) <- paste("CORUM:",names(corum))
names(itzhak2017) <- paste("Itzhak et al., 2017:",names(itzhak2017))

# Combine with mouse hallmark genes and some other larger datasets.
gene_lists <- c(gene_lists,
		mshallmark,
		#boyken2013presynapse,
		#takamori2006SV,
		#synsysnet,
		lopitDCpredictions,
		#peng2004psd,
		corum,
		itzhak2017)

# Examine size of pathways.
#knitr::kable(sapply(gene_lists,length))
#message(paste("\nAll pathway sizes:"))
#knitr::kable(sapply(gene_lists,length))

# Remove lists with less than 3 proteins.
drop <- which(sapply(gene_lists,length) < 3)
gene_lists <- gene_lists[-drop]

# Loop to perform GSE for each pathway.
message("\nPerforming GSE analysis for all modules:")
results <- list()
pbar <- txtProgressBar(max=length(gene_lists),style=3)
for (experiment in names(gene_lists)){
	# Get pathway specific genes.
	pathway_genes <- gene_lists[[experiment]]
	# Background is union of all network genes and pathway genes.
	background <- unique(c(all_entrez,pathway_genes))
	# Loop to perform hypergeometric test for enrichment.
	results_list <- list()
	for (i in c(1:length(module_entrez))) {
		results_list[[i]] <- hyperTest(pathway_genes,
					       module_entrez[[i]],
					       background)
	}
	names(results_list) <- paste0("M",names(modules))
	# Collect results in a data.table.
	hyper_dt <- as.data.table(do.call(rbind,results_list),
				  keep.rownames="Module")
	# Adjust p-values.
	hyper_dt$FDR <- p.adjust(hyper_dt$"P-value",method="BH")
	hyper_dt$P.adjust <- p.adjust(hyper_dt$"P-value",method="bonferroni")
	# Add module size annotation.
	sizes <- sapply(module_entrez,length)
	hyper_dt <- tibble::add_column(hyper_dt,"Module Size"=sizes,
				       .after="Module")
	# Add pathway annotation.
	hyper_dt <- tibble::add_column(hyper_dt,Pathway=experiment,
				       .after="Module Size")
	# Add total number of pathway genes.
	hyper_dt <- tibble::add_column(hyper_dt,
				       "Total Pathway Genes"=length(pathway_genes),
				       .after="Pathway")
	hyper_dt <- tibble::add_column(hyper_dt,
				       "N Pathway Genes"=sum(pathway_genes %in% all_entrez),
				       .after="Total Pathway Genes")
	# Number of pathway genes in each module.
	n <- sapply(module_entrez,function(x) sum(pathway_genes %in% x))
	hyper_dt <- tibble::add_column(hyper_dt,
				       "n Pathway Genes in Module"=n,
				       .after="N Pathway Genes")
	# Sort by fold enrichment.
	hyper_dt <- hyper_dt %>% arrange(desc(`Fold enrichment`))
	# Add gene names.
	genes <- getPPIs::getIDs(gene_lists[[experiment]],from="entrez",to="symbol",species="mouse")
	hyper_dt$"Pathway Genes" <- sapply(sapply(module_entrez,function(x) genes[names(genes) %in% x]),paste,collapse=";")
	# Return the results.
	results[[experiment]] <- hyper_dt
	setTxtProgressBar(pbar,value=match(experiment,names(gene_lists)))
}
close(pbar)

# Collect the results in a single data.table.
dt <- bind_rows(results)
sig_dt <- dt %>% filter(P.adjust < BF_alpha)

# Status:
n_mods <- length(unique(sig_dt$Module))
message(paste("\nNumber of modules with something interesting going on:",
	      n_mods))

# How many sig modules have been annotated?
nsig_sig <- sum(sig_modules %in% sig_dt$Module)
message(paste("\nNumber of significantly DA modules with",
	      "something interesting going on:", nsig_sig))

# Pretty print -- all modules.
temp_dt <- sig_dt
temp_dt$Pathway <- substr(temp_dt$Pathway,1,25)
knitr::kable(temp_dt %>% arrange(as.numeric(gsub("M","",Module))))

# Save the data.
myfile <- file.path(rdatdir,"GSEA_Results.csv")
fwrite(sig_dt,myfile)
