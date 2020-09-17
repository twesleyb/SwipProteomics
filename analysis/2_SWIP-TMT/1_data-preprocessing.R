#!/usr/bin/env Rscript

#' ---
#' title: Swip TMT Proteomics
#' description: Preprocessing and statistical analysis of Swip TMT proteomics.
#' authors: Tyler W A Bradshaw
#' ---

## INPUT data is in root/data/TMT_.zip.
# The zipped directory contains sample meta data and raw peptide data from PD:
zip_file = "TMT.zip" 
input_meta = "TMT-samples.csv"
input_data = "TMT-raw-peptide.csv"

## OPTIONS:
fig_height = 5.0 # Default height of figures (in).
fig_width = 5.0 # Default width of figures (in).

oldham_threshold = 2.5 # Sample connectivity threshold for detecting outliers.
# NOTE: No sample outliers were detected.

FDR_alpha = 0.1 # FDR threshold for differential abundance.
BF_alpha = 0.05 # Bonferroni threshold for differential abundance.
logFC_threshold = c(lwr=log2(0.8), upr=log2(1.2)) # FC threshold.

set.seed = as.numeric(Sys.time()) # seed for random operations.

## OUTPUT saved in root/.../tables:
# * Swip_TMT_Results.xlsx - an excel spreadsheet that contains:
#     - Sample meta data.
#     - Gene/protein identifiers.
#     - Raw peptide data.
#     - Final normalized data (log2 transformed).
#     - Statistical results for intrafraction contrasts (F4-10).
#     - Statistical results for WT v MUT contrast.

## OUTPUT for R package in root/data.
# Key datasets are saved in the data directory as read-only objects.
# * tmt_protein.rda
# NOTE: can be loaded in R with data(tmt_protein).

## OUTPUT for downstream analysis in root/rdata/
# Input files/temporary files that are passed to other scripts including 
# the Python Leidenalg clustering script are saved in root/rdata.

## Order of data processing operations:
# * Load the raw peptide intensity data from PD 2.2.
# * Initial sample loading normalization.
# * Impute missing peptide values with KNN algorithm (k=10).
# * Examine QC peptide reproducibility - remove QC outliers.
# * Summarize to protein level by summing all peptides for a protein.
# * Protein-level sample loading normalization.
# * IRS normalization - equalize protein measurements from different peptides.
# * Sample pool normalization - assumes that mean of pooled QC technical 
#   replicates is equal to mean of all biological replicates.
# * Protein level filtering - remove proteins identified by a single peptide;
#   remove proteins with too many missing values; remove proteins that are not 
#   reproducible (exhibit high inter-experimental variablility across the 3x 
#   biological replicates).
# * Assess intra-fraction differential abundance with a general linear model.
#   GLM: [Protein Abundance] ~ 0 + Fraction.Genotype.
# * Assess WT v MUT differential abundance across all fractions with a 
#   GLM: [Protein Abundance] ~ Fraction + Genotype.

## NOTE:
# The dependencies for this project were managed within a conda virtual 
# environment into which R and renv, an R package for managing R libraries,
# were installed. All additional R dependencies were installed using 
# renv. The anRichment package was installed from source (see installation
# script in root/).

#---------------------------------------------------------------------
## Misc function - getrd().
#---------------------------------------------------------------------

## NOTE: The functions used in this script are not robust. They were written
## to work with the input arguments that they are provided, and in many cases 
## will not perform as expected if passed different arguments. I attempted to
## keep the data in a tidy-ish format throughout. This decision makes some
## operations like plotting easier, but makes other operations like
## normalization more cumbersome and computationally costly.

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

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. 

# Load renv -- use renv::load NOT activate!
rootdir <- getrd()
renv::load(rootdir,quiet=TRUE) 

# Load required packages and functions.
suppressPackageStartupMessages({
	library(dplyr) # For manipulating data.
	library(data.table) # For working with tables.
	library(getPPIs) # For PPIs and mapping gene identifiers.
})

# Load project specific functions and data.
suppressWarnings({ devtools::load_all() })

# Project directories:
datadir <- file.path(rootdir, "data") # Key pieces of data saved as rda.
rdatdir <- file.path(rootdir, "rdata") # Temporary data files.
downdir <- file.path(rootdir, "downloads") # Misc downloads/temporary files.
tabsdir <- file.path(rootdir, "tables") # Output tables saved as excel files.

# Create project output directories if necessary.
# FIXME: 
if (!dir.exists(datadir)) { dir.create(datadir) }
if (!dir.exists(downdir)) { dir.create(downdir) }
if (!dir.exists(rdatdir)) { dir.create(rdatdir) }
if (!dir.exists(tabsdir)) { dir.create(tabsdir) }

#---------------------------------------------------------------------
## Load the raw data and sample info.
#---------------------------------------------------------------------

# Extract the raw TMT data from zipped file.
myfile <- file.path(datadir,zip_file)
unzip(myfile, exdir=downdir) # unzip into root/downloads/

# Read TMT.csv data into R with data.table::fread.
myfile <- file.path(downdir,tools::file_path_sans_ext(zip_file),input_data)
peptides <- fread(myfile)

# Load sample information.
myfile <- file.path(downdir,tools::file_path_sans_ext(zip_file),input_meta)
samples <- fread(myfile)

# Format Cfg force column -- this is the Force in g's used to obtain 
# the subcellular fraction.
samples$"Cfg Force (xg)" <- formatC(samples$"Cfg Force (xg)",big.mark=",")

rm(list="myfile") # remove temporary variables.

#---------------------------------------------------------------------
## Map all Uniprot accession numbers to stable entrez IDs.
#---------------------------------------------------------------------

# Map Uniprot IDs to Entrez using MGI batch query function.
# Map Entrez IDs to gene symbols using getPPIs::getIDs().
message("\nCreating gene identifier map.")

# First, remove any non-mouse proteins from the data.
peptides <- peptides %>% filter(grepl("OS=Mus musculus",Description))

# Remove these Immunoglobin proteins:
ig_prots <- c("P01631","P01646","P01665","P01680","P01746","P01750",
	      "P01786","P01864","P01878","P03975","P06330","P03987")
peptides <- peptides %>% filter(Accession %notin% ig_prots)

# Collect all Uniprot IDs.
uniprot <- unique(peptides$Accession)

# Map Uniprot IDs to Entrez using online MGI batch query function.
# This takes a couple minutes because the function currently
# downloads the MGI data each time the function is called.
entrez <- mgi_batch_query(uniprot,quiet=TRUE)
names(entrez) <- uniprot

# Map any remaining missing IDs by hand.
message("Mapping missing IDs by hand.\n")
missing <- entrez[is.na(entrez)]
mapped_by_hand <- c(P05214=22144, P0CG14=214987, P84244=15078)
entrez[names(mapped_by_hand)] <- mapped_by_hand

# Check: Have we successfully mapped all Uniprot IDs to Entrez?
check <- sum(is.na(entrez)) == 0
if (!check) { stop("Unable to map all UniprotIDs to Entrez.") }

# Map entrez ids to gene symbols using twesleyb/getPPIs.
symbols <- getPPIs::getIDs(entrez,from="entrez",to="symbol",species="mouse")

# Create gene identifier mapping data.table.
gene_map <- data.table(uniprot = names(entrez),
                       entrez = entrez,
	               symbol = symbols)
gene_map$id <- paste(gene_map$symbol,gene_map$uniprot,sep="|")

rm(list=c("ig_prots","uniprot","entrez","missing",
	  "mapped_by_hand", "check","symbols"))

#---------------------------------------------------------------------
## Tidy-up the input data from Proteome Discover.
#---------------------------------------------------------------------

# Convert PD df into tidy df. 
# Samples should contain the following columns:
# Treatment, Channel, Sample, Experiment
message("\nLoading raw data from Proteome Discover (PD2.2).")
cols <- colnames(peptides)[!grepl("Abundance",colnames(peptides))]
tidy_peptide <- tidyProt(peptides,id.vars=cols)

# Annotate tidy data with additional meta data from samples.
tidy_peptide <- left_join(tidy_peptide,samples,by="Sample")

# Summary of peptide/protein quantification:
message(paste("\nSummary of initial peptide/protein quantification",
	      "after removing contaiminants:"))
n_samples <- length(unique(tidy_peptide$Sample))
n_proteins <- length(unique(tidy_peptide$Accession))
n_peptides <- length(unique(tidy_peptide$Sequence))

df <- data.frame("Samples"=as.character(n_samples),
		 "Proteins"=formatC(n_proteins,big.mark=","),
		 "Peptides"=formatC(n_peptides,big.mark=","))
knitr::kable(df)

rm(list=c("cols","df","n_samples","n_proteins","n_peptides"))

#---------------------------------------------------------------------
## Perform sample loading normalization.
#---------------------------------------------------------------------

# Perform sample normalization. Normalization is done for each 
# experiment independently (group by Experiment:Sample).
# NOTE: Grouping by Experiment, Channel does't work because 
# Sample TMT channels (e.g. 126N were used in different expirements).
message("\nPerforming sample loading normalization.")
sl_peptide <- normSL(tidy_peptide, groupBy=c("Experiment","Sample"))

# Check the data before SL:
message("\nPeptide data before sample loading normalization:")
check_SL(tidy_peptide)

# Check the data after SL:
message("\nPeptide data after sample loading normalization:")
check_SL(sl_peptide) # Equal within an experiment.

#---------------------------------------------------------------------
## Impute missing peptide values.
#---------------------------------------------------------------------

# Impute missing peptide values with k-nearest neighbors (KNN) algorithm.
# * Missing QC values will not be imputed.
# * Peptides (rows) with more than 50% missingness will not be imputed.
# Values in these rows are masked (replaced with NA).
# NOTE: KNN imputing is done for each experimental group seperately.
message("\nImputing a small number of missing peptide values.")
imputed_peptide <- imputeKNNpep(sl_peptide, groupBy="Experiment",
				samples_to_ignore="SPQC",quiet=FALSE) 

#---------------------------------------------------------------------
## Examine reproducibility of QC measurements.
#---------------------------------------------------------------------
# Assess reproducibility of QC measurements and remove QC samples
# that are irreproducible.
# This strategy was adapted from Ping et al., 2018 (pmid: 29533394).
# For each experiment, the ratio of QC measurments is calculated. 
# These ratios are then binned based on average Intensity into
# 5 bins. For each bin, measuremnts that are outside 
# +/- 4x standard deviations from the bin's mean are removed. 

message("\nRemoving peptides with irreproducible QC measurements.")
filt_peptide <- filtQC(imputed_peptide,controls="SPQC",quiet=FALSE)

#---------------------------------------------------------------------
## Summarize to protein level.
#---------------------------------------------------------------------

# Protein intensity is just the sum of all peptides for that protein.
message("\nSummarizing proteins as the sum of their peptides.")
proteins <- summarize_prot(filt_peptide)

# Perform SL normalization across all experiments (grouped by Sample).
message("\nPerforming sample loading normalization between experiments.")
sl_protein <- normSL(proteins, groupBy="Sample")

# Column sums are equal across all experiments:
message("\nProtein data after sample loading normalization:")
check_SL(sl_protein)

#---------------------------------------------------------------------
## Perform IRS Normalization.
#---------------------------------------------------------------------

# Equalize QC measurements between experiments. Adjusts protein 
# measurements of biological replicates simultaneously.
# Accounts for protein quantification by different peptides in 
# each experiment.
message(paste("\nStandardizing protein measurements between",
	      "experiments by IRS normalization."))
irs_protein <- normIRS(sl_protein,controls="SPQC",robust=TRUE)

# Check, protein-wise QC measurements are now equal in a random protein.
# Generate list of QC proteins, and check the average of the three Exps.
qc_proteins <- irs_protein %>% filter(Treatment == "SPQC") %>% 
	group_by(Accession,Experiment) %>% 
	dplyr::summarize(Treatment = unique(Treatment),
		  "Mean(Intensity)" = mean(Intensity,na.rm=TRUE),
		  .groups = "drop") %>% group_by(Accession) %>% group_split()
message("\nIntra-experimental QC means are equal after IRS normalization:")
knitr::kable(sample(qc_proteins,1))

#---------------------------------------------------------------------
## Perform Sample Pool Normalization.
#---------------------------------------------------------------------

# QC samples were generated by pooling all biological replicates.
# Normalize mean of all QC samples to be equal to the mean of all 
# biological replicates.
message(paste("\nAccounting for experimental batch effect by",
	      "performing sample pool normalization."))
spn_protein <- normSP(irs_protein,pool=c("Control","Mutant"))

# Check, the mean of sample pool and biological replicates are now equal.
message(paste("\nThe mean of biological replicates and",
	      "pooled QC samples are now equal:"))
protein_list <- spn_protein %>% group_by(Accession)  %>% group_split()

# Drop any proteins with NA.
idx <- sapply(protein_list,function(x) any(is.na(x)))
protein_list <- protein_list[!idx]

# Get a random protein's data as an example.
protein <- sample(protein_list,1)[[1]]
protein$Sample_Pool <- protein$Treatment == "SPQC"
df <- protein %>% group_by(Sample_Pool) %>% 
	dplyr::summarize(Accession = unique(Accession),
		  Treatement=paste(unique(Treatment),collapse=" + "),
		  "Mean(Intensity)"=mean(Intensity), .groups="drop")
knitr::kable(df)

rm(list=c("df","protein_list","idx","protein"))

#---------------------------------------------------------------------
## Protein level filtering.
#---------------------------------------------------------------------

# Remove proteins that:
# * Were identified by a single peptide.
# * Contain too many (>50%) missing values.
# * Contain any missing QC values.
# * Have outlier protein measurements: 
#   mean(logRatio(replicates)) outside +/- nSD * mean(binned logRatio))
message(paste("\nFiltering proteins, this may take several minutes."))
filt_protein <- filtProt(spn_protein,
			 controls="SPQC",nbins=5,nSD=4,summary=TRUE)

# At this point there are no remaining missing values.
check <- sum(is.na(filt_protein$Intensity)) == 0
if (!check) { stop("Why are there missing values!") }

rm(list="check")

#---------------------------------------------------------------------
## Protein level imputing.
#---------------------------------------------------------------------

# If necessary, protein level imputing can be performed.
# Proteins with more than 50% missing values are ignored.

message(paste("\nImputing missing protein values."))
imputed_protein <- imputeKNNprot(filt_protein,ignore="SPQC")

#---------------------------------------------------------------------
## Identify Sample outliers.
#---------------------------------------------------------------------

# Calculate Oldham's normalized sample connectivity (zK) in order to 
# identify outlier samples.
# This approach was adapted from Oldham et al., 2012 (pmid: 22691535).
zK <- sampleConnectivity(filt_protein)
outlier_samples <- c(names(zK)[zK < -1 * oldham_threshold],
		     names(zK)[zK > +1 * oldham_threshold])

# There are no sample outliers.
check <- length(outlier_samples) == 0
if (check) { 
	message("\nNo sample level outliers were identified.")
} else { 
	stop("Why are there outlier samples?") 
}

#---------------------------------------------------------------------
## Intrafraction protein differential abundance.
#---------------------------------------------------------------------
# Use EdgeR glm to evaluate differential protein abundance for 
# WT v Mutant comparisons within a fraction.
# NOTE: Model = ~ Fraction + Treatment | Treatment = Genotype.Fraction
message(paste("\nEvaluating intra-fraction protein differential abundance."))

# Do the stastical comparisons with EdgeR glmFit and glmFTest.
comparisons <- "Genotype.Fraction"
glm_results <- glmDA(filt_protein,comparisons,samples,
		 samples_to_ignore="SPQC")

# Collect the results.
dge <- glm_results$dge
results_list <- glm_results$stats
results_list <- results_list[c("F4","F5","F6","F7","F8","F9","F10")]

# Annotate final normalized protein with sample meta data.
tmt_protein <- left_join(filt_protein,samples,
			 by=intersect(colnames(filt_protein),colnames(samples)))

# Summary of DA proteins:
message(paste0("\nSummary of differentially abundant proteins ",
	      "in each subceulluar fraction (FDR < ",FDR_alpha,"):"))
knitr::kable(t(sapply(results_list,function(x) sum(x$FDR < FDR_alpha))))

# Total number of unique DA proteins:
all_sig <- lapply(results_list, function(x) x$Accession[x$FDR < FDR_alpha])
n_sig <- length(unique(unlist(all_sig)))
message(paste("\nTotal number of unique DA proteins:",n_sig))

# Proteins that are commonly dysregulated:
combined_results <- rbindlist(results_list,idcol="Fraction")
idx <- match(combined_results$Accession,gene_map$uniprot)
combined_results$Gene <- gene_map$symbol[idx] 
df <- combined_results %>% group_by(Accession) %>% 
	dplyr::summarize(Gene = unique(Gene),
			 nSig = sum(FDR<0.1),
			 nFractions = length(FDR), 
			 .groups = "drop") %>% 
	filter(nSig==nFractions)

# Status:
message(paste("\nProteins that are differentially abundant in all fractions:\n",
	      paste(df$Gene,collapse=", ")))

rm(list=c("df","all_sig","n_sig","comparisons"))

#----------------------------------------------------------------------
## Assess differential abundance between WT vs Mutant groups.
#----------------------------------------------------------------------
# Use EdgeR glm to evaluate differential protein abundance.
# Compare WT v Mutant within a fraction.
# NOTE: Model: ~ Fraction + Genotype; We are interested in WT v MUT changes.
message(paste("\nEvaluating protein differential abundance between",
	      "WT and Mutant groups."))

# Do statistical comparisons with EdgeR.
comparisons <- "Genotype.Fraction"
samples_to_ignore = "SPQC"
alt_glm_results <- glmDA2(tmt_protein, comparisons, samples, samples_to_ignore)

# Extract key results.
alt_results <- alt_glm_results$stats

# Adjust pvalues using Bonferroni method.
PAdjust <- p.adjust(alt_results$PValue,method="bonferroni")
alt_results <- tibble::add_column(alt_results,PAdjust,.after="FDR")

# Summary of DA proteins with any logFC.
#bounds <- logFC_threshold
sig1 <- alt_results$FDR < 0.05
sig2 <- alt_results$PAdjust < BF_alpha 
#updown <- alt_results$logFC < bounds['lwr'] | alt_results$logFC < bounds['upr']
message("\nSummary of WT v MUT DA proteins:")
df <- data.table("FDR < 0.05"=sum(sig1),
		 "BF < 0.05"=sum(sig2))
knitr::kable(df)

# Column names are Adjusted.NAME
idy <- colnames(alt_results) %notin% c("Accession","Entrez","Symbol")
colnames(alt_results)[idy] <- paste0("Adjusted.",
				     colnames(alt_results)[idy])

#---------------------------------------------------------------------
## Calculate Protein abundance adjusted for fraction differences.
#---------------------------------------------------------------------

# NOTE: these adjusted values are not used for modeling, only 
# plotting purposes.

# Calculate protein abundance, adjusted for fraction differences.
# FIXME: Is CPM step needed? Log is necessary, but CPM norm?
logCPM <- edgeR::cpm(glm_results$dge, log=TRUE)

# Remove effect of fraction.
dm <- limma::removeBatchEffect(logCPM,
			       batch=dge$samples$fraction,
			       design=model.matrix(~treatment,data=dge$samples))

# Collect results.
dt <- as.data.table(dm,keep.rownames="Accession") %>% 
		reshape2::melt(id.var="Accession",
			       variable.name="Sample",
			       value.name="Adjusted.Intensity")
dt$Sample <- as.character(dt$Sample)

# Combine Adjusted protein with additional sample meta data.
tmp_prot <- tmt_protein %>% filter(Treatment != "SPQC") %>% 
	as.data.table()
adjusted_prot <- left_join(tmp_prot,dt,
			   by=intersect(colnames(tmp_prot),colnames(dt)))

# Add stats to adjusted protein data.
adjusted_prot <- left_join(adjusted_prot,alt_results,
			   by=intersect(colnames(adjusted_prot),
					colnames(alt_results)))

# Unlog the data.
adjusted_prot$Adjusted.Intensity <- 2^adjusted_prot$Adjusted.Intensity

#---------------------------------------------------------------------
## Combine final normalized TMT data and stats.
#---------------------------------------------------------------------

# Collect stats.
temp_dt <- bind_rows(results_list,".id"="Fraction")

# Merge tmt_protein data and glm stats.
tmt_protein <- left_join(tmt_protein,temp_dt,
			 by=intersect(colnames(tmt_protein),colnames(temp_dt)))

# Remove QC data.
tmt_protein <- tmt_protein %>% filter(Treatment != "SPQC") %>% 
	as.data.table()

# Combine tmt_protein with adjusted data and stats.
idy <- intersect(colnames(tmt_protein),colnames(adjusted_prot))
tmt_protein <- left_join(tmt_protein,adjusted_prot,by=idy)

# Annotate with gene identifiers.
idx <- match(tmt_protein$Accession,gene_map$uniprot)
tmt_protein$Entrez <- gene_map$entrez[idx]
tmt_protein$Symbol <- gene_map$symbol[idx]

# Organize the columns -- select columns of interest.
tmt_protein <- tmt_protein %>% dplyr::select(Experiment, Sample, Channel,
				      Treatment, Genotype, Fraction,
				      "Cfg Force (xg)",
				      "Date of Brain Fractionation",
				      "Mouse ID", Sex, DOB, "Age (mo)",
				      Accession, Entrez, Symbol, Peptides,
				      Intensity, Adjusted.Intensity, 
				      logFC, Adjusted.logFC,
				      PercentWT, Adjusted.PercentWT,
				      F, Adjusted.F,
				      PValue, Adjusted.PValue,
				      FDR, Adjusted.FDR)

rm(list=c("idx","idy"))

#---------------------------------------------------------------------
## Save the TMT data and stats as a single excel workbook
#---------------------------------------------------------------------
## Save TMT data as a single excel document.
## Create an excel workbook with the following sheets:
# * Samples
# * Raw Peptide
# * Norm Protein
# * Statistical results:
#     - Seperate sheet for each intra-fraction comparison.
#     - Include contrast specific data.

# Add normalized protein data to statistical results.
norm_df <- tmt_protein %>% as.data.table() %>%
	dcast(Accession ~ Sample,value.var="Intensity") 
norm_dm <- norm_df %>% as.matrix(rownames="Accession")

# Loop to add normalized protein data to glm statistical results.
message("\nSaving TMT data and statistical results.")
for (i in 1:length(results_list)){
	# Get the relevant data.
	df <- results_list[[i]]
	namen <- names(results_list)[i]
	subsamples <- samples$Sample[grepl(namen,samples$Fraction)]
	idy <- match(subsamples,colnames(norm_dm))
	dm <- norm_dm[,idy]
	# Sort the data by Exp.Sample
	ids <- sapply(sapply(strsplit(colnames(dm),"\\."),"[",c(6,7),
			     simplify=FALSE), paste,collapse=".")
	names(ids) <- colnames(dm)
	ids <- ids[order(ids)]
	dm <- dm[,names(ids)]
	# Coerce dm to dt and merge with stats df.
	dt <- as.data.table(dm,keep.rownames="Accession")
	dt_out <- left_join(df,dt,by="Accession") %>% as.data.table
	results_list[[i]] <- dt_out
} # Ends loop.

# Update data with Mean and SEM of WT and Mutant groups.
for (namen in names(results_list)) {
	df <- results_list[[namen]]
	idy <- grep("Abundance",colnames(df))
	dm <- df %>% dplyr::select(Accession, all_of(idy)) %>% 
		as.matrix(rownames="Accession")
	idy <- grep("Control",colnames(dm))
	WT_means <- apply(dm,1,function(x) log2(mean(x[idy])))
	WT_SEM <- apply(dm,1,function(x) log2(sd(x[idy])))/WT_means
	idy <- grep("Mutant",colnames(dm))
	MUT_means <- apply(dm,1,function(x) log2(mean(x[idy])))
	MUT_SEM <- apply(dm,1,function(x) log2(sd(x[idy])))/MUT_means
	df <- tibble::add_column(df,'WT Mean' = WT_means, .after = "Accession")
	df <- tibble::add_column(df,'WT SEM' = WT_SEM, .after = "WT Mean")
	df <- tibble::add_column(df,'MUT Mean' = MUT_means, .after = "WT SEM")
	df <- tibble::add_column(df,'MUT SEM' = MUT_SEM, .after = "MUT Mean")
	results_list[[namen]] <- as.data.table(df)
}

# Add adjusted protein values to alt stats.
df <- tmt_protein %>% as.data.table() %>%
	dcast(Accession ~ Sample, value.var="Adjusted.Intensity")
alt_results <- left_join(alt_results,df,by="Accession")
colnames(alt_results) <- gsub("Abundance","Adjusted.Abundance",
			      colnames(alt_results))

# Sort the columns.
idy <- grepl("Abundance",colnames(alt_results))
col_names <- gsub("Adjusted\\.","",colnames(alt_results)[idy])
idx <- match(col_names,samples$Sample)
e <- factor(samples$Experiment[idx],levels=c("Exp1","Exp2","Exp3"))
f <- factor(samples$Fraction[idx],levels=c("F4","F5","F6","F7","F8","F9","F10"))
g <- factor(samples$Genotype[idx],levels=c("WT","MUT"))
names(col_names) <- paste(e,g,f,sep=".")
col_order <- as.character(interaction(e,g,f))
col_order <- c(col_order[grepl("WT",col_order)],
	       col_order[grepl("MUT",col_order)])
sorted_cols <- paste("Adjusted",as.character(col_names[col_order]),sep=".")
all_cols <- c(colnames(alt_results)[colnames(alt_results) %notin% sorted_cols],
	      sorted_cols)
alt_results <- alt_results %>% dplyr::select(all_of(all_cols)) %>% 
	as.data.table()

# Calculate group mean and SEM.
df <- alt_results
idy <- grep("Abundance",colnames(df))
dm <- df %>% dplyr::select(Accession, all_of(idy)) %>% 
	as.matrix(rownames="Accession")
idy <- grep("Control",colnames(dm))
WT_means <- apply(dm,1,function(x) log2(mean(x[idy])))
WT_SEM <- apply(dm,1,function(x) log2(sd(x[idy])))/WT_means
idy <- grep("Mutant",colnames(dm))
MUT_means <- apply(dm,1,function(x) log2(mean(x[idy])))
MUT_SEM <- apply(dm,1,function(x) log2(sd(x[idy])))/MUT_means
df <- tibble::add_column(df,'WT Mean' = WT_means, .after = "Accession")
df <- tibble::add_column(df,'WT SEM' = WT_SEM, .after = "WT Mean")
df <- tibble::add_column(df,'MUT Mean' = MUT_means, .after = "WT SEM")
df <- tibble::add_column(df,'MUT SEM' = MUT_SEM, .after = "MUT Mean")
alt_results <- as.data.table(df)

# Add intra-fraction comparisons to list of results.
results_list[["WT v MUT"]] <- alt_results

# Annotate each df with gene identfifiers.
final_results <- list()
for (i in c(1:length(results_list))) {
	tmp_df <- results_list[[i]]
	idx <- match(tmp_df$Accession,gene_map$uniprot)
	Entrez <- gene_map$entrez[idx]
	Symbol <-  gene_map$symbol[idx]
	tmp_df <- tibble::add_column(tmp_df,Entrez,.after="Accession")
	tmp_df <- tibble::add_column(tmp_df,Symbol,.after="Entrez")
	final_results[[i]] <- tmp_df
}

# Save as excel workboook.
names(final_results) <- paste(names(results_list),"Results")
final_results <- c(list("Samples" = samples),
		   list("Raw Peptide" = peptides),
		   list("Norm Protein" = norm_df), final_results)
myfile <- file.path(tabsdir,"S2_Swip_TMT_Protein_GLM_Results.xlsx")
write_excel(final_results,myfile,rowNames=FALSE)

#---------------------------------------------------------------------
## Save output for downstream analysis.
#---------------------------------------------------------------------
# Save key results.

# Save raw data -- tidy_peptide.
myfile <- file.path(rdatdir,"tidy_peptide.csv")
fwrite(tidy_peptide,myfile)

# Save raw peptide data for analysis with MSstats or other workflows.
myfile <- file.path(rdatdir,"tidy_peptide.rda")
save(tidy_peptide,file=myfile,version=2)

# Save tmt_protein
myfile <- file.path(rdatdir,"tmt_protein.csv")
fwrite(tmt_protein,myfile)

# Save statistical results.
myfile <- file.path(rdatdir,"glm_results.RData")
saveRDS(glm_results,myfile)

# Save gene map.
myfile <- file.path(datadir,"gene_map.rda")
save(gene_map,file=myfile,version=2)

# Save tidy_protein (final normalized protein in tidy format) as rda object. 
myfile <- file.path(datadir,"swip_tmt.rda")
swip_tmt <- tmt_protein; save(swip_tmt,file=myfile,version=2)

message("\nDone!")
