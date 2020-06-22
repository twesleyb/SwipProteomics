#!/usr/bin/env Rscript

#' ---
#' title: Swip Proteomics
#' description: Preprocessing of Swip TMT proteomics data.
#' authors: Tyler W Bradshaw
#' ---

## NOTE: The functions used in this script are not robust. They were written
## to work with the input arguments that they are provided, and in many cases 
## will not perform as expected if passed different arguments.

## INPUT data in root/data/:
zip_file = "TMT.zip" # Zipped directory contains:
input_data = "TMT-raw-peptide.csv"
input_samples = "TMT-samples.csv"

## OPTIONS:
save_plots = TRUE # Should plots be saved in root/figs?
FDR_alpha = 0.1 # FDR threshold for differential abundance.
fold_change_delta = 0.2 # Fold change threshold.
sample_connectivity_threshold = 2.5 # Threshold for sample level outliers.
fig_height = 5.0 # Default height of figures.
fig_width = 5.0 # Default width of figures.

## OUTPUT saved in root/tables:
# * Swip_TMT_Results.xlsx

## OUTPUT for R package in root/data.
# * tmt_protein.rda

## Order of data processing operations:
# * Load the data from PD.
# * Initial sample loading normalization -- done within an experiment.
# * Impute missing peptide values with KNN algorithm (k=10).
# * Examine QC peptide reproducibility -- remove outliers.
# * Summarize to protein level by summing all peptides for a protein.
# * Protein-level sample loading normalization -- done across all experiments.

# * IRS normalization -- equalizes protein measurements made from different
#   peptides.
# * Sample pool normalization -- assumes that mean of QC technical replicates is
#   equal to mean of biological replicates..

# * Impute missing protein values with KNN (k=10).
# * Protein level filtering -- remove proteins identified by a single peptide;
#   remove proteins with too many missing values; remove proteins that are not 
#   reproducible (exhibit high inter-experimental variablility across the 3x 
#   biological replicates.
# * Asses differential abundance with a general linear model 
#   ([Abundance] ~ 0 + groups) as implemented by the Edge R package 
#   and its functions gmQLFfit() and glmQLFTest().

#---------------------------------------------------------------------
## Prepare the workspace.
#---------------------------------------------------------------------
# Prepare the R workspace for the analysis. 

# Load renv -- use renv::load NOT activate!
rootdir <- getrd()
renv::load(rootdir,quiet=TRUE) # NOTE: getrd is a f(x) in .Rprofile.

# Load required packages and functions.
suppressPackageStartupMessages({
	library(dplyr) # For manipulating data.
	library(ggplot2) # For making plots.
	library(data.table) # For working with tables.
})

# Load project specific functions and data.
suppressWarnings({ devtools::load_all() })

# Project directories:
datadir <- file.path(rootdir, "data") # Key pieces of data saved as rda.
fontdir <- file.path(rootdir, "fonts") # Arial font for plots.
rdatdir <- file.path(rootdir, "rdata") # Temporary data files.
tabsdir <- file.path(rootdir, "tables") # Output tables.
downdir <- file.path(rootdir, "downloads") # Misc stuff.
figsdir <- file.path(rootdir, "figs","TMT") # Output figures.

# Create output directories if necessary.
if (!dir.exists(datadir)) { dir.create(datadir) }
if (!dir.exists(downdir)) { dir.create(downdir) }
if (!dir.exists(rdatdir)) { dir.create(rdatdir) }
if (!dir.exists(tabsdir)) { dir.create(tabsdir) }
if (!dir.exists(figsdir)) { dir.create(figsdir, recursive = TRUE) }

# Set global plotting settings.
ggtheme()
set_font("Arial", font_path = fontdir)

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
myfile <- file.path(downdir,tools::file_path_sans_ext(zip_file),input_samples)
samples <- fread(myfile)

# Format Cfg force column -- this is the Force in g's used to obtain 
# the subcellular fraction.
samples$"Cfg Force (xg)" <- formatC(samples$"Cfg Force (xg)",big.mark=",")

rm(list="myfile")

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

# Collect all uniprot IDs.
uniprot <- unique(peptides$Accession)

# Map Uniprot IDs to entrez using online MGI batch query function.
# This takes a couple minutes because the function currently
# downloads the MGI data each time the function is called.
entrez <- mgi_batch_query(uniprot,quiet=TRUE)
names(entrez) <- uniprot

# Map any remaining missing ids by hand.
message("Mapping missing IDs by hand.\n")
missing <- entrez[is.na(entrez)]
mapped_by_hand <- c(P05214=22144, P0CG14=214987, P84244=15078)
entrez[names(mapped_by_hand)] <- mapped_by_hand

# Check: we have successfully mapped all uniprot ids.
check <- sum(is.na(entrez)) == 0
if (!check) { stop("Unable to map all Uniprot IDs to stable gene identifiers!") }

# Map entrez ids to gene symbols.
symbols <- getPPIs::getIDs(entrez,from="entrez",to="symbol",species="mouse")

# Create mapping data.table.
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
message("\nLoading raw data from Proteome Discover.")
cols <- colnames(peptides)[!grepl("Abundance",colnames(peptides))]
tidy_peptide <- tidyProt(peptides,id.vars=cols)

# Annotate tidy data with additional meta data from samples.
tidy_peptide <- left_join(tidy_peptide,samples,by="Sample")

rm(list="cols")

#---------------------------------------------------------------------
## Perform sample loading normalization.
#---------------------------------------------------------------------

# Perform sample normalization. Normalization is done for each 
# experiment independently (group by Experiment:Sample).
#  NOTE: Grouping by Experiment, Channel does't work because 
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
# NOTE: KNN imputing is done on each experimental group seperately.
message("\nImputing a small number of missing peptide values.")
imputed_peptide <- imputeKNNpep(sl_peptide, groupBy="Experiment",
				samples_to_ignore="SPQC",quiet=FALSE) 

# Generate missing value density plots.
tp <- sl_peptide
df <- tp %>% filter(Treatment != "QC") %>% filter(Experiment == "Exp1") %>%
	group_by(Accession,Sequence) %>%
	dplyr::summarize(Mean_Intensity = mean(Intensity,na.rm=TRUE), 
			      Contains_Missing = any(is.na(Intensity)),
			      .groups = "drop") %>%
	filter(!is.na(Mean_Intensity)) # Remove peptides with all missing values.
plot <- ggplot(df, aes(x=log2(Mean_Intensity),
		       fill=Contains_Missing,colour=Contains_Missing)) +
	geom_density(alpha=0.1,size=1) + ggtitle("Experiment 1")

# NOTE: From this plot, the distributions overlap. Missing values are
# missing at random or even missing completely at random.

if (save_plots) {
	myfile <- file.path(figsdir,"Peptide_MV_Density.pdf")
	ggsave(myfile,plot,height=fig_height,width=fig_width)
}

rm(list=c("tp","df","myfile","plot"))

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

# Perform SL normalization across all experiments (groupBy Sample).
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

# Generate list of QC proteins.
qc_proteins <- irs_protein %>% filter(Treatment == "SPQC") %>% 
	group_by(Accession,Experiment) %>% 
	dplyr::summarize(Treatment = unique(Treatment),
		  "Mean(Intensity)" = mean(Intensity,na.rm=TRUE),
		  .groups = "drop") %>% group_by(Accession) %>% group_split()

# Check, protein-wise QC measurements are now equal.
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
# * were identified by a single peptide.
# * contain too many missing values.
# * contain any missing QC values.
# * have outlier measurements.
message(paste("\nFiltering proteins, this may take several minutes."))
filt_protein <- filtProt(spn_protein,
			 controls="SPQC",nbins=5,nSD=4,summary=TRUE)

# There should be no missing values at this point.
check <- sum(is.na(filt_protein$Intensity)) == 0
if (!check) { stop("Why are there missing values!") }

rm(list="check")

#---------------------------------------------------------------------
## Protein level imputing.
#---------------------------------------------------------------------

# If necessary, protein level imputing can be performed.

# Impute missing values.
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
outlier_samples <- c(names(zK)[zK < -sample_connectivity_threshold],
		     names(zK)[zK > sample_connectivity_threshold])

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

# Use EdgeR glm to evaluate differential protein abundance.
# Compare WT v Mutant within a fraction.
message(paste("\nEvaluating intra-fraction protein differential abundance."))
# NOTE: Model = ~ Fraction + Treatment | Treatment = Genotype.Fraction
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

# Total number of unique DA proteeins:
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
# NOTE: model = "~ Fraction + Treatment" We are interested in 'Treatment'.
message(paste("\nEvaluating protein differential abundance between",
	      "WT and Mutant groups."))
comparisons <- "Genotype.Fraction"
samples_to_ignore = "SPQC"
alt_glm_results <- glmDA2(tmt_protein, comparisons, samples, samples_to_ignore)

# Extract key glm results.
alt_results <- alt_glm_results$stats

# Summary of DA proteins with logFC cutoff:
d <- fold_change_delta
sig <- alt_results$FDR < FDR_alpha 
DA <- alt_results$logFC < log2(1-d) | alt_results$logFC > log2(1+d)
message(paste("\nFor WT v Mutant contrast, there are..."))
message(paste("Total number of unique, differentially abundant proteins:",
	      sum(DA & sig)))
message(paste0("...Percent Change +/- ", round(100*fold_change_delta,2),"%."))
message(paste("...FDR <",FDR_alpha))

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

# Add adjusted protein values to tmt_protein.
# Rename Intensity column.
idy <- which(colnames(adjusted_prot) == "Intensity")
colnames(adjusted_prot)[idy] <- "Adjusted.Intensity"

# Add stats to adjusted protein data.
idy <- colnames(alt_results) %notin% c("Accession","Entrez","Symbol")
colnames(alt_results)[idy] <- paste0("Adjusted.",
				     colnames(alt_results)[idy])
adjusted_prot <- left_join(adjusted_prot,alt_results,
			   by=intersect(colnames(adjusted_prot),
					colnames(alt_results)))
adjusted_prot$Adjusted.Intensity <- 2^adjusted_prot$Adjusted.Intensity

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
# * Input - raw peptides
# * Output - normalized protein
# * Statistical results:
#     - Seperate sheet for each fraction/comparison.
#     - Include contrast specific data.

# Add normalized protein data 
norm_dm <- tmt_protein %>% as.data.table() %>%
	dcast(Accession ~ Sample,value.var="Intensity") %>%
	as.matrix(rownames="Accession")

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

# Add adjusted protein values to alt stats.
df <- adjusted_prot %>% as.data.table() %>%
	dcast(Accession ~ Sample, value.var="Adjusted.Intensity")
alt_results <- left_join(alt_results,df,by="Accession")
colnames(alt_results) <- gsub("Abundance","Adjusted.Abundance",
			      colnames(alt_results))

# Combine intra-fraction comparisons and wt v mutant comparisons.
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
final_results <- c(list("Samples" = samples), final_results)
myfile <- file.path(tabsdir,"Swip_TMT_Protein_GLM_Results.xlsx")
write_excel(final_results,myfile,rowNames=FALSE)

rm(list=c("myfile","idy"))

#---------------------------------------------------------------------
## Save output for downstream analysis.
#---------------------------------------------------------------------
# Save key results.

# Save raw data -- tidy_peptide.
myfile <- file.path(rdatdir,"tidy_peptide.csv")
fwrite(tidy_peptide,myfile)

# Save tmt_protein
myfile <- file.path(rdatdir,"tmt_protein.csv")
fwrite(tmt_protein,myfile)

# Save statistical results.
myfile <- file.path(rdatdir,"glm_results.RData")
saveRDS(glm_results,myfile)

# Save to select protein statistics to file.
myfile <- file.path(rdatdir,"Select_Protein_Stats.csv")
fwrite(ttest_dt,myfile)

# Save gene map.
myfile <- file.path(datadir,"gene_map.rda")
save(gene_map,file=myfile,version=2)

# Save tidy_protein as rda object. 
myfile <- file.path(datadir,"tmt_protein.rda")
save(tmt_protein,file=myfile,version=2)
