#!/usr/bin/env Rscript

#' ---
#' title: Swip Proteomics
#' description: Preprocessing of Swip TMT proteomics data.
#' authors: Tyler W Bradshaw
#' ---

## NOTE: The functions used in this script are not robust. They were written
## to work with the input arguments that they are provided, and in many cases 
## will not perform as expected if passed different arguments.

## Input data in root/input/:
zip_file = "TMT.zip"
input_data = "TMT-raw-peptide.csv"
input_samples = "TMT-samples.csv"

## Analysis Options:
save_plots = TRUE # Should plots be saved in root/figs?
FDR_alpha = 0.1 # FDR threshold for differential abundance.
fold_change_delta = 0.2 # Fold change threshold.
sample_connectivity_threshold = 2.5 # Threshold for sample level outliers.
fig_height = 5.0 # Default height of figures.
fig_width = 5.0 # Default width of figures.

## Input in root/data:
# * TMT.zip/TMT-samples.csv     - sample meta data.
# * TMT.zip/TMT-raw-peptide.csv - raw peptide data from PD.

## Output in root/tables:
# * Swip_TMT_Results.xlsx

## Output for R package in root/data.
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

start <- Sys.time()
message(paste("Starting analysis at:", start))

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
downdir <- file.path(rootdir, "downloads") # Misc trash.
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

# Format cfg force column -- this is the Force in g's used to obtain 
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
results <- glmDA(filt_protein,comparisons,samples,
		 samples_to_ignore="SPQC")

# Collect the results.
fit <- results$fit
qlf <- results$qlf
dge <- results$dge
glm_results <- results$stats
glm_results <- glm_results[c("F4","F5","F6","F7","F8","F9","F10")]

# Annotate final normalized protein with sample meta data.
tidy_protein <- left_join(filt_protein,samples,by=c("Experiment","Sample",
						     "Channel","Treatment"))

# Summary of DA proteins:
message(paste0("\nSummary of differentially abundant proteins ",
	      "in each subceulluar fraction (FDR < ",FDR_alpha,"):"))
knitr::kable(t(sapply(glm_results,function(x) sum(x$FDR < FDR_alpha))))

# Total number of unique DA proteeins:
# FIXME: NUMBER  is not right.
all_sig <- lapply(glm_results, function(x) x$Accession[x$FDR < FDR_alpha])
n_sig <- length(unique(unlist(all_sig)))
message(paste("\nTotal number of unique DA proteins:",n_sig))

# Proteins that are commonly dysregulated:
combined_results <- rbindlist(glm_results,idcol="Fraction")
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
## Alternative statistical model: WT vs Mutant differential abundance.
#----------------------------------------------------------------------

# Use EdgeR glm to evaluate differential protein abundance.
# Compare WT v Mutant within a fraction.
# NOTE: model = "~ Fraction + Treatment" We are interested in WT versus Mutant.
message(paste("\nEvaluating protein differential abundance between",
	      "WT and Mutant groups."))
comparisons <- "Genotype.Fraction"
samples_to_ignore = "SPQC"
alt_glm <- glmDA2(tidy_protein, comparisons, samples, samples_to_ignore)

# Extract key glm results.
alt_glm_results <- alt_glm$stats
fit <- alt_glm$fit
qlf <- alt_glm$qlf
dge <- alt_glm$dge

# Annotate results with gene Symbols.
idx <- match(alt_glm_results$Accession,gene_map$uniprot)
Symbol <- gene_map$symbol[idx]
alt_glm_results <- tibble::add_column(alt_glm_results, 
				      Symbol, .after="Accession")

# Summary of DA proteins with logFC cutoff:
d <- fold_change_delta
sig <- alt_glm_results$FDR < FDR_alpha 
DA <- alt_glm_results$logFC < log2(1-d) | alt_glm_results$logFC > log2(1+d)
message(paste("\nFor WT v Mutant contrast, there are..."))
message(paste("Total number of unique, differentially abundant proteins:",
	      sum(DA & sig)))
message(paste0("...Percent Change +/- ", round(100*fold_change_delta,2),"%."))
message(paste("...FDR <",FDR_alpha))

#----------------------------------------------------------------------
## Plot commonly dysregulated prots, adjusted for fraction differences.
#----------------------------------------------------------------------

# Calculate protein abundance, adjusted for fraction differences.
logCPM <- edgeR::cpm(dge, log=TRUE)

# Remove effect of freaction.
dm <- limma::removeBatchEffect(logCPM,batch=dge$samples$fraction,
				   design=model.matrix(~treatment,
						       data=dge$samples))

# Collect results.
dt <- as.data.table(dm,keep.rownames="Accession") %>% 
	reshape2::melt(id.var="Accession",
		       variable.name="Sample",
		       value.name="Intensity")
dt$Sample <- as.character(dt$Sample)

# Combine with additional meta data.
idy <- which(colnames(tidy_protein)=="Intensity")
temp_prot <- as.data.table(tidy_protein)[,-..idy] %>% filter(Treatment != "SPQC")
adjusted_prot <- left_join(temp_prot,dt, by=c("Accession","Sample"))

# Wash proteins + control.
prots <- c("Washc4","Washc1","Washc2","Washc5","Tubb4a")
names(prots) <- gene_map$uniprot[match(prots,gene_map$symbol)]

# Combine adjusted data and stats, subset to keep proteins of interest.
df <- adjusted_prot %>% filter(Accession %in% names(prots))

# Colors for WT and mutant groups.
colors = c(Control="#47b2a4",Mutant="#B86FAD")

# Labels will simply be WT and Mutant.
xlabels <- rep(c("WT","Mutant"),times=length(prots))

# Order of the factors.
factor_order <- paste(rep(names(prots),each=2), c("Control","Mutant"),sep=".")
df$Accession.Treatment <- as.character(interaction(df$Accession,df$Treatment))
df$Accession.Treatment <- factor(df$Accession.Treatment,levels=factor_order)

# Generate a plot.
plot <- ggplot(df, aes(x=Accession.Treatment, y=Intensity,
		       fill=Treatment)) + 
	geom_boxplot() + geom_point(aes(fill=Treatment,shape=Treatment)) +
	theme(axis.text.x=element_text(angle=45))
plot <- plot + scale_fill_manual(name="Genotype",values=colors)
plot <- plot + theme(legend.position="none")
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(panel.border=element_rect(colour="black",fill="NA",size=1))
plot <- plot + theme(axis.title.x = element_blank())
plot <- plot + scale_x_discrete(labels=xlabels)
plot <- plot + ylab("log2(Adjusted Intensity)")

# Add some lines to break up the data.
plot <- plot + geom_vline(xintercept=seq(2.5,length(proteins)*2,by=2),
			  linetype="dotted",size=0.5)

# Perform t-tests.
prot_list <- df %>% group_by(Accession) %>% group_split()
data_ttests <- lapply(prot_list,function(subdf) {
			      x <- subdf$Intensity
			      y <- subdf$Treatment
			      result <- t.test(x~y,paired=FALSE,
					       alternative="greater")
			      ttest_dt <- as.data.table(t(unlist(result)))
			      return(ttest_dt)
			  })
names(data_ttests) <- names(prots)
ttest_dt <- bind_rows(data_ttests,.id="Accession")
ttest_dt$p.adjust <- p.adjust(ttest_dt$p.value,method="bonferroni")

# Annotate the plot with stats.
stats <- df %>% filter(Treatment=="Control") %>% 
	group_by(Accession.Treatment) %>% 
	dplyr::summarize(ypos = 1.02*max(Intensity),.groups="drop")
stats$Accession.Treatment <- as.character(stats$Accession.Treatment)
stats$Accession <- sapply(strsplit(stats$Accession.Treatment,"\\."),"[",1)
stats <- stats %>% dplyr::select(Accession.Treatment,Accession,ypos)
stats$xpos <- seq(1.5,by=2,length.out=5)
stats$symbol <- "ns"
stats <- left_join(stats,ttest_dt,by="Accession")
stats$symbol[stats$p.adjust < 0.05] <- "*"
stats$symbol[stats$p.adjust < 0.005] <- "**"
stats$symbol[stats$p.adjust < 0.0005] <- "***"

# Add significance stars.
plot <- plot + 
	annotate("text",x=stats$xpos,y=stats$ypos,label=stats$symbol,size=4)

# Annotate with protein names.
symbols <- prots 
build <- ggplot_build(plot)
ymax <- build$layout$panel_params[[1]][["y.range"]][2]
plot <- plot + annotate("text",x=seq(1.5,length(prots)*2,by=2),
			y=ymax,label=symbols,size=5)

# Save.
if (save_plots) {
	myfile <- file.path(figsdir,
			    "Select_Proteins_Adjusted_Abundance.pdf")
	ggsave(myfile,plot, height = fig_height, width = fig_width)
}

rm(list=c("myfile","plot","ymax","symbols","build","stats","prot_list",
	  "data_ttests","xlabels","prots"))

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

# Sort the statistical results.
glm_results <- glm_results[c("F4","F5","F6","F7","F8","F9","F10")]

# Data to be added to stats.
norm_protein <- tidy_protein %>% as.data.table() %>%
	dcast(Accession ~ Sample,value.var="Intensity") %>%
	as.matrix(rownames="Accession")

# Loop to add normalized protein data to glm statistical results.
message("\nSaving TMT data and statistical results.")
for (i in 1:length(glm_results)){
	df <- glm_results[[i]]
	namen <- names(glm_results)[i]
	subsamples <- samples$Sample[grepl(namen,samples$Fraction)]
	idy <- match(subsamples,colnames(norm_protein))
	dm <- norm_protein[,idy]
	# Sort the data by Exp.Sample
	ids <- sapply(sapply(strsplit(colnames(dm),"\\."),"[",c(6,7),
			     simplify=FALSE), paste,collapse=".")
	names(ids) <- colnames(dm)
	ids <- ids[order(ids)]
	dm <- dm[,names(ids)]
	# Coerce dm to dt and merge with stats df.
	dt <- as.data.table(dm,keep.rownames="Accession")
	dt_out <- left_join(df,dt,by="Accession") %>% as.data.table
	glm_results[[i]] <- dt_out
} # Ends loop.

# Map uniprot to entrez IDs.
uniprot <- glm_results[[1]][["Accession"]]
idx <- match(uniprot,gene_map$uniprot)
entrez <- gene_map$entrez[idx]
names(entrez) <- uniprot

# Map missing ids by hand.
is_missing <- is.na(entrez)
missing_entrez <- entrez[is_missing]
mapped_by_hand <- c("P00848"=17705,
		    "P62806"=326619,
		    "Q64525"=319189,
		    "P06330"=111507, 
		    "Q5PR69"=320827, 
		    "Q80WG5"=241296,
		    "P03975"=15598, 
		    "Q80TK0"=100503185)
entrez[names(mapped_by_hand)] <- mapped_by_hand

# Check that we have successfully mapped all uniprot ids to entrez.
check <- sum(is.na(entrez))
if (check != 0) { stop("There are unmapped uniprot IDs.") }

# Map entrez to gene symbols.
idx <- match(entrez,gene_map$entrez)
symbols <- gene_map$symbol[idx]
names(symbols) <- entrez

# Loop to add Entrez IDs and Gene Symbols to data.
# Also fix column titles at the same time. 
for (i in 1:length(glm_results)){
	df <- glm_results[[i]]
	colnames(df) <- gsub("\\."," ",colnames(df))
	idy <- grepl("Abundance",colnames(df))
	colnames(df)[idy] <- paste("Log2",colnames(df)[idy])
	df <- tibble::add_column(df,"Entrez"=entrez[df$Accession],
				 .after="Accession")
	df <- tibble::add_column(df,"Gene"=symbols[as.character(df$Entrez)],
				 .after="Entrez")
	glm_results[[i]] <- df
}

## Save as excel workboook.
## FIXME: Add select protein ttest stats.
names(glm_results) <- paste(names(glm_results),"Results")
results <- c(list("Samples" = samples,
		"Normalized Protein" = tidy_protein),
	     glm_results)
myfile <- file.path(tabsdir,"Swip_TMT_Results.xlsx")
write_excel(results,myfile,rowNames=FALSE)

rm(list=c("myfile","idy","is_missing","entrez"))

#---------------------------------------------------------------------
## Combine final normalized TMT data and stats.
#---------------------------------------------------------------------

# Collect stats.
idx <- names(results)[grep("F[0-9]{1,2}",names(results))]
temp_list <- results[idx]
names(temp_list) <- gsub(" Results","",names(temp_list))
temp_list <- lapply(temp_list,function(x) x[,c(1:8)])
temp_dt <- bind_rows(temp_list,.id="Fraction")

# Merge protein data and stats.
tmt_protein <- left_join(tidy_protein,temp_dt,
			 by=intersect(colnames(tidy_protein),colnames(temp_dt)))

# Remove QC data.
tmt_protein <- tmt_protein %>% filter(Treatment != "SPQC") %>% 
	as.data.table()

# Add adjusted protein values to tmt_protein.
# Rename Intensity column.
idy <- which(colnames(adjusted_prot) == "Intensity")
colnames(adjusted_prot)[idy] <- "Adjusted.Intensity"

# Add stats to adjusted protein data.
colnames(alt_glm_results)[-c(1,2)] <- paste0("Adjusted.",
					     colnames(alt_glm_results)[-c(1,2)])
adjusted_prot <- left_join(adjusted_prot,alt_glm_results,by="Accession")

# tmt_protein with adjusted data and stats.
idy <- intersect(colnames(tmt_protein),colnames(adjusted_prot))
tmt_protein <- left_join(tmt_protein,adjusted_prot,by=idy)

# Select columns  of interest.
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

rm(list=c("idx","temp_list","temp_dt","idy"))

#---------------------------------------------------------------------
## Save output for downstream analysis.
#---------------------------------------------------------------------
# Save key results.

# Save raw data -- tidy_peptide.
myfile <- file.path(rdatdir,"tidy_peptide.csv")
fwrite(tidy_peptide,myfile)

# Save tidy_protein.
myfile <- file.path(rdatdir,"tidy_protein.csv")
fwrite(tidy_protein,myfile)

# Save statistical results.
myfile <- file.path(rdatdir,"glm_results.RData")
saveRDS(glm_results,myfile)

# Save to select protein statistics to file.
myfile <- file.path(rdatdir,"Select_Protein_Stats.csv")
fwrite(ttest_dt,myfile)

## Output saved in root/data:

# Save gene map.
myfile <- file.path(datadir,"gene_map.rda")
save(gene_map,file=myfile,version=2)

# Save tidy_protein as rda object. 
myfile <- file.path(datadir,"tmt_protein.rda")
save(tmt_protein,file=myfile,version=2)

# Done!
end <- Sys.time()
message(paste("\nCompleted analysis at:",end))
