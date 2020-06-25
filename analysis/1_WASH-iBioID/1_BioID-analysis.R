#!/usr/bin/env Rscript

# Analysing WASH (Washc1) iBioID proteomics.

## User parameters to change:
FDR_alpha = 0.1 # FDR significance threshold for protein enrichment.
enrichment_threshold = log2(3.0) # enrichment threshold.

## Input data in root/data
zipfile = "BioID.zip"
datafile = "BioID_raw_protein.csv" # In BioID.zip/

## Output in root/tables:
# * WASH_BioID_Results.xlsx

## Output in root/data:
# * wash_interactome.rda # The WASH iBioID proteome.

#-------------------------------------------------------------------------------
## Prepare the workspace.
#-------------------------------------------------------------------------------

# Activate renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Imports
suppressPackageStartupMessages({
	library(dplyr) # For manipulating the data.
	library(edgeR) # For statitical comparisons.
	library(geneLists) # For a list of mito proteins.
	library(getPPIs) # For mapping gene identifiers.
	library(data.table) # For working with data.tables.
})

# Load any additional functions in root/R.
suppressMessages({ devtools::load_all() })

# Project directories.
datadir <- file.path(root,"data")
rdatdir <- file.path(root,"rdata")
tabsdir <- file.path(root,"tables")
downdir <- file.path(root,"downloads")

# Create rdata dir if it doesn't exist.
if (!dir.exists(rdatdir)){ dir.create(rdatdir) }
if (!dir.exists(tabsdir)){ dir.create(tabsdir) }
if (!dir.exists(downdir)){ dir.create(downdir) }

#-------------------------------------------------------------------------------
## Load the raw data.
#-------------------------------------------------------------------------------

# Extract the raw data from zipped file.
myfile <- file.path(datadir,zipfile)
unzip(myfile) # unzip 

# Read into R.
myfile <- file.path(getwd(),tools::file_path_sans_ext(zipfile),datafile)
raw_prot <- fread(myfile)

# Clean-up.
myfile <- file.path(downdir,tools::file_path_sans_ext(zipfile))
unlink(myfile,recursive=TRUE)
unlink("./BioID", recursive=TRUE)

# Tidy-up the data.
message("\nLoading raw Swip BioID protein data.")
tidy_prot <- tidyProt(raw_prot,species="Mus musculus",
		      id.vars=c("Accession","Description","Peptides"))

# Load mitochondrial protein list from twesleyb/geneLists.
data(list=geneLists("mito"))
mito_entrez <- unlist(mitocarta2,use.names=FALSE)
mito_prot <- getIDs(mito_entrez,from="entrez",to="uniprot",species="mouse")

# Status.
nMito <- sum(mito_prot %in% tidy_prot$Accession)
message(paste(nMito,"mitochondrial proteins will be removed as contaminants."))

# Remove mitochondrial proteins as contaminants.
tidy_prot <- tidy_prot %>%  filter(Accession %notin% mito_prot)

# Summary:
nProt <- length(unique(tidy_prot$Accession))
message(paste0("\nTotal number of proteins quantified: ",
	      formatC(nProt,big.mark=","),"."))

#-------------------------------------------------------------------------------
## Sample loading normalization. 
#-------------------------------------------------------------------------------

# Sample loading normalization:
message("\nPerforming sample loading normalization.")
SL_prot <- normSL(tidy_prot,groupBy="Sample")

# Check, column sums should now be equal.
message("Total intensity sums are equal after sample loading normalization:")
df <- SL_prot %>% group_by(Sample) %>% 
	summarize("Total Intensity"=sum(Intensity,na.rm=TRUE))
knitr::kable(df)

#-------------------------------------------------------------------------------
## Sample pool normalization.
#-------------------------------------------------------------------------------

# Perform normalization to QC samples.
message("\nPerforming sample pool normalization used pooled QC samples.")
SPN_prot <- normSP(SL_prot,pool="QC")

#-------------------------------------------------------------------------------
## Protein level filtering.
#-------------------------------------------------------------------------------

# Status.
message("\nFiltering proteins...")

# Remove one hit wonders.
one_hit_wonders <- unique(SPN_prot$Accession[SL_prot$Peptides == 1])
n_ohw <- length(one_hit_wonders)
filt_prot <- SPN_prot %>% filter(Peptides>1)

# Status.
message(paste("... Number of one-hit-wonders:",n_ohw))

# Remove proteins with unreliable QC (more than 1 missing values).
df <- filt_prot %>% filter(grepl("QC",Sample)) %>% group_by(Accession) %>%
	summarize(N=length(Intensity),n_missing=sum(is.na(Intensity)))
out <- unique(df$Accession[df$n_missing>0])
filt_prot <- filt_prot %>% filter(Accession %notin% out)

# Status.
n_out <- length(out)
message(paste("... Number of proteins with missing QC data:",n_out))

# Remove proteins with more than 50% missingness as these cannot be imputed.
df <- filt_prot %>% filter(!grepl("QC",Sample)) %>% group_by(Accession) %>%
	summarize(N=length(Intensity),n_missing=sum(is.na(Intensity)))
out <- unique(df$Accession[df$n_missing>0.5*df$N])
filt_prot <- filt_prot %>% filter(Accession %notin% out)

# Status.
n_out <- length(out)
message(paste("... Number of proteins with too many missing values:",n_out))

# Insure that we are still working with a data.table.
filt_prot <- as.data.table(filt_prot)

# Status.
prots <- unique(filt_prot$Accession)
n_prot <- length(prots)
message(paste0("\nFinal number of quantifiable proteins: ",
	      formatC(n_prot,big.mark=","),"."))

#-------------------------------------------------------------------------------
## Impute missing values. 
#-------------------------------------------------------------------------------

# Proteins with missing values are less abundant than those without. 
# This is evidence that missing values are MNAR and can be imputed with
# the KNN algorith.
message("\nImputing missing protein values using the KNN algorithm (k=10).")
imp_prot <- imputeKNNprot(filt_prot,quiet=FALSE)

#-------------------------------------------------------------------------------
## Statistical testing -- EdgeR Exact Test. 
#-------------------------------------------------------------------------------

# Status.
message("\nEvaluating statistical enrichment with EdgeR's Exact Test.")

# Cast the data into a matrix to be passed to edgeR.
dm <- imp_prot %>% dcast(Accession ~ Sample,value.var="Intensity") %>%
	as.matrix(rownames=TRUE)

# Remove QC data.
qc_cols <- grepl("QC",colnames(dm))
dm_filt <- dm[,!qc_cols]

# Create DGEList object with mapping to genotype (group).
genotype <- factor(c(rep("WASH",3),rep("Control",3)))
dge <- DGEList(counts = dm_filt, group = genotype)

# Create design matrix.  
design <- model.matrix(~0+genotype, data=dge$samples)
colnames(design)[c(1,2)] <- levels(dge$samples$group)

# Perform TMM normalization.
dge <- calcNormFactors(dge,method="TMM")

# Estimate dispersion.
dge <- estimateDisp(dge,design,robust=TRUE)

# Extract normalized data.
norm_prot <- as.data.table(log2(dge$counts),keep.rownames="Accession")

# Add QC data back.
dt_qc <- as.data.table(log2(dm[,qc_cols]),keep.rownames="Accession")
norm_prot <- left_join(norm_prot, dt_qc,by="Accession")

# Perform exactTest.
data_ET <- edgeR::exactTest(dge, pair = c("Control", "WASH"))

# Call topTags to add FDR. Keep the data the same order by sort.by="none".
data_TT <- topTags(data_ET, n = Inf, sort.by = "none")

# Extract the statistical results from the topTags object. 
stats <- as.data.table(data_TT$table,keep.rownames="Accession")

# Merge with count normalized protein data.
results <- left_join(stats,norm_prot,by="Accession")

# Categorize candidates by enrichment and FDR.
up <- results$logFC > enrichment_threshold
results$candidate <- "no"
results$candidate[which(up & results$FDR <= 0.10 & results$FDR > 0.05)] <- "low" 
results$candidate[which(up & results$FDR <= 0.05 & results$FDR > 0.01)] <- "med"
results$candidate[which(up & results$FDR <= 0.01)] <- "high"
results$candidate <- factor(results$candidate, 
			    levels = c("high", "med", "low", "no"))

# Sort by P-value and logFC.
results <- results %>% arrange(PValue,logFC)
results <- results %>% arrange(candidate)

# Convert logCPM column to percent control.
results$logCPM <- 100*(2^results$logFC)
idy <- which(colnames(results)=="logCPM")
colnames(results)[idy] <- "Percent Control (%)"

# Summary:
sig <- results$FDR < FDR_alpha
up <- results$logFC > enrichment_threshold
nsig <- sum(sig & up)
message(paste0("\nNumber of significantly enriched proteins ",
	      "(log2FC > ",round(enrichment_threshold,2),
	      "; FDR < ",FDR_alpha,"): "), nsig,".")

# Map uniprot ids to entrez ids using mgi batch query.
# NOTE: this takes a little time as the function downloads the data from MGI.
message("\nMapping Uniprot IDs to stable Entrez IDs and gene symbols.")
entrez <- mgi_batch_query(results$Accession)

# Map any missing ids by hand.
is_missing <- is.na(entrez)
n_missing <- sum(is_missing)
message(paste("Mapping",n_missing,"missing gene identifiers by hand."))
mapped_by_hand <- c("P10853"=319180,"P62806"=326619,"P02301"=625328)
entrez[names(mapped_by_hand)] <- mapped_by_hand
check <- all(!is.na(entrez))
if (!check) { stop("Unable to map all uniprot to stable entrez ids.") }

# Map entrez to gene symbols.
symbols <- getIDs(entrez,from="entrez",to="symbol",species="mouse")
check <- all(!is.na(symbols))
if (!check) { stop("Unable to map all entrez ids to gene symbols.") }

# Add identifiers to data table.
results <- tibble::add_column(results,"Entrez"=entrez,.after="Accession")
results <- tibble::add_column(results,"Gene"=symbols,.after="Entrez")

# Create list of results, 
results_list <- list("Raw Protein" = tidy_prot, "BioID Results" = results)

# Add the mitochondrial proteins that were removed.
df <- raw_prot %>% filter(Accession %in% mito_prot) %>% 
	dplyr::select(Accession)
results_list[["Mitochondrial Contaiminants"]] <- df

# Write to file.
message("\nSaving results.")
myfile <- file.path(tabsdir,"WASH_BioID_Results.xlsx")
write_excel(results_list,myfile)

# Save results as rdata for downstream analysis.
myfile <- file.path(rdatdir,"WASH_BioID_Results.RData")
saveRDS(results,myfile)

# Save final results for R package in root/data.
wash_interactome <- results %>% filter(candidate != "no")
myfile <- file.path(datadir, "wash_interactome.rda")
save(wash_interactome,file=myfile,version=2)
