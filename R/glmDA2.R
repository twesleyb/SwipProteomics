glmDA2 <- function(tp, comparisons, samples, samples_to_ignore) {

  # glm Differential Abundance.
  # Imports.
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(edgeR)
    library(tibble)
  })

  # Cast tp into data matrix for EdgeR.
  tp_in <- tp <- as.data.table(tp)
  dm <- tp %>%
    dcast(Accession ~ Sample, value.var = "Intensity") %>%
    as.matrix(rownames = TRUE)

  # Create dge object.
  dge <- DGEList(counts = dm)

  # Perform TMM normalization.
  dge <- calcNormFactors(dge)

  # Drop samples to ignore from dge object.
  # We dont want to utilize QC when estimating dispersion or
  # performing statistical testing.
  drop <- grepl(samples_to_ignore, rownames(dge$samples))
  dge$samples <- dge$samples[!drop, ]
  drop <- grepl(samples_to_ignore, colnames(dge$counts))
  dge$counts <- dge$counts[, !drop]

  # Create sample groupings given contrasts of interest.
  subsamples <- samples %>% filter(Treatment != samples_to_ignore)
  contrasts <- unlist(strsplit(comparisons, "\\."))
  subsamples$Group <- apply(subsamples[, contrasts], 1, paste, collapse = ".")
  groups <- subsamples$Group
  names(groups) <- subsamples$Sample
  idx <- rownames(dge$samples)
  dge$samples$group <- as.factor(groups[idx])
  dge$samples$fraction <- sapply(strsplit(groups, "\\."), "[", 2)[idx]
  dge$samples$treatment <- sapply(strsplit(groups, "\\."), "[", 1)[idx]
  dge$samples$treatment <- as.factor(dge$samples$treatment)
  levels(dge$samples$treatment) <- c("WT", "MUT")

  # Create a design matrix for GLM.
  design <- model.matrix(~ fraction + treatment, data = dge$samples)
  design[, "treatmentMUT"] <- abs(design[, "treatmentMUT"] - 1) # Flip sign

  # Estimate dispersion.
  dge <- estimateDisp(dge, design, robust = TRUE)

  # Fit a general linear model.
  fit <- glmQLFit(dge, design, robust = TRUE)

  # Call glmQLFTest() to evaluate differences in contrasts.
  qlf <- glmQLFTest(fit)

  # Determine number of significant results with decideTests().
  summary_table <- summary(decideTests(qlf))

  # Call topTags to add FDR. Gather tabularized results.
  glm_results <- topTags(qlf, n = Inf, sort.by = "none")$table

  # Insure first column is Accession.
  Accession <- rownames(glm_results)
  glm_results <- add_column(glm_results, Accession, .before = 1)
  rownames(glm_results) <- NULL

  # Add percent WT and sort by pvalue.
  glm_results$logCPM <- 2^glm_results$logFC
  idy <- grep("logCPM", colnames(glm_results))
  colnames(glm_results)[idy] <- "PercentWT"
  glm_results <- glm_results[order(glm_results$PValue, decreasing = FALSE), ]

  # Return list of normalized data and results.
  return(list("stats" = glm_results, "fit" = fit, "dge" = dge, "qlf" = qlf))
}
