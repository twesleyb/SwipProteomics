glmDA <- function(tp,comparisons,samples,samples_to_ignore){

	# glm Differential Abundance.
	# Imports.
	suppressPackageStartupMessages({
		library(data.table)
		library(dplyr)
		library(edgeR)
		library(tibble)
	})

	# Cast tp into data matrix for EdgeR.
	tp_in <- tp
	dm <- tp %>%
		dcast(Accession ~ Sample, value.var="Intensity") %>%
		as.matrix(rownames=TRUE)

	# Create dge object.
	dge <- DGEList(counts=dm)

	# Perform TMM normalization.
	dge <- calcNormFactors(dge)

	# Drop samples to ignore from dge object.
	# We dont want to utilize QC when estimating dispersion or
	# performing statistical testing.
	drop <- grepl(samples_to_ignore,rownames(dge$samples))
	dge$samples <- dge$samples[!drop,]
	drop <- grepl(samples_to_ignore,colnames(dge$counts))
	dge$counts <- dge$counts[,!drop]

	# Create sample groupings given contrasts of interest.
	subsamples <- samples %>% filter(Treatment != samples_to_ignore)
	contrasts <- unlist(strsplit(comparisons,"\\."))
	subsamples$Group <- apply(subsamples[,contrasts],1,paste,collapse=".")
	groups <- subsamples$Group
	names(groups) <- subsamples$Sample
	idx <- rownames(dge$samples)
	dge$samples$group <- as.factor(groups[idx])
	dge$samples$fraction <- sapply(strsplit(groups,"\\."),"[",2)[idx]
	dge$samples$treatment <- sapply(strsplit(groups,"\\."),"[",1)[idx]

	# Create a design matrix for GLM.
	design <- model.matrix(~ 0 + group, data = dge$samples)
	colnames(design) <- levels(dge$samples$group)

	# Estimate dispersion.
	dge <- estimateDisp(dge, design, robust = TRUE)
	#save(dge,file="dge.rda",version=2)

	# Fit a general linear model.
	fit <- glmQLFit(dge, design, robust = TRUE)

	# Generate contrasts list.
	g <- as.character(groups)
	idx <- sapply(strsplit(g,"\\."),"[",2)
	contrast_list <- split(g,idx)
	contrast_list <- lapply(contrast_list,unique)
	contrast_list <- lapply(contrast_list,function(x) x[order(x)])
	contrast_list <- lapply(contrast_list,function(x) {
					 paste(x,collapse="-") })

	# Loop to generate contrasts for edgeR::glm.
	for (i in 1:length(contrast_list)) {
		contrast <- contrast_list[[i]]
		cmd <- paste0("makeContrasts(",contrast,", levels = design)")
		contrast_list[[i]] <- eval(parse(text=cmd))
	}

	# Call glmQLFTest() to evaluate differences in contrasts.
	contrasts <- colnames(design)
	qlf <- lapply(contrast_list,function(x) glmQLFTest(fit, contrast = x))

	# Determine number of significant results with decideTests().
	summary_tables <- lapply(qlf, function(x) summary(decideTests(x)))

	# Call topTags to add FDR. Gather tabularized results.
	glm_results <- lapply(qlf, function(x) {
				      topTags(x, n = Inf, sort.by = "none")$table
			  })

	# Insure first column is Accession.
	glm_results <- lapply(glm_results,function(x) {
				      Accession <- rownames(x)
				      x <- add_column(x,Accession,.before=1)
				      rownames(x) <- NULL
				      return(x)
			  })

	# Add percent WT and sort by pvalue.
	glm_results <- lapply(glm_results,function(x) {
				      x$logCPM <- 2^x$logFC
				      idy <- grep("logCPM",colnames(x))
				      colnames(x)[idy] <- "PercentWT"
				      x <- x[order(x$PValue,decreasing=FALSE),]
				      return(x)
			  })

	# Return list of normalized data and results.
	return(list("stats"=glm_results,"dge"=dge,"qlf"=qlf,"fit"=fit))
}
