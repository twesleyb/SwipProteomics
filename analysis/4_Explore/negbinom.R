#!/usr/bin/env Rscript

# title:
# author: tyler w a bradshaw

getrd <- function(here=getwd(), dpat= ".git") {
	# Get the repository's root directory.
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

root <- getrd()
renv::load(root)

suppressPackageStartupMessages({
	library(data.table)
	library(dplyr)
	library(edgeR)
	library(tibble)
})

data(swip_tmt) # tmt_protein
	
# Cast tidy protein into data matrix for EdgeR.
dm <- tmt_protein %>%  filter(Treatment != "SQPQC") %>%
	dcast(Accession ~ Sample, value.var="Intensity") %>% 
	as.matrix(rownames=TRUE)

# create dge object
dge <- DGEList(counts=dm)

# perform tmm normalization
dge <- calcNormFactors(dge)

# update edgeR samples with covariates
idx <- match(rownames(dge$samples),tmt_protein$Sample)
dge$samples$Treatment <- tmt_protein$Treatment[idx]
dge$samples$Fraction <- tmt_protein$Fraction[idx]
dge$samples$group <- interaction(dge$samples$Treatment,dge$samples$Fraction)

# create a design matrix for glm
design <- model.matrix(~ 0 + group, data = dge$samples)
colnames(design) <- levels(dge$samples$group)

# estimate dispersion
dge <- estimateDisp(dge, design, robust = TRUE)

# plot:
myfile <- file.path(root,"figs","Misc","plotBCV.pdf")
pdf(file=myfile,onefile=TRUE)
plotBCV(dge)
dev.off()

# fit a general linear model
fit <- glmQLFit(dge, design, robust = TRUE)

# plot
myfile <- file.path(root,"figs","Misc","plotQLDisp.pdf")
pdf(file=myfile,onefile=TRUE)
plotQLDisp(fit) # doesnt look like manual
dev.off()

# generate contrasts list
g <- as.character(groups)
idx <- sapply(strsplit(g,"\\."),"[",2)
contrast_list <- split(g,idx)
contrast_list <- lapply(contrast_list,unique)
contrast_list <- lapply(contrast_list,function(x) x[order(x)])
contrast_list <- lapply(contrast_list,function(x) { paste(x,collapse="-") })

# Loop to generate contrasts for edgeR::glm
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


# plot
for (i in c(1:length(pvals))) {
	namen <- names(pvals)[i]
	myfile <- file.path(root,"figs","Misc",paste0(namen,"_pval_hist.pdf"))
	pdf(file=myfile,onefile=TRUE)
	hist(pvals[[i]],breaks=20)
	dev.off()
}
