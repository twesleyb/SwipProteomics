#' module_analysis
#' @export
module_analysis <- function(data,partition,groups,nThreads){

	# This is a nasty bit of code that does the statistical comparisons of
	# module summary abundance with Kruskal-Wallis and Dunnett's 
	# tests. This is done in several steps:
	# 1. Calculate module summary abundance and other summary stats.
	# 2. Define a function to do the KW + Dunnett's tests.
	# 3. For every module, do 2.
	# 4. Merge the statisical results from 3.
	# 5. Combine with module stats from 1.
	# Imports.
	suppressPackageStartupMessages({
		library(dplyr)
		library(WGCNA)
		library(parallel)
		library(doParallel)
		library(data.table)
	})

	## 1. Calculate Module summary abundance and other summary stats.
	# Calculate module summary abundance (module eigenproteins).
	data_ME <- moduleEigengenes(data,colors=partition,
				    excludeGrey=TRUE,softPower=1,impute=FALSE)
	# Split partition into modules.
	modules <- split(partition,partition)
	names(modules) <- paste0("M",names(modules))
	k <- sum(names(modules) != "M0")
	n <- sapply(modules,length)
	# Extract ME data as a matrix.
	dm_ME <- as.matrix(data_ME$eigengenes)

	# Tidy-up the ME data.
	dt_ME <- melt(as.data.table(dm_ME,keep.rownames="Sample"),
		      id.vars="Sample",variable.name="Module",
	              value.name="ME")
	dt_ME$Module <- gsub("ME","",dt_ME$Module)
	# Add fraction label.
	Fraction <- sapply(strsplit(groups[dt_ME$Sample],"\\."),"[",1)
	dt_ME <- tibble::add_column(dt_ME,Fraction,.after="Module")
	# Extract module percent variance explained (PVE; ~quality).
	pve <- unlist(data_ME$varExplained)
	names(pve) <- gsub("X","M",names(pve))
	# Add PVE to data.table.
	dt_ME$PVE <- pve[paste0("M",dt_ME$Module)]
	# Add module size (n) to data.table.
	dt_ME$nProteins <- n[paste0("M",dt_ME$Module)]
	# Add total number of modules (k) to data.table.
	dt_ME$kModules <- k
	# Split ME matrix into list of column vectors.
	list_ME <- lapply(seq(ncol(dm_ME)),function(x) dm_ME[,x])
	names(list_ME) <- colnames(dm_ME)

	## 2. Define a function that does the statistical comparisons.
	KWDtest <- function(module){
		suppressPackageStartupMessages({
			library(dplyr)
			library(data.table)
		})
		x <- list_ME[[module]]
		g <- factor(groups[names(x)])
		# Perform KW test.
		kw_results <- kruskal.test(x ~ g)
		# Clean-up KW results.
		dt_KW <- as.data.table(t(do.call(rbind,kw_results)))
		dt_KW$statistic <- as.numeric(dt_KW$statistic)
		dt_KW$parameter <- as.numeric(dt_KW$parameter)
		dt_KW$p.value <- as.numeric(dt_KW$p.value)
		colnames(dt_KW)[3] <- "pval"
		dt_KW$method <- NULL
		dt_KW$data.name <- NULL
		colnames(dt_KW) <- paste("KW",colnames(dt_KW),sep=".")
		dt_KW <- tibble::add_column(dt_KW,Module=module,.before=1)
		# Perform Dunnett test with comparison to WT (control) fractions.
		all_samples <- unique(as.character(g))
		controls <- all_samples[grepl("WT",all_samples)]
		dtest_results <- DescTools::DunnettTest(x ~ g, control=controls)
		dtest_results$data.name <- NULL # Remove unncessary object.
		# Loop to extract comparisons between identical fractions from 
		# dtest_results. We will focus on intra-fraction comparisions 
		# e.g., F4.MUT versus F4.WT.
		dtest_list <- list()
		for (control in names(dtest_results)){
			df <- dtest_results[[control]]
			contrasts <- strsplit(gsub("\\.MUT|\\.WT","",
						     rownames(df)),"-")
			idx <- sapply(contrasts,function(x) {
					      do.call(identical,as.list(x)) })
			dtest_list[[control]] <- as.data.frame(subset(df,idx))
		}
		# Bind together the results of the loop.
		dt_DT <- do.call(rbind,dtest_list)
		# Clean-up the Dunnett's test data.table.
		colnames(dt_DT) <- paste("DT",colnames(dt_DT),sep=".")
		dt_DT$Module <- module
		dt_DT$Fraction <- gsub("\\.WT","",rownames(dt_DT))
		# Combine KW and DT results.
		KWDT_results <- left_join(dt_KW,dt_DT,by="Module")
		# Clean-up combined results.
		Fraction <- KWDT_results$Fraction
		KWDT_results$Fraction <- NULL
		KWDT_results <- tibble::add_column(KWDT_results,
						   Fraction,.after="Module")
		# Return KW + DT results.
		return(KWDT_results)
	}

	## 3. For every module, perform Kruskal-Wallis and Dunnett's tests
	## to evaluate changes in module summary abundance.
	# Utilize foreach and dopar to execute loop in parallel.
	if (nThreads > 1) {
		# Parallel loop.
		message(paste("Analyzing changes in module summary abundance",
			      "utilizing", nThreads,"parallel processors."))
		nodes <- makeCluster(nThreads)
		registerDoParallel(nodes)
		results <- foreach(i=seq(1,length(list_ME))) %dopar% { 
			KWDtest(i) }
		names(results) <- paste0("M",c(1:length(results)))
		suppressWarnings({ stopCluster(nodes) })
	} else {
		# Serial loop.
		# Add pbar!
		results <- list()
		for(i in seq(1,length(list_ME))) { 
			results[[i]] <- KWDtest(i) }
	}

	## 4. Merge the statistical results for all modules.
	dt_KWDT <- do.call(rbind,results)

	# 5. Merge module stats with statistical results.
	dt_KWDT$Module <- as.character(dt_KWDT$Module)
	results_combined <- left_join(dt_ME,dt_KWDT,by=c("Module","Fraction"))

	# 6. Finally, correct pvalues for k module comparisons. 
	# Bonferroni KW pvalue adjustment.
	KW.padj <- k * results_combined$KW.pval
	KW.padj[KW.padj > 1.0] <- 1
	results_combined <- tibble::add_column(results_combined,
					       KW.padj,.after="KW.pval")

	# Return the combined results.
	results_combined <- as.data.table(results_combined)
	return(results_combined)

} # EOF.
