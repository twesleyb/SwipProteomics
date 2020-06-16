ggplotPCA <- function(tp){
	# Create group names.
	groups <- paste(sapply(strsplit(tp$Sample,"\\."),"[",7),
			sapply(strsplit(tp$Sample,"\\."),"[",6))
	names(groups) <- tp$Sample
	# Perform PCA.
	dt <- tp %>% dcast(Accession ~ Sample,value.var="Intensity")
	dm <- as.matrix(dt,rownames=TRUE)
	dm[is.na(dm)] <- 1
	pca <- prcomp(t(log2(dm)))[["x"]][,1:2]
	# Generate plot.
	df <- as.data.table(pca)
	rownames(df) <- rownames(pca)
	rownames(df) <- groups[rownames(df)]
	plot <- ggplot(df,aes(x=PC1,y=PC2)) + geom_text(aes(label=rownames(df)))
	return(plot)
}
