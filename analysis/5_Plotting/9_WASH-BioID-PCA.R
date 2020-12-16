#!/usr/bin/env Rscript

#renv
devtools::load_all()

library(SwipProteomics)

data(wash_bioid)

data(bioid_gene_map)


dm <- wash_bioid %>% 
	reshape2::dcast(Accession ~ Sample, value.var = "Intensity") %>%
	as.data.table() %>% 
	as.matrix(rownames="Accession")


# PCA
pca <- prcomp(t(subdm[!is_neg,]))
pca_summary <- as.data.frame(t(summary(pca)$importance))
idx <- order(pca_summary[["Proportion of Variance"]],decreasing=TRUE)
pca_summary <- pca_summary[idx,]
top2_pc <- head(pca_summary[["Proportion of Variance"]],2)
names(top2_pc) <- head(rownames(pca_summary),2)

# Plot axis labels:
x_label <- paste0(names(top2_pc)[1],
		  " (PVE: ",round(100*top2_pc[1],2)," %)")
y_label <- paste0(names(top2_pc)[2],
		  " (PVE: ",round(100*top2_pc[2],2)," %)")

# Collect data for plotting.
df <- as.data.frame(pca$x[,names(top2_pc)])
colnames(df) <- c("x","y")

# annotate with group info
geno <- sapply(strsplit(rownames(df),"\\."),"[",2) # <- 
frac <- sapply(strsplit(rownames(df),"\\."),"[",3) # <-
df$group <- interaction(geno,frac)


## ---- generate the plot

plot <- ggplot(df, aes(x,y,color=group)) + geom_point(size=4)
plot <- plot + xlab(x_label)
plot <- plot + ylab(y_label)
plot <- plot + theme(axis.title.x = element_text(color = "black")) 
plot <- plot + theme(axis.title.x = element_text(size = 11))
plot <- plot + theme(axis.title.x = element_text(face = "bold"))
plot <- plot + theme(axis.title.y = element_text(color = "black")) 
plot <- plot + theme(axis.title.y = element_text(size = 11))
plot <- plot + theme(axis.title.y = element_text(face = "bold"))
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(panel.border = element_rect(fill=NA))
plot <- plot + scale_x_continuous(expand = c(0,0))
plot <- plot + scale_y_continuous(expand = c(0,0))


## ---- Save to file

myfile <- file.path(figsdir,"sample_PCA.pdf")
ggsave(myfile, plot, width=fig_width,height=fig_height)
message("saved: ", myfile)

##############################################################################


mungeDEP <- function(data_prot, gene_map) {

	#data_prot <- wash_bioid

	stopifnot(!any(data_prot$Abundance < 0))
	stopifnot(!any(is.na(data_prot$Abundance)))

	# cast to matrix, annotate with ID column (uniprot ids)
	df <- data_prot %>% 
		reshape2::dcast(Accession ~ Sample, value.var = "Intensity") %>%
		as.data.table() %>% 
		as.matrix(rownames="Accession") %>% 
		as.data.table(keep.rownames="ID") %>%
		mutate(name = gene_map$symbol[match(ID, gene_map$uniprot)])

	# insure name (gene symbols) are unique
	df$name <- make.unique(df$name, sep="-")

	## create experimental design matrix

	# samples are colnames
	samples <- colnames(df)[grepl("WASH|Control|QC",colnames(df))]

	# extract exp condition/group from sample name
	groups <- sapply(strsplit(samples,"\\ "),"[",2)

	# exp design df for DEP
	exp_design <- data.table(label = samples,
				 condition = groups,
				 replicate = make.unique(groups))


se <- make_se(df, columns = samples,  exp_design)


# save as pdf
myfile <- file.path(root,"figs","Samples","WASH_BioID-PCA.pdf")
pdf(file=myfile, onefile=TRUE)

DEP::plot_pca(se)

dev.off()

