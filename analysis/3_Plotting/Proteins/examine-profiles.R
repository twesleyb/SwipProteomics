#!/usr/bin/env Rscript

root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all()

data(swip)
data(partition)
data(gene_map)
data(msstats_prot)
data(msstats_results)

myfile <- file.path(root,"rdata","adjm.rda")
load(myfile) # adjm

myfile <- file.path(root,"rdata","ne_adjm.rda")
load(myfile) # ne_adjm

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(data.table)
})


getMedoid <- function(adjm, k = NULL, h = NULL, method = "complete") {
  # Find representative branches among groups in a heirarchy
  # The medoid of the group is the branch that is closest to
  # all branches in its group.
  hc <- hclust(as.dist(1 - adjm), method)
  hc_part <- cutree(hc, k, h)
  groups <- split(names(hc_part), hc_part)
  medoids <- sapply(seq(groups), function(x) {
	       idy = idx = groups[[x]]
	       col_sums <- apply(adjm[idx,idy], 2, sum)
	       medoid <- names(which(col_sums == min(col_sums)))
  })
  return(medoids)
}


prots <- getMedoid(adjm,k=5)

plots <- lapply(prots, plot_profile)
names(plots) <- prots

subadjm = adjm[prots,prots]
subadjm[lower.tri(subadjm,diag=TRUE)] <- NA
df = reshape2::melt(subadjm,na.rm=TRUE)
colnames(df) <- c("ProtA","ProtB","weight")
df <- df %>% arrange(desc(weight))

head(df)

protA <- names(plots)[1]
protB <- names(plots)[2]

corProt <- function(protA,protB,method="pearson") {
	# requires msstats_prot
	#library(dplyr)
	#library(data.table)
	df <- msstats_prot %>% filter(Protein == protA | Protein == protB)
	dm <- df %>% reshape2::dcast(Protein ~ Mixture + Channel, 
				     value.var = "Abundance") %>%
	        as.data.table() %>%
	        as.matrix(rownames="Protein")
	if (method %in% c("pearson", "kendall", "spearman")) {
		rho <- cor(t(dm), method = method)
	} else if (method == "bicor") {
		rho <- WGCNA::bicor(t(dm))
	} else {
		stop("method should be one of: ",
		    "c('pearson', 'kendall', 'spearman', 'bicor').")
	}
	return(rho)
}

corProt("E9Q8I9","P41216",method = "spearman")

corProt("E9Q8I9","P41216",method = "kendall")

corProt("E9Q8I9","P41216",method = "pearson")

corProt("E9Q8I9","P41216",method = "bicor")

plots[["E9Q8I9"]]

plots[["P41216"]] # these two are very different, but bicor is .2

plots[[1]]
plots[[2]] # these two are more similar

protA <- names(plots)[1]
protB <- names(plots)[2]

corProt(protA,protB,method = "spearman") # .60

corProt(protA,protB,method = "kendall") # .46

corProt(protA,protB,method = "pearson") # 0.92

corProt(protA,protB,method="bicor") # here bicor is -.12

# wow, this is pretty profound -- need to illustrate this example
