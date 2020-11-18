#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: generate some gof statistics for protein-wise models

## options

# threshold for defining proteins with poor fit
r2_threshold = 0.70


## prepare the env
root = "~/projects/SwipProteomics"
renv::load(root)

# library(SwipProteomics)
devtools::load_all(root)

# load the data
data(swip)
data(gene_map)
data(msstats_prot)
data(pd_annotation)


## imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(doParallel) 
})


## ---- prepare the data for variancePartition 

# cast the data into a matrix
form <- formula(Protein ~ Mixture + Genotype + BioFraction)
dm <- msstats_prot %>%
	reshape2::dcast(form, value.var= "Abundance") %>% 
	as.data.table() %>% as.matrix(rownames="Protein")

# we cannot work with missingness for variancePartition
idx <- apply(dm,1,function(x) any(is.na(x)))
if (sum(idx)>0) {
	warning("Removing ", sum(idx), " proteins with any missing values.")
	dm <- dm[!idx,]
}

stopifnot(!any(is.na(dm)))

# munge to create sample info from the dcast formula
info <- as.data.table(do.call(rbind,strsplit(colnames(dm),"_")))
colnames(info) <- strsplit(as.character(form)[3]," \\+ ")[[1]]


## ---- variancePartition analysis 

# calculate protein-wise variance explained by major covariates in the
# mixed-model:
form0 <- formula(~ (1|Mixture) + (1|Genotype) + (1|BioFraction))

# register parallel backend
n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(n_cores)

prot_varpart <- variancePartition::fitExtractVarPartModel(dm, form0, info)


## parse the response 

# collect results
df <- as.data.frame(prot_varpart)
class(df) <- "data.frame"

varpart_df <- as.data.table(df,keep.rownames="Protein")

# annotate with gene ids
idx <- match(varpart_df$Protein,gene_map$uniprot)
varpart_df$Symbol <- gene_map$symbol[idx]
varpart_df$Entrez <- gene_map$entrez[idx]

# sort cols and rows
varpart_df <- varpart_df %>%
	select(Protein,Symbol,Entrez,Mixture,Genotype,BioFraction,Residuals) %>%
	arrange(desc(Genotype))

# examine the top results
varpart_df %>% head() %>% knitr::kable()


## ---- fit all protein models 
# calculate Nakagawa coefficient of determination

# evaluate gof of all protein-wise models
# use proteins from varpart -- all proteins (-any with missing vals)
proteins <- unique(varpart_df$Protein)

# register parallel backend
n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(n_cores)

# REF: Nakagawa and Schielzeth 2013 and 2017
message("\nEvaluating Nakagawa goodness-of-fit, refitting modules...")

gof_list <- foreach (protein = proteins) %dopar% {
	# fit protein-wise model
	fx0 <- Abundance ~ 1 + Condition + (1 | Mixture)
	fm <- suppressMessages({
		try({
		  lmerTest::lmer(fx0, msstats_prot %>% filter(Protein==protein))
		}) 
	})
	# if error, return null
	if (inherits(fm,"try-error")) { 
             return(NULL)
	} else {
		# else, evaluate goodness-of-fit
		r2 <- setNames(as.numeric(r.squaredGLMM.merMod(fm)),
				       c("R2.fixef","R2.total"))
		# return coefficient of determination
		return(r2)
	}
} #EOL
names(gof_list) <- proteins

# collect results
idx <- !sapply(gof_list,is.null)
if (any(!idx)) { 
	warning("There were problems fitting ",sum(idx)," models.")
}

gof_df <- as.data.table(do.call(rbind,gof_list[idx]),keep.rownames="Protein")

# combine varpart -- the precent variance explained for each covariate and
# gof_df -- the overall R2 for total and fixed effects (which is the percent
# variance explained by our models).
protein_gof <- left_join(varpart_df,gof_df,by="Protein")

# examine results
protein_gof %>% head() %>% knitr::kable()

# highly spatially cohesive proteins
protein_gof %>% arrange(desc(BioFraction)) %>% head() %>% knitr::kable()


## ---- identify proteins with poor fit 

# use overall R2 -- the percent of variance explained by the model as a natural
# description of the overall quality of the fit
message("\nR2 threshold: ", r2_threshold)

# Summary
total <- length(unique(protein_gof$Protein))
out <- sum(protein_gof$R2.total < r2_threshold)
percent <- round(out/total,3)
final <- total-out
cbind(r2_threshold, out, percent, total,final) %>% knitr::kable()


# proteins with R2 less than threshold are poor_prots
poor_prots <- protein_gof$Protein[protein_gof$R2.total < r2_threshold]
message("\nNumber of proteins with poor fit: ", 
	formatC(length(poor_prots),big.mark=","))

# wash prots
message("\nWASHC* protein goodness-of-fit statistics:")
protein_gof %>% filter(Protein %in% mapID("Washc*")) %>% knitr::kable()


## ---- save results

# save character vector of poor_prots
myfile <- file.path(root,"data","poor_prots.rda")
save(poor_prots,file=myfile,version=2)
  
# save protein_gof data.table as rda --> used to annotate plots
myfile <- file.path(root,"data","protein_gof.rda")
save(protein_gof, file=myfile, version=2)
  
# TODO: save as excel table!
myfile <- file.path(root,"tables","Protein_GOF.csv")
fwrite(protein_gof,file=myfile)
