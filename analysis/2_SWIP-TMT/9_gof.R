#!/usr/bin/env Rscript

# prepare the env
root = "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

# load the data
data(fx0) # protein model
data(fx1) # module model
data(swip)
data(gene_map)
data(partition)
data(msstats_prot)
data(pd_annotation)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(variancePartition)
})


## prepare the data -----------------------------------------------------------

# anno msstats_prot with Module membership
msstats_prot <- msstats_prot %>% filter(Protein %in% names(partition)) %>% 
	mutate(Module = paste0("M",partition[Protein]))

# munge to split Condition (Geno.BioFrac) into geno annotation
samples <- pd_annotation
samples$Genotype <- sapply(strsplit(samples$Condition,"\\."),"[",1)

# munge to create mouse Subject identifier
samples$Subject <- as.numeric(interaction(samples$Genotype,samples$Experiment))

# munge to annotate msstats_prot with Subject (1-6)
msstats_prot$Subject <- as.numeric(interaction(msstats_prot$Genotype,msstats_prot$Mixture))

# cast the data into a matrix.
fx <- formula(Protein ~ Mixture + Genotype + BioFraction)
dm <- msstats_prot %>%
	reshape2::dcast(fx, value.var= "Abundance") %>% 
	as.data.table() %>% na.omit() %>% as.data.table() %>% as.matrix(rownames="Protein")

# munge to create sample info from the dcast fx
info <- as.data.table(do.call(rbind,strsplit(colnames(dm),"_")))
colnames(info) <- strsplit(as.character(fx)[3]," \\+ ")[[1]]

# calculate variance explained by major covariates
form <- formula(~ (1|Mixture) + (1|Genotype) + (1|BioFraction))
prot_varpart <- fitExtractVarPartModel(dm, form, info)


# collect results
df <- as.data.frame(prot_varpart)
varpart_df <- as.data.table(df,keep.rownames="Protein")

# annotate with gene ids
idx <- match(varpart_df$Protein,gene_map$uniprot)
varpart_df$Symbol <- gene_map$symbol[idx]
varpart_df$Entrez <- gene_map$entrez[idx]
# sort cols
varpart_df <- varpart_df %>%
	select(Protein,Symbol,Entrez,Mixture,Genotype,BioFraction,Residuals) %>%
	arrange(desc(Genotype))


## fit all protein models -----------------------------------------------------


# evaluate gof of all protein-wise models
proteins <- names(partition)

# register parallel backend
n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(n_cores)

gof_list <- foreach (protein = proteins) %dopar% {
	# fit model
	fx0 <- Abundance ~ 0 + Condition + (1 | Mixture)
	fm <- suppressMessages({
		try({
		  lmerTest::lmer(fx0, msstats_prot %>% filter(Protein==protein))
		}) 
	})
	if (inherits(fm,"try-error")) { 
             return(NULL)
	} else {
		# evaluate goodness-of-fit
		r2 <- setNames(as.numeric(r.squaredGLMM.merMod(fm)),
				       c("R2.fixef","R2.total"))
		return(r2)
	}
} #EOL
names(gof_list) <- proteins

# collect results
idx <- !sapply(gof_list,is.null)
gof_df <- as.data.table(do.call(rbind,gof_list[idx]),keep.rownames="Protein")

# combine varpart -- the precent variance explained for each covariate and
# gof_df -- the overall R2 for total and fixed effects
results_df <- left_join(varpart_df,gof_df,by="Protein")


## save results -----------------------------------------------------------------

# save as rda
protein_gof <- results_df
myfile <- file.path(root,"data","protein_gof.rda")
save(protein_gof, file=myfile, version=2)

# save the data
fwrite(protein_gof,file.path(root,"rdata","protein_gof.csv"))


## fit all module models -------------------------------------------------------

modules = split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))

moduleGOF <- function(module,msstats_prot){
  form1 <- Abundance ~ (1|Mixture) + (1|Genotype) + (1|BioFraction) + (1|Protein)
  fm1 <- lmer(form1,msstats_prot %>% filter(Module==module))
  vpart <- calcVarPart(fm1)
  form2 <- Abundance ~ 0 + Genotype:BioFraction + (1|Mixture) + (1|Protein)
  fm2 <- lmer(form2,msstats_prot %>% filter(Module==module))
  r2 <- setNames(as.numeric(r.squaredGLMM.merMod(fm2)),
		 nm=c("R2.fixef","R2.total"))
  rho <- c(vpart,r2)
  return(rho)
}

# loop to evaluate gof
results_list <- list()
pbar <- txtProgressBar(max=length(modules),style=3)
for (module in names(modules)) {
	gof <- tryCatch(expr = { moduleGOF(module,msstats_prot) },
			error = function(e) {}, # return null if error or 
			warning = function(w) {}) # warning
	results_list[[module]] <- gof
	setTxtProgressBar(pbar,value=match(module,names(modules)))
}
close(pbar)

# drop null results
idx <- sapply(results_list,is.null)
message("There were problems fitting ", sum(idx), " models.")

# collect results
df <- as.data.table(do.call(rbind,results_list[idx]),keep.rownames="Module")
df <- df %>% arrange(desc(Genotype))

# save the data
fwrite(df,file.path(root,"rdata","module_gof.csv"))

# save as rda
module_gof <- df
myfile <- file.path(root,"data","module_gof.rda")
save(module_gof, file=myfile, version=2)
