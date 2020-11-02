#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: generate some gof statistics for protein and module-level models

# prepare the env
root = "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

# load the data
data(fx0) # protein model
data(fx1) # module model
data(swip)
data(gene_map)
data(sigprots)
data(partition)
data(msstats_prot)
data(pd_annotation)


# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(doParallel)
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

## variancePartition -----------------------------------------------------------

# protein stuff

# calculate protein-wise variance explained by major covariates
form <- formula(~ (1|Mixture) + (1|Genotype) + (1|BioFraction))
prot_varpart <- variancePartition::fitExtractVarPartModel(dm, form, info)


## parse the response ----------------------------------------------------------

# protein stuff

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


# examine the top results
varpart_df %>% head() %>% knitr::kable()
# A1bg is interesting -- overall sig effect, but R2total is modest at 0.7
# this protein would have been flagged for poor fit and discarded before
# clustering


## fit all protein models -----------------------------------------------------

# protein-level stuff

# evaluate gof of all protein-wise models
proteins <- names(partition)

# register parallel backend
n_cores <- parallel::detectCores() - 1
doParallel::registerDoParallel(n_cores)

# FIXME: double check name, add ref
message("\nEvaluating Nakagawa goodness-of-fit, refitting modules...")

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
#sum(!idx) == 0
gof_df <- as.data.table(do.call(rbind,gof_list[idx]),keep.rownames="Protein")

# combine varpart -- the precent variance explained for each covariate and
# gof_df -- the overall R2 for total and fixed effects (which is the percent
# variance explained by our models).
results_df <- left_join(varpart_df,gof_df,by="Protein")


## identify proteins with poor fit ---------------------------------------------

# protein level stuff

loss <- function(R2_threshold,results_df,as_char=FALSE) {
	# a function that determines the number of proteins remove at a given r2
	idx <- results_df$R2.total< R2_threshold
	poor_prots <- results_df$Protein[idx]
	if (as_char) {
	  return(poor_prots)
	} else {
	  return(sum(idx))
	}
}

# how many removed at various thresh
thresh_range <- seq(0,1,length.out=100) # 12.78% loss
mylist = sapply(thresh_range,loss,results_df,as_char=TRUE)
x = sapply(mylist,function(x) sum(x %in% sigprots))
r2_threshold <- thresh_range[head(which(x>100),1)]

# proteins with poor fit
poor_prots = loss(r2_threshold,results_df,as_char=TRUE)
message("\nnumber of proteins with poor fits: ", length(poor_prots))

## fit all module models -------------------------------------------------------

# analysis of intra-module variance

modules = split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))

# a function the evaluates gof of modules with variancePartition
moduleGOF <- function(module,msstats_prot){
  ## Variance partition expects all factors to be modeled as mixed effects
  form1 <- Abundance ~ (1|Mixture) + (1|Genotype) + (1|BioFraction) + (1|Protein)
  fm1 <- lmer(form1,msstats_prot %>% filter(Module==module))
  vpart <- variancePartition::calcVarPart(fm1)
  form2 <- Abundance ~ 0 + Genotype:BioFraction + (1|Mixture) + (1|Protein)
  fm2 <- lmer(form2,msstats_prot %>% filter(Module==module))
  r2 <- setNames(as.numeric(r.squaredGLMM.merMod(fm2)),
		 nm=c("R2.fixef","R2.total"))
  rho <- c(vpart,r2)
  # fit 1 is the model for variancePartition; fit2 is the fit used for stats
  rho[["fit1.isSingular"]] = lme4::isSingular(fm1)
  rho[["fit2.isSingular"]] = lme4::isSingular(fm2)
  return(rho) # gof stats
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

# lot o boundary -- seems to be coming from variancePartition side of things
# lme4 doesnt like to fit the variance partition formula... that means that 
# protein does explain much var within module -- this might be good or bad... 
# OR mixture more likely 

# drop null results
idx <- sapply(results_list,is.null)
message("There were problems fitting ", sum(idx), " models.")

# collect results
df <- as.data.table(do.call(rbind,results_list[!idx]),keep.rownames="Module")
df <- df %>% arrange(desc(Genotype))
module_gof <- df


#------------------------------------------------------------------------------
# above we calculated the R2 for fixed and total affects, but we also are
# curious to know how much variation is explained by our clustering

# how much variation is explained by partition?
form0 <- formula(Abundance ~ (1|Mixture) + (1|Protein)) # what is the percent variance explained by each protein -- Genotype:BioFraction
form1 <- formula(Abundance ~ (1|Mixture) + (1|Module)) # how close do we get with modules -- if every protein is its own cluster then should be the same.
fm0 <- lmer(form0, data=msstats_prot %>% filter(Module != "M0"))
fm1 <- lmer(form1, data=msstats_prot %>% filter(Module != "M0"))
rho0 <- calcVarPart(fm0)
rho1 <- calcVarPart(fm1)

#knitr::kable(t(round(rho0,3))) # protein explains .647 of variance

#knitr::kable(t(round(rho1,3))) # module explains .123 of variance 

x = rho0["Protein"] # observed (experimental maximum)

y = rho1["Module"] # captured by module 

y/x # we capture about 20% of the variance with our partition?

# this criterion is maximized when every protein is its own module
# it would be interesting to see how this criterior looks across resolutions
# i.e. with cpm method 

#rho["Module"] # this is the quantity we want to maximize!
# 0.121 or 11 percent
# NOTE decrases to 0.08 with without recursive split


## save results -----------------------------------------------------------------

# save as rda

#myfile <- file.path(root,"data","poor_prots.rda")
#save(poor_prots,file=myfile,version=2)

## save protein-level results

# save as rda
protein_gof <- results_df
myfile <- file.path(root,"data","protein_gof.rda")
save(protein_gof, file=myfile, version=2)

# save the data
fwrite(protein_gof,file.path(root,"rdata","protein_gof.csv"))


# save the data
fwrite(df,file.path(root,"rdata","module_gof.csv"))

# save as rda
module_gof <- df
myfile <- file.path(root,"data","module_gof.rda")
save(module_gof, file=myfile, version=2)
