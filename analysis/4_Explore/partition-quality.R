#!/usr/bin/env/Rscript

root <- getwd()
renv::load(root)
devtools::load_all(root)

# load the data
data(msstats_prot)

library(dplyr)
library(data.table)

# load the partitions
myfile <- file.path(root,"rdata","ne_cpm_partition.csv")
part_df <- fread(myfile,drop=1)

# parse to list
part_list <- unlist(apply(part_df,1,list),recursive=FALSE)

# need to remove small modules!
res <- list()
for (i in c(1:100)) {
  print(i)
  part <- part_list[[i]] + 1
  part[part %in% which(table(part) < 5)] <- 0
  if (length(unique(part)) == 1) {
	  res[[i]] <- NA
  } else {
    filt_part <- part[part!=0]
    if (length(unique(filt_part)) == 1) {
	  res[[i]] <- NA
    } else {
    df <- msstats_prot %>% filter(Protein %in% names(filt_part)) %>% mutate(Module = paste0("M",filt_part[Protein]))
    fx <- "Abundance ~ (1|Module)"
    vp <- variancePartition::calcVarPart(lmerTest::lmer(fx, df))
    res[[i]] <- vp
  }
}
}
