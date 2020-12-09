
root = "~/projects/SwipProteomics"
renv::load(root)

data(pd_psm)
data(swip)

library(dplyr)
library(data.table)

# pd_psm
idy = grepl("Abundance",colnames(pd_psm))
df = pd_psm %>%
	reshape2::melt(id.var=colnames(pd_psm)[!idy], 
		       variable.name="Channel", value.name="Abundance")
df$Channel <- as.character(df$Channel)
df$Channel <- gsub("Abundance\\.\\.","",df$Channel)
df$Protein <- df$"Master.Protein.Accessions"

subdf = df %>% subset(Protein == swip)

colnames(subdf)
subdf %>% group_by(Protein, Charge, Spectrum.File) %>% summarize(Intensity = sum(Abundance))



ID <- apply(df %>% select(all_of(colnames(df)[!val_cols])),1,paste)

	grep("Abundance",colnames(df)),
	grep("Master.Protein.Accessions",colnames(df)))

