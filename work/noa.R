#!/usr/bin/env Rscript

renv

data(lopit_prots)
data(swip_partition)

path = "ER"
df <- data.frame(Protein = names(partition))
df[[path]] <- df$Protein %in% lopit_prots[[path]]
fwrite(df,"noa.csv")
