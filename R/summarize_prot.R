# Summarize to protein level.
summarize_prot <- function(tp){
	# Imports.
	suppressPackageStartupMessages({
		library(data.table)
		library(dplyr)
	})
	# Sum to protein level.
	tp$Intensity[is.na(tp$Intensity)] <- 0
	tp <- ungroup(tp)
	proteins <- tp %>% group_by(Experiment,Sample,Channel,
				    Treatment,Accession) %>%
		summarize(Peptides = length(Intensity),
			  Intensity = sum(Intensity,na.rm=TRUE))
	proteins$Intensity[proteins$Intensity==0] <- NA
	return(proteins)
}
