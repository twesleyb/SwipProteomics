tidyProt <- function(raw_data,id.vars,species=NULL,
		     samples=NULL,summary=FALSE){
	#----------------------------------------------------------------------
	# title: tidyProt
	# description: a function to tidy proteomics data.
	#----------------------------------------------------------------------

	# misc function, lquote -----------------------------------------------
	lquote <- function(string, single = TRUE) {
	  # Wrap a string in single or double quotes.
	  single_quote <- "'"
	  double_quote <- "\""
	  if (single) {
	    # Single quote.
	    return(paste0(single_quote, string, single_quote))
	  } else {
	    # Double quote.
	    return(paste0(double_quote, string, double_quote))
	  }
	}

	# Imports --------------------------------------------------------------
	suppressPackageStartupMessages({
		library(data.table)
		#library(TBmiscr)
		library(tibble)
		library(dplyr)
	})

	# Collect relevant columns.
	#colNames <- colnames(PDdata)
	#colIDs <- c(grep(accessionCol,colNames),
	#	    grep(sequenceCol,colNames),
	#	    grep(modCol,colNames),
	#	    grep(dataCol,colNames))

	# Melt and fix column names.
	#dt <- PDdata %>% dplyr::select(colIDs) %>%  as.data.table()
	dt <- as.data.table(raw_data) %>%
		melt(id.vars = id.vars,
		     variable.name="Sample",
		     value.name="Intensity",
		     variable.factor=FALSE) # Don't coerce to factor.
	# Remove proteins that do not coorespond to species of interest.
	idx <- grepl(paste0("OS=",species),dt$Description)
	if (!all(idx)) {
		n_out <- length(unique(dt$Accession[!idx]))
		msg <- paste(n_out,"proteins are not from", lquote(species),
			     "and will be removed.")
		warning(msg,call.=FALSE)
		dt <- dt %>% filter(grepl(paste0("OS=",species),Description))
	}
	# Insure zeros are NA.
	dt$Intensity[dt$Intensity==0] <- NA
	#df <- dt %>% melt(id.vars = c(accessionCol,sequenceCol,modCol))
	#colnames(df)[c(4,5)] <- c("Sample","Intensity")
	# Add experiment, channel, and treatment annotations.
	#idx <- match(df$Sample,samples$Sample)
	#df <- tibble::add_column(df,Experiment=samples$Experiment[idx],.after=1)
	#df <- tibble::add_column(df,Channel=samples$Channel[idx],.after=1)
	#df <- tibble::add_column(df,Treatment=samples$Treatment[idx],.after=1)
	# Reorder columns.
	#df <- df %>%
	#	dplyr::select(c("Experiment","Sample","Channel","Treatment",
	#		 "Accession","Sequence","Modifications","Intensity"))
	#tp <- as.data.table(df)
	# Status.
	#nProt <- formatC(length(unique(dt$Accession)),big.mark=",")
	#nPep <- formatC(length(unique(tp$Sequence)),big.mark=",")
	#if (summary){
	#	message(paste("Number of unique peptides identified:", nPep))
	#	message(paste("Number of unique proteins identified:", nProt))
	#}
	# Return tidy df.
	return(as.data.table(dt))
}
