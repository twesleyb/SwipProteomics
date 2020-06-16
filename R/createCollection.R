# Create anRichment gene reference collection object.
createRefCollection <- function(gene_list, description, 
			     data_source, organism) {
	# Imports.
	suppressPackageStartupMessages({
		library(anRichment)
	})
	# Check the input.
	if (is.null(names(gene_list))) { 
		stop("Please provide a named list of genes.")
	}
	# Loop to create geneSets. 
	geneSets <- list()
	for (i in c(1:length(gene_list))) {
		geneSets[[i]] <- newGeneSet(geneEntrez = gene_list[[i]],
					    geneEvidence = "IC",
					    geneSource = "Custom Gene List",
					    ID = names(gene_list)[i],
					    name = names(gene_list)[i],
					    description = description,
					    source = data_source,
					    organism = organism,
					    internalClassification="myGenes",
					    groups = "myGroup",
					    lastModified = Sys.Date())
	} # Ends loop.
	# Annotate with group.
	myGroup <- newGroup(name = "myGroup", 
			    description = description,
			    source = data_source)
	# Create gene reference collection.
	refCollection <- newCollection(geneSets, list(myGroup))
  return(refCollection)
}
