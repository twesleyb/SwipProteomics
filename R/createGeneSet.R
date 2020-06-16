# Create geneSet object.
createGeneSet <- function(genes, pathway_name, 
			     description, data_source, organism = "mouse") {
  suppressPackageStartupMessages({
	library(anRichment)
  })
  geneSet <- newGeneSet(
    geneEntrez = genes,
    geneEvidence = "IEA",
    geneSource = "Custom Gene List",
    ID = pathway_name, # diseaseId
    name = pathway_name, # Shortened disease name
    description = description,
    source = data_source,
    organism = organism,
    internalClassification = "myGeneCollection",
    groups = "myGroup",
    lastModified = Sys.Date()
  )
  return(geneSet)
}
