#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
## Installation of the anRichment library.
#------------------------------------------------------------------------------
# This script should serve as a guide to install the Horvarth lab's 
# anRichment package which is useful for Gene Ontology analysis.

# Adapted from: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/

# Run this script in project root for correct installation!
renv::load()

# First install some additional dependencies. 
dependencies <- c(
		  "TxDb.Hsapiens.UCSC.hg19.knownGene",
		  "TxDb.Mmusculus.UCSC.mm10.knownGene",
		  "XML",
		  "AnnotationDbi",
		  "GO.db",
		  "org.Hs.eg.db",
		  "org.Mm.eg.db",
		  "WGCNA")

# A Function to install dependencies. 
rip <- function(package, method = "utils", ...) {
	# Install a R package. Supports the following methods:
	#     utils::install.packages()
	#     BiocManager::install()
	#     devtools::install_github()
	#     source - installs the package from Cran provided its source url, 
	#              this method depends upon the Linux bash utility, rip..
	# If method is source, parse the package name from its url.
	if (method == "source") {
		url <- package
		package <- strsplit(strsplit(url,"/")[[1]][6],"_")[[1]][1] 
	}
	# Insure that the package is not already installed.
	if (requireNamespace(package, quietly = TRUE)) {
		message(paste(package,"is already installed!"))
	} else if (method == "BiocManager") {
		BiocManager::install(package, ...)
	} else if (method == "utils") {
		utils::install.packages(package, ...)
	} else if (method == "devtools") {
		devtools::install_github(package, ...)
	} else if (method == "source") {
		cmd <- paste("rip", url, ...)
		system(cmd)
	} else stop("problem installing package")
}

sapply(dependencies,function(x) rip(x, method = "BiocManager"))

# Download AnRichment source code.
#urls <- c("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/anRichmentMethods_0.90-1.tar.gz",
#	  "https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/anRichment_1.01-2.tar.gz")
# Updated urls:
urls <- c("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/anRichmentMethods_0.91-94.tar.gz",
	  "https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/anRichment_1.10-1.tar.gz")

download.file(urls[1], basename(urls[1]))
download.file(urls[2], basename(urls[2]))

untar(basename(urls[1]))
untar(basename(urls[2]))

# Install from source code. 
dir <- getwd()
install.packages(paste(dir,"anRichmentMethods", sep="/"), 
		 repos = NULL, type = "source")
install.packages(paste(dir,"anRichment", sep="/"), 
		 repos = NULL, type = "source")

# Remove temporary files. 
#unlink("anRichment_1.01-2.tar.gz")
#unlink("anRichmentMethods_0.90-1.tar.gz")
#unlink("anRichmentMethods", recursive = TRUE)
