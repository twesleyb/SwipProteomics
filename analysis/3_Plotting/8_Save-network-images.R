#!/usr/bin/env Rscript

#' ---
#' title:
#' description:
#' authors: Tyler W Bradshaw
#' ---

## Analysis options:

#--------------------------------------------------------------------
## Set-up the workspace.
#--------------------------------------------------------------------

start <- Sys.time()
message(paste("Starting analysis at:",start))

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Global imports.
suppressPackageStartupMessages({
  library(RCy3)
  library(dplyr)
  library(igraph)
  library(data.table)
})

# Functions.
suppressWarnings({ devtools::load_all() })

# Project directories.
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Output directory for cytoscape networks.
netwdir <- file.path(root,"networks")
if (!dir.exists(netwdir)) {
	dir.create(netwdir)
}

#---------------------------------------------------------------------
## Combine modules images (SVG) as a single pdf.
#---------------------------------------------------------------------
# Convert to svg to pdf.

# Directory for grouped pdfs.
output_dir <- file.path(figsdir,"Final-Markers")
dir.create(output_dir,recursive=TRUE)

message("Combining markers and saving as single pdf...")

for (module in names(marker_modules)){
# Get subset of proteins.
prots <- paste0(gsub("\\|","_",marker_modules[[module]]),".pdf")
check <- all(prots %in% names(all_plots))
if (!check) { stop("We are missing some plots!") } 
# Create directory for grouped plots.
to_dir <- file.path(output_dir,gsub(" ","_",module))
dir.create(to_dir)
# Remove any existing plots.
unlink(list.files(to_dir,full.names=TRUE))
# Copy to new directory.
namen <- names(all_plots[prots])
response <- file.copy(from=all_plots[prots],
		      to=file.path(to_dir,namen))
if (!all(response)) { stop("Problem saving pdfs.") }

# Combined into a single pdf with ghostscript cli utility.
# Create gs command.
myfile <- file.path(to_dir,paste0(gsub(" ","_",module),".pdf"))
cmd <- c("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite",
	 "-sOutputFile=OUTPUT.pdf", "INPUT.pdf")
cmd <- gsub("OUTPUT.pdf",myfile,cmd)
cmd <- gsub("INPUT.pdf",file.path(to_dir,"*.pdf"),cmd)
cmd <- paste(cmd,collapse=" ")

# Execute system command
response <- system(cmd,intern=TRUE)

# Remove input pdfs.
unlink(file.path(to_dir,namen))

# Done!
end <- Sys.time()
message(paste("\nCompleted analysis at:",end))
