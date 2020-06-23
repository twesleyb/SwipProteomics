#!/usr/bin/env Rscript

## Options:
niter = 500

# Load renv.
root <- getrd()
renv::load(root,quiet=TRUE)

# Load additional functions.
devtools::load_all()

# Loop to generate random colors on color.mind.
all_colors <- list()
pbar <- txtProgressBar(max=niter, style=3)
for (i in c(1:niter)) {
	# Generate colors with color.mind.api.
	script <- file.path(root,"bin","color-mind.sh")
	colors <- rjson::fromJSON(system(script,intern=TRUE))[[1]]
	# Collect the colors as hex codes.
	hex_colors <- sapply(colors, function(x) {
				     rgb(x[1],x[2],x[3], maxColorValue= 255) })
	all_colors[[i]] <- unique(hex_colors)
	setTxtProgressBar(pbar,value=i)
}
colormind <- unique(unlist(all_colors))

message(paste("\nCollected",length(colormind),"colors!"))

# Save the colors.
myfile <- file.path(root,"data","colormind.rda")
save(colormind,file=myfile, version=2)

## Examine the colors.
#library(ggplot2)
#plots <- list()
#for (color in unique_colors) {
#	p <- ggplot() + theme(panel.background = element_rect(fill = color))
#	plots[[color]] <- p + theme(panel.background = element_rect(fill=color))
#}
