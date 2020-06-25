#!/usr/bin/env bash

# SVG's are saved from Cytoscape.
# Convert SVG to PDF and then combine.

# First, convert each file.svg to pdf.
# The density argument is to preserve quality of the SVG image.
# https://stackoverflow.com/questions/21828182/imagemagick-svg-to-pdf-conversion-image-quality-is-bad
FILES="$(ls ~/projects/SwipProteomics/figs/Networks/SVG/*.svg)"
for file in $FILES
do
	echo "Converting "${file%.*}" to pdf."
	convert -density 343 $file ""${file%.*}".pdf"
done
