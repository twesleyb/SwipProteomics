#!/usr/bin/env bash

# SVG's are saved from Cytoscape.
# Each SVG can be converted to a high quality tiff with svg2tiff.
# If all SVGs are combined as a single pdf, its too big. Split every 1000
# proteins.

# First, convert each file.svg to pdf.
# The density argument is to preserve quality of the SVG image.
# https://stackoverflow.com/questions/21828182/imagemagick-svg-to-pdf-conversion-image-quality-is-bad
FILES="$(ls ~/projects/SwipProteomics/figs/Networks/SVG/*.svg)"
for file in $FILES
do
	echo "Converting "${file%.*}" to pdf."
	convert -density 343 $file ""${file%.*}".pdf"
done

# Move all pdfs to temporary PDF directory.
mkdir ~/projects/SwipProteomics/figs/Networks/PDF
mv ~/projects/SwipProteomics/figs/Networks/SVG/*.pdf ~/projects/SwipProteomics/figs/Networks/PDF/

# Combine PDFs in groups of 100.
OUT="~/projects/SwipProteomics/figs/Networks/001-100_Module_Graphs.pdf"
FILES="$(ls ~/projects/SwipProteomics/figs/Networks/SVG/*.pdf | head -100)"
mv "$FILES" ~/projects/SwipProteomics/figs/Networks/PDF/
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="$OUT" ~/projects/SwipProteomics/figs/Networks/PDF/*.pdf
rm "~/projects/SwipProteomics/figs/Networks/PDF/*.pdf"
