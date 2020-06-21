#!/usr/bin/env bash

# SVG's are saved from Cytoscape.
# Each SVG can be converted to a high quality tiff with svg2tiff.

# Convert to pdf.
# https://stackoverflow.com/questions/21828182/imagemagick-svg-to-pdf-conversion-image-quality-is-bad
FILES="$(ls ../../figs/Networks/*.svg)"
for file in $FILES
do
	echo "Converting "${file%.*}" to pdf."
	convert -density 343 $file ""${file%.*}".pdf"
done

# Combine pdfs.
echo "See this script for command to combine pdfs."

gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="$output" "$input"

# Remove temporary pdfs.
#rm ../../fig/Networks/*.pdf

# Move output file to final resting place.
mv $output ../../figs/Networks/

