#!/usr/bin/env bash

# Combine PDFs.
OUTPUT="Module_Graphs.pdf"
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="$OUTPUT" ./*.pdf
