#!/usr/bin/env bash

# R -- write a color to file as rbg with 4 unknowns.
# R -- execute script
# bash -- gets the colors.
# R -- parse the colors. Get a new seed?
# Repeat...

# Get 5 fresh colors.
INPUT='{"input":["N","N","N","N","N"],"model":"default"}'
COLORS="$(curl -q 'http://colormind.io/api/' --data-binary "$INPUT")"

# Pass input as text.
curl -q 'http://colormind.io/api/' --data-binary "@input.txt"
