#!/usr/bin/env bash

# generate latex model with R::stargazer
./lmer-latex.R &> model.tex && \
	pdflatex model.tex && \
	# clean-up
	rm model.tex && \
	rm model.log && \
	rm model.aux 
