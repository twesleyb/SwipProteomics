#!/usr/bin/env bash

# run the analysis
./0_PD-data-preprocess.R && \
	./1_MSstatsTMT-analysis.R && \
	./2_protein-gof.R && \
	./3_generate-networks.R && \
	./4_leidenalg-clustering.py && \
	./5_post-leidenalg.R && \
	./6_module-preservation.R && \
	./7_lmerTest-module-analysis.R && \
	./8_module-gof.R && \
	./9_module-GSEA.R 
