#!/usr/bin/env bash

echo "Clustering protein covariation network..." && \
	./1_generate-network.R && \
	./2_leidenalg-clustering.py && \
	./3_post-leidenalg.R
