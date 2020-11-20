#!/bin/bash

./1_generate-network.R && \
	./2_leidenalg-clustering.py && \
	./3_post-leidenalg.R && \
	./4_module-preservation.R
