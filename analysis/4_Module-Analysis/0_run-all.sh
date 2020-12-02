#!/usr/bin/env bash

./1_module-lmerTest-analysis.R && \
	./2-module-variancePartition.R && \
	./3_module-gsea.R
