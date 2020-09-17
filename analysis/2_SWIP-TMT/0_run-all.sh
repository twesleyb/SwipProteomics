#!/usr/bin/env bash
: '
title: SynaptopathyProteomics/1_Main-Analysis/0_run-all.sh
description: execute this script to run the analysis
author: tyler w bradshaw

* 0_run-all.sh
* 1_data-preprocessing.R
* 2_network-generation.R
* 3_leidenalg-clustering.py
* 4_module-preservation.R
* 5_Module-GLM-analysis.R
* 6_Module-GSEA.R
'


## INPUT -------------------------------------------------------------

# preprocess the both the 'Cortex' and 'Striatum' data
CORTEX='Cortex'
STRIATUM='Striatum'

# use the 'Surprise' METHOD to optimize partitions
METHOD="Surprise"

# input adjm in root/rdata
CORTEX_ADJM="~/projects/SynaptopathyProteomics/rdata/cortex_ne_adjm.csv"
STRIATUM_ADJM="~/projects/SynaptopathyProteomics/rdata/striatum_ne_adjm.csv"
COX_PPI_ADJM="~/projects/SynaptopathyProteomics/rdata/cortex_ppi_adjm.csv"
STR_PPI_ADJM="~/projects/SynaptopathyProteomics/rdata/striatum_ppi_adjm.csv"

# output
PARTITION="~/projects/SynaptopathyProteomics/rdata/multiplex_partition.csv"


## FUNCTIONS ---------------------------------------------------------

spin() {
	# title: spin
	# description: a simple progres spinner
	# author: William Pursell
	# reference: https://stackoverflow.com/questions/12498304/
	pid=$!
	spin='-\|/'
	i=0
	while kill -0 $pid 2>/dev/null
	do
		i=$(( (i+1) %4 ))
		printf "\r${spin:$i:1}"
		sleep 0.1 # NOTE: can be adjusted
	done
}


## MAIN ------------------------------------------------------------------------
# Run each analysis script, passing Cortex or Striatum as input.
# Redirect stdout and stderr from each script to an indiviual report. Tee is
# used so that stderr and stdout are also returned by this script.
# Use `spin` to monitor progress.
# NOTE: append to file with `tee -a`.

# exit when any command fails
set -e

# remove existing reports
read -p "Removing existing reports and run the analysis? [Y]es/[N]o:" -n 1 -r && echo
if [[ $REPLY =~ ^[Yy] ]]; then
	rm -f *.report
else
	exit 0
fi

# STEP 1
echo "Preprocessing 'Cortex' data." && \
	# run script, return stderr & stdout, tee to report, spin to monitor progress
	./1_*.R "Cortex" 2>&1 | tee .1_data-preprocessing.report & spin && \
	echo "Preprocessing 'Striatum' data." && \
	./1_*.R "Striatum" 2>&1 | tee -a .1_data-preprocessing.report & spin

# STEP 2
echo "Generating multiplex partition of the 'Cortex' and 'Striatum' graphs." && \
./2_*.py "$CORTEX_ADJM" "$STRIATUM_ADJM" "$COX_PPI_ADJM" "$STR_PPI_ADJM" \
        -m "Surprise" "Surprise" "Surprise" "Surprise"\
        --output "$PARTITION" 2>&1 | tee .2_multi-leiden.report & spin

# STEP 3
echo "Enforcing 'Cortex' module-self preservation." && \
	./3_*.R "$CORTEX" 2>&1  | tee .3_module-preservation.report & spin && \
	echo "Enforcing 'Striatum' module-self preservation." && \
	./3_*.R "$STRIATUM" 2>&1 | tee -a .3_module-preservation.report & spin

# STEP 4
echo "Performing module-level analysis of 'Cortex' modules." && \
	./4_*.R "$CORTEX" 2>&1 | tee .4_module-glm-analysis.report & spin && \
	echo "Performing module-level analysis of 'Striatum' modules." && \
	./4_*.R "$STRIATUM" & 2>&1 | tee -a .4_module-glm-analysis.report & spin

# STEP 5
echo "Analyzing 'Cortex' modules for gene list enrichment." && \
	./5_*.R "$CORTEX" 2>&1 | tee .5_module-GSEA.report & spin && \
	echo "Analyzing 'Striatum' modules for gene list enrichment." && \
	./5_*.R "$STRIATUM" 2>&1 | tee -a .5_module-GSEA.report & spin && \
	echo "Completed analysis."
# EOF
