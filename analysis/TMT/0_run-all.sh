#!/usr/bin/env bash
# 0_run-analysis.sh - execute this script to run the analysis.

# Output log files:
REPORT="Swip-TMT.report"

# Define a progress spinner function.
# From William Pursell: https://stackoverflow.com/questions/12498304/
spin() {
	pid=$! # Process ID of previously executed command.
	spin='-\|/'
	i=0
	while kill -0 $pid 2>/dev/null
	do
		i=$(( (i+1) %4 ))
		printf "\r${spin:$i:1}"
		sleep 0.1 # Can be adjusted.
	done
}

# Remove any existing reports.
rm -f *.report

# STEP 1.
echo "Processing TMT data."
./1_data-preprocessing.R &> "$REPORT" & spin   

# STEP 2.
echo "Generating protein networks."
./2_network-generation.R &>> "$REPORT" & spin   

# STEP 3.
echo "Clustering the protein co-variation network."
./3_leidenalg-clustering.py &>> "$REPORT" & spin

# STEP 4.
echo "Enforcing module preservation."
./4_module-preservation.R &>> "$REPORT" & spin

# STEP 5.
echo "Analyzing modules for changes in abundance."
./5_module-diff-abundance.R &>> "$REPORT" & spin

# STEP 6.
echo "Analyzing modules for GO enrichment."
./6_module-go-analysis.R &>> "$REPORT" & spin

# STEP 7.
echo "Analyzing modules for WASH BioID enrichment."
./7_module-wash-enrichment.R &>> "$REPORT" & spin

# STEP 8.
echo "Analyzing modules for NDD gene enrichment."
./8_module-ndd-enrichment.R &>> "$REPORT" & spin

# STEP 9.
#echo "Clustering protein network using Leiden algorithm."
#9_create-cytoscape-graphs.R &>> "$REPORT" & spin

