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
	echo -e "\n"
}

# Remove any existing reports.
rm -f *.report

# STEP 1.
echo "Processing TMT data."
./1_*.R &> "$REPORT" & spin   

# STEP 2.
echo "Generating protein networks."
./2_*.R &>> "$REPORT" & spin   

# STEP 3.
echo "Clustering the protein co-variation network."
./3_*.py &>> "$REPORT" & spin

# STEP 4.
echo "Enforcing module preservation."
./4_*.R &>> "$REPORT" & spin

# STEP 5.
echo "Analyzing modules for changes in abundance."
./5_*.R &>> "$REPORT" & spin

# STEP 6.
echo "Analyzing modules for GO enrichment."
./6_*.R &>> "$REPORT" & spin

# STEP 7.
echo "Analyzing modules for enrichment of WASH BioID proteins."
./7_*.R &>> "$REPORT" & spin

# STEP 8.
echo "Analyzing modules for enrichment of NDD-associated genes."
./8_*.R &>> "$REPORT" & spin

# STEP 9.
echo "Generating Cytoscape networks."
#9_*.R &>> "$REPORT" & spin

