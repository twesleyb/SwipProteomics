#!/usr/bin/env bash
# Execute this script to run the analysis.

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

# Run the analysis 
./1_BioID-analysis.R &> iBioID-Proteomics.report & spin
./2_PPI-network.R &>> iBioID-Proteomics.report & spin
