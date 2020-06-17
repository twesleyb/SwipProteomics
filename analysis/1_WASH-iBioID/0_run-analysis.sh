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
	echo -e "\n"
}

# Remove any existing reports.
rm -f *.report

# Run the analysis 
echo "Processing iBioID proteomics data."
./1_*.R &> iBioID-Proteomics.report & spin

echo "Building PPI network."
./2_*.R &>> iBioID-Proteomics.report & spin

cat ./iBioID-Proteomics.report
