#!/usr/bin/env bash
# 0_run-analysis.sh - execute this script to run the analysis.

# Output log files:
REPORT="Plotting.report"

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
echo "Generating protein, module, and community color assignments."
./1_*.R &> "$REPORT" & spin   

# STEP 2.
echo "Plotting protein abundance across all fractions."
./2_*.R &>> "$REPORT" & spin   

# STEP 3.
echo "Plotting protein PCA with module color annotations."
./3_*.R &>> "$REPORT" & spin

# STEP 4.
echo "Plotting all proteins from a module together."
./4_*.R &>> "$REPORT" & spin

# STEP 5.
echo "Plotting adjusted module intensity data to compare genotypes."
./5_*.R &>> "$REPORT" & spin

# STEP 6.
echo "Creating Cytoscape graphs."
./6_*.R &>> "$REPORT" & spin

cat "$REPORT"
