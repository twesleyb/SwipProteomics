#!/usr/bin/env bash

# Run the analysis.
./1_BioID-analysis.R &> BioID.report
./2_PPI-network.R &> PPI.report
