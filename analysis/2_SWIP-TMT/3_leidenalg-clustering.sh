#!/usr/bin/env bash

## INPUTs:
root="$HOME/projects/SwipProteomics" # project root
METHOD="Surprise" # optimization method

ADJM0="$root/rdata/ne_adjm.csv" # enhanced adjacency matrix
ADJM1="$root/rdata/ppi_adjm.csv" # ppi adjacency matrix

OUT0="$root/rdata/leidenalg_partition.csv" # network partition
OUT1="$root/rdata/multi_leidenalg_partition.csv" 

# [1] only analyze the co-variation network with surprise
$root/Py/mleiden.py "$ADJM0" -m $METHOD -o $OUT0 --recursive True

# [2] optimize the multiplex partition of the covariation and PPI graphs
#$root/Py/mleiden.py "$ADJM0" "$ADJM1" -m $METHOD $METHOD -o $OUT1
