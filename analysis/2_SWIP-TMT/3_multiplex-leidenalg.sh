#!/usr/bin/env bash

## INPUTs:
root="$HOME/projects/SwipProteomics"

METHOD="Surprise"

ADJM0="$root/rdata/ne_adjm.csv"
ADJM1="$root/rdata/ppi_adjm.csv"

OUT0="$root/rdata/leidenalg_partition.csv"
OUT1="$root/rdata/multi_leidenalg_partition.csv"

# [1] run the script: just enhanced adjm
$root/Py/mleiden.py "$ADJM0" -m $METHOD -o $OUT0

# [2] run the script: multi-leiden (combined ne adjm and PPI adjm)
$root/Py/mleiden.py "$ADJM0" "$ADJM1" -m $METHOD $METHOD -o $OUT1
