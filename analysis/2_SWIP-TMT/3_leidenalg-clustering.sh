#!/usr/bin/env bash
: '
title: 3_leidenalg-clustering.sh

usage: ./3_leidenalg-clustering.sh

description: this executable script performs clustering of the protein covariation 
graph using the Leidenalg python library. This script just passes the args to a
command line function mleiden.py that performs clustering of input adjacency 
matrices given as a csv files and specified methods. 

two appraoahces are specified below: 

1) perform clustering of the co-variation
graph using the Surprise quality statistic. Split any module with > 100 nodes
again recursively.

2) Perform clustering of the co-varation graph and PPI graph independently.
combine these partitions as a single, multiplex partition which optimizes the
quality of both graphs.
'

## INPUTs:
root="$HOME/projects/SwipProteomics" # project root
METHOD="Surprise" # optimization method

ADJM0="$root/rdata/ne_adjm.csv" # enhanced adjacency matrix
ADJM1="$root/rdata/ppi_adjm.csv" # ppi adjacency matrix

OUT0="$root/rdata/leidenalg_partition.csv" # network partition
OUT1="$root/rdata/multi_leidenalg_partition.csv" 

# [1] analyze the co-variation network with surprise + recursive
# FIXME: Modularity statistic needs to be updated
$root/Py/mleiden.py "$ADJM0" -m $METHOD -o $OUT0 --recursive 1

# [2] optimize the multiplex partition of the covariation and PPI graphs
# FIXME: I don't think multiplex is working after addition of recursive option
#$root/Py/mleiden.py "$ADJM0" "$ADJM1" -m $METHOD $METHOD -o $OUT1
