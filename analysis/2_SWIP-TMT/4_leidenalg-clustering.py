#!/usr/bin/env python3

'''
title: SwipProteomics
description: Leidenalg Clustering of the Enhanced Protein Covariation Network
authors: Tyler W A Bradshaw
'''

## Inputs ---------------------------------------------------------------------

## Project root:
root = "~/projects/SwipProteomics"

## Parameters for multiresolution methods:
## NOTE: some parameters are not used if they are not required by the set 
## optimization method
rmin = 0 # Min resolution for multi-resolution methods.
rmax = 1 # Max resolution for multi-resolution methods.
nsteps = 100 # Number of steps to take between rmin and rmax.
max_size = 100 # Maximum allowable size of a module.

## General optimization methods:
output_name = 'cpm' # Prefix out output partition, saved as .csv.
optimization_method = 'CPM'
n_iterations = -1  # Not the number of recursive iterations, but the number
# of optimization iterations.

## Recursive option:
recursive = False # If module_size > max_size, then cluster recursively.
recursive_method = 'Surprise'

## Input data:
# Input adjacency matrix should be in root/rdata/
#adjm_file = 'ne_adjm.csv'
adjm_file = 'adjm.csv' # only cpm is applicable bc we have negative and positive
# edges

## Output:
# Saved in root/rdata/
# [output_name]_partitions.csv


## Prepare the workspace ------------------------------------------------------

import os
import sys
import glob
from sys import stderr
from os.path import dirname

import numpy as np
from numpy import linspace
from importlib import import_module
from progressbar import ProgressBar
from leidenalg import Optimiser, find_partition

from igraph import Graph
from pandas import read_csv, DataFrame

# project Directories:
rdatdir = os.path.join(root,"rdata")
funcdir = os.path.join(root,"Py")

# Load user defined functions.
sys.path.append(root)
from myfun import * # try putting SwipProtomics in bashrc's python path
#from Py.myfun import *


## Leidenalg qualitity metrics ------------------------------------------------

# Leidenalg supports the following optimization methods:
methods = {
        # Modularity
        "Modularity": {'partition_type' : 'ModularityVertexPartition',
            'weights' : True, 'signed' : False,
            'resolution_parameter' : None, 'n_iterations' : n_iterations,
            'multi_resolution' : False },
        # Surprise
        "Surprise": {'partition_type' : 'SurpriseVertexPartition',
            'weights' : True, 'signed' : False,
            'resolution_parameter' : None, 'n_iterations' : n_iterations,
            'multi_resolution' : False },
        # RBConfiguration
        "RBConfiguration": {'partition_type' : 'RBConfigurationVertexPartition',
            'weights' : True, 'signed' : False,
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps},
            'n_iterations' : n_iterations,
            'multi_resolution' : True },
        # RBER
        "RBER": {'partition_type' : 'RBERVertexPartition',
            'weights' : True, 'signed' : False,
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps},
            'n_iterations' : n_iterations,
            'multi_resolution' : True },
        # CPM
        "CPM": {'partition_type' : 'CPMVertexPartition',
            'weights' : True, 'signed' : True,
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps},
            'n_iterations' : n_iterations,
            'multi_resolution' : True },
        # Significance
        # FIXME: Significance method doesn't seem to be working.
        "Significance":
        {'partition_type' : 'SignificanceVertexPartition',
            'weights': None, 'signed' : False,
            'resolution_parameter' : None,
            'n_iterations' : n_iterations,
            'multi_resolution' : False }
        }

# Get method specific parameters for clustering.
parameters = methods.get(optimization_method)
method = parameters.get('partition_type')

# Status report.
print("Performing Leidenalg clustering utilizing the {}".format(method),
        "method to find optimal partition(s).", file=stderr)


## Load input adjacency matrix and create an igraph object --------------------

# Load graph adjacency matrix.
myfile = os.path.join(rdatdir,adjm_file)
adjm = read_csv(myfile, header = 0, index_col = 0)

# Create igraph graph object and add to parameters dictionary.
# Note, this can take several minutes.
if parameters.get('weights') is not None:
    # Create a weighted graph.
    g = graph_from_adjm(adjm,weighted=True,signed=parameters.pop('signed'))
    parameters['weights'] = 'weight'
    parameters['graph'] = g
else:
    # Create an unweighted graph.
    g = graph_from_adjm(adjm,weighted=False,signed=parameters.pop('signed'))
    parameters['graph'] = g
#EIS

# the input graph
print(g.summary()) # NOTE: UNW is UNdirected and Weighted!


## Community detection with the Leiden algorithm -----------------------------

# FIXME: really we need two functions:
# * single-resolution (+/- recursive)
# * multi-resolution methods

# Update partition type parameter.
# Dynamically load the partition_type class.
# This is the method to be used to optimize the clustering.
parameters['partition_type'] = getattr(import_module('leidenalg'), method)

# Remove any None type parameters.
out = [key for key in parameters if parameters.get(key) is None]
for key in out: del parameters[key]


## Perform Leidenalg module detection -----------------------------------------

## MULTIRESOLUTION METHODS

if parameters.pop('multi_resolution') is True:
  # collect clustering params
  profile = list()
  pbar = ProgressBar()
  res_range = linspace(**parameters.get('resolution_parameter'))
  #n_iter = parameters.pop(n_iterations)
  # for res in pbar(linespace(**p.pop('resolution_range'))):
  for resolution in pbar(res_range):
        parameters['resolution_parameter'] = resolution
        partition = find_partition(**parameters)
        optimiser = Optimiser()
        diff = optimiser.optimise_partition(partition,parameters['n_iterations'])
        profile.append(partition)
  # EOL

 # else:

 # ## SINGLE RESOLUTION METHODS
 # # Single resolution methods: first iteration if recursive
 # profile = list()
 # partition = find_partition(**parameters)
 # optimiser = Optimiser()
 # diff = optimiser.optimise_partition(partition,n_iterations=-1)
 # profile.append(partition)
 # if not recursive:
 #   print("... Final partition: " + partition.summary() + ".", file=stderr)
 # # Recursively split modules that are too big.
 # if recursive:
 #    print("... Initial partition: " + partition.summary() + ".", file=stderr)
 #    # Update optimization method.
 #    method = methods.get(recursive_method).get('partition_type')
 #    if type(method) is str:
 #       parameters['partition_type'] = getattr(import_module('leidenalg'),
 #               method)
 #    elif type(method) == 'type':
 #            parameters['partition_type'] = method
 #    # Initial module membership.
 #    initial_membership = partition.membership
 #    subgraphs = partition.subgraphs()
 #    too_big = [subg.vcount() > max_size for subg in subgraphs]
 #    n_big = sum(too_big)
 #    msg = "\nSplitting {} modules that contain more than {} nodes."
 #    print(msg.format(n_big,max_size),file=stderr)
 #    while any(too_big):
 #        # Perform clustering for any subgraphs that are too big.
 #        idx = [i for i, too_big in enumerate(too_big) if too_big]
 #        parameters['graph'] = subgraphs.pop(idx[0])
 #        part = find_partition(**parameters)
 #        optimiser = Optimiser()
 #        diff = optimiser.optimise_partition(part,n_iterations=-1)
 #        # Add to list.
 #        subgraphs.extend(part.subgraphs())
 #        too_big = [subg.vcount() > max_size for subg in subgraphs]
 #    #EOL to split modules
 #    # Collect subgraph membership as a single partition.
 #    nodes = [subg.vs['name'] for subg in subgraphs]
 #    parts = [dict(zip(n,[i]*len(n))) for i, n in enumerate(nodes)]
 #    new_part = {k: v for d in parts for k, v in d.items()}
 #    # Set membership of initial graph.
 #    membership = [new_part.get(node) for node in partition.graph.vs['name']]
 #    partition.set_membership(membership)
 #    # Replace partition in profile list.
 #    profile[0] = partition
 #    print("... Final partition: " + partition.summary() + ".", file=stderr)
# EIS


## Save Leidenalg clustering results ------------------------------------------

# collect matrix in which each row is a partition (nrow = nRes)
df = DataFrame(columns = profile[0].graph.vs['name'])
for i in range(len(profile)):
    df.loc[i] = profile[i].membership
#EOL

# save the data as csv
myfile = os.path.join(rdatdir, output_name + "_partition.csv")
df.to_csv(myfile)
