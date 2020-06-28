#!/usr/bin/env python3
'''
title: Leidenalg Clustering 
description: clustering the protein network with Leidenalg + Surprise 
authors: Tyler W A Bradshaw
'''

## Parameters: 
rmin = 0 # Min resolution for multi-resolution methods.
rmax = 1 # Max resolution for multi-resolution methods.
nsteps = 100 # Number of steps to take between rmin and rmax.
max_size = 100 # Maximum allowable size of a module.
recursive = False # If module_size > max_size, then cluster recursively.
n_iterations = -1 # Num of Lalg iter. If -1, then repeat until no improvement.
optimization_method = 'Modularity' # Optimization method.
# Methods: Modularity Surprise CPM RBConfiguration RBER Significance

## Input data:
# Input adjacency matrix should be in root/rdata/
adjm_file = 'filt_ne_adjm.csv'

## Output:
# Saved in root/rdata/
# [output_name]_partitions.csv
output_name = 'FiltNE' # Prefix out output partition, saved as .csv.

#------------------------------------------------------------------------------
## Prepare the workspace.
#------------------------------------------------------------------------------

# Imports.
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

# Directories.
here = os.getcwd()
root = dirname(dirname(here))
rdatdir = os.path.join(root,"rdata")
funcdir = os.path.join(root,"Py")

# Load user defined functions.
sys.path.append(root)
from Py.myfun import *

# Leidenalg supports the following optimization methods:
methods = {
        # 1. Modularity
        "Modularity": {'partition_type' : 'ModularityVertexPartition', 
            'weights' : True, 'signed' : False,
            'resolution_parameter' : None, 'n_iterations' : n_iterations},
        # 2. Surprise
        "Surprise": {'partition_type' : 'SurpriseVertexPartition', 
            'weights' : True, 'signed' : False,
            'resolution_parameter' : None, 'n_iterations' : n_iterations},
        # 3. RBConfiguration
        "RBConfiguration": {'partition_type' : 'RBConfigurationVertexPartition', 
            'weights' : True, 'signed' : False,
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps},
            'n_iterations' : n_iterations},
        # 4. RBER
        "RBER": {'partition_type' : 'RBERVertexPartition', 
            'weights' : True, 'signed' : False,
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps},
            'n_iterations' : n_iterations},
        # 5. CPM
        "CPM": {'partition_type' : 'CPMVertexPartition', 
            'weights' : True, 'signed' : True,
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps},
            'n_iterations' : n_iterations},
        # 6. Significance
        # FIXME: Significance method doesn't seem to be working.
        "Significance": 
        {'partition_type' : 'SignificanceVertexPartition', 
            'weights': None, 'signed' : False,
            'resolution_parameter' : None,
            'n_iterations' : n_iterations}
        }

# Get method specific parameters for clustering.
parameters = methods.get(optimization_method)
method = parameters.get('partition_type')

# Status report.
print("Performing Leidenalg clustering utilizing the {}".format(method),
        "method to find optimal partition(s)...", file=stderr)

#---------------------------------------------------------------------
## Load input adjacency matrix and create an igraph object.
#---------------------------------------------------------------------

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

#--------------------------------------------------------------------
## Community detection with the Leiden algorithm.
#--------------------------------------------------------------------

# Update partition type parameter.
# Dynamically load the partition_type class. 
# This is the method to be used for optimizing the clustering.
parameters['partition_type'] = getattr(import_module('leidenalg'),method)

# Remove any None type parameters.
out = [key for key in parameters if parameters.get(key) is None]
for key in out: del parameters[key]

# Perform Leidenalg module detection. 
if parameters.get('resolution_parameter') is None:
    # Single resolution methods.
    profile = list()
    partition = find_partition(**parameters)
    optimiser = Optimiser()
    diff = optimiser.optimise_partition(partition,n_iterations=-1)
    profile.append(partition)
    if not recursive: 
        print("... Initial partition: " + partition.summary() + ".", file=stderr)
    # Recursively split modules that are too big.
    if recursive:
        initial_membership = partition.membership
        subgraphs = partition.subgraphs()
        too_big = [subg.vcount() > max_size for subg in subgraphs]
        n_big = sum(too_big)
        print("Recursively spliting {} large module(s)...".format(n_big),
                file=stderr)
        print("... Initial partition: " + partition.summary() + ".",file=stderr)
        while any(too_big):
            # Perform clustering for any subgraphs that are too big.
            idx = [i for i, too_big in enumerate(too_big) if too_big] 
            parameters['graph'] = subgraphs.pop(idx[0])
            part = find_partition(**parameters)
            optimiser = Optimiser()
            diff = optimiser.optimise_partition(part,n_iterations=-1)
            # Add to list.
            subgraphs.extend(part.subgraphs())
            too_big = [subg.vcount() > max_size for subg in subgraphs]
        # Collect subgraph membership as a single partition.
        nodes = [subg.vs['name'] for subg in subgraphs]
        parts = [dict(zip(n,[i]*len(n))) for i, n in enumerate(nodes)]
        new_part = {k: v for d in parts for k, v in d.items()}
        # Set membership of initial graph.
        membership = [new_part.get(node) for node in partition.graph.vs['name']]
        partition.set_membership(membership)
        # Replace partition in profile list.
        profile[0] = partition
        print("... Final partition: " + partition.summary() + ".", file=stderr)
else:
    # Loop to perform multi-resolution clustering methods.
    pbar = ProgressBar()
    profile = list()
    resolution_range = linspace(**parameters.get('resolution_parameter'))
    for resolution in pbar(resolution_range):
        # Update resolution parameter.
        parameters['resolution_parameter'] = resolution
        partition = find_partition(**parameters)
        optimiser = Optimiser()
        diff = optimiser.optimise_partition(partition,n_iterations=-1)
        profile.append(partition)
        # Ends loop.
    print("... Final partition: " + partition.summary() + ".", file=stderr)
# Ends If/else.

#------------------------------------------------------------------------------
## Save Leidenalg clustering results.
#------------------------------------------------------------------------------

if recursive:
    # Save initial partition.
    df = DataFrame(columns = profile[0].graph.vs['name'])
    df.loc['Membership'] = initial_membership
    myfile = os.path.join(rdatdir, output_name + "_initial_partition.csv")
    df.to_csv(myfile)

# Collect partition results and save as csv. 
if len(profile) == 1:
    # Single resolution profile:
    results = {
            'Modularity' : [partition.modularity for partition in profile],
            'Membership' : [partition.membership for partition in profile],
            'Summary'    : [partition.summary() for partition in profile]}
else: 
    # Multi-resolution profile:
    results = {
        'Modularity' : [partition.modularity for partition in profile],
        'Membership' : [partition.membership for partition in profile],
        'Summary'    : [partition.summary() for partition in profile],
        'Resolution' : [partition.resolution_parameter for partition in profile]}
# Ends if/else

# Save cluster membership vectors.
myfile = os.path.join(rdatdir, output_name + "_partition.csv")
df = DataFrame(results['Membership'])
df.columns = profile[0].graph.vs['name']
df.to_csv(myfile)