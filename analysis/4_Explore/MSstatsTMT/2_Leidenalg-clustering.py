#!/usr/bin/env python3
'''
title: Leidenalg Clustering 
description: clustering the protein network with Leidenalg + Surprise 
authors: Tyler WA Bradshaw
'''

## General optimization settings:
max_size = 100 # Maximum allowable size of a module.
optimization_method = 'Surprise'
n_iterations = -1 # The number of optimization iters 
# If -1, then the algorithm is run until no improvement further improvement.

## Recursive option:
recursive = False # If module_size > max_size, then cluster recursively.
recursive_method = 'Surprise'

## Input data:
# Input adjacency matrix should be in root/rdata/
adjm_file = "ne_adjm.csv"

## Output:
# Saved in root/rdata/
# [output_name]_partitions.csv

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
        # Modularity
        "Modularity": {'partition_type' : 'ModularityVertexPartition', 
            'weights' : True, 'signed' : False,
            'resolution_parameter' : None, 'n_iterations' : n_iterations},
        # Surprise
        "Surprise": {'partition_type' : 'SurpriseVertexPartition', 
            'weights' : True, 'signed' : False,
            'resolution_parameter' : None, 'n_iterations' : n_iterations}
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
# This is the method to be used to optimize the clustering.
parameters['partition_type'] = getattr(import_module('leidenalg'),method)

# Remove any None type parameters.
out = [key for key in parameters if parameters.get(key) is None]
for key in out: del parameters[key]

# Perform Leidenalg module detection. 
profile = list()
partition = find_partition(**parameters)
optimiser = Optimiser()
diff = optimiser.optimise_partition(partition,n_iterations=-1)
profile.append(partition)

# The intial partition:
initial_membership = partition.membership
subgraphs = partition.subgraphs()
too_big = [subg.vcount() > max_size for subg in subgraphs]
n_big = sum(too_big)

# Status.
msg="Using the {} method to split {} modules that contain more than {} nodes."
print(msg.format(recursive_method,n_big,max_size),file=stderr)
print("... Initial partition: " + partition.summary() + ".", file=stderr)

# Update optimization method.
method = methods.get(recursive_method).get('partition_type')
if type(method) is str:
    parameters['partition_type'] = getattr(import_module('leidenalg'),method)
elif type(method) == 'type':
    parameters['partition_type'] = method

# while loop to recursively split subgraphs that are too big.
while any(too_big):
    n_iter = 0
    idx = [i for i, too_big in enumerate(too_big) if too_big] 
    parameters['graph'] = subgraphs.pop(idx[0])
    part = find_partition(**parameters)
    optimiser = Optimiser()
    diff = optimiser.optimise_partition(part,n_iterations=-1)
    # Add to list.
    subgraphs.extend(part.subgraphs())
    too_big = [subg.vcount() > max_size for subg in subgraphs]
    n_iter += 1

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
