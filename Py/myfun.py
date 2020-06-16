#!/usr/bin/env Python3

#--------------------------------------------------------------------
## xstr
#--------------------------------------------------------------------

def xstr(s):
    ''' Convert NoneType to blank ('') string.'''
    if s is None:
        return ''
    else:
        return str(s)
# EOF

#--------------------------------------------------------------------
## contains
#--------------------------------------------------------------------

def contains(mylist,value,return_index=False):
    ''' Check if list contains a value. 
    Like list.index(value) but returns False if the provided list
    does not contain value. 
    '''
    list_as_dict = dict(zip(mylist,range(len(mylist))))
    index = list_as_dict.get(value)
    if not return_index: 
        return type(index) is int # True if in list.
    else:
        return index
# EOF

#--------------------------------------------------------------------
## filterModules
#--------------------------------------------------------------------

import numpy as np 

def filterModules(partition,min_size=5,unassigned=0):
    """ Set modules with size less than minimum size to 0. """
    membership = np.array(partition.membership)
    sizes = np.array(partition.sizes())
    remove = np.array(range(len(sizes)))[sizes < min_size]
    out = [node in remove for node in membership]
    membership[out] = unassigned
    partition.set_membership(membership)
    m = np.array(partition.membership)
    unassigned = sum(m==unassigned)/len(m)
    return partition
# EOF

#--------------------------------------------------------------------
## graph_from_adjm
#--------------------------------------------------------------------

import numpy as np
from igraph import Graph
from pandas import DataFrame

def graph_from_adjm(adjm,weighted=True,signed=True):
    if not signed: adjm = abs(adjm)
    edges = adjm.stack().reset_index()
    edges.columns = ['nodeA','nodeB','weight']
    edges = edges[edges.weight != 0]
    edge_tuples = list(zip(edges.nodeA,edges.nodeB,edges.weight))
    if weighted: g = Graph.TupleList(edge_tuples,weights=True)
    if not weighted: g = Graph.TupleList(edge_tuples,weights=False)
    return g

#--------------------------------------------------------------------
## add_method
#--------------------------------------------------------------------

from functools import wraps 

def add_method(cls):
    ''' Add a method to an existing class.
    Example:
    @add_method(A)
    def foo():
        print('hello world!')
    '''
    def decorator(func):
        @wraps(func) 
        def wrapper(self, *args, **kwargs): 
            return func(*args, **kwargs)
        setattr(cls, func.__name__, wrapper)
        # Note we are not binding func, but wrapper which accepts self but does exactly the same as func
        return func # returning func means func can still be used normally
    return decorator

#--------------------------------------------------------------------
## apply_best_threshold
#--------------------------------------------------------------------

import igraph
import numpy as np
from igraph import Graph
from pandas import DataFrame

def apply_best_threshold(graph,num=20):
    ''' Function to apply a threshold to a graph such that it is a 
    single component.'''
    #FIXME: depends upon myfun!
    # Remove multiple edges.
    graph = graph.simplify(multiple=False) # has names.
    nodes = graph.vs['name']
    # Get adjm and get good start and stop for search space.
    adjm = np.array(graph.get_adjacency(attribute="weight").data)
    start = min(adjm.max(axis=1))
    stop = max(adjm.max(axis=1))
    # Loop to check if graph is single component after thresholding.
    i = 0
    is_connected = True
    space = np.linspace(start,stop,num)
    adjm = np.array(graph.get_adjacency(attribute="weight").data)
    while is_connected:
        threshold = space[i]
        mask = adjm > threshold
        df = DataFrame(adjm * mask, index=nodes,columns=nodes)
        edges = df.stack().reset_index()
        edges.columns = ['nodeA','nodeB','weight']
        edges = edges[edges.weight != 0]
        edge_tuples = list(zip(edges.nodeA,edges.nodeB,edges.weight))
        subg_filt = Graph.TupleList(edge_tuples,weights=True)
        components = subg_filt.components()
        n_components = len(set(components.membership))
        is_connected = n_components == 1
        if is_connected: i +=1
    # Ends while loop.
    # Get best threshold.
    if i == 0:
        threshold = space[i]
    else:
        threshold = space[i-1]
    # Apply best threshold.
    adjm = np.array(graph.get_adjacency(attribute="weight").data)
    mask = adjm > threshold
    df = DataFrame(adjm * mask, index=nodes,columns=nodes)
    edges = df.stack().reset_index()
    edges.columns = ['nodeA','nodeB','weight']
    edge_tuples = list(zip(edges.nodeA,edges.nodeB,edges.weight))
    g_filt = Graph.TupleList(edge_tuples,weights=True)
    return g_filt
# Ends function.

#--------------------------------------------------------------------
## Cluster MCL
#--------------------------------------------------------------------

import os
import igraph
import subprocess
import numpy as np
from pandas import DataFrame
from importlib import import_module
from leidenalg import VertexPartition

def clusterMCL(graph, inflation=1.2, weight='weight',
        quality_measure='Modularity', ncores =8, quiet=True):
    ''' A function to perform MCL clustering.'''
    quality_measures = {
            # Modularity
            "Modularity": 'ModularityVertexPartition', 
            # Surprise
            "Surprise": 'SurpriseVertexPartition', 
            # RBConfiguration
            "RBConfiguration": 'RBConfigurationVertexPartition', 
            # RBER
            "RBER": 'RBERVertexPartition', 
            # CPM
            "CPM": 'CPMVertexPartition', 
            # Significance
            "Significance": 'SignificanceVertexPartition'}
    edges = graph.get_edgelist()
    nodes = [graph.vs[edge]['name'] for edge in edges]
    weights = graph.es[weight]
    edge_list = [node + [edge] for node,edge in zip(nodes,weights)]
    DataFrame(edge_list).to_csv(".tempnet.csv",sep="\t",header=False,index=False)
    # Execute MCL. Pipe stdout back into python.
    cmd = ["mcl", ".tempnet.csv","--abc", # input file.
            "-I", str(inflation), # inflation parameter.
            "-te", str(ncores), # number of cores.
            "-o","-"] # output file.
    if quiet:
        # Suppress stderr by piping to devnull.
        process = subprocess.Popen(cmd,stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL)
    else :
        process = subprocess.Popen(cmd,stdout=subprocess.PIPE)
    # Parse the output.
    out = process.communicate()
    os.remove(".tempnet.csv")
    modules = list(out)[0].decode('utf-8').split("\n") # decode
    modules = [module.split("\t") for module in modules]
    modules = [module for module in modules if module != ['']] # remove empty
    # Create partition.
    k = len(modules)
    module_size = [len(module) for module in modules]
    module_names = list(range(1,k+1))
    partition = dict(zip(sum(modules, []), 
        np.repeat(module_names, module_size, axis=0)))
    clusters = graph.clusters()
    membership = [partition.get(protein) for protein in graph.vs['name']]
    # Replace NoneType with 0.
    membership = [0 if m is None else m for m in membership]
    # Set membership.
    mcl_clusters = VertexPartition.CPMVertexPartition(graph)
    # Dynamically load partition type.
    PartitionType = getattr(import_module('leidenalg'),
        quality_measures.get(quality_measure))
    mcl_clusters = PartitionType(graph)
    mcl_clusters.set_membership(membership)
    mcl_clusters.recalculate_modularity()
    mcl_clusters.renumber_communities()
    return(mcl_clusters)
# Done.

#--------------------------------------------------------------------
## clusterMaxMCL
#--------------------------------------------------------------------

def clusterMaxMCL(graph,inflation):
    ''' A wrapper around clusterMCL to find best inflation 
    parameter.'''
    mcl_clusters = list()
    Q = list()
    for i in inflation:
        result = clusterMCL(graph,inflation=i)
        mcl_clusters.append(result)
        Q.append(result.recalculate_modularity())
    # Done.
    idx = Q.index(max(Q))
    best_Q = Q[idx]
    best_i = inflation[idx]
    best_clusters = mcl_clusters[idx]
    #print("Best Inflation : {}".format(best_i) +
    #        "\nBest Modularity: {}".format(best_Q))
    return(best_clusters)
# Done.

#--------------------------------------------------------------------
## thresholdGraph
#--------------------------------------------------------------------

import numpy as np
from pandas import DataFrame

# Second attempt at thresholding function. Removes edges until graph
# splits into multiple components. 
def thresholdGraph(graph,start=None,step_size=0.05,weight='weight',quiet=True):
    # Get graph's adjacency matrix.
    graph = graph.simplify(multiple=False)
    A = DataFrame(graph.get_adjacency(attribute='weight').data)
    # Determine a good starting place.
    if start is None: 
        start = max(A.min(axis=1))
    if not quiet:
        print("Starting search at edge weight: {}".format(start))
    cutoff = start
    not_connected = 0
    while not_connected < 1:
        # Apply a threshold.
        mask = A > cutoff
        Ax = A * mask
        # Diagonal of Degree matrix.
        D = np.diag(Ax.sum(axis=1))
        # Laplacian.
        L = D - Ax
        # EigenValues.
        evals = np.linalg.eig(L)[0] 
        # If EigenValue = 0 then unnconnected.
        not_connected = sum(evals==0)
        if not_connected < 1: cutoff += step_size
    # Ends loop.
    # Apply best cutoff.
    graph.es.select(weight_lt=cutoff).delete()
    # Check.
    L = np.matrix(graph.laplacian('weight'))
    evals = np.linalg.eig(L)[0] 
    # Status report.
    if not quiet:
        print("Final Edge weight cutoff        : {}".format(cutoff))
        print("Number of unconnected components: {}".format(sum(evals==0)))
    # Return thresholded graph.
    return(graph)
# Ends function.

#--------------------------------------------------------------------
## getrd
#--------------------------------------------------------------------

from os import getcwd, listdir
from os.path import isdir, join, dirname

def getrd(here=getcwd(),dpattern=".git",max_tries=5):
    ''' Find a git repositories root directory.'''
    # Initial params
    tries = 0
    root = False
    here = getcwd()
    # While loop to try and find project root.
    while not root and tries < max_tries:
        onlydirs = [d for d in listdir(here) if isdir(join(here, d))]
        root = dpattern in onlydirs
        if not root: 
            here = dirname(here)
            tries += 1
    # Return root directory.
    root_dir = here
    return(root_dir)
# EOF
