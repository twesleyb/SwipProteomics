#!/usr/bin/env python3

'''
title: multi-leiden.py

description:

    Cluster input graph(s) with the Leiden algorithm.

    Cluster an adjacency matrix using the Leiden algorithm and one of
    seven quality metrics for optimization. A wrapper around leidenalg.

    This is an executable script that takes command line inputs.

author: tyler w a bradshaw

usage:
    multi-leiden [adjms] [methods] [OPTIONS]

input:
    * adjms: str, filepath to one or more input adjacency matrices as csv files
    * methods: str, Leiden alg quality metric for optimization of the graph(s)
output:
    * [output_name]_partition.csv
options:
    * n_iterations: int, the number of optimization iterations
'''


## imports --------------------------------------------------------------------

import sys
import numpy
import pickle
import igraph
import pandas
import argparse
import importlib
import leidenalg
import __main__ as main


## functions ------------------------------------------------------------------

def interactive():
    ''' check if python is being run interactively '''
    # @import __main__ as main
    # if executable, then main will have attr '__file__'
    is_interactive = not hasattr(main,'__file__')
    return is_interactive
#EOF


def pickle_args(args,args_file='.args.pickle'):
    ''' pickle argument dictionary returned by argparse '''
    # @import pickle
    with open(args_file,'wb') as stream:
        pickle.dump(args,stream,protocol=pickle.HIGHEST_PROTOCOL)
#EOF


def check_input(args):
    ''' check that an optimization method was specified for each input graph '''
    # @import sys
    if len(args.get('adjms')) != len(args.get('methods')):
        sys.exit('''An optimization method must be specified for each adjm.''')
    #EIS
#EOF


def check_methods(args,methods):
    ''' check that specified input methods match methods dictionary keys'''
    # @import sys
    is_valid_method=[method in methods.keys() for method in args.get('methods')]
    if not all(is_valid_method):
        sys.exit('''Specify an optimization method for each graph.
        Methods: 'Modularity','Surprise','RBConfiguration','RBER','CPM','Significance'
        ''')
    #EIS
#EOF


def load_args(args_file='.args.pickle'):
    ''' load pickled args '''
    # @import pickle
    with open (args_file, 'rb') as stream:
        args = pickle.load(stream)
    return args
#EOF


def parse_args():
    ''' Parse command line input with argparse.ArgumentParser '''
    # @import argparse
    ap = argparse.ArgumentParser('''
    Cluster an input adjacency matrix with the Leiden algorithm.
    ''')
    # required:
    ap.add_argument('adjms', type=str, nargs="+",
            help='path to input adjacency matrices as csv files')
    ap.add_argument('-m','--methods', type=str, nargs="+",
            required = True,
            help='the optimization method for clustering')
    # options:
    ap.add_argument('-n','--niter',type = int, default = -1,
            help = 'the number of optimization iterations')
    #ap.add_argument('-r','--resolution', nargs='+',
    #        default=[1], type=float, help = '''the resolution
    #        parameter for multiresolution methods''')
    ap.add_argument('-o', '--output', type=str,
            default='partition.csv', help =  'output filename')
    # collect and return arg dictionary
    args=vars(ap.parse_args())
    return args
#EOF


def graph_from_adjm(adjm, subset=None, weighted=True, signed=True):
    ''' Melts an input adjacency matrix given as a pandas dataframe
    into an edge list and then creates an igraph.Graph object using
    the Graph.TupleList function. While there are many ways to make
    a graph in igraph, construction of the igraph is done this way
    because it works reproducibly for me.
    '''
    # @import igraph
    if not signed:
    # If graph is unsigned (signed=F), then all edges must be +
        adjm = abs(adjm)
    # melt adjm into edge list dataframe
    edges = adjm.stack().reset_index()
    edges.columns = ['nodeA','nodeB','weight']
    # insure zero-value edges are removed
    edges = edges[edges.weight != 0]
    # if list of nodes to 'subset' is passed, subset the df by
    # collecting rows in which both 'nodeA' and 'nodeB' are in
    # 'subset'; combining nodeA, nodeB, and edge weight as a tuple
    if subset is not None:
        idx = edges['nodeA'].isin(subset) & \
                edges['nodeB'].isin(subset)
        edges = edges[idx]
    edge_tuples = list(zip(edges.nodeA,edges.nodeB,edges.weight))
    # construct the graph:
    if weighted:
        # create a weighted graph
        g = igraph.Graph.TupleList(edge_tuples,weights=True)
    else:
        # create an unweighted graph
        g = igraph.Graph.TupleList(edge_tuples,weights=False)
    # return an igraph object
    return g
#EOF


def load_LA_class(submodule,module='leidenalg'):
    ''' Load a Leidenalg partition class. '''
    # @import importlib
    # @import leidenalg
    partition_class = getattr(importlib.import_module(module),submodule)
    return partition_class
#EOF


def get_clustering_parameters(args,methods):
    '''
    get clustering parameters for an input graph and its optimization method
    returns: a list of clustering parameters
    NOTE: dict.copy is important to avoid objects referencing each other
          in the event the same method is used to optimize both graphs!
    '''
    params = [methods.get(m).copy() for m in args.get('methods')]
    return params
#EOF


def subset_adjm(adjm,subset):
    ''' subset a symmetric adjm given as a pandas df '''
    # @import pandas
    idx = adjm.index.isin(subset)
    sub_adjm = adjm[adjm.index[idx]].iloc[idx] # rows and cols
    return sub_adjm
# EOF


def is_connected(adjm):
    ''' retuns series of connected nodes in adjm as a pandas df '''
    # @import pandas
    unconnected = adjm.sum(axis=1) == 0
    nodes = adjm.columns[~unconnected]
    return nodes # connected nodes
# EOF


## define leidenalg methods dictionary -----------------------------------------
# NOTE: input 'methods' should match the methods dictionary keys
# NOTE: utilizes load_LA_class() to dynamically load the partition class

# Leidenalg supports the following clustering methods:
methods = {
        'Modularity' : { # optimization of modularity for weighted graphs
            'partition_type' : load_LA_class('ModularityVertexPartition'),
            'weighted' : True, # The graph can be weighted
            'signed' : False, # The graph can NOT be signed (contain -negative edges)
            'multi_resolution' : False}, # Single resolution only
        'Surprise' : { # optimization of surprise for weighted graphs
            'partition_type' : load_LA_class('SurpriseVertexPartition'),
            'weighted' : True, # The graph can be weighted
            'signed' : False, # The graph can NOT be signed
            'multi_resolution' : False}, # Single resolution only
        'RBConfiguration' : { # optimization of RBC for weighted, multiresolution graphs
            'partition_type' : load_LA_class('RBConfigurationVertexPartition'),
            'weighted' : True, # The graph can be weighted
            'signed' : False, # The graph can NOT be signed
            'multi_resolution' : True}, # multi-resolution
        'RBER' : { # optimization of RBER for weighted, multiresolution graphs
            'partition_type' : load_LA_class('RBERVertexPartition'),
            'weighted' : True, # The graph can be weighted
            'signed' : False, # The graph cannot be signed.
            'multi_resolution' : True}, # multi-resolution
        'CPM' : { # optimization of cpm for weighted, multiresolution graphs
            'partition_type' : load_LA_class('CPMVertexPartition'),
            'weighted' : True, # The graph cannot be weighted.
            'signed' : True, # The graph can be signed.
            'multi_resolution' : True}, # multi-resolution
        'Significance' : { # optimization of significance for unweighted graphs
            'partition_type' :  load_LA_class('SignificanceVertexPartition'),
            'weighted' :  False, # The graph cannot be weighted
            'signed' : False, # The graph cannot be signed
            'multi_resolution' : False} # single resolution only
        }


## parse input, get clustering params -----------------------------------------
#save args: pickle_args(args)
#interactive: args=load_args()

# parse input
args = parse_args()

# perform some checks on the input
check_input(args) # user should specify an optimization method for each graph
check_methods(args,methods) # user specified methods should match methods.keys()

# Get clustering parameters for each input graph
# each object in the 'clustering_params' list is a dictionary containing the
# clustering parameters for a given graph that will be passed to leidenalg
params = get_clustering_parameters(args,methods)

## load all adjacency matrices from file --------------------------------------
# NOTE: expect a csv with both a header row and an index column

adjms = [ pandas.read_csv(adjm,header=0,index_col=0) for adjm in args['adjms'] ]

## insure input adjms match ---------------------------------------------------
# if there are multiple graphs, only the union of their nodes will be analyzed,
# subset the graphs keeping overlapping nodes:

n = [ adjm.columns for adjm in adjms]
nodes = list(set(n[0]).intersection(*n))
adjms = [ subset_adjm(adjm,subset=nodes) for adjm in adjms ]

assert all([a.shape[0] == len(nodes) for a in adjms])

## identify  unconnected nodes ------------------------------------------------
# if nodes are unconnected, then this creates problems
# ensure that only connected nodes are kept, by finding all connected nodes

k = [is_connected(a).tolist() for a in adjms]
keep = list(set(k[0]).intersection(*k)) # finds union of nodes in list k

## build graphs ---------------------------------------------------------------
# loop to build graphs for each set of parameters:
# NOTE: this step can be time consuming

for i in range(len(adjms)):
    print("\nBuilding input graph: {}".format(i+1),file=sys.stderr)
    p = params[i]
    adjm = adjms[i]
    g = graph_from_adjm(adjm,subset=keep,weighted=p['weighted'],signed=p['signed'])
    # update graph building parameters
    p.update({'graph' : g })
    p.update({'weights' : 'weight'}) # Critical: la must knows about weights!
    del p['weighted']
    del p['signed']
#EOL


## perform clustering of the individual graphs --------------------------------

parts_list = list() # add partition for each graph to parts_list
for i in range(len(params)):

    msg = "\nUsing {} to find an optimal partition in graph {}."
    print(msg.format(args['methods'][i],i+1),file=sys.stderr)
    p = params[i].copy()

    if p.pop('multi_resolution') is True:
        # multiresolutiion methods:
        # NOTE: currently only supports analyzing a single resolution
        p.update({'resolution_parameter' : 1})
        partition = leidenalg.find_partition(**p)
        optimiser = leidenalg.Optimiser()
        diff = optimiser.optimise_partition(partition,n_iterations=args['niter'])
    else:

        # single resolution methods:
        partition = leidenalg.find_partition(**p)
        optimiser = leidenalg.Optimiser()
        diff = optimiser.optimise_partition(partition,n_iterations=args['niter'])

        partition.summary()

    #EIS
    # status report:
    msg = 'Initial partition: {}'
    print(msg.format(partition.summary()),file=sys.stderr)
    print('Modularity: {}'.format(partition.modularity))
    print('Quality: {}'.format(partition.quality()))
    parts_list.append(partition)
#EOL


## Optimize Multiplex partition ------------------------------------------------

# Given the two partitions, optimize multiplex.
# INPUT GRAPHS MUST BE DEFINED ON THE SAME VERTICES.
if len(parts_list) > 1:
    print("\nOptimizing multiplex partition.")
    optimiser = leidenalg.Optimiser()
    diff = optimiser.optimise_partition_multiplex(
            parts_list,
            layer_weights=[1 for i in range(len(params))], # all == 1
            n_iterations=args['niter'])
    # The input partitions will be updated.
    partition = parts_list[0]
    print("Final multiplex partition: {}".format(partition.summary()))
    print('Modularity: {}'.format(partition.modularity))
    print('Quality: {}'.format(partition.quality()))
    # Save final partition.
    membership = partition.membership
    df = pandas.DataFrame([membership])
    g0 = params[0]['graph']
    df.columns = g0.vs['name']
    df.to_csv(args['output'])
#EIS
