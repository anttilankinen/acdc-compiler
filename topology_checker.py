#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
topology checker for acdc networks
Created on Wed Jan 29 10:16:36 2020

@author: antti
"""
import numpy as np

def find_all_paths(graph, start, end, mode = 'OUT'):
    """
    find all simple paths between two nodes
    args:
        graph - igraph Graph
        start - int, index of start node
        end - int, index of end node
        mode - str, either 'OUT' or 'IN'. If 'IN', the roles of the nodes
            are reversed
    returns:
        all_paths - nested list
    """
    def find_all_paths_aux(adjlist, start, end, path):
        path = path + [start]
        if start == end:
            return [path]
        paths = []
        for node in adjlist[start] - set(path):
            paths.extend(find_all_paths_aux(adjlist, node, end,
                                            path))
        return paths
    adjlist = [set(graph.neighbors(node, mode = mode)) \
        for node in range(graph.vcount())]
    all_paths = []
    start = start if type(start) is list else [start]
    end = end if type(end) is list else [end]
    for s in start:
        for e in end:
            all_paths.extend(find_all_paths_aux(adjlist, s, e, []))
    return all_paths

def is_valid(A, g):
    """
    check if a graph is a valid ACDC network
    args:
        g - iGraph graph
    returns:
        bool, validity of the graph
        None / str, potential error message
    """
    if np.sum(np.abs(A.diagonal())) > 0:
        return False, 'Self-connections are not allowed.'
    for i in range(g.vcount()-1):
        for j in range(i+1, g.vcount()):
            paths_out = find_all_paths(g, i, j, mode='OUT')
            paths_in = find_all_paths(g, i, j, mode='IN')
            for path in paths_out + paths_in:
                edges = [(p, p+1) for p in path[:-1]]
                for e in edges[1:-1]:
                    if A[e[0], e[1]] == -1:
                        return False, 'Deactivation reaction in the middle '+ \
                    'of a cascade.'
            if len(paths_out) > 0:
                # some paths exists this way
                lengths_out = [len(k)-1 for k in paths_out]
                if len(lengths_out) > 1:
                    # a FFL exists
                    lengths_out_invalid_diff = \
                        any([(k-min(lengths_out)) % 2 for k in lengths_out])
                    if min(lengths_out) < 2:
                        return False, 'Short feedforward loop in network.'
                    elif lengths_out_invalid_diff:
                        return False, 'Invalid branch structure in network.'
            if len(paths_in) > 0:
                # some paths exists this way
                lengths_in = [len(k)-1 for k in paths_in]
                if len(lengths_in) > 1:
                    # a FFL exists
                    lengths_in_invalid_diff = \
                        any([(k-min(lengths_in)) % 2 for k in lengths_in])
                    if min(lengths_in) < 2:
                        return False, 'Short feedforward loop in network.'
                    elif lengths_in_invalid_diff:
                        return False, 'Invalid branch structure in network.'
            if len(paths_in) > 0 and len(paths_out) > 0:
                # a FBL exists
                loop_lengths = [k + l for k in lengths_out for l in lengths_in]
                if any(k % 2 for k in loop_lengths):
                    return False, 'Odd cycle in network.'
                elif any(k < 6 for k in loop_lengths):
                    return False, 'Short cycle in network.'
            
    return True, None
            
            
        





