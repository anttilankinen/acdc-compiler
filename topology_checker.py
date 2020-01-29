#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
topology checker for acdc networks
Created on Wed Jan 29 10:16:36 2020

@author: antti
"""


def find_all_paths(graph, start, end, mode = 'OUT', maxlen = None):
    def find_all_paths_aux(adjlist, start, end, path, maxlen = None):
        path = path + [start]
        if start == end:
            return [path]
        paths = []
        if maxlen is None or len(path) <= maxlen:
            for node in adjlist[start] - set(path):
                paths.extend(find_all_paths_aux(adjlist, node, end,
                                                path, maxlen))
        return paths
    adjlist = [set(graph.neighbors(node, mode = mode)) \
        for node in range(graph.vcount())]
    all_paths = []
    start = start if type(start) is list else [start]
    end = end if type(end) is list else [end]
    for s in start:
        for e in end:
            all_paths.extend(find_all_paths_aux(adjlist, s, e, [], maxlen))
    return all_paths

def is_valid(g):
    # check if a graph g is a valid acdc network
    for i in range(g.vcount()-1):
        for j in range(i+1, g.vcount()):
            paths_out = find_all_paths(g, i, j, mode='OUT')
            paths_in = find_all_paths(g, i, j, mode='IN')
            if len(paths_out) > 0:
                # some paths exists this way
                lengths_out = [len(k)-1 for k in paths_out]
                if len(lengths_out) > 1:
                    # a FFL exists
                    lengths_out_invalid_diff = \
                        any([(k-min(lengths_out)) % 2 for k in lengths_out])
                    if min(lengths_out) < 2 or lengths_out_invalid_diff:
                        # FFl shorter branch is too short or some longer branch
                        # has invalid length
                        return False
            if len(paths_in) > 0:
                # some paths exists this way
                lengths_in = [len(k)-1 for k in paths_in]
                if len(lengths_in) > 1:
                    # a FFL exists
                    lengths_in_invalid_diff = \
                        any([(k-min(lengths_in)) % 2 for k in lengths_in])
                    if min(lengths_in) < 2 or lengths_in_invalid_diff:
                        # FFl shorter branch is too short or some longer branch
                        # has invalid length
                        return False
            if len(paths_in) > 0 and len(paths_out) > 0:
                # a FBL exists
                loop_lengths = [k + l for k in lengths_out for l in lengths_in]
                if any(k % 2 for k in loop_lengths) or \
                    any(k < 6 for k in loop_lengths):
                    return False
            
    return True
            
            
        





