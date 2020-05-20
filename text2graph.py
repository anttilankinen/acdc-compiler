#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 17:41:46 2020

@author: antti
"""
import numpy as np
import igraph

def text2graph(filename):
    """
    create an iGraph graph based on a text file
    args:
        filename (str) - name of input file
    returns:
        A - numpy NDArray
        g - iGraph graph
        names - list
    """
    f = open(filename, 'r').readlines()
    # filter newlines out
    f = [line[:-1] if '\n' in line else line for line in f]
    # number of species in the network
    names = list(set([line.split()[i] for line in f for i in (0, -1)]))
    #names = list(set([line[i] for line in f for i in (-1, 0)]))
    names.sort()
    # create the adjacency matrix
    A = np.zeros([len(names), len(names)])
    for i, line in enumerate(f):
        if '->' in line: # activation reaction, edge weight 1
            A[names.index(line.split()[0]), names.index(line.split()[-1])] = 1
        elif '-|' in line: # repression reaction, edge weight -1
            A[names.index(line.split()[0]), names.index(line.split()[-1])] = -1
        else:
            raise ValueError('No relation symbol found on line %d of %s' %
                             (i, filename))
    # create graph
    g = igraph.Graph.Adjacency((A != 0).tolist())
    g.es['weight'] = A[A.nonzero()]   
    return A, g, names
    
    