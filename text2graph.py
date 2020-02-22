#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 17:41:46 2020

@author: antti
"""
import numpy as np
import igraph

def text2graph(filename):
    f = open(filename, 'r').readlines()
    # filter newlines out
    f = [line[:-1] if '\n' in line else line for line in f]
    # figure out how many 
    names = list(set([line[i] for line in f for i in (-1, 0)]))
    
    A = np.zeros([len(names), len(names)])
    for i, line in enumerate(f):
        if '->' in line:
            A[names.index(line[0]), names.index(line[-1])] = 1
        elif '-|' in line:
            A[names.index(line[0]), names.index(line[-1])] = -1
        else:
            raise ValueError('No relation symbol found on line %d of %s' %
                             (i, filename))
    g = igraph.Graph.Adjacency((A != 0).tolist())
    g.es['weight'] = A[A.nonzero()]   
    return A, g, names
    
    