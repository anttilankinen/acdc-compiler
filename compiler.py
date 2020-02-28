#!/usr/bin/env python
import argparse
import os
from subprocess import call
from text2graph import text2graph
from topology_checker import is_valid
from domain_enumerator import make_domains
from make_nupack_script import design_script, defect_script

def get_args():
    parser = argparse.ArgumentParser(description='Check validity of an ACDC ' +
                                     'network and compute corresponding ' +
                                     'DNA sequences, if possible.')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Name of output file')
    parser.add_argument('-n', '--nupack', action='store_true')
    parser.add_argument('-c', '--central', action='store_true')
    parser.add_argument('-s', '--stop', type=int, default=5)
    return parser.parse_args()


#args = get_args()
A, g, names = text2graph('testgraph2.txt') # adjacency matrix, graph, names of nodes
result, error = is_valid(g)
if result:
    # YAY!
    species = make_domains(A, g, names, central_mismatch=True)
    design_script(species, central_mismatch=True, stop=5)
    if True:
        if os.path.isfile('design_0.npo'):
            os.remove('design_0.npo')
        call(['multitubedesign', 'design.np'])
        if os.path.isfile('design_0.npo'):
            domains = [d if not '*' in d else d[:-1] for s in species 
               for d in s.state_strand + s.id_strand + s.active_state_strand]
            domains = list(set(domains))
            print('Success! Created a network of %d ' % len(names) +
                  'species using %d domains.' % len(domains))
            f = open('design_0.npo', 'r')
            npo = f.readlines()
            f.close()
            print('Defect: %s' % npo[-5].split()[-1])
            
            
else:
    print('Error: %s. Aborting' % error) 
    


