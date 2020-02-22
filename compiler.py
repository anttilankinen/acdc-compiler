#!/usr/bin/env python
import argparse
from subprocess import call
from text2graph import text2graph
from topology_checker import is_valid
from domain_enumerator import make_domains
from make_nupack_script import make_nupack_script

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


args = get_args()
A, g, names = text2graph(args.input) # adjacency matrix, graph, names of nodes
result, error = is_valid(g)
if result:
    # YAY!
    species = make_domains(A, g, names, central_mismatch=args.central)
    make_nupack_script(species, central_mismatch=args.central, stop=args.stop)
    if args.nupack:
        call(['multitubedesign', 'script'])
else:
    print('Error: %s. Aborting' % error) 
    


