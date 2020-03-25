#!/usr/bin/env python
import argparse
import os
from subprocess import call
from text2graph import text2graph
from topology_checker import is_valid
from domain_enumerator import make_domains
from make_nupack_script import design_script

def get_args():
    parser = argparse.ArgumentParser(description='Check validity of an ACDC ' +
                                     'network and compute corresponding ' +
                                     'DNA sequences, if possible.')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Name of output file')
    parser.add_argument('-n', '--nupack', action='store_true')
    parser.add_argument('-c', '--central', action='store_true')
    parser.add_argument('-d', '--defect', type=int, default=5)
    return parser.parse_args()


args = get_args()
print('ACDC Compiler')
print('Input file:', args.input)
print('Target defect (%):', args.defect)
print('Central domain mismatches:', args.central)
print('Run NUPACK:', args.nupack)
print('Checking if graph is a valid ACDC graph...', end='', flush=True)
A, g, names = text2graph(args.input) # adjacency matrix, graph, names of nodes
result, error = is_valid(g)
if error is None:
    print(' Done')
else:
    print(' Error:', error, 'Exiting')
if result:
    # YAY!
    print('Enumerating domains...', end='', flush=True)
    species = make_domains(A, g, names, central_mismatch=args.central)
    print(' Done')
    print('Creating NUPACK script...', end='', flush=True)
    design_script(species, central_mismatch=args.central, stop=args.defect)
    print(' Done')
    if args.nupack:
        if os.path.isfile('design_0.npo'):
            os.remove('design_0.npo')
        print('Running NUPACK...', end='', flush=True)
        call(['multitubedesign', 'design.np'])
        if os.path.isfile('design_0.npo'):
            print(' Done')
            domains = [d if not '*' in d else d[:-1] for s in species 
               for d in s.state_strand + s.id_strand + s.active_state_strand]
            domains = list(set(domains))
            f = open('design_0.npo', 'r')
            npo = f.readlines()
            defect = float(npo[-5].split()[-1])
            f.close()
            if defect < args.defect / 100:
                print('Success! Created a network of %d ' % len(names) +
                  'species using %d domains.' % len(domains))
            else:
                print('Maximum number of iterations reached, did not reach ' +
                      'defect objective. Created a network of %d ' %
                      len(names) + 'species using %d domains.' % len(domains))
            print('Defect: %s' % npo[-5].split()[-1])
            
            
else:
    print('Error: %s. Aborting' % error) 
    


