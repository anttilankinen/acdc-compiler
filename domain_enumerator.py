#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 14:05:57 2020

@author: antti
"""
import string
import numpy as np
import igraph

class Species:

    def initialise_strands(self):
        
        self.active_state_strand = ['0'] * 5
        self.state_strand = ['0'] * 5
        self.id_strand = ['0'] * 5
    
    def fill_domains(self, domain_counter):
        # fills currently empty domains with new domains
        # first need to figure which of the strands contains complementary
        # domains
        if domain_counter == 0:
            id_compl = True
        else:
            id_compl = any(['*' in d for d in self.id_strand])
            
        # first fill out stuff based on neighbours
        inherit_active = sum([d != '0' for d in self.active_state_strand]) > \
        sum([d != '0' for d in self.state_strand])
        print(inherit_active)
        
        if inherit_active:
            for i, d in enumerate(self.active_state_strand):
                if d == '0':
                    self.active_state_strand[i] = \
                    label(domain_counter, not(id_compl))
                    domain_counter += 1
            if len(self.activate_u) == 0 and len(self.repress_u) == 0:
                self.state_strand = self.active_state_strand
            else:
                self.state_strand[1:] = self.active_state_strand[1:]
                self.state_strand[2] = \
                swap_central(self.active_state_strand[2])
                self.state_strand[0] = label(domain_counter, not(id_compl))
                domain_counter += 1
            if not '2' in self.state_strand[4]:
                if id_compl:
                    self.state_strand[4] += '2'
                else:
                    self.state_strand[4] = \
                    self.state_strand[4][:-1] + '2*'
            
                    
        else:
            for i, d in enumerate(self.state_strand):
                if d == '0':
                    self.state_strand[i] = \
                    label(domain_counter, not(id_compl))
                    domain_counter += 1
                    if i == 2 and domain_counter == 3:
                        self.state_strand[2] = 'c1'
            if len(self.activate_u) == 0 and len(self.repress_u) == 0:
                self.active_state_strand = self.state_strand
            else:
                self.active_state_strand[1:] = self.state_strand[1:]
                self.active_state_strand[2] = \
                swap_central(self.state_strand[2])
                self.active_state_strand[0] = \
                label(domain_counter, not(id_compl))
                domain_counter += 1      
            if not '2' in self.active_state_strand[4]:
                if id_compl:
                    self.active_state_strand[4] += '2'
                else:
                    self.active_state_strand[4] = \
                    self.active_state_strand[4][:-1] + '2*'
                
            
            
        #for i, d in enumerate(self.state_strand):
         #   if d == '0':
          #      self.state_strand[i] =  label(domain_counter, not(id_compl))
                #if i == 2:
                #    self.state_strand[2] += '1'
           #     domain_counter += 1"""
        self.id_strand[1:4] = complement(self.state_strand[3:0:-1])
        self.id_strand[2] = 'c*' if id_compl else 'c'
        
        #if self.id_strand[2] == '0':
        #   self.id_strand[2] = 'c2*'
        for i in [0,4]:
            if self.id_strand[i] == '0':
                self.id_strand[i] = label(domain_counter, id_compl)
                domain_counter += 1
                
        return domain_counter
    
    def fill_complementary_domains(self):
        # todo:
        # add checks for active state strand
        """
            check which of the neighbours have some domains already and fill
            own domains accordingly 
        """
        # first downstream neighbours
        for s in self.activate_d:
            if not any([d == '0' for d in s.state_strand[4:1:-1]]):
                self.active_state_strand[0:3] = \
                complement(s.state_strand[4:1:-1],
                           switch_central=False)
                self.id_strand[4] = complement(s.id_strand[0],
                              switch_central=False)
                break
        for s in self.repress_d:
            if not any([d == '0' for d in s.active_state_strand[4:1:-1]]):
                self.active_state_strand[0:3] = \
                complement(s.active_state_strand[4:1:-1],
                           switch_central=False)
                self.id_strand[4] = complement(s.id_strand[0],
                              switch_central=False)
                break

        # then upstream neighbours
        for s in self.activate_u:
            if not any([d == '0' for d in s.active_state_strand[0:3]]):
                self.state_strand[4:1:-1] = \
                complement(s.active_state_strand[0:3],
                           switch_central=False)
                self.id_strand[0] = complement(s.id_strand[4],
                              switch_central=False)
                break
        for s in self.repress_u:
            if not any([d == '0' for d in s.active_state_strand[0:3]]):
                self.active_state_strand[4:1:-1] = \
                complement(s.active_state_strand[0:3],
                           switch_central=False)
                self.id_strand[0] = complement(s.id_strand[4],
                              switch_central=False)
                break

        
        
        
        
def label(n, compl=False):
    if n == 0:
        return 'a' + '*' if compl else 'a'
    digits = []
    
    while n:
        digits.append(int(n % 26))
        n //= 26
    if len(digits) > 1:
        digits[-1] -= 1
        
    letters = ''    
    for d in digits[::-1]:
        letters += string.ascii_lowercase[d]
    return letters + '*' if compl else letters

def complement(d, switch_central=False):
    if type(d) == str:
        if d == '0':
            return d
        if 'c' in d:
            if switch_central:
                return central_complement(d)
            
        if '*' in d:
            return d[:-1]
        else:
            return d + '*'
    elif type(d) == list:
        return [complement(c, switch_central=switch_central) for c in d]
        
def central_complement(d):
    if '*' in d:
        d = d[:-2] + str(int(5 - int(d[-2]))) + d[-1]
        return d[:-1]
    else:
        d = d[:-1] + str(int(5 - int(d[-1])))
        return d + '*'
        
def switch_state_central_domain(d):        
    return complement(central_complement(d), switch_central=False)
A = np.array([[0, 1,0 ],[0,0, 1], [0, 0, 0]])
g = igraph.Graph.Adjacency((A != 0).tolist())

# Add edge weights and node labels.
g.es['weight'] = A[A.nonzero()]    

def create_species(A):
    # create species from adjacency matrix A
    species_list = [Species() for i in range(A.shape[0])]
    for i, s in enumerate(species_list):
        s.activate_u = [species_list[j] for j in \
                        np.argwhere(A[:,i] == 1).reshape(-1).tolist()]
        s.repress_u = [species_list[j] for j in \
                       np.argwhere(A[:,i] == -1).reshape(-1).tolist()]
        s.activate_d = [species_list[j] for j in \
                        np.argwhere(A[i,:] == 1).reshape(-1).tolist()]
        s.repress_d = [species_list[j] for j in \
                       np.argwhere(A[i,:] == -1).reshape(-1).tolist()]
        s.initialise_strands()
        
    return species_list
        
def smallest_common(l1, l2):
    # smallest element of sorted list l1 that is contained in l2
    for k in l1:
        if k in l2:
            return k
        
def circle_forward(d):
    if len(d) == 3: # contains "*"
        if int(d[1]) < 3:
            return 'c' + str(int(d[1]) + 1)
        else:
            return 'c1'
    else:
        if int(d[1]) < 3:
            return 'c' + str(int(d[1]) + 1) + '*'
        else:
            return 'c1*'

def circle_backward(d):
    if len(d) == 3: # contains "*"
        if int(d[1]) > 1:
            return 'c' + str(int(d[1]) - 1)
        else:
            return 'c3'
    else:
        if int(d[1]) > 1:
            return 'c' + str(int(d[1]) - 1) + '*'
        else:
            return 'c3*'
        
def swap_central(d):
    print(d)
    if len(d) == 3:
        return d[0] + str(3 - int(d[1])) + d[2]
    else:
        return d[0] + str(3 - int(d[1]))
    
def enumerate_domains(species_list, g):
    domain_counter = 0
    not_enumerated = np.arange(1, len(species_list))
    domain_counter += species_list[0].fill_domains(domain_counter)
    neighbourhood = g.neighbors(0)

    for _ in range(len(not_enumerated)):
        i = smallest_common(not_enumerated, neighbourhood)
        species_list[i].fill_complementary_domains()
        domain_counter = species_list[i].fill_domains(domain_counter)
        not_enumerated = np.delete(not_enumerated,
                                   np.argwhere(not_enumerated == i))
        neighbourhood.extend(g.neighbors(i))
        neighbourhood = np.unique(neighbourhood).tolist()
        neighbourhood.sort()
    
    return species_list

def unique(l1, l2):
    unique_list = []
    for item in l1:
        if not item in unique_list:
            unique_list.append(item)
    for item in l2:
        if not item in unique_list:
            unique_list.append(item)
    return unique_list
    

def make_chains(species_list):
    chains = []
    for s in species_list:
        print('\n')
        print(chains)
        # if any of downstream repress active state strand is in any chain:
        # add self.active_state_strand to that chain
        # if any downstream activate passive state strand is in any chain:
        # add self.active_state_strand to that chain
        # if any upstream repressor active state strand is in any chain:
        # add self.active_state_strand to that chain
        as_neighbours = [k.active_state_strand
                                   for k in
                                   s.repress_d + s.repress_u] + \
                                  [k.state_strand for k in s.activate_d]
        chains = [c + [s.active_state_strand] if
                  any([a in c for a in as_neighbours]) else c
                  for c in chains]
        print()
        print(chains)
        
        if not any([s.active_state_strand in c for c in chains]):
            chains.append([s.active_state_strand])
        print()
        print(chains)
        # if any upstream activator active state strand is in any chain:
        # add self.state_strand to that chain
        ps_neighbours = [k.active_state_strand
                               for k in s.activate_u]
        
        chains = [c + [s.state_strand] if
                  any([p in c for p in ps_neighbours]) else c
                  for c in chains]
        print()
        print(chains)
        if not any([s.state_strand in c for c in chains]):
            chains.append([s.state_strand])
        print()
        print(chains)


        new_chains = []
        merged = []
        for i in range(len(chains)):
            for j in range(i+1, len(chains)):
                if len(chains[i]) + len(chains[j]) > \
                len(unique(chains[i],chains[j])) and not j in merged:
                    # overlap between chains, merge them
                    chains[i] = unique(chains[i], chains[j])
                    merged.append(j)
            new_chains.append(chains[i])
        chains = new_chains
        
    return chains
        
                    
        
species = create_species(A)
added_domains = enumerate_domains(species, g)

#for s in species:
#    print(s.state_strand)
#    print(s.id_strand[::-1])
#    print(s.active_state_strand)
#    print('')
    
chains = make_chains(species)
    
    
            
        