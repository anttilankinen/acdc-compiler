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
        if len(self.activate_u) > 0 or len(self.repress_u) > 0:
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
        if hasattr(self, 'active_state_strand'):
            if self.active_state_strand[0] == '0':
                self.active_state_strand[0] = label(domain_counter,
                                                    not(id_compl))
            if self.active_state_strand[1] != '0':
                self.state_strand[1] = self.active_state_strand[1]
            elif self.state_strand[1] != '0':
                self.active_state_strand[1] = self.state_strand[1]
        
            if self.active_state_strand[2] != '0':
                self.state_strand[2] = \
                switch_state_central_domain(self.active_state_strand[2]) 
            elif self.state_strand[2] != '0':
                self.active_state_strand[2] = \
                switch_state_central_domain(self.state_strand[2]) 
            if self.active_state_strand[3] != '0':
                self.state_strand[3] = self.active_state_strand[3]  
            elif self.state_strand[3] != '0':
                self.active_state_strand[3] = self.state_strand[3]
            print(self.active_state_strand)
            print(self.state_strand)
            if self.active_state_strand[4] != '0':
                if '*' in self.active_state_strand[4]:
                    self.state_strand[4] = self.active_state_strand[4][:-2] +\
                    self.active_state_strand[4][-1]
                else:
                    self.state_strand[4] = self.active_state_strand[4][:-1]
            elif self.state_strand[4] != '0':
                if '*' in self.state_strand[4]:
                    self.active_state_strand[4] = self.state_strand[4][:-1] +\
                    '2*'
                else:
                    self.active_state_strand[4] = self.state_strand[4] + '2'
                
            
            
        for i, d in enumerate(self.state_strand):
            if d == '0':
                self.state_strand[i] =  label(domain_counter, not(id_compl))
                if i == 2:
                    self.state_strand[2] += '1'
                domain_counter += 1
        self.id_strand[1] = complement(self.state_strand[3])
        self.id_strand[3] = complement(self.state_strand[1])
        
        if self.id_strand[2] == '0':
            self.id_strand[2] = 'c2*'
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
                try:
                    self.active_state_strand[0:3] = \
                    complement(s.state_strand[4:1:-1],
                               switch_central=False)
                    self.id_strand[4:1:-1] = complement(s.id_strand[0:3],
                                  switch_central=True)
                    break
                except:
                    self.state_strand[0:3] = \
                    complement(s.state_strand[4:1:-1],
                               switch_central=False)
                    self.id_strand[4:1:-1] = complement(s.id_strand[0:3],
                                  switch_central=True)
                    break
        for s in self.repress_d:
            if not any([d == '0' for d in s.active_state_strand[4:1:-1]]):
                try:
                    self.active_state_strand[0:3] = \
                    complement(s.active_state_strand[4:1:-1],
                               switch_central=False)
                    self.id_strand[4:1:-1] = complement(s.id_strand[0:3],
                                  switch_central=True)
                    break
                except:
                    self.state_strand[0:3] = \
                    complement(s.active_state_strand[4:1:-1],
                               switch_central=False)
                    self.id_strand[4:1:-1] = complement(s.id_strand[0:3],
                                  switch_central=True)
        # then upstream neighbours
        for s in self.activate_u:
            try:
                if not any([d == '0' for d in s.active_state_strand[0:3]]):
                    self.state_strand[4:1:-1] = \
                        complement(s.active_state_strand[0:3],
                                   switch_central=False)
                    self.id_strand[0:3] = complement(s.id_strand[4:1:-1],
                                  switch_central=True)
                    break
            except:
                if not any([d == '0' for d in s.state_strand[0:3]]):
                    self.state_strand[4:1:-1] = \
                        complement(s.state_strand[0:3],
                                   switch_central=False)
                    self.id_strand[0:3] = complement(s.id_strand[4:1:-1],
                                  switch_central=True)
                    break
        for s in self.repress_u:
            try:
                if not any([d == '0' for d in s.active_state_strand[0:3]]):
                    self.active_state_strand[4:1:-1] = \
                        complement(s.active_state_strand[0:3],
                                   switch_central=False)
                    self.id_strand[0:3] = complement(s.id_strand[4:1:-1],
                                  switch_central=True)
                    break
            except:
                if not any([d == '0' for d in s.state_strand[0:3]]):
                    self.active_state_strand[4:1:-1] = \
                        complement(s.state_strand[0:3],
                                   switch_central=False)
                    self.id_strand[0:3] = complement(s.id_strand[4:1:-1],
                                  switch_central=True)
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
g = igraph.Graph.Adjacency((A > 0).tolist())

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
        
    
def enumerate_domains(species_list):
    domain_counter = 0
    not_enumerated = np.arange(1, len(species_list))
    domain_counter += species_list[0].fill_domains(domain_counter)
    neighbourhood = g.neighbors(0)

    for _ in range(len(not_enumerated)):
        i = smallest_common(not_enumerated, neighbourhood)
        print('Now filling species %d' % i)
        species_list[i].fill_complementary_domains()
        domain_counter = species_list[i].fill_domains(domain_counter)
        not_enumerated = np.delete(not_enumerated,
                                   np.argwhere(not_enumerated == i))
        neighbourhood.extend(g.neighbors(i))
        neighbourhood = np.unique(neighbourhood).tolist()
        neighbourhood.sort()
    
    return species_list
        
species = create_species(A)
ss1 = species[1].state_strand
added_domains = enumerate_domains(species)
    
    
            
        