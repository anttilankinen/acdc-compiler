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
    
    def __init__(self, name):
        self.name = name

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
        if len(self.activate_u) == 0 and len(self.repress_u) == 0:
                self.state_strand = self.active_state_strand
        
        """if inherit_active:
            for i, d in enumerate(self.active_state_strand):
                if d == '0':
                    self.active_state_strand[i] = \
                    label(domain_counter, not(id_compl))
                    domain_counter += 1
            if len(self.activate_u) == 0 and len(self.repress_u) == 0:
                self.state_strand = self.active_state_strand
            else:
                self.state_strand[1:] = self.active_state_strand[1:]
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
                    #if i == 2 and domain_counter == 3:
                    #    self.state_strand[2] = 'c1'
            if len(self.activate_u) == 0 and len(self.repress_u) == 0:
                self.active_state_strand = self.state_strand
            else:
                self.active_state_strand[1:] = self.state_strand[1:]
                self.active_state_strand[0] = \
                label(domain_counter, not(id_compl))
                domain_counter += 1      
            if not '2' in self.active_state_strand[4]:
                if id_compl:
                    self.active_state_strand[4] += '2'
                else:
                    self.active_state_strand[4] = \
                    self.active_state_strand[4][:-1] + '2*'"""
                    
        for i in range(len(self.state_strand)):
            if i == 0 and self.state_strand[0] == '0':
                self.state_strand[0] = \
                    label(domain_counter, not(id_compl))
                domain_counter += 1
            if i in [1, 2, 3] and self.active_state_strand[i] != '0' and \
            self.state_strand[i] == '0':
                self.state_strand[i] = self.active_state_strand[i]
            elif i in [1, 2, 3] and self.state_strand[i] != '0' and \
            self.active_state_strand[i] == '0':
                self.active_state_strand[i] = self.state_strand[i]
            elif i in [1, 2, 3] and self.state_strand[i] == '0' and \
            self.active_state_strand[i] == '0':
                self.state_strand[i] = \
                label(domain_counter, not(id_compl))
                self.active_state_strand[i] = self.state_strand[i]
                domain_counter += 1
            if i == 4 and self.state_strand[4] != '0' and \
            self.active_state_strand[4] == '0':
                if '2*' in self.state_strand[4]:
                    self.active_state_strand[4] = \
                    self.state_strand[4][:-2] + '*'
                elif '2' in self.state_strand[4]:
                    self.active_state_strand[4] = \
                    self.state_strand[4][:-1]
                elif '*' in self.state_strand[4]:
                    self.active_state_strand[4] = \
                    self.state_strand[4][:-1] + '2*'
                else:
                    self.active_state_strand[4] = \
                    self.state_strand[4] + '2'
            elif i == 4 and self.active_state_strand[4] != '0' and \
            self.state_strand[4] == '0':
                if '2*' in self.active_state_strand[4]:
                    self.state_strand[4] = \
                    self.active_state_strand[4][:-2] + '*'
                elif '2' in self.active_state_strand[4]:
                    self.state_strand[4] = \
                    self.active_state_strand[4][:-1]
                elif '*' in self.active_state_strand[4]:
                    self.state_strand[4] = \
                    self.active_state_strand[4][:-1] + '2*'
                else:
                    self.state_strand[4] = \
                    self.active_state_strand[4] + '2'
            elif i == 4  and self.active_state_strand[4] == '0' and \
            self.state_strand[4] == '0':
                if id_compl:
                    self.active_state_strand[4] = \
                    label(domain_counter) + '2'
                    self.state_strand[4] = \
                    label(domain_counter)
                    domain_counter += 1
                else:
                    self.active_state_strand[4] = \
                    label(domain_counter) + '2*'
                    self.state_strand[4] = \
                    label(domain_counter) + '*'
                    domain_counter += 1
                    
            
            
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
        if self.active_state_strand[0] == '0':
            # do this only here to make sure domain names run smoothly
                self.active_state_strand[0] = \
                    label(domain_counter, not(id_compl))
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
                complement(s.state_strand[4:1:-1])
                self.id_strand[4] = complement(s.id_strand[0])
                break
        for s in self.repress_d:
            if not any([d == '0' for d in s.active_state_strand[4:1:-1]]):
                self.active_state_strand[0:3] = \
                complement(s.active_state_strand[4:1:-1])
                self.id_strand[4] = complement(s.id_strand[0])
                break

        # then upstream neighbours
        for s in self.activate_u:
            if not any([d == '0' for d in s.active_state_strand[0:3]]):
                self.state_strand[4:1:-1] = \
                complement(s.active_state_strand[0:3])
                self.id_strand[0] = complement(s.id_strand[4])
                break
        for s in self.repress_u:
            if not any([d == '0' for d in s.active_state_strand[0:3]]):
                self.active_state_strand[4:1:-1] = \
                complement(s.active_state_strand[0:3])
                self.id_strand[0] = complement(s.id_strand[4])
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

def complement(d):
    if type(d) == str:
        if d == '0':
            return d
            
        if '*' in d:
            return d[:-1]
        else:
            return d + '*'
    elif type(d) == list:
        return [complement(c) for c in d]

def create_species(A, names):
    # create species from adjacency matrix A
    species_list = [Species(names[i]) for i in range(len(names))]
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
    members = []
    for l, s in enumerate(species_list):
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
        members = [members[i] + [l] if
                  any([a in chains[i] for a in as_neighbours]) else members[i]
                  for i in range(len(members))]
        
        if not any([s.active_state_strand in c for c in chains]):
            chains.append([s.active_state_strand])
            members.append([l])
        # if any upstream activator active state strand is in any chain:
        # add self.state_strand to that chain
        ps_neighbours = [k.active_state_strand
                               for k in s.activate_u]
        
        chains = [c + [s.state_strand] if
                  any([p in c for p in ps_neighbours]) else c
                  for c in chains]
        members = [members[i] + [l] if
                  any([p in chains[i] for p in ps_neighbours]) else members[i]
                  for i in range(len(members))]
        if not any([s.state_strand in c for c in chains]):
            chains.append([s.state_strand])
            members.append([l])


        new_chains = []
        new_members = []
        merged = []
        for i in range(len(chains)):
            for j in range(i+1, len(chains)):
                if len(chains[i]) + len(chains[j]) > \
                len(unique(chains[i],chains[j])) and not j in merged:
                    # overlap between chains, merge them
                    chains[i] = unique(chains[i], chains[j])
                    members[i] = unique(members[i], members[j])
                    merged.append(j)
            new_chains.append(chains[i])
            new_members.append(members[i])
        chains = new_chains
        members = new_members
        
    return chains, members

def move_to_centre(l, i):
    # move i to centre of list l
    l[l.index(i)] = l[1]
    l[1] = i
    return l

def suitable_domains(l):
    # return all numbers 1,2,3 that aren't on list l
    return [i for i in [2, 3, 4] if i not in l]


def allocate_central_domains(chains, members, species):
    differ_pairs = [[i for i in range(len(members)) if j in members[i]]
    for j in range(len(species))]
    differ_pairs = [pair for pair in differ_pairs if len(pair) == 2]
    # these have to be made such that the one that the one in focus
    # is in the middle
    if len(differ_pairs) == 1:
        # there's only two chains, so just set them to be different and you're
        # good
        for i in range(2):
            for strand in chains[i]:
                if '*' in strand[2]:
                    strand[2] = 'c' + str(i+2) + '*'
                else:
                    strand[2] = 'c' + str(i+2)
    else:
        triplets = []
        for c in range(len(chains)):
            for i in range(len(differ_pairs)):
                for j in range(i+1, len(differ_pairs)):
                    if c in differ_pairs[i] and c in differ_pairs[j]:
                        triplets.append(
                                move_to_centre(
                                        unique(
                                                differ_pairs[i],
                                                differ_pairs[j]), c))
        triplets = [t for t in triplets if len(t) == 3]
                    
        for triplet in triplets:
            triplet_central_domains = [chains[chain][0][2]
            for chain in triplet]
            # ok
            allocated_numbers = [-1] * 3
            for i in range(3):
                t = triplet_central_domains[i]
                if len(t) > 1:
                    if t[1] != '*':
                        allocated_numbers[i] = int(t[1])
            new_domains = suitable_domains([a for a in allocated_numbers
                                            if a != -1])
            if len(new_domains):
                # not all three 
                if any([a == -1 for a in allocated_numbers]):
                    # something's empty, fill in with new domains
                    for d in new_domains:
                        try:
                            allocated_numbers[allocated_numbers.index(-1)] = d
                        except:
                            pass
        
                if allocated_numbers[1] == allocated_numbers[0] or \
                allocated_numbers[1] == allocated_numbers[2] and \
                allocated_numbers[1] != -1:
                    allocated_numbers[1] = new_domains[0]
        
            # then set the chains accordingly
            for i in range(3):
                for strand in chains[triplet[i]]:
                    if '*' in strand[2]:
                        strand[2] = 'c' + str(allocated_numbers[i]) + '*'
                    else:
                        strand[2] = 'c' + str(allocated_numbers[i])
                        
def make_domains(A, g, names, central_mismatch=False):
    species = create_species(A, names)
    enumerate_domains(species, g)
    if central_mismatch:
        chains, members = make_chains(species)
        allocate_central_domains(chains, members, species)
    return species
    
if __name__ == '__main__':      
    A = np.array([[0, 1, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1],
                  [0, 0, 0, 0]])
 
    names = ['K', 'X', 'Y', 'Z']
    species = create_species(A, names)
    g = igraph.Graph.Adjacency((A != 0).tolist())
    g.es['weight'] = A[A.nonzero()]   
    enumerate_domains(species, g)   
    chains, members = make_chains(species)
    allocate_central_domains(chains, members, species)

    for s in species:
        print('Species: ', s.name)
        print('State: ', s.state_strand)
        print('ID: ', s.id_strand[::-1])
        print('Active: ', s.active_state_strand)
        print('')


    
    
            
        