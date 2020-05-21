#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 14:05:57 2020

@author: antti
"""
import string
from itertools import combinations
import numpy as np

class Species:
    """
    class describing an ACDC species
    has the following attributes:
        name - str
        id_strand - list
        state_strand - list
        active_state_strand - list
        activate_d - list, set of species which this species activates
        repress_d - list, set of species which this species represses
        activate_u - list, set of species which activate this species
        repress_u - list, set of species which repress this species
    """
    
    def __init__(self, name):
        self.name = name

    def initialise_strands(self):
        # all strands are initially the same with no assigned domains
        self.active_state_strand = ['0'] * 5
        self.state_strand = ['0'] * 5
        self.id_strand = ['0'] * 5
    
    def fill_domains(self, domain_counter):
        """
        fill currently empty domains with new domains
        these domains do not require any complementarity with currently
        existing domains in other species
        args:
            self
            domain_counter - int, incremental labelling for new domains
        returns:
            domain_counter - int
        """

        if domain_counter == 0:
            id_compl = True # id strand's domains will contain '*'
        else:
            id_compl = any(['*' in d for d in self.id_strand])
            
            
        if len(self.activate_u) == 0 and len(self.repress_u) == 0:
            # no upstream species, so passive and active state strand
            # of this species aren't required in this network
            self.state_strand = self.active_state_strand
    
        # fill domains on the two state strands
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
            # downstream interface
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
        # id strand

        self.id_strand[1] = complement(self.state_strand[3]) 
        self.id_strand[2] = 'c*' if id_compl else 'c'
        # mismatch on inner toehold
        self.id_strand[3] = inner_toehold_mismatch_complement(
                self.state_strand[1])
        
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
        """
        check which of the neighbours have some domains already and fill
        own domains accordingly 
        args:
            self
        """
        # downstream neighbours
        for s in self.activate_d:
            # from this set, have to inherit their state strand
            # downstream interface to self's active state strand upstream
            # interface
            if not any([d == '0' for d in s.state_strand[4:1:-1]]):
                self.active_state_strand[0:3] = \
                complement(s.state_strand[4:1:-1])
                self.id_strand[4] = complement(s.id_strand[0])
                break
        for s in self.repress_d:
            # from this set, have to inherit their active state strand
            # downstream interface to self's active state strand upstream
            # interface
            if not any([d == '0' for d in s.active_state_strand[4:1:-1]]):
                self.active_state_strand[0:3] = \
                complement(s.active_state_strand[4:1:-1])
                self.id_strand[4] = complement(s.id_strand[0])
                break

        # then upstream neighbours
        for s in self.activate_u:
            # from this set, have to inherit their active state strand
            # upstream interface to self's state strand downstream interface
            if not any([d == '0' for d in s.active_state_strand[0:3]]):
                self.state_strand[4:1:-1] = \
                complement(s.active_state_strand[0:3])
                self.id_strand[0] = complement(s.id_strand[4])
                break
        for s in self.repress_u:
            # from this set, have to inherit their active state strand
            # upstream interface to self's active state strand downstream
            # interface
            if not any([d == '0' for d in s.active_state_strand[0:3]]):
                self.active_state_strand[4:1:-1] = \
                complement(s.active_state_strand[0:3])
                self.id_strand[0] = complement(s.id_strand[4])
                break

        
        
        
        
def label(n, compl=False):
    """
    transform a number to a base-26 representation and compute corresponding
    alphabets
    args:
        n - int
        compl - bool
    returns:
        str
    """
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
    """
    return the complement of a domain
    args:
        d - str or list
    returns:
        str or list
    """
    if type(d) == str:
        if d == '0':
            return d
            
        if '*' in d:
            return d[:-1]
        else:
            return d + '*'
    elif type(d) == list:
        return [complement(c) for c in d]
    
def inner_toehold_mismatch_complement(d):
    """
    return the complementary variant of an upstream interface inner toehold
    args:
        d - str
    returns:
        str
    """
    if '*' in d:
        return d[:-1] + '2'
    else:
        return d + '2*'
    

def create_species(A, names):
    """
    create a set of Species objects and compute their relations in the ACDC
    network
    args:
        A - numpy NDArray
        names - list
    returns:
        species_list - list
    """
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
    """
    assign domains to each species given the whole network structure
    args:
        species_list - list
        g - iGraph graph
    """
    
    domain_counter = 0 # keep track of which domain names have been used
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
    # unique list made from items in two lists
    unique_list = []
    for item in l1:
        if not item in unique_list:
            unique_list.append(item)
    for item in l2:
        if not item in unique_list:
            unique_list.append(item)
    return unique_list
    

def make_chains(species_list):
    """
    decompose the set of strands in the whole network to chains such that the
    within one chain, the central domain variants are the same for all strands
    args:
        species_list - list
    returns:
        chains - list
        members - list
    """
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

        # each strand can only belong to one chain, so merge chains
        # correspondingly
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
    # move item i to centre of list l of length 3
    l[l.index(i)] = l[1]
    l[1] = i
    return l

def suitable_domains(l):
    # return all numbers 2,3,4 that aren't on list l
    return [i for i in [2, 3, 4] if i not in l]


def make_pairs(differ_sets):
    """
    make pairs of chains which have to have different central domain variants
    args:
        differ_sets - list
    returns:
        differ_pairs - list
    """
    differ_pairs = []
    for s in differ_sets:
        pairs = combinations(s, 2)
        for p in pairs:
            differ_pairs.append(p)
    return differ_pairs

def allocate_central_domains(chains, members, species):
    """
    assign central domain variants to each strand
    args:
        chains - list
        members - list
        species - list
    """
    # make pairs of chains that have to have different central domain variants
    differ_sets= [[i for i in range(len(members)) if j in members[i]]
    for j in range(len(species))]
    differ_pairs = make_pairs(differ_sets)
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
        # merge these pairs into triplets ordered triplets
        # such that the central member of the triplet has to have different
        # central domain variants with both of the outer members.
        # however, the two outer members do not necessarily have to have
        # different variants
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
        # assign variants to each triplet
        for triplet in triplets:
            triplet_central_domains = [chains[chain][0][2]
            for chain in triplet]
            
            allocated_numbers = [-1] * 3
            for i in range(3):
                t = triplet_central_domains[i]
                if len(t) > 1:
                    if t[1] != '*':
                        allocated_numbers[i] = int(t[1])
            new_domains = suitable_domains([a for a in allocated_numbers
                                            if a != -1])
            if len(new_domains):
                # not all three have been assigned a domain
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
        
            # then set the variants accordingly
            for i in range(3):
                for strand in chains[triplet[i]]:
                    if '*' in strand[2]:
                        strand[2] = 'c' + str(allocated_numbers[i]) + '*'
                    else:
                        strand[2] = 'c' + str(allocated_numbers[i])
                        
                        
def make_domains(A, g, names, central_mismatch=False):
    # put everything together
    species = create_species(A, names)
    enumerate_domains(species, g)
    """
    in testing
    if central_mismatch:
        # assign central domain variants
        chains, members = make_chains(species)
        allocate_central_domains(chains, members, species)
        """

    return species

    
    
            
        