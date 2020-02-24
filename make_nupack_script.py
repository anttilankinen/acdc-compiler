#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 18:14:05 2020

@author: antti
"""
from time import localtime
def design_script(species_list, central_mismatch, stop):
    # no mismatches here yet
    f = open('design.np', 'w')
    lt = localtime()
    f.write('# Created %d.%d.%d %d.%d.%d #\n' % (lt[2], lt[1], lt[0],
            lt[3], lt[4], lt[5]))
    f.write('material = dna\n')
    f.write('temperature[C] = 25.0\n')
    f.write('trials = 10\n')
    if central_mismatch:
        species_format1 = '.' * 5 + '(' * 10 + '.' + '(' * 22 + '.' * 5 + \
        '+' + '.' * 5 + ')' * 22 + '.' + ')' * 10 + '.' * 5
        species_format2 = '.' * 5 + '(' * 16 + '.' + '(' * 16 + '.' * 5 + \
        '+' + '.' * 5 + ')' * 16 + '.' + ')' * 16 + '.' * 5
        species_format3 = '.' * 5 + '(' * 22 + '.' + '(' * 10 + '.' * 5 + \
        '+' + '.' * 5 + ')' * 10 + '.' + ')' * 22 + '.' * 5
        # 1,2
        fuel_format1 = '.' * 10 + '(' * 5 + '.' + '(' * 5 + '.' + \
        '(' *  16 + '.' + '(' * 4 + '+' + ')' * 4 + '.' + ')' * 16 + \
        '.' + ')' * 5 + '.' + ')' * 5 + '.' * 10
        # 1,3
        fuel_format2 = '.' * 10 + '(' * 5 + '.' + '(' * 11 + '.' + \
        '(' *  10 + '.' + '(' * 4 + '+' + ')' * 4 + '.' + ')' * 10 + \
        '.' + ')' * 11 + '.' + ')' * 5 + '.' * 10
        # 2,3
        fuel_format3 = '.' * 10 + '(' * 11 + '.' + '(' * 5 + '.' + \
        '(' *  10 + '.' + '(' * 4 + '+' + ')' * 4 + '.' + ')' * 10 + \
        '.' + ')' * 5 + '.' + ')' * 11 + '.' * 10
        intermediate_product_format = \
        '.' * 10 + '(' * 33 + '+' + ')' * 33 + '.' * 10
    else:
        species_format = '.' * 5 + '(' * 25 + '.' * 5 + '+' +\
        '.' * 5 + ')' * 25 + '.' * 5
        fuel_format = '.' * 10 + '(' * 20 + '.' + '(' * 4 + '+' +\
        ')' * 4 + '.' + ')' * 20 + '.' * 10
        intermediate_product_format = \
        '.' * 10 + '(' * 25 + '+' + ')' * 25 + '.' * 10
        
            
    f.write('\n# domains #\n\n')
    domains = [d if not '*' in d else d[:-1] for s in species_list 
               for d in s.state_strand + s.id_strand + s.active_state_strand]
    domains = list(set(domains))
    domains.sort()
    for d in domains:
        if d == 'c':
            if central_mismatch:
                f.write('domain c = N5CN5CN5CN5\n')
                #f.write('domain c = N23\n')
            else:
                f.write('domain c = N15\n')
        elif d == 'c2':
            f.write('domain c2 = N5GN5CN5CN5\n')
            #f.write('domain c2 = N23\n')
        elif d == 'c3':
            f.write('domain c3 = N5CN5GN5CN5\n')
            #f.write('domain c3 = N23\n')
        elif d == 'c4':
            f.write('domain c4 = N5CN5CN5GN5\n')
            #f.write('domain c4 = N23\n')
        else:
            if False:#'2' in d and d[0] in domains:
               # force toehold C-C mismatch 
                for s in species_list:
                    if d in s.active_state_strand + s.state_strand:
                       # figure out which of the neighbours
                       # has the other toehold
                        for s2 in s.activate_u + s.repress_u + \
                        s.activate_d + s.repress_d:
                            if d[0] + '*' in s2.active_state_strand + \
                            s2.state_strand:
                               # we have our pair of species
                                break
                        # now we need to figure out what
                        # is their relation
                        if s2 in s.repress_u + s.activate_u:
                            # s2 is upstream of s
                            print(s.name, s2.name, d, 1)
                            f.write('domain %s = GN4\n' % d[0])
                            f.write('domain %s = CN4\n' % d)
                            break
                        else:
                            # s2 is downstream of s
                            print(s.name, s2.name, d, 2)
                            f.write('domain %s = N4G\n' % d[0])
                            f.write('domain %s = N4C\n' % d)
                            break
                    elif d + '*' in s.active_state_strand + s.state_strand:
                        for s2 in s.activate_u + s.repress_u + \
                        s.activate_d + s.repress_d:
                            if d[0] in s2.active_state_strand + \
                            s2.state_strand:
                                # we have our pair of species
                                break
                        # now we need to figure out what
                        # is their relation
                        if s2 in s.repress_u + s.activate_u:
                            # s2 is upstream of s
                            print(s.name, s2.name, d, 3)
                            f.write('domain %s = N4C\n' % d[0])
                            f.write('domain %s = N4G\n' % d)
                            break
                        else:
                            # s2 is downstream of s
                            print(s.name, s2.name, d, 4)
                            f.write('domain %s = CN4\n' % d[0])
                            f.write('domain %s = GN4\n' % d)
                            break
                    elif d[0] in s.active_state_strand + s.state_strand:
                        for s2 in s.activate_u + s.repress_u + \
                        s.activate_d + s.repress_d:
                            if d + '*' in s2.active_state_strand + \
                            s2.state_strand:
                                # we have our pair of species
                                break
                        # now we need to figure out what
                        # is their relation
                        if s2 in s.repress_u + s.activate_u:
                            # s2 is upstream of s
                            print(s.name, s2.name, d, 5)
                            f.write('domain %s = CN4\n' % d[0])
                            f.write('domain %s = GN4\n' % d)
                            break
                        else:
                            # s2 is downstream of s
                            print(s.name, s2.name, d, 6)
                            f.write('domain %s = N4C\n' % d[0])
                            f.write('domain %s = N4G\n' % d)
                            break
                    elif d[0] + '*' in s.active_state_strand + s.state_strand:
                        for s2 in s.activate_u + s.repress_u + \
                        s.activate_d + s.repress_d:
                            if d in s2.active_state_strand + \
                            s2.state_strand:
                                # we have our pair of species
                                break
                        # now we need to figure out what
                        # is their relation
                        if s2 in s.repress_u + s.activate_u:
                            # s2 is upstream of s
                            print(s.name, s2.name, d, 7)
                            f.write('domain %s = N4G\n' % d[0])
                            f.write('domain %s = N4C\n' % d)
                            break
                        else:
                            # s2 is downstream of s
                            print(s.name, s2.name, d, 8)
                            f.write('domain %s = GN4\n' % d[0])
                            f.write('domain %s = CN4\n' % d)
                            break
                           
                    
                # force toehold C-C mismatch
            elif '2' in d and d[0] in domains:
                for s in species_list:
                    if d[0] in s.active_state_strand[0]:
                        if s.active_state_strand[0] == d[0]:
                            # no invert, d[0] C
                            f.write('domain %s = N4C\n' % d[0])
                            f.write('domain %s = N4G\n' % d)
                            break
                        elif s.active_state_strand[0] == d:
                            # no invert, d[0] G
                            f.write('domain %s = N4G\n' % d[0])
                            f.write('domain %s = N4C\n' % d)
                            break
                        elif s.active_state_strand[0] == d[0] + '*':
                            # invert, d[0] G
                            f.write('domain %s = GN4\n' % d[0])
                            f.write('domain %s = CN4\n' % d)
                            break
                        elif s.active_state_strand[0] == d + '*':
                            f.write('domain %s = CN4\n' % d[0])
                            f.write('domain %s = GN4\n' % d)
                            break
                        break
                    elif d[0] in s.active_state_strand[4]:
                        if s.active_state_strand[4] == d[0]:
                            # no invert, d[0] C
                            f.write('domain %s = CN4\n' % d[0])
                            f.write('domain %s = GN4\n' % d)
                            break
                        elif s.active_state_strand[4] == d:
                            # no invert, d[0] G
                            f.write('domain %s = GN4\n' % d[0])
                            f.write('domain %s = CN4\n' % d)
                            break
                        elif s.active_state_strand[4] == d[0] + '*':
                            # invert, d[0] G
                            f.write('domain %s = N4G\n' % d[0])
                            f.write('domain %s = N4C\n' % d)
                            break
                        elif s.active_state_strand[4] == d + '*':
                            f.write('domain %s = N4C\n' % d[0])
                            f.write('domain %s = N4G\n' % d)
                            break
                        break
            elif not d + '2' in domains:
                f.write('domain %s = N5\n' % d)
               
               
                 
               
                
                       
                 
            #   f.write('similarity ' + d[0] + '3 = N5\n')
             #  f.write(d[0] + '.similarity[frac] '+ d[0] +
               #        '3 = [0.80, 0.80]\n')
              # f.write(d + '.similarity[frac] ' + d[0] +
                #       '3 = [1.00, 1.00]\n')

    
    f.write('\n# strands # \n\n')
    strand_counter = 1
    for s in species_list:
        # make the strands to structures string here
        if len(s.activate_u + s.repress_u):
            strand = list_to_string(s.state_strand)
            f.write('strand std%d = %s\n' % (strand_counter, strand))
            s.state_strand_name = 'std' + str(strand_counter)
            strand_counter += 1
        strand = list_to_string(s.id_strand)
        f.write('strand std%d = %s\n' % (strand_counter, strand))
        s.id_strand_name = 'std' + str(strand_counter)
        strand_counter += 1
        strand = list_to_string(s.active_state_strand)
        f.write('strand std%d = %s\n' % (strand_counter, strand))
        s.active_state_strand_name = 'std' + str(strand_counter)
        strand_counter += 1
        
        
        #active_strand_string = 's' + str(strand_counter - 1) + ' s' + \
        #str(strand_counter - 2)
        #if len(s.activate_u + s.repress_u):
        #    strand_string = 'std' + str(strand_counter - 3) + \
        #    ' std' + str(strand_counter - 2)
        #    strand_strings.append(strand_string)
        #strand_strings.append(active_strand_string) + str(strand_counter)
            
    f.write('\n# complexes #\n\n')

    for s in species_list:
        if len(s.activate_u + s.repress_u):
            f.write('complex %s = ' % (s.name) + s.state_strand_name +
                    ' ' + s.id_strand_name + '\n')
            f.write('complex %sact = ' % (s.name) +
                    s.active_state_strand_name + ' ' + s.id_strand_name + '\n')
        else:
            f.write('complex %s = ' % (s.name) + s.active_state_strand_name + 
                    ' ' + s.id_strand_name + '\n')
        for downstream_species in s.activate_d:
            f.write('complex WASTE_' + s.name + downstream_species.name * 2 +
                    'act = ' +
                    downstream_species.state_strand_name +
                    ' ' + s.active_state_strand_name + '\n')
            f.write('complex ' + s.name + downstream_species.name + ' = ' +
                    s.id_strand_name + ' ' +
                    downstream_species.id_strand_name + '\n')
            f.write('complex FUEL_' + s.name + downstream_species.name * 2 +
                    'act = ' +
                    downstream_species.active_state_strand_name +
                    ' ' + s.active_state_strand_name + '\n')
        for downstream_species in s.repress_d:
            f.write('complex WASTE_' + s.name + downstream_species.name +
                    'act' + downstream_species.name + ' = ' +
                    downstream_species.active_state_strand_name +
                    ' ' + s.active_state_strand_name + '\n')
            f.write('complex ' + s.name + downstream_species.name + ' = ' +
                    s.id_strand_name + ' ' +
                    downstream_species.id_strand_name + '\n')
            f.write('complex FUEL_' + s.name + downstream_species.name +
                    'act' + downstream_species.name + ' = ' +
                    downstream_species.state_strand_name + ' ' +
                    s.active_state_strand_name + '\n')
    
    f.write('\n# structures of complexes #\n\n')
    # species
    for s in species_list:
        if len(s.activate_u + s.repress_u):
            if central_mismatch:
                if s.state_strand[2] in ['c2', 'c4*']: 
                    f.write('%s.structure = %s\n' % (s.name, species_format1))
                elif s.state_strand[2] in ['c3', 'c3*']:
                    f.write('%s.structure = %s\n' % (s.name, species_format2))
                elif s.state_strand[2] in ['c4', 'c2*']:
                    f.write('%s.structure = %s\n' % (s.name, species_format3))
                else:
                    raise ValueError('Incorrect central domain')
                if s.active_state_strand[2] in ['c2', 'c4*']: 
                    f.write('%s.structure = %s\n' %
                            (s.name + 'act', species_format1))
                elif s.active_state_strand[2] in ['c3', 'c3*']:
                    f.write('%s.structure = %s\n' %
                            (s.name + 'act', species_format2))
                elif s.active_state_strand[2] in ['c4', 'c2*']:
                    f.write('%s.structure = %s\n' %
                            (s.name + 'act', species_format3))
                else:
                    raise ValueError('Incorrect central domain')
            else:
                f.write('%s.structure = %s\n' % (s.name, species_format))
                f.write('%s.structure = %s\n' % (s.name + 'act',
                                                 species_format))
        else:
            if central_mismatch:
                if s.active_state_strand[2] in ['c2', 'c4*']: 
                    f.write('%s.structure = %s\n' %
                            (s.name, species_format1))
                elif s.active_state_strand[2] in ['c3', 'c3*']:
                    f.write('%s.structure = %s\n' %
                            (s.name, species_format2))
                elif s.active_state_strand[2] in ['c4', 'c2*']:
                    f.write('%s.structure = %s\n' %
                            (s.name, species_format3))
                else:
                    raise ValueError('Incorrect central domain')
            else:
                f.write('structure %s = %s\n' % (s.name, species_format))
    # fuels, intermediates, wastes  
    for s in species_list:
        for downstream_species in s.activate_d:
            f.write('%s.structure = %s\n' %
                    (s.name + downstream_species.name,
                     intermediate_product_format))
            
            # catalyst state strand is on the bottom
            if central_mismatch:
                if downstream_species.active_state_strand[2] == 'c2' \
                and s.active_state_strand[2] == 'c3*' \
                or downstream_species.active_state_strand[2] == 'c3' \
                and s.active_state_strand[2] == 'c2*'\
                or downstream_species.active_state_strand[2] == 'c3*' \
                and s.active_state_strand[2] == 'c4' \
                or downstream_species.active_state_strand[2] == 'c4*' \
                and s.active_state_strand[2] == 'c3':
                    f.write('%s.structure = %s\n' %
                        ('FUEL_' + s.name + downstream_species.name * 2 +
                         'act', fuel_format1))
                elif downstream_species.active_state_strand[2][0:2] == 'c4' \
                and s.active_state_strand[2][0:2] == 'c2' \
                or downstream_species.active_state_strand[2][0:2] == 'c2' \
                and s.active_state_strand[2][0:2] == 'c4':
                    f.write('%s.structure = %s\n' %
                        ('FUEL_' + s.name + downstream_species.name * 2 +
                         'act', fuel_format2))
                elif downstream_species.active_state_strand[2] == 'c3*' \
                and s.active_state_strand[2] == 'c2' \
                or downstream_species.active_state_strand[2] == 'c2*' \
                and s.active_state_strand[2] == 'c3'\
                or downstream_species.active_state_strand[2] == 'c3' \
                and s.active_state_strand[2] == 'c4*' \
                or downstream_species.active_state_strand[2] == 'c4' \
                and s.active_state_strand[2] == 'c3*':
                    f.write('%s.structure = %s\n' %
                        ('FUEL_' + s.name + downstream_species.name * 2 +
                         'act', fuel_format3))
            else:
                f.write('%s.structure = %s\n' %
                        ('FUEL_' + s.name + downstream_species.name * 2 +
                         'act', fuel_format))
            f.write('%s.structure = %s\n' %
                    ('WASTE_' + s.name + downstream_species.name * 2 +
                     'act',
                     intermediate_product_format))
        for downstream_species in s.repress_d:
            f.write('%s.structure = %s\n' %
                    (s.name + downstream_species.name,
                     intermediate_product_format))
            if central_mismatch:
                if downstream_species.state_strand[2] == 'c2' \
                and s.active_state_strand[2] == 'c3*' \
                or downstream_species.state_strand[2] == 'c3' \
                and s.active_state_strand[2] == 'c2*'\
                or downstream_species.state_strand[2] == 'c3*' \
                and s.active_state_strand[2] == 'c4' \
                or downstream_species.state_strand[2] == 'c4*' \
                and s.active_state_strand[2] == 'c3':
                    f.write('%s.structure = %s\n' %
                        ('FUEL_' + s.name + downstream_species.name + 'act' +
                         downstream_species.name, fuel_format1))
                elif downstream_species.state_strand[2][0:2] == 'c4' \
                and s.active_state_strand[2][0:2] == 'c2' \
                or downstream_species.state_strand[2][0:2] == 'c2' \
                and s.active_state_strand[2][0:2] == 'c4':
                    f.write('%s.structure = %s\n' %
                        ('FUEL_' + s.name + downstream_species.name + 'act' +
                         downstream_species.name, fuel_format2))
                elif downstream_species.state_strand[2] == 'c3*' \
                and s.active_state_strand[2] == 'c2' \
                or downstream_species.state_strand[2] == 'c2*' \
                and s.active_state_strand[2] == 'c3'\
                or downstream_species.state_strand[2] == 'c3' \
                and s.active_state_strand[2] == 'c4*' \
                or downstream_species.state_strand[2] == 'c4' \
                and s.active_state_strand[2] == 'c3*':
                    f.write('%s.structure = %s\n' %
                        ('FUEL_' + s.name + downstream_species.name + 'act' +
                         downstream_species.name, fuel_format3))
            else:
                f.write('%s.structure = %s\n' %
                        ('FUEL_' + s.name + downstream_species.name + 'act' +
                         downstream_species.name, fuel_format))
            f.write('%s.structure = %s\n' %
                    ('WASTE_' + s.name + downstream_species.name + 'act' +
                     downstream_species.name,
                     intermediate_product_format))
            
            
    f.write('\n# prevent sequence patterns #\n\n')
    f.write('prevent = AAAA, CCCC, GGGG, UUUU, MMMMMMM, KKKKKK, ' +
            'WWWWWW, SSSSSS, RRRRRR, YYYYYY\n')
    f.write('\n# max defect #\n\n')
    f.write('stop = %.2f\n' % (stop / 100))
    f.close()

            
def list_to_string(l):
    out = ''
    for item in l:
        out += item + ' '
    return out

def defect_script(species_list, central_mismatch=False):
    # perform mismatches and re-analyse system
    f = open('design.np', 'r')
    np = f.readlines()
    f.close()
    f = open('design_0.npo', 'r')
    npo = f.readlines()
    f.close()
    
    f = open('defect.np', 'w')
    
    # find where we start the domain sequences in npo file
    domains = [d if not '*' in d else d[:-1] for s in species_list 
               for d in s.state_strand + s.id_strand + s.active_state_strand]
    domains_start = ['domains:' in line for line in npo].index(True)
    domain_line = domains_start + 1
    for line in np:
        if not 'domain ' in line:
            f.write(line)
        else:
            # pick sequence from the 
            d = npo[domain_line].split()[0]
            sequence = npo[domain_line].split()[2]
            if d == 'c' and central_mismatch:
                sequence[5] = 'C'
                sequence[11] = 'C'
                sequence[17] = 'C'
            if d == 'c2' and central_mismatch:
                sequence[5] = 'G'
                sequence[11] = 'C'
                sequence[17] = 'C'
            if d == 'c3' and central_mismatch:
                sequence[5] = 'C'
                sequence[11] = 'G'
                sequence[17] = 'C'
            if d == 'c4' and central_mismatch:
                sequence[5] = 'C'
                sequence[11] = 'C'
                sequence[17] = 'G'
            if '2' in d and d[0] in domains:
                for s in species_list:
                    if d in s.active_state_strand + s.state_strand:
                        # figure out which of the neighbours
                        # has the other toehold
                        for s2 in s.activate_u + s.repress_u + \
                        s.activate_d + s.repress_d:
                            if d[0] + '*' in s2.active_state_strand + \
                            s2.state_strand:
                                # we have our pair of species
                                break
                        # now we need to figure out what
                        # is their relation
                        if s2 in s.repress_u + s.activate_u:
                            # s2 is upstream of s
                            sequence = 'C' + sequence[1:]
                            break
                        else:
                            # s2 is downstream of s
                            sequence = sequence[:4] + 'C'
                            break
                    elif d + '*' in s.active_state_strand + s.state_strand:
                        for s2 in s.activate_u + s.repress_u + \
                        s.activate_d + s.repress_d:
                            if d[0] in s2.active_state_strand + \
                            s2.state_strand:
                                # we have our pair of species
                                break
                        # now we need to figure out what
                        # is their relation
                        if s2 in s.repress_u + s.activate_u:
                            # s2 is upstream of s
                            sequence = sequence[:4] + 'G'
                            break
                        else:
                            # s2 is downstream of s
                            sequence = 'G' + sequence[1:]
                            break
                    elif d[0] in s.active_state_strand + s.state_strand:
                        for s2 in s.activate_u + s.repress_u + \
                        s.activate_d + s.repress_d:
                            if d + '*' in s2.active_state_strand + \
                            s2.state_strand:
                                # we have our pair of species
                                break
                        # now we need to figure out what
                        # is their relation
                        if s2 in s.repress_u + s.activate_u:
                            # s2 is upstream of s
                            sequence = 'G' + sequence[1:]
                            break
                        else:
                            # s2 is downstream of s
                            sequence = sequence[:4] + 'G'
                            break
                    elif d[0] + '*' in s.active_state_strand + s.state_strand:
                        for s2 in s.activate_u + s.repress_u + \
                        s.activate_d + s.repress_d:
                            if d in s2.active_state_strand + \
                            s2.state_strand:
                                # we have our pair of species
                                break
                        # now we need to figure out what
                        # is their relation
                        if s2 in s.repress_u + s.activate_u:
                            # s2 is upstream of s
                            sequence = sequence[:4] + 'C'
                            break
                        else:
                            # s2 is downstream of s
                            sequence = 'C' + sequence[1:]
                            break
            if d + '2' in domains:
                for s in species_list:
                    if d + '2' in s.active_state_strand + s.state_strand:
                        # figure out which of the neighbours
                        # has the other toehold
                        for s2 in s.activate_u + s.repress_u + \
                        s.activate_d + s.repress_d:
                            if d + '*' in s2.active_state_strand + \
                            s2.state_strand:
                                # we have our pair of species
                                break
                        # now we need to figure out what
                        # is their relation
                        if s2 in s.repress_u + s.activate_u:
                            # s2 is upstream of s
                            sequence = 'G' + sequence[1:]
                            break
                        else:
                            # s2 is downstream of s
                            sequence = sequence[:4] + 'G'
                            break
                    elif d + '2*' in s.active_state_strand + s.state_strand:
                        for s2 in s.activate_u + s.repress_u + \
                        s.activate_d + s.repress_d:
                            if d in s2.active_state_strand + \
                            s2.state_strand:
                                # we have our pair of species
                                break
                        # now we need to figure out what
                        # is their relation
                        if s2 in s.repress_u + s.activate_u:
                            # s2 is upstream of s
                            sequence = sequence[:4] + 'C'
                            break
                        else:
                            # s2 is downstream of s
                            sequence = 'C' + sequence[1:]
                            break
                    elif d in s.active_state_strand + s.state_strand:
                        for s2 in s.activate_u + s.repress_u + \
                        s.activate_d + s.repress_d:
                            if d + '2*' in s2.active_state_strand + \
                            s2.state_strand:
                                # we have our pair of species
                                break
                        # now we need to figure out what
                        # is their relation
                        if s2 in s.repress_u + s.activate_u:
                            # s2 is upstream of s
                            sequence = 'C' + sequence[1:]
                            break
                        else:
                            # s2 is downstream of s
                            sequence = sequence[:4] + 'C'
                            break
                    elif d + '*' in s.active_state_strand + s.state_strand:
                        for s2 in s.activate_u + s.repress_u + \
                        s.activate_d + s.repress_d:
                            if d + '2' in s2.active_state_strand + \
                            s2.state_strand:
                                # we have our pair of species
                                break
                        # now we need to figure out what
                        # is their relation
                        if s2 in s.repress_u + s.activate_u:
                            # s2 is upstream of s
                            sequence = sequence[:4] + 'G'
                            break
                        else:
                            # s2 is downstream of s
                            sequence = 'G' + sequence[1:]
                            break
            f.write('domain %s = %s\n' % (d, sequence))
            domain_line += 2 # npo file defines complements as well
    f.close()
    
    
    
    
        
        
        
        
        
        
        