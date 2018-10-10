#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 15:52:21 2018

@author: xiangyin
"""

# create codon table
with open('RNA_codon_table_1.txt') as handle:
    RNA_codon_table = handle.read().splitlines()

RNA_codon_dict = dict()

for item in RNA_codon_table:
    key = item[:3]
    value = item[-1]
    RNA_codon_dict[key] = value
    
def reverse_complement(text):
    tmp_dict = {'A':'T', 'T':'A', 'C':'G','G':'C'}
    tmp_list = list()
    for item in text:
        tmp_list += tmp_dict[item]
    tmp_list.reverse()
    tmp_str = ''
    for item in tmp_list:
        tmp_str += item
        
   
    return tmp_str

def dna_2_rna(text):
    tmp_dict = {'A':'A', 'T':'U', 'C':'C','G':'G'}
    tmp_str = ''
    for item in text:
        tmp_str += tmp_dict[item]
    return tmp_str

def rna_2_dna(text):
    tmp_dict = {'A':'A', 'U':'T', 'C':'C','G':'G'}
    tmp_str = ''
    for item in text:
        tmp_str += tmp_dict[item]
    return tmp_str

def get_dna_peptide(input_dna, peptide):
    
    output_dna_list  = []
    step = len(peptide) * 3
    for i in range(len(input_dna) // step):
        tmp_text = ''
        start = i * step
        stop = (i + 1) * step
        print(start)
        tmp_dna = input_dna[start : stop]

        #tmp_dna = 'AAGGAAGTATTTGAGCCTCATTATTAC'
        tmp_rna = dna_2_rna(tmp_dna)
        for j in range(len(peptide)):
            j_start = j * 3
            j_stop = (j + 1) * 3

            key = tmp_rna[j_start : j_stop]
            tmp_text += RNA_codon_dict[key]
        #print(i)    
        #print(tmp_rna)
        #print(tmp_text)
        if tmp_text == peptide:
            tmp_dna = rna_2_dna(tmp_rna)
            output_dna_list.append(tmp_dna)
            
    return output_dna_list

# read input data     
with open('input.txt') as handle:
    input_data = handle.read().splitlines()



input_dna = input_data[0]  
rc_input_dna = reverse_complement(input_dna)
peptide = input_data[1]



dna_list = []
rc_dna_list = []

for i in range(3):
    dna_list.append(input_dna[i:])
    rc_dna_list.append(rc_input_dna[i:])


#print(get_rna_peptide(rc_rna_list, peptide))

tmp_list = []
for item in dna_list:
    for dna_item in get_dna_peptide(item, peptide):
        #item = reverse_complement(item)
        tmp_list.append(dna_item)
#print(get_dna_peptide(dna_list[2], peptide))        
for item in rc_dna_list:
    for dna_item in get_dna_peptide(item, peptide):
        #print(item)
        dna_item = reverse_complement(dna_item)
        tmp_list.append(dna_item)
        #print(item)
print(tmp_list)        
#print(tmp_list)
text = ''       
for item in tmp_list:
    text += item
    text += '\n'
with open('output.txt', 'w') as handle:
    handle.write(text)

    
#input_rna = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'

"""
output_amino_acid  = []
for i in range(len(input_rna) // 3):
    start = i * 3
    stop = i * 3 + 3
    key = input_rna[start : stop]
    output_amino_acid.append(RNA_codon_dict[key])

output_text = ''
for item in output_amino_acid:
    output_text += item
print(output_text)

with open('output.txt', 'w') as handle:
    handle.write(output_text)
    """