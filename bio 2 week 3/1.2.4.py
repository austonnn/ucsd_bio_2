#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 14:37:13 2018

@author: xiangyin
"""

with open('RNA_codon_table_1.txt') as handle:
    RNA_codon_table = handle.read().splitlines()

RNA_codon_dict = dict()

for item in RNA_codon_table:
    key = item[:3]
    value = item[-1]
    RNA_codon_dict[key] = value
    
    
with open('dataset_96_4.txt') as handle:
    input_rna = handle.read()
  
#input_rna = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'

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
    