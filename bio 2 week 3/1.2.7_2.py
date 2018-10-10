#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 23:23:17 2018

@author: xiangyin
"""




with open('RNA_codon_table_1.txt') as handle:
    RNA_codon_table = handle.read().splitlines()

RNA_codon_dict = dict()

for item in RNA_codon_table:
    key = item[:3]
    value = item[-1]
    RNA_codon_dict[key] = value
    
rv_RNA_codon_dict = dict()

for item in RNA_codon_table:
    val = item[:3]
    key = item[-1]
    if key not in rv_RNA_codon_dict.keys():
        rv_RNA_codon_dict[key] = []       
    rv_RNA_codon_dict[key].append(val)
    
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

def check_match(tmp_dna, peptide, hpattern):
    tmp_rna = dna_2_rna(tmp_dna)
    tmp_text = ''
    for j in range(len(peptide)):
            j_start = j * 3
            j_stop = (j + 1) * 3

            key = tmp_rna[j_start : j_stop]
            tmp_text += RNA_codon_dict[key]
    if hash(tmp_text) == hpattern:
        if tmp_text == peptide:
            return True
    else:
        return False


    
    
def get_dna_patterns(peptide):
    rna_patterns = ['']
    for item in peptide:
        tmp_list = []
        for rv_RNA in rv_RNA_codon_dict[item]:
            for rna_seqs in rna_patterns:
                tmp_list.append(rna_seqs + rv_RNA)
        rna_patterns = tmp_list
    dna_patterns = []
    for item in rna_patterns:
        dna_patterns.append(rna_2_dna(item))
        dna_patterns.append(reverse_complement(rna_2_dna(item)))
        
    return set(dna_patterns)

#print(get_dna_patterns('N'))    

def get_dna_peptide(input_dna, peptide):
    # 虽然我不知道微生物怎么做的，但是通过读取dna片段，
    # 然后检查它的互补片段是否满足要求是比较简洁的思路。
    # 比生成dna片段，它的互补片段，检查完了再转化回去要清晰。
    # 更多的优化包括，把整个序列转化好，然后检查位置，然后用位置信息去调取原来的片段。
    # 因为这里有大量的重复计算，在反复把序列转化为蛋白序列的过程中。
    dna_patterns = get_dna_patterns(peptide)
    output_dna_list  = []
    step = len(peptide) * 3
    for i in range((len(input_dna) - step) // 3 + 1):
        start = i * 3
        stop = i * 3 + step
        #print(start)
        tmp_dna = input_dna[start : stop]
        #print(tmp_dna)

        #tmp_dna = 'AAGGAAGTATTTGAGCCTCATTATTAC'
        if tmp_dna in dna_patterns:
            output_dna_list.append(tmp_dna)


            
    return output_dna_list


# read input data     
with open('input.txt') as handle:
    input_data = handle.read().splitlines()
    
#peptide = 'KEVFEPHYY'
#for item in input_data:
#    print('#')
#    print(check_match(item, peptide))
#    print(check_match(reverse_complement(item), peptide))

input_dna = input_data[0]  
peptide = input_data[1]


with open('Bacillus_brevis.txt') as handle:
    input_dna = handle.read().replace('\n', '')
peptide = 'VKLFPWFNQY'
peptide = 'QIQVLEG'
dna_list = []
tmp_list = []

for i in range(3):
    dna_list.append(input_dna[i:])
#dna_list = [input_dna[:]]
for item in dna_list:
    for dna_item in get_dna_peptide(item, peptide):
        #item = reverse_complement(item)
        tmp_list.append(dna_item)
    
text = ''       
for item in tmp_list:
    text += item
    text += '\n'
with open('output.txt', 'w') as handle:
    handle.write(text)

    