#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 10:16:13 2018

@author: xiangyin
"""

# A Branch-and-Bound Algorithm for Cyclopeptide Sequencing



aminoAcid = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 
             'N', 'D', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']

aminoAcidMass = [ 57, 71, 87, 97, 99, 101, 103, 113,  
                 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

aminoAcid_dict = dict(zip(aminoAcid, aminoAcidMass))


def Mass(Peptide):
    tmp_num = 0
    for item in Peptide:
        tmp_num += aminoAcid_dict[item]
    return tmp_num
    
    
    
def CycloSpectrum(Peptide):
    PrefixMass = [0]
    for i in range(len(Peptide)):
        tmp_num = PrefixMass[i] + aminoAcid_dict[Peptide[i]]
        PrefixMass.append(tmp_num)
    peptideMass = PrefixMass[len(Peptide)]
    CyclicSpectrum = [0]
    for i in range(len(Peptide)):
        for j in range(i + 1, len(Peptide) + 1):
            tmp_mass = PrefixMass[j] - PrefixMass[i]
            CyclicSpectrum.append(tmp_mass)
            if i > 0 and j < len(Peptide):
                CyclicSpectrum.append(peptideMass - tmp_mass)
    CyclicSpectrum.sort()
    return CyclicSpectrum

def LinearSpectrum(Peptide):
    PrefixMass = [0]
    for i in range(len(Peptide)):
        tmp_num = PrefixMass[i] + aminoAcid_dict[Peptide[i]]
        PrefixMass.append(tmp_num)
    LinearSpectrum = [0]
    for i in range(len(Peptide)):
        for j in range(i + 1, len(Peptide) + 1):
            LinearSpectrum.append(PrefixMass[j] - PrefixMass[i])
    LinearSpectrum.sort()
    return LinearSpectrum

def Expand(Peptides):
    tmp_list = []
    for peptide in Peptides:
        for item in aminoAcid:
            #tmp_text
            tmp_text = peptide + item
            #peptide += item
            tmp_list.append(tmp_text)    
    return tmp_list

def Expand_word(Peptides, word):
    tmp_list = []
    for peptide in Peptides:
        for item in word:
            #tmp_text
            tmp_text = peptide + item
            #peptide += item
            tmp_list.append(tmp_text)    
    return tmp_list


def linear_check_consistent(Peptide, Spectrum):
    #print('#')
    #print(Peptide)
    tmp_list = LinearSpectrum(Peptide)
    #print(tmp_list)
    for item in tmp_list:
        if str(item) not in Spectrum:
            #print(item)
            #print('a')
            return False       
    return True



def cyclo_check_consistent(Peptide, Spectrum):
    #print('#')
    #print(Peptide)
    tmp_list = CycloSpectrum(Peptide)
    #print(tmp_list)
    for item in tmp_list:
        if str(item) not in Spectrum:
            #print(item)
            #print('a')
            return False       
    return True

def CyclopeptideSequencing(Spectrum):
    # Peptides ← a set containing only the empty peptide
    Peptides = ['']

    ans_list = []
    word = []
    while Peptides:
        if len(word) == 0:
            # 做一个小的优化，只把可能出现的基础字母放进去后续的循环中。
            Peptides = Expand(Peptides)
            tmp_Peptides = list(Peptides)
            for Peptide in Peptides:
                if not linear_check_consistent(Peptide, Spectrum):
                    tmp_Peptides.remove(Peptide)
            Peptides = tmp_Peptides
            word = Peptides
            #print(word)
            
        Peptides = Expand_word(Peptides, word)
        #print(Peptides)
        #print(Peptides)
        tmp_Peptides = list(Peptides)
        for Peptide in Peptides:
            #print(Peptide)
            #print(Mass(Peptide))
            #print(Spectrum[-1])
            if Mass(Peptide) == int(Spectrum[-1]):

                # 完整比较用 cyclo_check
                if cyclo_check_consistent(Peptide, Spectrum):
                    #print(Peptide)
                    ans_list.append(Peptide)
                tmp_Peptides.remove(Peptide)
            # 一般性的比较用 linear_check
            elif not linear_check_consistent(Peptide, Spectrum):
                tmp_Peptides.remove(Peptide)
        Peptides = tmp_Peptides
        #print(Peptides)
        
    return ans_list
                 
                 
with open('dataset_100_6.txt') as handle:
    Spectrum = handle.read().split()
#Peptides = set([''])
#Spectrum = [0, 113]
ans = CyclopeptideSequencing(Spectrum)
#Peptide = "NCGND"
#ans = not check_consistent(Peptide, Spectrum)
text = ''
for item in ans:
    for i in item:
        text += str(aminoAcid_dict[i])
        text += '-'
    text = text[:-1]
    text += ' '
    
with open('output.txt', 'w') as handle:
    handle.write(text)
    
#print(text)
#Spectrum = [0, 113, 128, 186, 241, 299, 314, 427]
#CyclopeptideSequencing(Spectrum)