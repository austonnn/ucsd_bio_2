#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 13:31:09 2018

@author: xiangyin
"""

aminoAcid = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 
             'N', 'D', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']

aminoAcidMass = [ 57, 71, 87, 97, 99, 101, 103, 113,  
                 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

aminoAcid_dict = dict(zip(aminoAcid, aminoAcidMass))

aminoAcid_dict['Q'] = 128
aminoAcid_dict['L'] = 113



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

def Score(Peptide, Spectrum):
    ans_num = 0
    list_peptide = CycloSpectrum(Peptide)

    for item in list_peptide:
        #print(type(item))
        tmp_Spectrum = list(Spectrum)
        #print(tmp_Spectrum)
        if str(item) in Spectrum:
            ans_num += 1
            tmp_Spectrum.remove(str(item))
            Spectrum = list(tmp_Spectrum)
    
    return ans_num


def LinearScore(Peptide, Spectrum):
    ans_num = 0
    list_peptide = LinearSpectrum(Peptide)

    for item in list_peptide:
        #print(type(item))
        tmp_Spectrum = list(Spectrum)
        #print(tmp_Spectrum)
        if str(item) in Spectrum:
            ans_num += 1
            tmp_Spectrum.remove(str(item))
            Spectrum = list(tmp_Spectrum)
    
    return ans_num
    

with open('input.txt') as handle:
    Peptide, Spectrum = handle.read().splitlines()
Spectrum = Spectrum.split()
ans = LinearScore(Peptide, Spectrum)
print(ans)