#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 15:52:31 2018

@author: xiangyin
"""
aminoAcid = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 
             'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']

aminoAcidMass = [ 57, 71, 87, 97, 99, 101, 103, 113, 113, 
                 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]


aminoAcid_dict = dict(zip(aminoAcid, aminoAcidMass))
# 1.7 CS: Generating the Theoretical Spectrum of a Peptide


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


def CyclicSpectrum(Peptide):
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

with open('dataset_98_4.txt') as handle:
    Peptide = handle.read().strip()
#Peptide = 'NQEL'
ans = CyclicSpectrum(Peptide)
text = ''
for item in ans:
    text += str(item)
    text += ' ' 

with open('output.txt', 'w') as handle:
    handle.write(text)