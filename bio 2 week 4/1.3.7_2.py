#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 22 11:17:57 2018

@author: xiangyin
"""

"""
aminoAcid = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 
             'N', 'D', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']

aminoAcidMass = [ 57, 71, 87, 97, 99, 101, 103, 113,  
                 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

aminoAcid_dict = dict(zip(aminoAcid, aminoAcidMass))

aminoAcid_dict['Q'] = 128
aminoAcid_dict['L'] = 113
"""

def extended_mass_table():
    aminoAcid_dict = {}
    for i in range(57,201):
        aminoAcid_dict[chr(i)] = i
    return aminoAcid_dict

aminoAcid_dict = extended_mass_table()

mass_dict = {'':0}
score_dict = {'':0}

def Mass(Peptide):
    if Peptide in mass_dict.keys():
        return mass_dict[Peptide]
    tmp_num = 0
    for item in Peptide:
        tmp_num += aminoAcid_dict[item]
    #tmp_num = mass_dict[Peptide[:-1]] + aminoAcid_dict[Peptide[-1]]
    mass_dict[Peptide] = tmp_num
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

def Score(Peptide, Spectrum):
    if Peptide in score_dict.keys():
        return score_dict[Peptide]
    ans_num = 0
    list_peptide = CycloSpectrum(Peptide)

    for item in list_peptide:
        #print(type(item))
        tmp_Spectrum = list(Spectrum)
        #print(tmp_Spectrum)
        if str(item) in tmp_Spectrum:
            ans_num += 1
            tmp_Spectrum.remove(str(item))
            Spectrum = list(tmp_Spectrum)
    score_dict[Peptide] = ans_num
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
"""    
def Expand(Leaderboard):
    tmp_list = []
    for Peptides in Leaderboard:
        tmp_list.extend(ExpandPeptides(Peptides))    
    return tmp_list
"""
def Expand(Leaderboard):
    tmp_list = []
    for Peptide in Leaderboard:
        for item in aminoAcid_dict.keys():
            tmp_Peptide  = Peptide + item
#            if tmp_Peptide not in Peptide_dict.keys():
#                tmp_score = Score(tmp_Peptide, Spectrum)
#                tmp_mass = Mass(tmp_Peptide)
#                Peptide_dict[tmp_Peptide] = [tmp_score, tmp_mass]

            tmp_list.append(tmp_Peptide) 
    return tmp_list

def Trim(Leaderboard, Spectrum, N):
    """
    Trim(Leaderboard, Spectrum, N, AminoAcid, AminoAcidMass)
    for j ← 1 to |Leaderboard|
        Peptide ← j-th peptide in Leaderboard
        LinearScores(j) ← LinearScore(Peptide, Spectrum)
    sort Leaderboard according to the decreasing order of scores in LinearScores
    sort LinearScores in decreasing order
    for j ← N + 1 to |Leaderboard|
        if LinearScores(j) < LinearScores(N)
            remove all peptides starting from the j-th peptide from Leaderboard
            return Leaderboard
    return Leaderboard
    """
    LinearScores_list = []
    ans_Leaderboard = []
    for j in range(len(Leaderboard)):
        Peptide = Leaderboard[j]
        LinearScores_list.append(LinearScore(Peptide, Spectrum))
        #LinearScores_list.append(Peptide_dict[Peptide][0])
#
    Leaderboard_list = list(zip(Leaderboard, LinearScores_list))
    Leaderboard_list = sorted(Leaderboard_list, key = lambda item:item[1], 
                      reverse = True)
    if N < len(Leaderboard):
        for j in range(N, len(Leaderboard)):
            if Leaderboard_list[j][1] < Leaderboard_list[N][1]:
                for item in Leaderboard_list[0 : j - 1]:
                    ans_Leaderboard.append(item[0])
                print('#')
                print(len(ans_Leaderboard))
                return ans_Leaderboard
                #return tmp_Leaderboard[0 : j - 1]
    for item in Leaderboard_list:
        ans_Leaderboard.append(item[0])
    print(len(ans_Leaderboard))
    return ans_Leaderboard



def LeaderboardCyclopeptideSequencing(Spectrum, N):
    #Leaderboard ← set containing only the empty peptide
    Leaderboard = ['']
    #Leaderboard = [['', 0, 0]]
    #LeaderPeptide ← empty peptide
    LeaderPeptide = ''
    LeaderPeptide_list = []
    #LeaderPeptide = ['', 0, 0]
    #while Leaderboard is non-empty
    while Leaderboard:
        #Leaderboard ← Expand(Leaderboard)
        Leaderboard = Expand(Leaderboard)
        #print(Leaderboard)
        #print(Leaderboard)
        #for each Peptide in Leaderboard
        tmp_Leaderboard = list(Leaderboard)
        for Peptide in Leaderboard:
            #if Mass(Peptide) = ParentMass(Spectrum)
            if Mass(Peptide) == int(Spectrum[-1]):
            #if Peptide_dict[Peptide][1] == int(Spectrum[-1]):
                
                if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum):
                #if Peptide_dict[Peptide][0] > Peptide_dict[LeaderPeptide][0]:
                    LeaderPeptide_list = []
                    LeaderPeptide_list.append(Peptide)
                    LeaderPeptide = Peptide
                elif Score(Peptide, Spectrum) == Score(LeaderPeptide, Spectrum):
                #if Peptide_dict[Peptide][0] > Peptide_dict[LeaderPeptide][0]:
                    #LeaderPeptide_list = []
                    LeaderPeptide_list.append(Peptide)
                    
            #else if Mass(Peptide) > ParentMass(Spectrum)
            elif Mass(Peptide) > int(Spectrum[-1]):
            #elif Peptide_dict[Peptide][1]  > int(Spectrum[-1]):
                #remove Peptide from Leaderboard
                tmp_Leaderboard.remove(Peptide)
        #Leaderboard ← Trim(Leaderboard, Spectrum, N)
        #print(Leaderboard)
        Leaderboard = tmp_Leaderboard
        Leaderboard = Trim(Leaderboard, Spectrum, N)
        #print(LeaderPeptide)
        #print(Leaderboard)

    #output LeaderPeptide
    return LeaderPeptide_list


def SpectralConvolution(Spectrum):
    ans_list = []
    tmp_list = [0]
    tmp_list.extend(Spectrum)
    for i in tmp_list:
        for j in Spectrum:
            if 57 <= (i - j) <= 200:
                ans_list.append(i - j)
            pass
    ans_list = sorted(ans_list)
    return ans_list

import collections

# imput data
with open('input.txt') as handle:
#with open('input.txt') as handle:
    M, N, Spectrum = handle.read().splitlines()
#A cyclic peptide LeaderPeptide with amino acids taken only 
#from the top M elements (and ties) of the convolution of
#Spectrum that fall between 57 and 200,
M = int(M)
#  the size of Leaderboard is restricted to the top N (and ties).
N = int(N)   

Spectrum = Spectrum.split()
tmp_list = []
for item in Spectrum:
    tmp_list.append(int(item))
Spectrum = sorted(tmp_list)





Conv = SpectralConvolution(Spectrum)

List = collections.Counter(Conv).most_common(M)
frequency = List[-1][1]

List1 = dict(collections.Counter(Conv))
ExtendedAlphabetMass = [k for k, v in List1.items() if v >= frequency]

ExtendedAlphabet = ''
for i in range(0,len(ExtendedAlphabetMass)):
    ExtendedAlphabet += chr(ExtendedAlphabetMass[i])

global AminoAcid
global AminoAcidMass
AminoAcid = list(ExtendedAlphabet)
AminoAcidMass = ExtendedAlphabetMass

AminoAcid_dict = dict(zip(AminoAcid, AminoAcidMass))



ans = LeaderboardCyclopeptideSequencing(Spectrum, N)


"""
tmp_ans = []
for key,val in score_dict.items():
    if val == score_dict[ans[0]]:
        tmp_ans.append(key)
"""

#print(top_M_acid_list)
text = ''
for item in ans:
    for i in item:
        text += str(aminoAcid_dict[i])
        text += '-'
        
    text = text[:-1]
    text += '\n'
text = text[:-1]
#print(text)

with open('output.txt', 'w') as handle:
    handle.write(text)
