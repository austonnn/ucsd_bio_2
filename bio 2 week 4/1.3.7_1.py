#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 22 04:02:47 2018

@author: xiangyin
"""
import collections

mass_dict = {'':0}
score_dict = {'':0}
linear_Score_dict = {'':0}

import timeit
#aminoAcid = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 
#             'N', 'D', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
#
#aminoAcidMass = [ 57, 71, 87, 97, 99, 101, 103, 113,  
#                 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]
#
#aminoAcid_dict = dict(zip(aminoAcid, aminoAcidMass))
#
#aminoAcid_dict['Q'] = 128
#aminoAcid_dict['L'] = 113

def Mass(Peptide):
    if Peptide in mass_dict.keys():
        return mass_dict[Peptide]
    tmp_num = 0
#    for item in Peptide:
#        tmp_num += aminoAcid_dict[item]
    tmp_num = mass_dict[Peptide[:-1]] + aminoAcid_dict[Peptide[-1]]
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
    tmp_Spectrum = list(Spectrum)
    for item in list_peptide:
        #print(type(item))
        #tmp_Spectrum = list(Spectrum)
        #print(tmp_Spectrum)
        if item in tmp_Spectrum:
            ans_num += 1
            tmp_Spectrum.remove(item)
            #Spectrum = list(tmp_Spectrum)
    score_dict[Peptide] = ans_num
    return ans_num


def LinearScore(Peptide, Spectrum):
    if Peptide in linear_Score_dict.keys():
        return linear_Score_dict[Peptide]
    ans_num = 0
    list_peptide = LinearSpectrum(Peptide)
    tmp_Spectrum = list(Spectrum)
    for item in list_peptide:
        #print(type(item))
        #tmp_Spectrum = list(Spectrum)
        #print(tmp_Spectrum)
        if item in tmp_Spectrum:
            ans_num += 1
            tmp_Spectrum.remove(item)
            #Spectrum = list(tmp_Spectrum)
    linear_Score_dict[Peptide] = ans_num
    #print(ans_num)
    return ans_num
    
#def Expand1(Leaderboard):
#    #print(len(Leaderboard))
#    tmp_list = []
#    for Peptides in Leaderboard:
#        print(len(Peptides))
#        tmp_list.extend(ExpandPeptides(Peptides))    
#    return tmp_list
#
#def ExpandPeptides(Peptides):
#    print(Peptides)
#    print(len(Peptides))
#    tmp_list = []
#    for peptide in Peptides:
#        #print(AminoAcid)
#        for item in aminoAcid_dict.keys():
#            #tmp_text
#            tmp_text = peptide + item
#            #peptide += item
#            tmp_list.append(tmp_text)  
#    
#    return tmp_list

def Expand(Leaderboard):
    tmp_list = []
    for Peptide in Leaderboard:
        # amino acids taken only from the top M elements (and ties) 
        # of the convolution of Spectrum that fall between 57 and 200
        #print(len(top_M_acid_list))
        for item in AminoAcid:
            tmp_Peptide  = Peptide + item
            tmp_list.append(tmp_Peptide) 
    
    return tmp_list

def Trim(Leaderboard, Spectrum, N):
    #print('##')
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

    Leaderboard_list = list(zip(Leaderboard, LinearScores_list))
    Leaderboard_list = sorted(Leaderboard_list, key = lambda item:item[1], 
                      reverse = True)
    #print(Leaderboard_list[0])
    if N < len(Leaderboard):
        for j in range(N, len(Leaderboard)):
            #
            if Leaderboard_list[j][1] < Leaderboard_list[N][1]:
                
                for item in Leaderboard_list[0 : j - 1]:
                    ans_Leaderboard.append(item[0])
                return ans_Leaderboard
    for item in Leaderboard_list:
        ans_Leaderboard.append(item[0])
    #print(len(ans_Leaderboard))
    return ans_Leaderboard

def AnotherTrim(Leaderboard, Spectrum, N):
    # ref http://bioinformaticsalgorithms.com/faqs/antibiotics.html#week3
    #print(len(Leaderboard))
    ScoreHistogram = []
    LinearScore_list = []
    tmp_Leaderboard = list(Leaderboard)
    for i in range(len(Spectrum) + 1):
        ScoreHistogram.append(0)

    for j in range(len(Leaderboard)):
        Peptide = Leaderboard[j]
        #tmp_LinearSpectrum  = LinearSpectrum(Peptide)
        #tmp_LinearScore = Score(Peptide, tmp_LinearSpectrum)
        
        tmp_LinearScore = LinearScore(Peptide, Spectrum)
        #print(tmp_LinearScore)
        LinearScore_list.append(tmp_LinearScore)
        ScoreHistogram[tmp_LinearScore] = ScoreHistogram[tmp_LinearScore] + 1
    #print(ScoreHistogram)
    tmp_sum = 0
    for i in range(len(Spectrum) + 1):
        tmp_sum =  tmp_sum + ScoreHistogram[i]
        if tmp_sum > len(Leaderboard) - N:
            ScoreThreshold =  i - 1
            print(ScoreThreshold)
            break
            #print(ScoreThreshold)
    for j in range(len(Leaderboard)):
        #print(Peptide)
        Peptide =  Leaderboard[j]
        #print(tmp_Leaderboard)
        if LinearScore_list[j] <= ScoreThreshold:
                tmp_Leaderboard.remove(Peptide)
    return tmp_Leaderboard

def LeaderboardCyclopeptideSequencing(Spectrum, N):
    #Leaderboard ← set containing only the empty peptide
    Leaderboard = ['']
    #LeaderPeptide ← empty peptide
    LeaderPeptide = ''
    LeaderPeptide_list = []
    #while Leaderboard is non-empty
    #print(len(Leaderboard))
    while Leaderboard:
        #Leaderboard ← Expand(Leaderboard)
        Leaderboard = Expand(Leaderboard)
        tmp_Leaderboard = list(Leaderboard)
        #print(len(Leaderboard))
        #for each Peptide in Leaderboard
        for Peptide in Leaderboard:
            #if Mass(Peptide) = ParentMass(Spectrum)
            if Mass(Peptide) == Spectrum[-1]:
                
                #if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum)
                if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum):
                    LeaderPeptide_list = []
                    LeaderPeptide_list.append(Peptide)
                    LeaderPeptide = Peptide
                elif Score(Peptide, Spectrum) == Score(LeaderPeptide, Spectrum):
                #if Peptide_dict[Peptide][0] > Peptide_dict[LeaderPeptide][0]:
                    #LeaderPeptide_list = []
                    LeaderPeptide_list.append(Peptide)
            #else if Mass(Peptide) > ParentMass(Spectrum)
            elif Mass(Peptide) > Spectrum[-1]:
                #remove Peptide from Leaderboard
                tmp_Leaderboard.remove(Peptide)
        #Leaderboard ← Trim(Leaderboard, Spectrum, N)
        #print(len(Leaderboard))
        #print(Leaderboard)
        Leaderboard = Trim(tmp_Leaderboard, Spectrum, N)
        print(len(Leaderboard))
    #output LeaderPeptide
    return LeaderPeptide_list


with open('dataset_104_8.txt') as handle:
#with open('dataset_104_7.txt') as handle:
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

def SpectralConvolution(Spectrum):
    ans_list = []
    tmp_list = [0]
    #Spectrum = sorted(Spectrum)
    tmp_list.extend(Spectrum)
    tmp_list = sorted(tmp_list)
    for i in tmp_list:
        for j in tmp_list:
            if 56 <(i - j)<201:
                ans_list.append(i - j)

    return ans_list

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
AminoAcid = ExtendedAlphabet
AminoAcidMass = ExtendedAlphabetMass

global aminoAcid_dict
aminoAcid_dict = dict(zip(AminoAcid, AminoAcidMass))

#ans = LeaderboardCyclopeptideSequencing(Spectrum, N)

start = timeit.default_timer()
ans = LeaderboardCyclopeptideSequencing(Spectrum, N)
stop = timeit.default_timer()
print (stop - start)


text = ''
for item in ans:
    for i in item:
        text += str(aminoAcid_dict[i])
        text += '-'
        
    text = text[:-1]
    text += '\n'
text = text[:-1]
#print(text)
#print(text)    
with open('output.txt', 'w') as handle:
    handle.write(text)


