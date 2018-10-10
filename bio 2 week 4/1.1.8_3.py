#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 15:04:19 2018

@author: xiangyin
"""
# currect =  current&correct answer

import timeit
aminoAcid = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 
             'N', 'D', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']

aminoAcidMass = [ 57, 71, 87, 97, 99, 101, 103, 113,  
                 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

aminoAcid_dict = dict(zip(aminoAcid, aminoAcidMass))

aminoAcid_dict['Q'] = 128
aminoAcid_dict['L'] = 113

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
    
def Expand(Leaderboard):
    tmp_list = []
    for Peptides in Leaderboard:
        tmp_list.extend(ExpandPeptides(Peptides))    
    return tmp_list

def ExpandPeptides(Peptides):
    tmp_list = []
    for peptide in Peptides:
        for item in aminoAcid:
            #tmp_text
            tmp_text = peptide + item
            #peptide += item
            tmp_list.append(tmp_text)    
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

    Leaderboard_list = list(zip(Leaderboard, LinearScores_list))
    tmp_list = sorted(Leaderboard_list, key = lambda item:item[1], 
                      reverse = True)
    for j in range(N, len(Leaderboard)):
        if tmp_list[j][1] < tmp_list[N][1]:
            for item in tmp_list[0 : j - 1]:
                ans_Leaderboard.append(item[0])
            return ans_Leaderboard
    for item in tmp_list:
        ans_Leaderboard.append(item[0])
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
    #while Leaderboard is non-empty
    while Leaderboard:
        #Leaderboard ← Expand(Leaderboard)
        Leaderboard = Expand(Leaderboard)
        tmp_Leaderboard = list(Leaderboard)
        #for each Peptide in Leaderboard
        for Peptide in Leaderboard:
            #if Mass(Peptide) = ParentMass(Spectrum)
            if Mass(Peptide) == Spectrum[-1]:
                
                #if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum)
                if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum):
                    LeaderPeptide = Peptide
            #else if Mass(Peptide) > ParentMass(Spectrum)
            elif Mass(Peptide) > Spectrum[-1]:
                #remove Peptide from Leaderboard
                tmp_Leaderboard.remove(Peptide)
        #Leaderboard ← Trim(Leaderboard, Spectrum, N)
        print(len(Leaderboard))
        Leaderboard = Trim(tmp_Leaderboard, Spectrum, N)
    #output LeaderPeptide
    return LeaderPeptide


with open('input.txt') as handle:
    N, Spectrum = handle.read().splitlines()
#Leaderboard = Leaderboard.split()
Spectrum = Spectrum.split()
N = int(N)

start = timeit.default_timer()
ans = LeaderboardCyclopeptideSequencing(Spectrum, N)
stop = timeit.default_timer()
print (stop - start)

print(ans)
text = ''
for item in ans:
    text += item
    text += ' '
with open('output.txt', 'w') as handle:
    handle.write(text)
    
    
    