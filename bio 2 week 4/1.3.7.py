#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 15:23:09 2018

@author: xiangyin
"""

"""
Spectral Convolution Problem: Compute the convolution of a spectrum.
     Input: A collection of integers Spectrum.
     Output: The list of elements in the convolution of Spectrum. If an element has multiplicity k, it should appear exactly k times;
    you may return the elements in any order.
"""

# import from 1.2
"""
aminoAcid = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 
             'N', 'D', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']

aminoAcidMass = [ 57, 71, 87, 97, 99, 101, 103, 113,  
                 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

aminoAcid_dict = dict(zip(aminoAcid, aminoAcidMass))

aminoAcid_dict['Q'] = 128
aminoAcid_dict['L'] = 113
"""
import timeit
#def extended_mass_table():
#    aminoAcid_dict = {}
#    for i in range(57,201):
#        aminoAcid_dict[chr(i)] = i
#    return aminoAcid_dict
#
#aminoAcid_dict = extended_mass_table()

mass_dict = {'':0}
score_dict = {'':0}
linear_Score_dict = {'':0}

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
    #print(len(Spectrum))
    tmp_Spectrum = list(Spectrum)
    for item in list_peptide:
        #print(type(item))
        #print(tmp_Spectrum)
        if item in tmp_Spectrum:
            #tmp_Spectrum = list(Spectrum)
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

        #print(tmp_Spectrum)
        if item in tmp_Spectrum:
            ans_num += 1
            tmp_Spectrum.remove(item)
    linear_Score_dict[Peptide] = ans_num
    return ans_num

   
#def Expand(Leaderboard):
#    tmp_list = []
#    for Peptides in Leaderboard:
#        tmp_list.extend(ExpandPeptides(Peptides))    
#    return tmp_list

def Expand(Leaderboard, top_M_acid_list):
    tmp_list = []
    for Peptide in Leaderboard:
        # amino acids taken only from the top M elements (and ties) 
        # of the convolution of Spectrum that fall between 57 and 200
        #print(len(top_M_acid_list))
        for item in top_M_acid_list:
            tmp_Peptide  = Peptide + item
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
        #print(sorted(LinearScores_list,reverse = True)[0])
    #print('LinearScores_list#')
    #print(LinearScores_list)
    Leaderboard_list = list(zip(Leaderboard, LinearScores_list))
    Leaderboard_list = sorted(Leaderboard_list, key = lambda item:item[1], 
                      reverse = True)
    #print(Leaderboard_list)
    #print(Leaderboard_list)
    if N < len(Leaderboard):
        #print(N)
        #print(len(Leaderboard))
        for j in range(N, len(Leaderboard)):
            #print(j)
            if Leaderboard_list[j][1] < Leaderboard_list[N][1]:


                for item in Leaderboard_list[0 : j-1]:
                    ans_Leaderboard.append(item[0])
                #print(len(ans_Leaderboard))
                #print('#')
                #print(len(ans_Leaderboard))
                return ans_Leaderboard
                #return tmp_Leaderboard[0 : j - 1]
    for item in Leaderboard_list:
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
            break

    for j in range(len(Leaderboard)):
        #print(Peptide)
        Peptide =  Leaderboard[j]
        #print(tmp_Leaderboard)
        if LinearScore_list[j] < ScoreThreshold:
                tmp_Leaderboard.remove(Peptide)
    return tmp_Leaderboard

def LeaderboardCyclopeptideSequencing(Spectrum, N, top_M_acid_list):
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
        Leaderboard = Expand(Leaderboard, top_M_acid_list)
        
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
                
                tmp_Leaderboard.remove(Peptide)      
            #else if Mass(Peptide) > ParentMass(Spectrum)
            elif Mass(Peptide) > int(Spectrum[-1]):
            #elif Peptide_dict[Peptide][1]  > int(Spectrum[-1]):
                #remove Peptide from Leaderboard
                tmp_Leaderboard.remove(Peptide)

        #Leaderboard ← Trim(Leaderboard, Spectrum, N)
        #print(Leaderboard)
        #Leaderboard = tmp_Leaderboard
        Leaderboard = Trim(tmp_Leaderboard, Spectrum, N)
        print(len(Leaderboard))
        #print(len(Leaderboard))
        #print(LeaderPeptide)
        #print(Leaderboard)

    #output LeaderPeptide
    return LeaderPeptide_list

# end import 1.2
def Convolution_Spectrum(Spectrum):
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

def Convolution_Spectrum_sorted(Spectrum):
    tmp_list = Convolution_Spectrum(Spectrum)
    ans_dict = {}
    for item in tmp_list:
        if item not in ans_dict.keys():
            ans_dict[item] = 0
        ans_dict[item] += 1
    ans_list = sorted(ans_dict.items(), key = lambda item:item[1], 
                      reverse = True)
    #print(ans_list)               
    return ans_list
    
def top_M_acid(sorted_acid_list, M):
    ans_list = []
#    print('#')
#    print(sorted_acid_list)
#    print(len(sorted_acid_list))
    if M < len(sorted_acid_list):
        for j in range(M, len(sorted_acid_list)):
            if sorted_acid_list[j][1] < sorted_acid_list[M][1]:
                for item in sorted_acid_list[0 : j]:

                    ans_list.append(chr(item[0]))
                return ans_list
                #return tmp_Leaderboard[0 : j - 1]
    for item in sorted_acid_list:
        ans_list.append(chr(item[0]))
    
    return ans_list


def ConvolutionCyclopeptideSequencing(Spectrum, N, M):
    sorted_acid_list = Convolution_Spectrum_sorted(Spectrum)
    #print(sorted_acid_list)
    top_M_acid_list = top_M_acid(sorted_acid_list, M)
    #print(len(top_M_acid_list))
    print_top_M(top_M_acid_list)
    ans = LeaderboardCyclopeptideSequencing(Spectrum, N, top_M_acid_list)
    
    
    return ans

def print_top_M(top_M_acid_list):
    text = ''
    for item in top_M_acid_list:
        text += str(aminoAcid_dict[item])
        text += ' '
    print('top_M_acid_list')
    print(text)
    return text


def get_score_mass(peptide_text):
    tmp_list = peptide_text.split('-')
    Peptide = ''
    for item in tmp_list:
        Peptide += chr(int(item))
          
    ans = Score(Peptide, Spectrum)
    return ans
        
    
    
    

#with open('input.txt') as handle:
with open('dataset_104_8.txt') as handle:
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

start = timeit.default_timer()
sorted_acid_list = Convolution_Spectrum_sorted(Spectrum)
#print(sorted_acid_list)
top_M_acid_list = top_M_acid(sorted_acid_list, M)

#print(len(top_M_acid_list))
#print_top_M(top_M_acid_list)
aminoAcid_dict = {}
for item in top_M_acid_list:
    aminoAcid_dict[item] = ord(item)
    
ans = LeaderboardCyclopeptideSequencing(Spectrum, N, top_M_acid_list)

#ans = ConvolutionCyclopeptideSequencing(Spectrum, N, M)
stop = timeit.default_timer()
print (stop - start)


tmp_ans = []
for key,val in score_dict.items():
    if val == score_dict[ans[0]]:
        tmp_ans.append(key)

#print(top_M_acid_list)
text = ''
for item in tmp_ans:
    for i in item:
        text += str(aminoAcid_dict[i])
        text += '-'
        
    text = text[:-1]
    text += '\n'
text = text[:-1]
#print(text)

with open('output.txt', 'w') as handle:
    handle.write(text)
  
#num = 0   
#for val in score_dict.values():
#    if val == 21:
#        num += 1
#print(num)
#        

def get_cyc_score(pep_num, Spectrum):
    pep_num = pep_num.split('-')
    Peptide = ''
    for item in pep_num:
        Peptide += chr(int(item))
    print('Peptide')
    print(Peptide)
    print('line score')
    print(LinearScore(Peptide, Spectrum))
    print('cyc score')
    print(Score(Peptide, Spectrum))

#print(get_cyc_score('99-71-137-57-72-57', Spectrum))
#Spectrum = [0,113,114,128,129,22]
    
