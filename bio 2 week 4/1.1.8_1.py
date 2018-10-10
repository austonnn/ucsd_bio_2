#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 15:04:19 2018

@author: xiangyin
"""

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


def extended_mass_table():
    aminoAcid_dict = {}
    for i in range(57,201):
        aminoAcid_dict[chr(i)] = i
    return aminoAcid_dict


aminoAcid_dict = extended_mass_table()

#mass_dict = dict(aminoAcid_dict)
#global mass_dict
#mass_dict[''] = 0
#mass_dict = {'':0}

score_dict = {'':0}

line_score_dict = {'':0}

class Peptide_cls:
    def __init__(self, name):
        self.name = name
        self.cyclo_score = None
        self.line_score = None
        self.mass = None
        self.name_len = len(name)
    def get_cyclo_score(self, Spectrum):
        if self.cyclo_score == None:
            self.cyclo_score = Score(self, Spectrum)
        return self.cyclo_score
    def get_line_score(self, Spectrum):
        if self.line_score  == None:
            self.line_score = LinearScore(self, Spectrum)
        return self.line_score
    def get_mass(self):
        if self.mass  == None:
            self.mass = Mass(self.name)
        return self.mass
        

def Mass(Peptide):
    #print(len(mass_dict))
#    if Peptide in mass_dict.keys():
#        return mass_dict[Peptide]
    tmp_num = 0
    for item in Peptide:
        tmp_num += aminoAcid_dict[item]
#    tmp_num = mass_dict[Peptide[:-1]] + aminoAcid_dict[Peptide[-1]]
#    mass_dict[Peptide] = tmp_num
    return tmp_num

def CycloSpectrum(Peptide):
    PrefixMass = [0]
    for i in range(Peptide.name_len):
        tmp_num = PrefixMass[i] + aminoAcid_dict[Peptide.name[i]]
        PrefixMass.append(tmp_num)
    #peptideMass = PrefixMass[len(Peptide.name)]
    peptideMass = Peptide.get_mass()
    CyclicSpectrum = [0]
    for i in range(Peptide.name_len):
        for j in range(i + 1, Peptide.name_len + 1):
            tmp_mass = PrefixMass[j] - PrefixMass[i]
            CyclicSpectrum.append(tmp_mass)
            if i > 0 and j < Peptide.name_len:
                CyclicSpectrum.append(peptideMass - tmp_mass)
    CyclicSpectrum.sort()
    return CyclicSpectrum

def LinearSpectrum(Peptide):
    PrefixMass = [0]
    for i in range(Peptide.name_len):
        tmp_num = PrefixMass[i] + aminoAcid_dict[Peptide.name[i]]
        PrefixMass.append(tmp_num)
    LinearSpectrum = [0]
    for i in range(Peptide.name_len):
        for j in range(i + 1, Peptide.name_len + 1):
            LinearSpectrum.append(PrefixMass[j] - PrefixMass[i])
    LinearSpectrum.sort()
    return LinearSpectrum

def Score1(Peptide, Spectrum):
    #if Peptide in score_dict.keys():
    #    return score_dict[Peptide]
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
    #score_dict[Peptide] = ans_num
    return ans_num

def Score(Peptide, Spectrum):
    #if Peptide in score_dict.keys():
    #    return score_dict[Peptide.name]
    ans_num = 0
    Cyclo_Spectrum_list = CycloSpectrum(Peptide)
    tmp_Spectrum = list(Spectrum)
    for item in Cyclo_Spectrum_list:
        #print(type(item))
        #tmp_Spectrum = list(Spectrum)
        #print(tmp_Spectrum)
        if str(item) in tmp_Spectrum:
            ans_num += 1
            tmp_Spectrum.remove(str(item))
            #Spectrum = list(tmp_Spectrum)
    #score_dict[Peptide.name] = ans_num
    return ans_num

def LinearScore1(Peptide, Spectrum, line_score_dict):
    if Peptide in line_score_dict.keys():
        return line_score_dict[Peptide]
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
    line_score_dict[Peptide.name] = ans_num
    return ans_num

def LinearScore(Peptide, Spectrum):
#    if Peptide in line_score_dict.keys():
#        return line_score_dict[Peptide]
    ans_num = 0
    list_peptide = LinearSpectrum(Peptide)
    tmp_Spectrum = list(Spectrum)
    for item in list_peptide:
        #print(type(item))
        
        #print(tmp_Spectrum)
        if str(item) in tmp_Spectrum:
            ans_num += 1
            tmp_Spectrum.remove(str(item))
            
#    line_score_dict[Peptide] = ans_num
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
            tmp_Peptide  = Peptide_cls(Peptide.name + item)
#            if tmp_Peptide not in Peptide_dict.keys():
#                tmp_score = Score(tmp_Peptide, Spectrum)
#                tmp_mass = Mass(tmp_Peptide)
#                Peptide_dict[tmp_Peptide] = [tmp_score, tmp_mass]

            tmp_list.append(tmp_Peptide) 
    return tmp_list

def Trim1(Leaderboard, Spectrum, N):
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
        LinearScores_list.append(LinearScore(Peptide, Spectrum, line_score_dict))
        #LinearScores_list.append(Peptide_dict[Peptide][0])

    Leaderboard_list = list(zip(Leaderboard, LinearScores_list))
    Leaderboard_list = sorted(Leaderboard_list, key = lambda item:item[1], 
                      reverse = True)
    
    if N < len(Leaderboard):
        for j in range(N, len(Leaderboard)):
            if Leaderboard_list[j][1] < Leaderboard_list[N][1]:
                for item in Leaderboard_list[0 : j - 1]:
                    ans_Leaderboard.append(item[0])
                print("leardboard")
                print(Leaderboard_list[N][1])
                print(Leaderboard_list[j][1])
                return ans_Leaderboard
                #return tmp_Leaderboard[0 : j - 1]
    for item in Leaderboard_list:
        ans_Leaderboard.append(item[0])
    return ans_Leaderboard


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
#    LinearScores_list = []
    #ans_Leaderboard = []
#    for j in range(len(Leaderboard)):
#        Peptide = Leaderboard[j]
#        LinearScores_list.append(LinearScore(Peptide, Spectrum, line_score_dict))
#        #LinearScores_list.append(Peptide_dict[Peptide][0])
#
#    Leaderboard_list = list(zip(Leaderboard, LinearScores_list))
    Leaderboard_list = sorted(Leaderboard, 
                              key = lambda item:item.get_line_score(Spectrum), 
                              reverse = True)
    
    if N < len(Leaderboard):
        for j in range(N, len(Leaderboard)):
            if Leaderboard_list[j].get_line_score(Spectrum) < Leaderboard_list[N].get_line_score(Spectrum):
#                for item in Leaderboard_list[0 : j - 1]:
#                    ans_Leaderboard.append(item)
#                return ans_Leaderboard
                return Leaderboard_list[0 : j - 1]
#    for item in Leaderboard_list:
#        ans_Leaderboard.append(item[0])
#    return ans_Leaderboard
    return Leaderboard_list

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

def remove_peptide(Leaderboard, Peptide):
    #Leaderboard = list(Leaderboard)
    for item in Leaderboard:
        if item.name == Peptide.name:
            Leaderboard.remove(item)
            return Leaderboard

def LeaderboardCyclopeptideSequencing(Spectrum, N, mass_dict):
    #Leaderboard ← set containing only the empty peptide
    Peptide = Peptide_cls('')
    Leaderboard = []
    Leaderboard.append(Peptide)
    #Leaderboard = [['', 0, 0]]
    #LeaderPeptide ← empty peptide
    LeaderPeptide = Peptide_cls('')
    LeaderPeptide_list = []
    #LeaderPeptide = ['', 0, 0]
    #while Leaderboard is non-empty
    while Leaderboard:
        #Leaderboard ← Expand(Leaderboard)
        Leaderboard = Expand(Leaderboard)
        #print(len(Leaderboard))
        #print(len(mass_dict))
        #print(Leaderboard)
        #print(Leaderboard)
        #for each Peptide in Leaderboard
        tmp_Leaderboard = list(Leaderboard)
        #tmp_dict = {'':0}
        for Peptide in Leaderboard:

            #tmp_dict[Peptide.name] = Peptide.get_mass()
            #if Mass(Peptide) = ParentMass(Spectrum)
            if Peptide.get_mass() < int(Spectrum[-1]) and Peptide.get_mass() > int(Spectrum[-1]) - 57:
                tmp_Leaderboard = remove_peptide(tmp_Leaderboard, Peptide)  
            elif Peptide.get_mass() == int(Spectrum[-1]):
            #if Peptide_dict[Peptide][1] == int(Spectrum[-1]):
                if Peptide.get_cyclo_score(Spectrum) > LeaderPeptide.get_cyclo_score(Spectrum):
                #if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum):
                #if Peptide_dict[Peptide][0] > Peptide_dict[LeaderPeptide][0]:
                    LeaderPeptide_list = []
                    LeaderPeptide_list.append(Peptide)
                    LeaderPeptide = Peptide
                elif Peptide.get_cyclo_score(Spectrum) == LeaderPeptide.get_cyclo_score(Spectrum):
                #if Peptide_dict[Peptide][0] > Peptide_dict[LeaderPeptide][0]:
                    #LeaderPeptide_list = []
                    LeaderPeptide_list.append(Peptide)
                tmp_Leaderboard = remove_peptide(tmp_Leaderboard, Peptide)
            #else if Mass(Peptide) > ParentMass(Spectrum)
            elif Peptide.get_mass() > int(Spectrum[-1]):
            #elif Peptide_dict[Peptide][1]  > int(Spectrum[-1]):
                #remove Peptide from Leaderboard
                tmp_Leaderboard = remove_peptide(tmp_Leaderboard, Peptide)
        #Leaderboard ← Trim(Leaderboard, Spectrum, N)
        #print(tmp_dict)
        #mass_dict = dict(tmp_dict)
        Leaderboard = tmp_Leaderboard
        Leaderboard = AnotherTrim(Leaderboard, Spectrum, N)
        print('#')
        print(len(Leaderboard))
        #print(LeaderPeptide)
        #print(Leaderboard)

    #output LeaderPeptide
    return LeaderPeptide_list


#with open('input.txt') as handle:
with open('dataset_102_8.txt') as handle:
    N, Spectrum = handle.read().splitlines()

Spectrum = Spectrum.split()
N = int(N)
ans = LeaderboardCyclopeptideSequencing(Spectrum, N, mass_dict)

#print(ans)


#text = ''
#for item in ans:
#    text += item
#    text += ' '
#with open('output.txt', 'w') as handle:
#    handle.write(text)

    
text = ''
for item in ans:
    for i in item.name:
        text += str(aminoAcid_dict[i])
        text += '-'
        
    text = text[:-1]
    text += '\n'
text = text[:-1]
#print(text)
#print(text)    
with open('output.txt', 'w') as handle:
    handle.write(text)
    
    
#print(score('VKLFPWFNQY'))  
#print(Score('VKLFPWFNQY', Spectrum)) 
#print(Mass('VYKNFWPFIK'))