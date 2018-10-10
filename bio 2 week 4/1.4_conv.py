#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 14:42:45 2018

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
import timeit


#aminoAcid = ['V', 'K', 'L', 'F', 'P', 'W', 'F', 'N', 'Q', 'Y']
#
#aminoAcidMass = [99, 128, 113, 147, 97, 186, 147, 114, 128, 163]

aminoAcid = ['V', 'K', 'L', 'F', 'P', 'W', 'N', 'Y']

aminoAcidMass = [99, 128, 113, 147, 97, 115, 114, 163]

aminoAcid_dict = dict(zip(aminoAcid, aminoAcidMass))





def extended_mass_table():
    aminoAcid_dict = {}
    for i in range(57,201):
        aminoAcid_dict[chr(i)] = i
    return aminoAcid_dict

#aminoAcid_dict_ext = extended_mass_table()
#aminoAcid_dict = extended_mass_table()

mass_dict = {'':0}
score_dict = {'':0}
line_score_dict = {'':0}
tb_score_dict = {'':0}
tb_line_score_dict = {'':0}

def Mass(Peptide):
    if Peptide in mass_dict.keys():
        return mass_dict[Peptide]
    tmp_num = 0
    for item in Peptide:
        tmp_num += aminoAcid_dict[item]
#    tmp_num = mass_dict[Peptide[:-1]] + aminoAcid_dict[Peptide[-1]]
    mass_dict[Peptide] = tmp_num
    return tmp_num

#print(Mass('NQYVKLFPW'))
#print(Mass('WFNQYVKL'))
#print(Mass('WFNQYVKL'))
#print(Mass('FPWFNQYVKL'))
#print(CycloSpectrum('FPWFNQYVKL'))

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
    #ans_num += tb_Score(Peptide, Spectrum)
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

def tb_Score(Peptide, Spectrum):
    if Peptide in tb_score_dict.keys():
        return tb_score_dict[Peptide]
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
            
    tb_score_dict[Peptide] = ans_num
    return ans_num



def LinearScore(Peptide, Spectrum):
    if Peptide in line_score_dict.keys():
        return line_score_dict[Peptide]
    ans_num = 0
    #ans_num += tb_LinearScore(Peptide, Spectrum_tb)
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
    line_score_dict[Peptide] = ans_num
    return ans_num

def tb_LinearScore(Peptide, Spectrum_tb):
    if Peptide in tb_line_score_dict.keys():
        return tb_line_score_dict[Peptide]
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
    tb_line_score_dict[Peptide] = ans_num
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



def LeaderboardCyclopeptideSequencing_old(Spectrum, N):
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
#            if Mass(Peptide) < int(Spectrum[-1]):
#                if Mass(Peptide) > int(Spectrum[-1]) - 57:
#                    tmp_Leaderboard.remove(Peptide)
            
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
                #tmp_Leaderboard.remove(Peptide)
                #tmp_Leaderboard.remove(Peptide)   
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
            #print('$$')
            #print(len(Peptide))
            #if Mass(Peptide) = ParentMass(Spectrum)
#            if Mass(Peptide) < int(Spectrum[-1]):
#                if Mass(Peptide) > int(Spectrum[-1]) - 57:
#                    tmp_Leaderboard.remove(Peptide)
            #if len(Peptide) == 10:
#            if Mass(Peptide) <= 1308:
#                tmp_Leaderboard.remove(Peptide)
            if len(Peptide) == 10:    
            #if Mass(Peptide) == 1283:#1322:#1283:
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
                #tmp_Leaderboard.remove(Peptide)
                tmp_Leaderboard.remove(Peptide)   
            #else if Mass(Peptide) > ParentMass(Spectrum)
            elif len(Peptide) > 10: 
            #elif Mass(Peptide) > 1283:#1322:#1283:
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

def Convolution_Spectrum(Spectrum):
    ans_list = []
    # if Spectrum does not include a 0 in it， 
    # should add one to the following list.
    # Be careful!
    tmp_list = []
    #Spectrum = sorted(Spectrum)
    for item in Spectrum:
        tmp_list.append(int(item))
    #tmp_list.extend(Spectrum)
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

def round_number(x):
    pass

with open('real_spectrum.txt') as handle:
#with open('Tyrocidine_B1_Spectrum_25_real.txt') as handle:
    Spectrum = handle.read().split()
#A cyclic peptide LeaderPeptide with amino acids taken only 
#from the top M elements (and ties) of the convolution of
#Spectrum that fall between 57 and 200,
#M = int(M)
#  the size of Leaderboard is restricted to the top N (and ties).
N = 1000



#Spectrum = Spectrum.split()
# next line, spectrum do not have 0 in it


#tmp_list = [0]
#因为我注意到完全没有200以下的，既然我们可以假设已知组成，那么也可以补充上这几个质量。
#tmp_list = [0, 97, 99, 113, 114, 128,  147,  163, 115]
tmp_list = [0]
#tmp_list.extend(aminoAcid_dict.values())
tmp_list1 = []
for item in Spectrum:
    tmp_num = float(item) - 1
    print('#')
    #tmp_list.append(tmp_num)
    #tmp_list.append(round(tmp_num))
    tmp_list.append(int(tmp_num))
    #tmp_list.append(round(tmp_num + 0.3))
    #print(tmp_num)
    #tmp_list.append(round(tmp_num - 0.3))
    #print(tmp_num)
#test_pep = 'NQYVKLFPFF'   
test_pep = 'NKYVKLFPFF' 

Spectrum = sorted(tmp_list)
#Spectrum.append(1283)
#Spectrum_tb = [0, 97, 99, 113, 114, 128, 128, 147, 147, 163, 186, 227, 241, 242, 244, 260, 261, 262, 283, 291, 333, 340, 357, 388, 389, 390, 390, 405, 430, 430, 447, 485, 487, 503, 504, 518, 543, 544, 552, 575, 577, 584, 631, 632, 650, 651, 671, 672, 690, 691, 738, 745, 747, 770, 778, 779, 804, 818, 819, 835, 837, 875, 892, 892, 917, 932, 932, 933, 934, 965, 982, 989, 1031, 1039, 1060, 1061, 1062, 1078, 1080, 1081, 1095, 1136, 1159, 1175, 1175, 1194, 1194, 1208, 1209, 1223, 1225, 1322]
#Spectrum.append(Spectrum_tb)

#Spectrum_test = CycloSpectrum(test_pep)

#start = timeit.default_timer()
#sorted_acid_list = Convolution_Spectrum_sorted(Spectrum)
##print(sorted_acid_list)
#M = 20
#top_M_acid_list = top_M_acid(sorted_acid_list, M)
##top_M_acid_list = [99, 128, 113, 147, 97, 186, 147, 114, 128, 163]
##print(len(top_M_acid_list))
##print_top_M(top_M_acid_list)
#aminoAcid_dict = {}
#for item in top_M_acid_list:
#    aminoAcid_dict[item] = ord(item)

#因为我注意到完全没有200以下的，既然我们可以假设已知组成，那么也可以补充上这几个质量。
#Spectrum.extend(aminoAcid_dict.values())
#amino_Acid_list = [0, 97, 99, 113, 114, 128,  147,  163, 115]
Spectrum.extend(aminoAcid_dict.values())
Spectrum = sorted(Spectrum)

start = timeit.default_timer()   
ans = LeaderboardCyclopeptideSequencing(Spectrum, N)

#ans = ConvolutionCyclopeptideSequencing(Spectrum, N, M)
stop = timeit.default_timer()
print (stop - start)
test_pep = 'NKYVKLFPFF'
b25 = 'PWFNKYVKLF'
#print(ans)


#text = ''
#for item in ans:
#    text += item
#    text += ' '
#with open('output.txt', 'w') as handle:
#    handle.write(text)

   
text = ''
for item in ans:
    for i in item:
        text += str(aminoAcid_dict[i])
        text += ' '
        
    text = text[:-1]
    text += '\n'
text = text[:-1]
print(text)
#print(text)    
with open('output.txt', 'w') as handle:
    handle.write(text)
    
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
    print(Score(Peptide, Spectrum))#print(score('VKLFPWFNQY'))  
#print(Score('VKLFPWFNQY', Spectrum)) 
#print(Mass('VYKNFWPFIK'))
for Peptide in ans:
    if Score(Peptide, Spectrum) != 82:
        print('#*#')
        print (Peptide)
  
print('score')
#print(Score('NKYVKLFPFF', Spectrum))
#print(Score('LFPWFNKYVK', Spectrum))
