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



#def extended_mass_table():
#    aminoAcid_dict = {}
#    for i in range(57,201):
#        aminoAcid_dict[chr(i)] = i
#    return aminoAcid_dict
#
#aminoAcid_dict = extended_mass_table()
#
#mass_dict = {'':0}
#score_dict = {'':0}

import sys
import timeit
import heapq
from itertools import product
from collections import Counter

#table = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]

table = []
for i in range(57, 201):
    table.append(i)
    


def read_file(input_file):
    f = open(input_file)
    N,Spectrum = [item.strip() for item in f.readlines()]
    f.close()
    return (int(N),map(int,Spectrum.split(' ')))

def score(theorectical_spectrum,experimental_spectrum):
    return len(list((Counter(theorectical_spectrum) & Counter(experimental_spectrum)).elements()))

def generate_subspectrums(peptide):
    l = len(peptide)
    looped = peptide + peptide
    return [0,sum(peptide)]+[sum(looped[start:start+length]) for start,length in product(range(0,l),range(1,l))]

def cut(Leaderboard,Spectrum,N):
    if len(Leaderboard) > N:
        results = []
        for Peptide in Leaderboard:
            try:
                Peptide_experimental_spectrum = generate_subspectrums(Peptide)
            except:
                Peptide = Peptide[0]+[Peptide[1]]
                Peptide_experimental_spectrum = generate_subspectrums(Peptide)
            results.append((Peptide,score(Spectrum,Peptide_experimental_spectrum)))
        tie = heapq.nlargest(N,results,key=lambda x: x[1])[-1][1]
        res = list(filter(lambda x: x[1]>=tie,results))
        return list(zip(*res))[0]
    else:
        return Leaderboard
            
def LeaderboardCyclopeptideSequencing(Spectrum,N):
    Leaderboard = [0]
    LeaderPeptide = []
    while Leaderboard != []:
        Leaderboard = [list(pt) for pt in product(Leaderboard,table)]
        for Peptide in Leaderboard:
            try:
                Peptide_experimental_spectrum = generate_subspectrums(Peptide)
            except:
                Leaderboard = [Peptide[0]+[Peptide[1]] if x == Peptide else x for x in Leaderboard]
                Peptide = Peptide[0]+[Peptide[1]]
                Peptide_experimental_spectrum = generate_subspectrums(Peptide)
            if max(Peptide_experimental_spectrum) == max(Spectrum):
                LeaderPeptide_experimental_spectrum = generate_subspectrums(LeaderPeptide)
                if score(Spectrum,Peptide_experimental_spectrum) > score(Spectrum,LeaderPeptide_experimental_spectrum):
                    LeaderPeptide = Peptide
            elif max(Peptide_experimental_spectrum) > max(Spectrum):
                Leaderboard.remove(Peptide)
        Leaderboard = cut(Leaderboard,Spectrum,N)
    return LeaderPeptide

def result(filename):
    #with open('input.txt') as handle:
    with open('dataset_102_8.txt') as handle:
        N, Spectrum = handle.read().splitlines()
        #
    Spectrum = Spectrum.split()
    tmp_spectrum = []
    for item in Spectrum:
        tmp_spectrum.append(int(item))
    Spectrum = tmp_spectrum 
    N = int(N)
    #print(Spectrum)
    results = LeaderboardCyclopeptideSequencing(Spectrum,N)
    return results[1:]

if __name__ == "__main__":

    start = timeit.default_timer()
    #results = result(sys.argv[-1])
    results = result('input.txt')
    print ('-'.join(map(str,results)))
    print ('')
    stop = timeit.default_timer()
    print (stop - start)
#
##with open('input.txt') as handle:
#with open('dataset_102_8.txt') as handle:
#    N, Spectrum = handle.read().splitlines()
#
#Spectrum = Spectrum.split()
#N = int(N)
#ans = LeaderboardCyclopeptideSequencing(Spectrum, N)
#
##print(ans)
#
#
##text = ''
##for item in ans:
##    text += item
##    text += ' '
##with open('output.txt', 'w') as handle:
##    handle.write(text)
#
#    
#text = ''
#for item in ans:
#    for i in item:
#        text += str(aminoAcid_dict[i])
#        text += '-'
#        
#    text = text[:-1]
#    text += '\n'
#text = text[:-1]
##print(text)
##print(text)    
#with open('output.txt', 'w') as handle:
#    handle.write(text)
#    
#    
##print(score('VKLFPWFNQY'))  
##print(Score('VKLFPWFNQY', Spectrum)) 
##print(Mass('VYKNFWPFIK'))