#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 11:57:07 2018

@author: xiangyin
"""

#Counting Peptides with Given Mass Problem: Compute the number of peptides of given mass.
#     Input: An integer m.
#     Output: The number of linear peptides having integer mass m.


aminoAcid = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 
             'N', 'D', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']

aminoAcidMass = [ 57, 71, 87, 97, 99, 101, 103, 113,  
                 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

aminoAcid_dict = dict(zip(aminoAcid, aminoAcidMass))


def get_seqs(total_mass):
    #tmp_num = 0
    #mass_dict = {}
    text_list = []

    for i in aminoAcid:
        tmp_text = ''
        tmp_text += i
        if total_mass - aminoAcid_dict[i] < 0:
            #print(tmp_text)
            break
        #print(total_mass - aminoAcid_dict[i])

        elif total_mass - aminoAcid_dict[i] > 0:
            for item in get_seqs(total_mass - aminoAcid_dict[i]):
                tmp_item = tmp_text + item
                text_list.append(tmp_item)
        elif total_mass - aminoAcid_dict[i] == 0:
            #text_list.append(tmp_text)
            #print(i)
            #print(tmp_num)
            text_list.append(tmp_text)
            return text_list
    return text_list


mass_dict = {}

def count_seqs(total_mass):
    tmp_num = 0
    #text_list = []

    for i in aminoAcid:
        if (total_mass - aminoAcid_dict[i]) in mass_dict.keys():
            tmp_num += mass_dict[(total_mass - aminoAcid_dict[i])]
        elif total_mass - aminoAcid_dict[i] < 0:
            #print(tmp_text)
            #mass_dict[(total_mass - aminoAcid_dict[i])] = 0
            break
        #print(total_mass - aminoAcid_dict[i])
        elif total_mass - aminoAcid_dict[i] == 0:
            #text_list.append(tmp_text)
            #print(i)
            #print(tmp_num)
            tmp_num += 1
            return tmp_num
        elif total_mass - aminoAcid_dict[i] > 0:
            tmp_num += count_seqs(total_mass - aminoAcid_dict[i])
    mass_dict[total_mass] = tmp_num
    return tmp_num
       
total_mass = 1334   
ans_num =  count_seqs(total_mass)
print(ans_num)