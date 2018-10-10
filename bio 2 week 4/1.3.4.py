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

def Convolution_Spectrum(Spectrum):
    ans_list = []
    tmp_list = [0]
    tmp_list.extend(Spectrum)
    for i in tmp_list:
        for j in Spectrum:
            if i > j:
                ans_list.append(i - j)
            pass
    ans_list = sorted(ans_list)
    return ans_list
    



with open('dataset_104_4.txt') as handle:
    Spectrum = handle.read() #.splitlines()
    
Spectrum = Spectrum.split()
tmp_list = []
for item in Spectrum:
    tmp_list.append(int(item))
Spectrum = tmp_list


ans = Convolution_Spectrum(Spectrum)


#print(ans)

text = ''
for item in ans:
    text += str(item)
    text += ' '
print(text)
with open('output.txt', 'w') as handle:
    handle.write(text)