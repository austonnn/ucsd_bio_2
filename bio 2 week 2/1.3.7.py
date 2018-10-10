# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 12:58:13 2017

@author: xiangyin
"""


k = 3
d = 2
input_text = 'TAATGCCATGGGATGTT'
def list_sorted_kmer(text, k, d):
    kmer_list = []
    for i in range(len(text) - k - d - k + 1):
        tmp_list = []
        tmp_list.append(text[i: i +k])
        tmp_list.append(text[i+k+d:i+k+d+k])
        kmer_list.append(tmp_list)
        #sort(kmer_list)
        kmer_list.sort()
    print(len(kmer_list))
    return kmer_list

text = ''
for item in list_sorted_kmer(input_text, k, d):
    text += '('
    text += str(item[0])
    text += '|'
    text += str(item[1])
    text += ')'
    
#print(list_sorted_kmer(input_text, k, d))
print(text)







