# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 20:07:50 2017

@author: xiangyin
"""

def string_spelled_gapped_patterns(gapped_patterns, k, d):
        #FirstPatterns ← the sequence of initial k-mers from GappedPatterns
        #SecondPatterns ← the sequence of terminal k-mers from GappedPatterns
        first_patterns = []
        second_patterns = []
        for item in gapped_patterns:
            first_patterns.append(item[:k])
            second_patterns.append(item[-k:])
        #PrefixString ← StringSpelledByGappedPatterns(FirstPatterns, k)
        #SuffixString ← StringSpelledByGappedPatterns(SecondPatterns, k)        
        prefix_string = string_spelled_patterns(first_patterns, k)
        suffix_string = string_spelled_patterns(second_patterns, k)
        for idx in range( k + d + 1, len(prefix_string)):
            #if the i-th symbol in PrefixString does not equal the (i - k - d)-th symbol in SuffixString
            if prefix_string[idx] != suffix_string[idx - k - d]:
                return "there is no string spelled by the gapped patterns"
        #return PrefixString concatenated with the last k + d symbols of SuffixString
        
        return prefix_string + suffix_string[-(k + d):]
    
def string_spelled_patterns(patterns, k):
    string = ''
    for item in patterns:
        string += item[0]
    string += item[1:]
        
    return string

    
file = open("dataset_6206_7.txt")
k, d = file.readline().strip().split()
k = int(k)
d = int(d)
gapped_patterns = file.read().splitlines()
#text = str(text[0])
file.close()

ans = string_spelled_gapped_patterns(gapped_patterns, k, d)
print(ans)

file = open('1.5.txt','w')
file.write(ans)
file.close()
