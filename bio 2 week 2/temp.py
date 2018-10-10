# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

def break_reads(patterns, k):
    tmp_pattern = []
    
    for item in patterns:
        item = item.split('|')        
        key = item[0][: -1] + '|' + item[1][: -1]
        tmp_pattern.append(key)
        value = item[0][1 :] + '|' + item[1][1 :]
        tmp_pattern.append(value)
    
    return tmp_pattern


"""
file = open("break_reads.txt")
#k, d = file.readline().strip().split()
k = 10 #int(k)
#d = 1000 #int(d)
patterns = file.read().splitlines()
#text = str(text[0])
file.close()
"""

tmp_list = ['aa','bb','aa']
tmp_set = set(tmp_list)
tmp_list1 = list(tmp_set)
print(tmp_list1)
