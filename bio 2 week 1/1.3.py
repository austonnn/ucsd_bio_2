# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 20:34:31 2017

@author: xiangyin
"""

def genome_path_to_string(genome_path):
    tmp_string = [] 
    tmp_string.append(genome_path[0])
    for text in genome_path[1:]:
        tmp_string.append(text[-1])
    ans = ''.join(tmp_string)
    return ans

"""

file = open("dataset_198_3.txt")
#k, t, n = file.readline().strip().split()
genome_path = file.read().splitlines()
file.close()


#print(genome_path_to_string(genome_path))
f = open('data.txt','w')
f.write(genome_path_to_string(genome_path))
f.close()
"""


def prefix(text):
    return text[:-1]
    
def suffix(text):
    return text[1:]

def overlap(patterns):
    graph = []    
    for item in patterns:
        tmp_list = list(patterns)
        tmp_list.remove(item)
        for tmp_item in tmp_list:
            if suffix(item) == prefix(tmp_item):
                graph.append([item, tmp_item])                
    return graph


file = open("dataset_198_10.txt")
#k, t, n = file.readline().strip().split()
patterns = file.read().splitlines()
file.close()


graph = overlap(patterns)
tmp_list = []
for item in graph:
    tmp_str = item[0] + ' -> ' + item[1]
    tmp_list.append(tmp_str)

txt = '\n'.join(tmp_list)
f = open('1.3.txt','w')
f.write(txt)
f.close()
