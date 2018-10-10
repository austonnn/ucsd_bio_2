# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 10:01:55 2017

@author: xiangyin
"""

def de_bruijn_pattern(patterns):
    graph_dict = dict()
    for item in patterns:
        
        key = item[: -1]
        value = item[1 :]
        
        if key in graph_dict.keys():
            graph_dict[key].append(value)
        else:
            graph_dict[key] = []
            graph_dict[key].append(value)
    
    
    return graph_dict


file = open("dataset_200_8.txt")
#k = file.readline().strip()
patterns = file.read().splitlines()
#text = str(text[0])
file.close()


graph_dict = de_bruijn_pattern(patterns)
tmp_list = []
for key, values in graph_dict.items():
    if len(values) == 0:
        pass
    values = sorted(values)
    tmp_str_values = ','.join(values)

    tmp_str = key + ' -> ' + tmp_str_values
    tmp_list.append(tmp_str)

txt = '\n'.join(tmp_list)
f = open('1.5.txt','w')
f.write(txt)
f.close()