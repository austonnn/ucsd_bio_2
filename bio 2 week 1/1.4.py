# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 10:23:02 2017

@author: xiangyin
"""

def de_bruijn_graph(text, k):
    k = int(k)
    graph_dict = dict()
    for index in range(0 , len(text) - k + 1):
        tmp_text = text[index : index + k]

        key = tmp_text[: -1]
        value = tmp_text[1 :]
        
        if key in graph_dict.keys():
            graph_dict[key].append(value)
        else:
            graph_dict[key] = []
            graph_dict[key].append(value)
    
    return graph_dict





file = open("dataset_199_6.txt")
k = file.readline().strip()
text = file.read().splitlines()
text = str(text[0])
file.close()



#print(de_bruijn_graph(text, k))

graph_dict = de_bruijn_graph(text, k)
tmp_list = []
for key, values in graph_dict.items():
    if len(values) == 0:
        pass
    values = sorted(values)
    tmp_str_values = ','.join(values)

    tmp_str = key + ' -> ' + tmp_str_values
    tmp_list.append(tmp_str)

txt = '\n'.join(tmp_list)
f = open('1.4.txt','w')
f.write(txt)
f.close()