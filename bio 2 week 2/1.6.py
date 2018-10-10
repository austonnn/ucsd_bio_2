# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 10:34:40 2017

@author: xiangyin
"""


def text_2_graph(text):
    graph_dict = dict()
    for item in text:
        key, value = item.strip().split(' -> ')
        graph_dict[key] = []
        # find out start
        #start_dict[key] = 1
        for val in value.split(','):
            graph_dict[key].append(val)
            # find out start
    return graph_dict




def get_in_out_dict(graph):
    in_out_dict = dict()
    for key, value in graph.items():
        for val in value:
            if key in in_out_dict.keys():
                in_out_dict[key][0] += 1
            if key not in in_out_dict.keys():
                in_out_dict[key] = [1, 0]
            if val in in_out_dict.keys():
                in_out_dict[val][1] += 1
            if val not in in_out_dict.keys():
                in_out_dict[val] = [0, 1]    
    
    return in_out_dict


def maximal_non_branching_paths(graph):
    paths = []
    in_out_dict =get_in_out_dict(graph)
    #print(in_out_dict )
    """
    Paths ← empty list
        for each node v in Graph
            if v is not a 1-in-1-out node
                if out(v) > 0
                    for each outgoing edge (v, w) from v
                        NonBranchingPath ← the path consisting of the single edge (v, w)
                        while w is a 1-in-1-out node
                            extend NonBranchingPath by the outgoing edge (w, u) from w 
                            w ← u
                        add NonBranchingPath to the set Paths
        for each isolated cycle Cycle in Graph
            add Cycle to Paths
        return Paths
    """
    cycle_list = []
    for key, val in in_out_dict.items():
        
        if val != [1, 1]:
            if val[0] > 0:
                for out_going_edge in graph[key]:
                    non_branching_path = []
                    non_branching_path.append([key, out_going_edge])
                    tmp_key = out_going_edge                   
                    while in_out_dict[tmp_key] == [1, 1]:
                        for tmp_out_edge in graph[tmp_key]:                       
                            non_branching_path.append([tmp_key, tmp_out_edge])
                        tmp_key = tmp_out_edge
                    paths.append(non_branching_path)
        #if val == [1, 1]:
        else:
            if key not in cycle_list:

                out_going_edge = graph[key][0]
                non_branching_path = []
                non_branching_path.append([key, out_going_edge])
                tmp_key = out_going_edge
    
                while in_out_dict[tmp_key] == [1, 1]:
                    tmp_out_edge = graph[tmp_key][0]
                    non_branching_path.append([tmp_key, tmp_out_edge])
                    tmp_key = tmp_out_edge
                    if tmp_key == key:
                        paths.append(non_branching_path)                        
                        for item in non_branching_path:
                            cycle_list.append(item[0])                            
                        break
        
                    
                        
    
    
    return paths



file = open("a.txt")
#k = file.readline().strip()
text = file.read().splitlines()

text1 = []
for item in text:
    item = item.strip()
    text1.append(item)


text = text1
#text = str(text[0])
file.close()


graph= text_2_graph(text)
#in_out_dict =in_out_graph(graph)
max_paths = maximal_non_branching_paths(graph)

# 输出结果
text = ''
for path in max_paths:
    for item in path:
        text += item[0] + ' -> '
    text += item[-1]    
    text += '\n'
#print(text)
f = open('1.6.txt','w')
f.write(text)
f.close()