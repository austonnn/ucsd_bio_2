# -*- coding: utf-8 -*-
"""
Created on Sat Dec 23 12:41:25 2017

@author: xiangyin
"""


import random

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


file = open("dataset_203_7.txt")
k = file.readline().strip()
text = file.read().splitlines()

text1 = []
for item in text:
    item = item.strip()
    text1.append(item)


text = text1
#text = str(text[0])
file.close()


graph_dict = de_bruijn_pattern(text)
tmp_list = []
for key, values in graph_dict.items():
    if len(values) == 0:
        pass
    values = sorted(values)
    tmp_str_values = ','.join(values)

    tmp_str = key + ' -> ' + tmp_str_values
    tmp_list.append(tmp_str)

text = '\n'.join(tmp_list)
text = tmp_list







def eulerian_cycle(graph):
    #form a cycle Cycle by randomly walking in Graph (don't visit the same edge twice!)
    ans_cycle = []
    unexplored_graph = dict(graph)
    #tmp_key = random.choice(list(unexplored_graph.keys()))
    start_point = random.choice(list(unexplored_graph.keys()))
    #start_point = tmp_key
    tmp_start_point = start_point
    ans_cycle.append(tmp_start_point)
    #print(tmp_start_point)
    tmp_start_points = []
    #cycle_list = []
    tmp_cycle = []
    while len(unexplored_graph) > 0:
        #print(unexplored_graph)
        #
        unexplored_graph, tmp_end_point, tmp_start_points = update_graph(unexplored_graph, tmp_start_point, tmp_start_points)        

        
        
        tmp_start_point = tmp_end_point
        ans_cycle.append(tmp_end_point)

        if tmp_end_point == start_point:
            #cycle_list.append(tmp_cycle)
            if len(tmp_start_points) > 0:
                tmp_start_point = tmp_start_points.pop()
                tmp_cycle = get_new_cycle(ans_cycle, tmp_start_point)
                start_point = tmp_cycle[0]
                #tmp_cycle.pop()
                ans_cycle = list(tmp_cycle)

                
            
        
        # got ans_cycle
        
        #tmp_cycle  append cycle from new start
        
        
    """
        
        while there are unexplored edges in Graph
            select a node newStart in Cycle with still unexplored edges
            form Cycle’ by traversing Cycle (starting at newStart) and then randomly walking 
            Cycle ← Cycle’
    """
    #tmp_cycle.append(tmp_cycle[0])
    return ans_cycle
    #return graph


"""
def text_2_graph(text):
    graph_dict = dict()
    for item in text:
        key, value = item.strip().split(' -> ')
        graph_dict[key] = []
        for val in value.split(','):
            graph_dict[key].append(val)    
    return graph_dict
    """

def update_graph(graph, start_point, start_points):
    #print(graph)
    #print(start_points)

    end_point = random.choice(graph[start_point])
    graph[start_point].remove(end_point)
    if start_point in start_points:
        start_points.remove(start_point)
    if len(graph[start_point]) > 0:
        start_points.append(start_point)    
    if len(graph[start_point]) == 0:
        #if start_point in start_points:
            #start_points.remove(start_point)
        graph.pop(start_point)
    #if end_point in start_points:
        #start_points.remove(end_point)
        
    
    
    return graph, end_point, start_points
"""
def get_new_cycle(cycle, new_start):
    tmp_list = []
    index = cycle.index(new_start)
    tmp_list = cycle[index : -1] + cycle[:index]
    tmp_list.append(new_start)
    return tmp_list   
    """

def text_2_graph(text):
    graph_dict = dict()
    start_dict = dict()
    for item in text:
        key, value = item.strip().split(' -> ')
        graph_dict[key] = []
        # find out start
        #start_dict[key] = 1
        for val in value.split(','):
            graph_dict[key].append(val)
            # find out start
            if key in start_dict.keys():
                start_dict[key] += 1
            if key not in start_dict.keys():
                start_dict[key] = 1
            if val in start_dict.keys():
                start_dict[val] -= 1
            if val not in start_dict.keys():
                start_dict[val] = -1
    #graph_dict[]
    start_point = ''
    end_point = ''
    for key, val in start_dict.items():
        if val == 1:
            start_point = key
        if val == -1:
            end_point = key
    if end_point not in graph_dict.keys():
        graph_dict[end_point] = []
        
    graph_dict[end_point].append(start_point)
    return graph_dict, start_point, end_point

graph, start_point, end_point = text_2_graph(text)
#circuit = []
#print(eulerian_cycle(graph))
#node_i = random.choice(list(graph.keys()))
#find_circuit(node_i)
#print(circuit)

#circuit.reverse()
#ans = circuit
ans = eulerian_cycle(graph)


def get_new_cycle(cycle, start_point, end_point):
    tmp_list = []
    for index in range(len(cycle)):
        if cycle[index] == end_point and cycle[index + 1] == start_point:
            break
    #index = cycle.index(new_start)
    tmp_list = cycle[index + 1 : -1] + cycle[:index]
    tmp_list.append(end_point)
    return tmp_list
ans = get_new_cycle(ans, start_point, end_point)
text = ''
for item in ans[:-1]:
    text += item + '->'
text += ans[-1]    


def genome_path_to_string(genome_path):
    tmp_string = [] 
    tmp_string.append(genome_path[0])
    for text in genome_path[1:]:
        tmp_string.append(text[-1])
    ans = ''.join(tmp_string)
    return ans


ans = genome_path_to_string(ans) 
print(ans)
f = open('1.2.7.txt','w')
f.write(ans)
f.close()


