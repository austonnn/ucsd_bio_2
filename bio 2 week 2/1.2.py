# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 13:55:42 2017

@author: xiangyin
"""


import random


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

file = open("dataset_203_2.txt")
#k = file.readline().strip()
text = file.read().splitlines()
#text = str(text[0])
file.close()

def text_2_graph(text):
    graph_dict = dict()
    for item in text:
        key, value = item.strip().split(' -> ')
        graph_dict[key] = []
        for val in value.split(','):
            graph_dict[key].append(val)    
    return graph_dict

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

def get_new_cycle(cycle, new_start):
    tmp_list = []
    index = cycle.index(new_start)
    tmp_list = cycle[index : -1] + cycle[:index]
    tmp_list.append(new_start)
    return tmp_list
        
    

graph = text_2_graph(text)
#print(eulerian_cycle(graph))
ans = eulerian_cycle(graph)
text = ''
for item in ans[:-1]:
    text += item + '->'
text += ans[-1]    
print(text)
f = open('1.2.txt','w')
f.write(text)
f.close()