#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 23 15:31:55 2018

@author: xiangyin
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 12:58:13 2017

@author: xiangyin
"""
import random

# 1.3中的程序

def break_reads(patterns):
    tmp_pattern = []
    
    for item in patterns:
        item = item.split('|')        
        key = item[0][: -1] + '|' + item[1][: -1]
        tmp_pattern.append(key)
        value = item[0][1 :] + '|' + item[1][1 :]
        tmp_pattern.append(value)
    
    return list(set(tmp_pattern))

def break_reads_ng(patterns):
    tmp_pattern = []
    
    for item in patterns:
        #item = item.split('|')        
        key = item[: -1] #+ '|' + item[1][: -1]
        tmp_pattern.append(key)
        value = item[1 :] #+ '|' + item[1][1 :]
        tmp_pattern.append(value)
    
    return list(set(tmp_pattern))

def de_bruijn_pattern(patterns):
    graph_dict = dict()
    for item in patterns:
        #print(patterns)
        item = item.split('|')        
        key = item[0][: -1] + '|' + item[1][: -1]
        value = item[0][1 :] + '|' + item[1][1 :]        
        if key not in graph_dict.keys():
            graph_dict[key] = []
        graph_dict[key].append(value)    
    return graph_dict

def pattern_2_graph(patterns):
    graph_dict = dict()
    for item in patterns:
        item = item.split('|')        
        key = item[0][: -1]  + '|' + item[1][: -1]
        value = item[0][1 :]  + '|' + item[1][1 :]        
        if key not in graph_dict.keys():
            graph_dict[key] = []
        graph_dict[key].append(value)    
    return graph_dict

def pattern_2_graph_ng(patterns):
    graph_dict = dict()
    for item in patterns:
        #item = item.split('|')        
        key = item[: -1] # + '|' + item[1][: -1]
        value = item[1 :] # + '|' + item[1][1 :]        
        if key not in graph_dict.keys():
            graph_dict[key] = []
        graph_dict[key].append(value)    
    return graph_dict

def balance_graph(graph):
    start = ''
    end = ''
    tmp_dict = {}
    for key, values in graph.items():
        for val in values:
            if key in tmp_dict.keys():
                tmp_dict[key] += 1
            if key not in tmp_dict.keys():
                tmp_dict[key] = 1
            if val in tmp_dict.keys():
                tmp_dict[val] -= 1
            if val not in tmp_dict.keys():
                tmp_dict[val] = -1
            
    for key, val in tmp_dict.items():
        if val == 1:
            start = key
        if val == -1:
            end = key
    if end not in graph.keys():
        graph[end] = []
    graph[end].append(start)
        
    
    return graph, start, end

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

def get_new_cycle(cycle, start_point, end_point):
    tmp_list = []
    for index in range(len(cycle)):
        if cycle[index] == end_point and cycle[index + 1] == start_point:
            break
    #index = cycle.index(new_start)
    tmp_list = cycle[index + 1 : -1] + cycle[:index]
    tmp_list.append(end_point)
    return tmp_list








def string_spelled_gapped_patterns(gapped_patterns, k, d):
        #FirstPatterns ← the sequence of initial k-mers from GappedPatterns
        #SecondPatterns ← the sequence of terminal k-mers from GappedPatterns
        k = k - 1
        d = d + 1
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



# 1.6  中的程序
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
                        #cycle_list.append(tmp_key)
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

file = open("reads.txt")
#k, d = file.readline().strip().split()
k = 120 #int(k)
d = 1000 #int(d)
gapped_patterns = file.read().splitlines()
#text = str(text[0])
file.close()


#ans = string_spelled_gapped_patterns(gapped_patterns, k, d)

#graph = de_bruijn_pattern(gapped_patterns)

t = 100
for i in range(t):
    gapped_patterns = break_reads(gapped_patterns)
print('kmer is %s'%len(gapped_patterns[0]))
graph = pattern_2_graph(gapped_patterns)

max_paths = maximal_non_branching_paths(graph)
print('num of max path is %s'%len(max_paths))
# 输出结果
text = ''
for path in max_paths:
    for item in path:
        text += item[0] + ' -> '
    text += item[-1]    
    text += '\n'
#print(text)
f = open('1.4.10.txt','w')
f.write(text)
f.close()
#balanced_graph, start, end = balance_graph(graph)
"""
gapped_cycle = eulerian_cycle(balanced_graph)
new_gapped_cycle = get_new_cycle(gapped_cycle, start, end)
ans = string_spelled_gapped_patterns(new_gapped_cycle, k, d)
#print(ans)

file = open('1.4.10.txt','w')
file.write(ans)
file.close()
"""