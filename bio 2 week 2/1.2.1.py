# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 20:43:34 2017

@author: xiangyin
"""
import random

circuit_pos = 0
def find_circuit(node_i):
    if len(graph[node_i]) == 0:
        #graph.pop(node_i)
        circuit.append(node_i)
        graph.pop(node_i)
        #circuit[circuit_pos] = node_i
        #circuit_pos += 1
    else:
        while node_i in graph.keys() and len(graph[node_i]) > 0:
            node_j = random.choice(graph[node_i])
            #delete_edges(node_j, node_i)
            graph[node_i].remove(node_j)
            #if len(graph[node_i]) == 0:
                #graph.pop(node_i)
            #print(graph)
            find_circuit(node_j)
            #tmp_circuit = find_circuit(node_j)
            #for item in tmp_circuit:
                #circuit.append(item)
        circuit.append(node_i)
        #circuit[circuit_pos] = node_i
        #circuit_pos += 1
    #return circuit


def eulerian_cycle(graph):
    global circuit
    circuit = []
    node_i = random.choice(list(graph.keys()))
    find_circuit(node_i)
    circuit.reverse()
    return circuit     
"""        
file = open("a.txt")
#k = file.readline().strip()
text = file.read().splitlines()
#text = str(text[0])
file.close()
"""

with open('a.txt', 'r') as read_handle:
    text = read_handle.read().splitlines()

def text_2_graph(text):
    graph_dict = dict()
    for item in text:
        key, value = item.strip().split(' -> ')
        graph_dict[key] = []
        for val in value.split(','):
            graph_dict[key].append(val)    
    return graph_dict

graph = text_2_graph(text)
#circuit = []
#print(eulerian_cycle(graph))
#node_i = random.choice(list(graph.keys()))
#find_circuit(node_i)
#print(circuit)

#circuit.reverse()
#ans = circuit
ans = eulerian_cycle(graph)
text = ''
for item in ans[:-1]:
    text += item + '->'
text += ans[-1]    
print(text)
f = open('1.2.txt','w')
f.write(text)
f.close()