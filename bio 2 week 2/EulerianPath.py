import random
    
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


def eulerian_cycle(graph, start_point):
    global circuit
    circuit = []
    node_i = start_point
    #node_i = random.choice(list(graph.keys()))
    find_circuit(node_i)
    circuit.reverse()
    return circuit     
        
#file = open("dataset_203_6.txt")
#k = file.readline().strip()
#text = file.read().splitlines()
#text = str(text[0])
#file.close()

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
ans = eulerian_cycle(graph, start_point)


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
print(text)