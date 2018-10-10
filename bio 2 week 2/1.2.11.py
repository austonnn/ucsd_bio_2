# -*- coding: utf-8 -*-
"""
Created on Mon Dec 25 16:08:37 2017

@author: xiangyin
"""

from itertools import product
k = 8
kmers = [ ''.join(x) for x in product('01', repeat=k) ]


import circular_string_problem


graph_dict = circular_string_problem.de_bruijn_pattern(kmers)

ans = circular_string_problem.eulerian_cycle(graph_dict)




ans = circular_string_problem.genome_path_to_string(ans) 

ans = ans[: len(ans) - k + 1]


f = open('1.2.11.txt','w')
f.write(ans)
f.close()