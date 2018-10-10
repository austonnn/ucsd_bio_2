# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 15:33:05 2017

@author: xiangyin
"""

def composition(text, k):
    tmp_list = []
    for i in range(len(text) - k + 1):
        #print(text[i : i + k])
        tmp_list.append(text[i : i + k])
    return tmp_list

file = open("dataset_197_3.txt")
k = file.readline().strip()#
text = file.read().strip()#
file.close()

k = int(k)


word_list = composition(text, k)
word_list = sorted(word_list)

#for item in word_list:
#    print(item)
        
#file_object = open('thefile.txt', 'w')
#file_object.writelines(word_list)
#file_object.close( )

txt = '\n'.join(word_list)
f = open('data.txt','w')
f.write(txt)
f.close()

























