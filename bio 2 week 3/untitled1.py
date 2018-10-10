#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 10:50:14 2018

@author: xiangyin
"""


def test(l1, l2):
    if not l1 or not l2: 
        print(not l1 or not l2)
        print(not l1)
        print(not l2)
        return l1 or l2

    
print(test([],[1,2]))

        