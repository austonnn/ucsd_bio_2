#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 09:12:58 2018

@author: xiangyin
"""

import matplotlib.pyplot as plt
import numpy as np


def firstn(n):
    num1, num2 = 0, 1
    count = 0
    while count < n:
        yield num1
        tmp_num = num1 + num2
        num1 = num2
        num2 = tmp_num
        count += 1
for item in firstn(10):
    print(item)

x = np.arange(0., 6., 0.1)
y = 3* x**4 - 8 * x**3 + 6 * x**2 - 12
plt.plot(x, y)
plt.show()


def y1(x):
    return (x -3) ** 2 / 60
tmp_sum = 0
for item in x:
    tmp_sum += y1(item)
print(tmp_sum)
    

x**3 / 18 - 9