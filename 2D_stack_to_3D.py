#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  3 06:50:51 2019

@author: keshavpatil
"""

import numpy as np

a = [[1,2,3],[4,5,6]]
b = [[7,8,9],[5,9,1]]
a = np.array(a)
b = np.array(b)
s = np.zeros(a.shape)
s = np.array([s])
for i in range(5):
    s = np.concatenate((s,[a]))

s = np.delete(s, (0), axis=0)
print(s)    
    
    