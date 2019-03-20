#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 13:49:30 2019
@author: brent
"""
import matplotlib.pyplot as plt

def readFloats(pathToFile):

    with open(pathToFile) as f:
        a = [float(x) for x in f.read().split()]
    return a

data = readFloats("./2T_18GHz.txt")

for i in range(len(data)):
    data[i] = data[i] *1.0E15

plt.figure(figsize=(20,10))    
num_bins = 50
plt.hist(data, num_bins, alpha=0.75)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.xlabel('Power (fW)', fontsize = 20)
plt.ylabel('Counts', fontsize = 20)
plt.title(r'Histogram', fontsize = 30)
plt.subplots_adjust(left=0.15)
plt.show()