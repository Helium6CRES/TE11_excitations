#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 13:49:30 2019
@author: brent
"""

import numpy as np
import matplotlib.pyplot as plt

def readFloats(pathToFile):

    with open(pathToFile) as f:
        a = [float(x) for x in f.read().split()]
    return a

data_lo = readFloats("./2T_18GHz.txt")
data_hi = readFloats("./2T_24GHz.txt")

print(len(data_lo))
print(len(data_hi))

for i in range(len(data_lo)):
    data_lo[i] = data_lo[i] *5.0E14

for i in range(len(data_hi)):
    data_hi[i] = data_hi[i] *5.0E14

plt.figure(figsize=(12,6))    
num_bins = 50
bins = np.linspace(0.0, 5, 100)
plt.hist(data_hi, bins, alpha=0.75, label='24 GHz, N = 18317')
plt.hist(data_lo, bins, alpha=0.6, facecolor = 'red', label='18 GHz, N = 11996')
plt.legend(loc='upper left', prop={'size': 20})
plt.tick_params(axis='both', which='major', labelsize=20)
plt.xlabel('Power (fW)', fontsize = 20)
plt.ylabel('Counts', fontsize = 20)
plt.title(r'2.0 Tesla Power Distribution', fontsize = 20)
plt.subplots_adjust(left=0.15)
plt.show()
