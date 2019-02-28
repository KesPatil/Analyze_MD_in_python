#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 14:29:34 2018

@author: keshavpatil
"""

# This code reads the fes.dat file generated by PLUMED 
# and plots the contour plot for free energy, if there are only two
# collective variables.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


with open('fes.dat') as f:
    
    lines = f.readlines()
    
    
    
#lines = lines[9:2659]

new_lines = []    
# Get rid of empty lines
for line in lines:
    # Strip whitespace, should leave nothing if empty line was just "\n"
    if not line.strip():
        continue
    # We got something, save it
    else:
        new_lines.append(line)
        
lines = new_lines[9:]
    
bin_size_active = 302 # This was set in the PLUMED simulation
bin_size_inactive = 312
x = []
y1 = [] # stores the inactive coordinates
y2 = [] # stores the active coordinates
z1 = []
z = np.zeros((1,bin_size_inactive))
for i in range(0,len(lines)):
    
    words = lines[i].split()
    int_lst = [float(x) for x in words]
   
    y1.append(int_lst[0]) #y1 takes in inactive rmsd coordinates
    y2.append(int_lst[1]) #y2 takes in active rmsd coordinates
    z1.append(int_lst[2])  #z takes the free energy
 
            
y1 = np.array(y1)
y2 = np.array(y2)

y1 = y1[0:bin_size_inactive]
y2 = y2[0:len(lines):bin_size_inactive]
j = 1


# making a 2-D matrix z
for i in range(0,len(z1),bin_size_inactive):
    z = np.r_[z,[z1[i:i+bin_size_inactive]]] 


# deleting the first row of all zeros in z
z = np.delete(z, (0), axis=0)


x,y = np.meshgrid(y1,y2)


#plt.contourf(x,y,z-z.min(),cmaps = plt.cm.rainbow)
z = z - z.min()

# masking the matrix z 

                             
contours = plt.contour(y1, y2, z, 10, colors='black',linewidths=0.1)
z =np.ma.masked_where((14 < z) & (z < 100), z)
#plt.clabel(contours, inline=True, fontsize=8)
plt.imshow(z, extent=[min(y1), max(y1), min(y2), max(y2)], origin='lower',
           cmap='jet')#,ylabel='kkk')
plt.clim(0,14)
#plt.colorbar()
#plt.axis(aspect='image');


plt.xlabel('RMSD to inactive',fontsize=20)
plt.ylabel('RMSD to active',fontsize=20)

#plt.xticks(np.arange(np.min(y1), np.max(y1), 10))
#plt.yticks(np.arange(np.min(y2), np.max(y2), 10))                                
plt.colorbar()
#plt.title('HILLS',fontsize=20)
plt.savefig('HILLS.png', bbox_inches = 'tight', dpi = 300)

plt.show()
