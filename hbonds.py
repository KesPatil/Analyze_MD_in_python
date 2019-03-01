#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 17:06:02 2019

@author: keshavpatil
"""

# This code computes the H-bond occupancy
# Here we are considering only the H-bonds from the 
# Has to be upgraded to use pandas

import MDAnalysis
import MDAnalysis.analysis.hbonds


var1 = '/Users/keshavpatil/Desktop/BRAF/100ns_unbiased/active_braf/unb100.xtc'
var2 = '/Users/keshavpatil/Desktop/BRAF/100ns_unbiased/active_braf/unb100.gro'
u = MDAnalysis.Universe(var2, var1)
hbonds = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u, 'protein', 'protein', distance=3.0, angle=120.0)
hbonds.run()


hbonds_occupancy=hbonds.count_by_type()
hbonds_occupancy_list=hbonds_occupancy.tolist()

# get the specific H-bonds:
# By specific I mean: the H-bonds whose either the 
# donor_atom or the acceptor_atom or both belong to 
# alpha-helix or the activation loop
# This is ofcourse dependent on the 
# type of kinase. The github repository which
# has the info on this for somewhere around hundreds of 
# kinases is: https://github.com/ejjordan/mutator/blob/master/all_kinase.csv
# Take this input from the user
print('Enter the threshold for the h-bond occupany')
threshold = input()
print('Enter the beginning residue number for alpha-helix of your kinase')
begin_alpha = input()
print('Enter the ending residue number for alpha-helix of your kinase')
end_alpha = input()
print('Enter the beginning residue number for catalytic loop of your kinase')
begin_catloop = input()
print('Enter the ending residue number for catalytic loop of your kinase')
end_catloop = input()

# accumulator stores the h-bond occupancy of all the residues that 
# you are concerned with.
accumulator  = 0
donor_resid = hbonds_occupancy.donor_resid
acceptor_resid = hbonds_occupancy.acceptor_resid

for i in range(0,len(donor_resid)):
    if (begin_alpha <= donor_resid[i] <= end_alpha or 
        begin_catloop <= donor_resid[i] <= end_alpha or
        begin_alpha <= acceptor_resid[i] <= end_alpha or 
        begin_catloop <= acceptor_resid[i] <= end_alpha and 
        hbonds_occupancy.frequency[i] > threshold):
        # then add the occupancy to the accumulator
        accumulator = accumulator + hbonds_occupancy.frequency[i]
    




