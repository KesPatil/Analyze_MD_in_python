#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 15:19:48 2019

@author: keshav
"""

import numpy as np
import time
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import XTC, GRO
import MDAnalysis.analysis.rms
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.analysis import align

start = time.time()
active_ref1 = '/run/media/keshav/DELENDA_backup/alk_meta/meta_walkers_meta_struct/active_wt/s03-sim/active_wt_alk_CA_meta_min.pdb'
inactive_ref1 =  '/run/media/keshav/DELENDA_backup/alk_meta/meta_walkers_meta_struct/active_wt/s03-sim/inactive_wt_alk_CA_meta_min.pdb'
active_ref = MDAnalysis.Universe(active_ref1)
inactive_ref = MDAnalysis.Universe(inactive_ref1)


x1 = 1  # RMSD distance to active 
x2 = 1  # RMSD distance to inactive

tol_x1 = 0.15 # RMSD tolerance
tol_x2 = 0.15 # RMSD tolerance

b = 0 # dummy variable that's required
var2=  '/run/media/keshav/DELENDA_backup/alk_meta/meta_walkers_meta_struct/active_wt/s03-sim/md.part0004.gro'
u1   = mda.Universe(var2)
spe = u1.select_atoms("protein")
pos_protein = spe.positions
s = np.zeros(pos_protein.shape)
stackk = np.zeros(pos_protein.shape)
stackk = np.array([stackk])





# Looping through all the xtc files in the four folder one by one

# First folder : active_wt

for i in range(4,13,1):
    if i < 10:
        i = '0'+str(i)
        
    var1 = '/run/media/keshav/DELENDA_backup/alk_meta/meta_walkers_meta_struct/active_wt/s03-sim/md.part00'+str(i)+'.xtc'
    var2=  '/run/media/keshav/DELENDA_backup/alk_meta/meta_walkers_meta_struct/active_wt/s03-sim/md.part0004.gro'
    u   = mda.Universe(var2, var1)
    a   = u.select_atoms("protein and name CA")
    a_protein = u.select_atoms("protein")
    active_CA = active_ref.select_atoms("protein and name CA")
    activeCA  = active_CA.positions
    inactive_CA = inactive_ref.select_atoms("protein and name CA")
    inactiveCA = inactive_CA.positions
    for j in range(0,len(u.trajectory)):
        u.trajectory[j]
        mobileCA = a.positions
        mobile_protein = a_protein.positions
        r_active = rmsd(mobileCA, activeCA, superposition=True)
        r_inactive = rmsd(mobileCA, inactiveCA, superposition=True)
        if r_active >= x1 - tol_x1 and r_active <= x1 + tol_x1 and r_inactive >= x2 - tol_x2 and r_inactive <= x2 + tol_x2:
            s = s + mobile_protein
            b = b + 1
            stackk = np.concatenate((stackk,[mobile_protein]))
        
            

# Second folder : active_wt_replicate

for i in range(4,14,1):
    if i < 10:
        i = '0'+str(i)
    var1 = '/run/media/keshav/DELENDA_backup/alk_meta/meta_walkers_meta_struct/active_wt_replicate/s03-sim/md.part00'+str(i)+'.xtc'
    var2=  '/run/media/keshav/DELENDA_backup/alk_meta/meta_walkers_meta_struct/active_wt_replicate/s03-sim/md.part0006.gro'
    u   = mda.Universe(var2, var1)
    a   = u.select_atoms("protein and name CA")
    a_protein = u.select_atoms("protein")
    active_CA = active_ref.select_atoms("protein and name CA")
    activeCA  = active_CA.positions
    inactive_CA = inactive_ref.select_atoms("protein and name CA")
    inactiveCA = inactive_CA.positions
    for j in range(0,len(u.trajectory)):
        u.trajectory[j]
        mobileCA = a.positions
        mobile_protein = a_protein.positions
        r_active = rmsd(mobileCA, activeCA, superposition=True)
        r_inactive = rmsd(mobileCA, inactiveCA, superposition=True)
        if r_active >= x1 - tol_x1 and r_active <= x1 + tol_x1 and r_inactive >= x2 - tol_x2 and r_inactive <= x2 + tol_x2:
            s = s + mobile_protein
            b = b + 1
            stackk = np.concatenate((stackk,[mobile_protein]))
            
            

# Third folder : inactive_wt

for i in range(4,15,1):
    if i < 10:
        i = '0'+str(i)
    var1 = '/run/media/keshav/DELENDA_backup/alk_meta/meta_walkers_meta_struct/inactive_wt/s03-sim/md.part00'+str(i)+'.xtc'
    var2=  '/run/media/keshav/DELENDA_backup/alk_meta/meta_walkers_meta_struct/inactive_wt/s03-sim/md.part0006.gro'
    u   = mda.Universe(var2, var1)
    a   = u.select_atoms("protein and name CA")
    a_protein = u.select_atoms("protein")
    active_CA = active_ref.select_atoms("protein and name CA")
    activeCA  = active_CA.positions
    inactive_CA = inactive_ref.select_atoms("protein and name CA")
    inactiveCA = inactive_CA.positions
    for j in range(0,len(u.trajectory)):
        u.trajectory[j]
        mobileCA = a.positions
        mobile_protein = a_protein.positions
        r_active = rmsd(mobileCA, activeCA, superposition=True)
        r_inactive = rmsd(mobileCA, inactiveCA, superposition=True)
        if r_active >= x1 - tol_x1 and r_active <= x1 + tol_x1 and r_inactive >= x2 - tol_x2 and r_inactive <= x2 + tol_x2: 
            s = s + mobile_protein
            b = b + 1
            stackk = np.concatenate((stackk,[mobile_protein]))
            
            
# Fourth folder : inactive_wt_replicate

for i in range(3,13,1):
    if i < 10:
        i = '0'+str(i)
    var1 = '/run/media/keshav/DELENDA_backup/alk_meta/meta_walkers_meta_struct/inactive_wt_replicate/s03-sim/md.part00'+str(i)+'.xtc'
    var2=  '/run/media/keshav/DELENDA_backup/alk_meta/meta_walkers_meta_struct/inactive_wt_replicate/s03-sim/md.part0006.gro'
    u   = mda.Universe(var2, var1)
    a   = u.select_atoms("protein and name CA")
    a_protein = u.select_atoms("protein")
    active_CA = active_ref.select_atoms("protein and name CA")
    activeCA  = active_CA.positions
    inactive_CA = inactive_ref.select_atoms("protein and name CA")
    inactiveCA = inactive_CA.positions
    for j in range(0,len(u.trajectory)):
        u.trajectory[j]
        mobileCA = a.positions
        mobile_protein = a_protein.positions
        r_active = rmsd(mobileCA, activeCA, superposition=True)
        r_inactive = rmsd(mobileCA, inactiveCA, superposition=True)
        if r_active >= x1 - tol_x1 and r_active <= x1 + tol_x1 and r_inactive >= x2 - tol_x2 and r_inactive <= x2 + tol_x2: 
            s = s + mobile_protein
            b = b + 1
            stackk = np.concatenate((stackk,[mobile_protein]))

s_avg = s/b # this stores the average coordinate file = this is fictitious
stackk = np.delete(stackk, (0), axis=0)
dist = (stackk - s_avg)**2

q = []
for i in range(0,b):
    nor = np.linalg.norm(dist[i][0][0])
    q.append(nor)

index_req = q.index(min(q))
stack_req = dist[index_req][0][0]

end = time.time()
print(end - start)