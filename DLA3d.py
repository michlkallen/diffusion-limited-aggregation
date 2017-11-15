# -*- coding: utf-8 -*-
"""
Created: 2017-10-20
Script to generate off-lattice diffusion-limited aggregates.
Currently only uses identical particle sizes, but could be extended.
Format is the same as DLA.py except addition of phi angle to properly sample
3D space (spherically random)
Tested up to 100 particles
"""

import numpy as np
from scipy.spatial.distance import cdist

r = 1.0 #min particle radius
l = 1.0 #step length
n = 100 #number of particles desired

centers = np.zeros((n,3)) #x,y,z
rMax = np.zeros((n,1))

space = np.zeros((n,1))

rCurrent = 2.5*r
gapCurrent = 2.5*r
gapInd = 0

for i in range(1,n):
    rRelease = rCurrent + 5.0*r #release radius
    theta = 2*np.pi*np.random.random()
    phi = np.arccos(2*np.random.random()-1)
    posCurrent = rRelease*np.array([[np.cos(theta)*np.sin(phi),
                                     np.sin(theta)*np.sin(phi),
                                     np.cos(phi)]]) #release pt  
    centers[i,:] = posCurrent #update particle center to current position
    
    space = cdist(centers[:i,],posCurrent)
    
    gapCurrent = np.amin(space[:i]) #calculate current min distance
    gapInd = np.argmin(space[:i]) #current min distance index
    
    while gapCurrent > (2*r):
        theta = 2*np.pi*np.random.random()
        phi = np.arccos(2*np.random.random()-1)
        step = l*np.array([[np.cos(theta)*np.sin(phi),
                                     np.sin(theta)*np.sin(phi),
                                     np.cos(phi)]])
        posCurrent += step
        rCurrent = cdist(posCurrent,centers[[0,]])
        
        #reset distance if too far outside of range (kill radius)
        if rCurrent > 2.0*rRelease:
            posCurrent = rRelease*np.array([[np.cos(theta)*np.sin(phi),
                                     np.sin(theta)*np.sin(phi),
                                     np.cos(phi)]])
            continue
        
        space = cdist(centers[:i,],posCurrent)
    
        gapCurrent = np.amin(space[:i]) #calculate new min distance
        gapInd = np.argmin(space[:i]) #new min distance index
        
        #is there intersection?
        if gapCurrent <= (2*r):
            print('Found intersection ',i)
            #end point of step
            p2 = posCurrent
            #step vector direction
            angle = l/r*step
            #prior step location
            p1 = p2 - angle
            #end point adjustment to correct location without overlap
            Asq = np.dot((p2[0,] - p1[0,]), (p2[0,] - p1[0,]))
            Bsq = np.dot((p2[0,] - centers[gapInd,]), (p2[0,] - centers[gapInd,]))
            AB = np.dot((p2[0,] - p1[0,]), (p2[0,] - centers[gapInd,]))
            alpha = np.roots([Asq, 2*AB, (Bsq - (2*r)**2)])
            
            #potential final adjusted center
            pf1 = p2 + alpha[0]*angle
            pf2 = p2 + alpha[1]*angle
            
            #select location closest to the prior step location
            if cdist(pf1,p1) > cdist(pf2,p1):
                pf = pf2
            else:
                pf = pf1
            
            centers[i] = pf
            rMax[i] = cdist(pf,centers[[0,]])
            rCurrent = np.amax(rMax)
            break

#save text file of center locations for external use
#np.savetxt('/centers.txt',centers,delimiter = ' ')