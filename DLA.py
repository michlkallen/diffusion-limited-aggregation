# -*- coding: utf-8 -*-
"""
Created: 2017-10-20
Script to generate off-lattice diffusion-limited aggregates.
Currently only uses identical particle sizes, but could be extended.
Tested up to 10k particles
"""
import numpy as np
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

r = 1.0 #min particle radius
l = 1.0 #step length
n = 200 #number of particles desired

#initialize arrays to hold center points, aggregate radius, dist btwn particles
centers = np.zeros((n,2)) #x,y
rMax = np.zeros((n,1)) #aggregate radius

space = np.zeros((n,1)) #dist btwn current particle and others

#initialize the starting position, starting gap and index
rCurrent = 2.5*r
gapCurrent = 2.5*r
gapInd = 0

for i in range(1,n):
    rRelease = rCurrent + 5.0*r #release radius
    theta = 2*np.pi*np.random.random()
    posCurrent = rRelease*np.array([[np.cos(theta),np.sin(theta)]]) #release pt  
    centers[i,:] = posCurrent #update particle center to current position
    
    space = cdist(centers[:i,],posCurrent)
    
    gapCurrent = np.amin(space[:i]) #calculate current min distance
    gapInd = np.argmin(space[:i]) #current min distance index
    
    while gapCurrent > (2*r):
        theta = 2*np.pi*np.random.random()
        step = l*np.array([[np.cos(theta), np.sin(theta)]])
        posCurrent += step
        rCurrent = cdist(posCurrent,centers[[0,]])
        
        #reset distance if too far outside of range (kill radius)
        if rCurrent > 2.0*rRelease:
            posCurrent = rRelease*np.array([[np.cos(theta),np.sin(theta)]])
            continue
        
        space = cdist(centers[:i,],posCurrent)
    
        gapCurrent = np.amin(space[:i]) #calculate new min distance
        gapInd = np.argmin(space[:i]) #new min distance index
        
        #is there intersection?
        if gapCurrent <= (2*r):
            print('Found intersection ', i)
            p2 = posCurrent #current step location
            angle = l/r*step #step vector direction
            p1 = p2 - angle #prior step location
            
            #end point adjustment to correct location without overlap
            Asq = np.dot((p2[0,] - p1[0,]), (p2[0,] - p1[0,]))
            Bsq = np.dot((p2[0,] - centers[gapInd,]), (p2[0,] - centers[gapInd,]))
            AB = np.dot((p2[0,] - p1[0,]), (p2[0,] - centers[gapInd,]))
            alpha = np.roots([Asq, 2*AB, (Bsq - (2*r)**2)])
            
            #two pontetial corrected center locations (from the 2 roots)
            pf1 = p2 + alpha[0]*angle
            pf2 = p2 + alpha[1]*angle
            
            #select location closest to the prior step location
            if cdist(pf1,p1) > cdist(pf2,p1):
                pf = pf2
            else:
                pf = pf1
            
            #update the various arrays
            centers[i] = pf
            rMax[i] = cdist(pf,centers[[0,]])
            rCurrent = np.amax(rMax)
            break

#save text file of centers for external plotting or for other scripts
#np.savetxt('/centers.txt',centers,delimiter = ' ')

#PLOT
plt.clf()

fig = plt.figure()
ax = fig.add_subplot(111)
patches = []

#draw circle patch at each particle center, adjust radius if you want overlap
for i in range(n):
    circle = Circle((centers[i]), radius = r)
    patches.append(circle)

p = PatchCollection(patches, color = 'black')
ax.add_collection(p)
plt.axis('off')
plt.axis('equal')
plt.show()