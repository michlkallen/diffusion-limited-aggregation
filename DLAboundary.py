# -*- coding: utf-8 -*-
"""
Created: 2017-10-31
Script to generate and draw boundary circle for DLA.py
"""

import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

#load the text file containing the centers
#in this case the file is located on the D: drive, change according to your loc
cntrs = np.loadtxt('D:/centers.txt')
n = int(np.size(cntrs)/2) #number of particles
r = 1.0 #particle size from DLA file

#determine the convex hull for the data set, use to find bounding points
hull = ConvexHull(cntrs)
hullSize = np.size(hull.simplices[:,0])

hullInd = np.unique(np.ravel(hull.simplices[:]))
hullPts = cntrs[hullInd]

#find distances from hull points to center
rHull = cdist(hullPts,np.zeros((1,2)))

#set first bounding circle point to max hull radius
p1 = hullPts[[np.argmax(rHull)]]

#set second bounding circle point to max dist between 1st and other hull pts
distP1 = cdist(hullPts,p1)
p2 = hullPts[[np.argmax(distP1)]]

#set third bounding circle point to max dist between midpoint and other hull pts
distPav = cdist(hullPts,(p1+p2)/2)
p3 = hullPts[[np.argmax(distPav)]]

#reorient point arrays for convenience
p1,p2,p3 = p1.ravel(), p2.ravel(), p3.ravel()

#find center of bounding circle
A = np.array([[2*(p1[0] - p2[0]), 2*(p1[1] - p2[1])],
               [2*(p1[0] - p3[0]), 2*(p1[1] - p3[1])]])
b = np.array([[(p1[0]**2 - p2[0]**2) + (p1[1]**2 - p2[1]**2)],
               [(p1[0]**2 - p3[0]**2) + (p1[1]**2 - p3[1]**2)]])

x = np.linalg.solve(A,b) #center

#bounding radius
rBound = np.sqrt((p1[0] - x[0])**2 + (p1[1] - x[1])**2)

#PLOT
plt.clf()

fig = plt.figure()
ax = fig.add_subplot(111)
patches = []

#plot boundary circle patch
bound = Circle((x[0],x[1]), radius = rBound + 1.0, fill = False)
ax.add_patch(bound)

#plot the particles
for i in range(n):
    circle = Circle((cntrs[i]), radius = r)
    patches.append(circle)

p = PatchCollection(patches, color = 'black')
ax.add_collection(p)
plt.axis('off')
plt.axis('equal')
plt.show()