# -*- coding: utf-8 -*-
"""
Connectivity diagram for aggregation files generated via DLA.py
"""

import types
import numpy as np
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
from matplotlib.backend_bases import GraphicsContextBase, RendererBase
from matplotlib.collections import LineCollection

#load the text file containing the centers
#in this case the file is located on the D: drive, change according to your loc
cntrs = np.loadtxt('D:/centers.txt')
n = int(np.size(cntrs)/2) #number of particles
r = 1.0 #particle size from DLA file

#calculate the distance between all particle centers
distMatrix = np.round(cdist(cntrs,cntrs),1)

#determine neighbor indices
ind = np.argwhere(distMatrix == 2.)

#the following correct the endcap style of the line segments to rounded
class GC(GraphicsContextBase):
    def __init__(self):
        super().__init__()
        self._capstyle = 'round'

def custom_new_gc(self):
    return GC()

RendererBase.new_gc = types.MethodType(custom_new_gc, RendererBase)

#PLOT
plt.clf()

fig = plt.figure()
ax = fig.add_subplot(111)
lines = []

for i in range(int(np.size(ind)/2)):
    seg = np.array([cntrs[ind[i,0]],
                    cntrs[ind[i,1]]])
    lines.append(seg)

p = LineCollection(lines, linestyle = 'solid',
                   linewidth = 0.2, color = 'black')
ax.add_collection(p)
plt.axis('off')
plt.axis('equal')
plt.show()