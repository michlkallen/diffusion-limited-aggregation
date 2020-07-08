# -*- coding: utf-8 -*-
"""
Script to generate off-lattice diffusion-limited aggregates.
Currently only uses identical particle sizes, but could be extended.
Particles are release in spherically random 3D space.
Tested up to 100 particles. Could view by generating OpenSCAD script, for example.

Functions
---------
aggregate3D : base function - use to generate a 3D aggregate with options for
    plotting and exporting the particle centers as a text file

Typical Method
--------------
>>> from dla_3d import *
>>> aggregate3D(50)  # generate 50 particle aggregate
"""
import numpy as np
from scipy.spatial.distance import cdist


# --------------------------------------------------------------------------------------
def aggregate3D(n, r=1.0, len=1.0):
    """
    Generates 3D diffusion-limited aggregate with 'n' particles.

    aggregate3D(n, r=1.0, len=1.0)

    Parameters
    ----------
    n : number of particles in the aggregate
    r : particle radius
    len : step length (recommend keeping the same as 'r')

    Example
    -------
    aggregate3D(200)
    """

    centers = np.zeros((n, 3))  # x, y, z
    r_max = np.zeros((n, 1))

    space = np.zeros((n, 1))

    r_current = 2.5*r
    gap_current = 2.5*r
    gap_idx = 0

    for i in range(1, n):
        r_release = r_current + 5.0*r  # release radius
        theta = 2*np.pi*np.random.random()
        phi = np.arccos(2*np.random.random() - 1)
        pos_current = r_release*np.array([[np.cos(theta)*np.sin(phi),
                                         np.sin(theta)*np.sin(phi),
                                         np.cos(phi)]])  # release pt
        centers[i, :] = pos_current  # update particle center to current position

        space = cdist(centers[:i, ], pos_current)

        gap_current = np.amin(space[:i])  # calculate current min distance
        gap_idx = np.argmin(space[:i])  # current min distance index

        while gap_current > (2*r):
            theta = 2*np.pi*np.random.random()
            phi = np.arccos(2*np.random.random() - 1)
            step = len*np.array([[np.cos(theta)*np.sin(phi),
                                np.sin(theta)*np.sin(phi),
                                np.cos(phi)]])
            pos_current += step
            r_current = cdist(pos_current, centers[[0, ]])

            # reset distance if too far outside of range (kill radius)
            if r_current > 2.0*r_release:
                pos_current = r_release*np.array([[np.cos(theta)*np.sin(phi),
                                                   np.sin(theta)*np.sin(phi),
                                                   np.cos(phi)]])
                continue

            space = cdist(centers[:i, ], pos_current)

            gap_current = np.amin(space[:i])  # calculate new min distance
            gap_idx = np.argmin(space[:i])  # new min distance index

            # is there intersection?
            if gap_current <= (2*r):
                print('Found intersection ', i)
                # end point of step
                p2 = pos_current
                # step vector direction
                angle = len/r*step
                # prior step location
                p1 = p2 - angle
                # end point adjustment to correct location without overlap
                Asq = np.dot((p2[0, ] - p1[0, ]),
                             (p2[0, ] - p1[0, ]))
                Bsq = np.dot((p2[0, ] - centers[gap_idx, ]),
                             (p2[0, ] - centers[gap_idx, ]))
                AB = np.dot((p2[0, ] - p1[0, ]),
                            (p2[0, ] - centers[gap_idx, ]))
                alpha = np.roots([Asq, 2*AB, (Bsq - (2*r)**2)])

                # potential final adjusted center
                pf1 = p2 + alpha[0]*angle
                pf2 = p2 + alpha[1]*angle

                # select location closest to the prior step location
                if cdist(pf1, p1) > cdist(pf2, p1):
                    pf = pf2
                else:
                    pf = pf1

                centers[i] = pf
                r_max[i] = cdist(pf, centers[[0, ]])
                r_current = np.amax(r_max)
                break

    # save text file of center locations for external use
    np.savetxt('3d_centers.txt', centers, delimiter=' ')
