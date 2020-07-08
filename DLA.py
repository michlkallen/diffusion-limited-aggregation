# -*- coding: utf-8 -*-
"""
Scripts to generate off-lattice diffusion-limited aggregates and related information.
Currently only uses identical particle sizes, but could be extended.
Tested up to 10k particles

Functions
---------
aggregate2D : base function - use to generate a 2D aggregate with options for
    plotting and exporting the particle centers as a text file (used in other functions)
bounded_aggregate : plot the aggregate and its bounding circle
    (uses text file generated from aggregate2D)
aggregate_conn : plot connectivity diagram for a given aggregate
    (uses text file generated from aggregate2D)

Typical Method
--------------
>>> from dla import *
>>> aggregate2D(1000)  # generate 1000 particle aggregate
>>> bounded_aggregate()  # plot aggregate + bounding circle for newly created aggregate
>>> aggregate_conn()  # plot connectivity diagram for newly created aggregate
"""
import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial import ConvexHull
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection, LineCollection


# --------------------------------------------------------------------------------------
def aggregate2D(n, r=1.0, len=1.0, output_text=True, plot=True):
    """
    Generates 2D diffusion-limited aggregate with 'n' particles.

    aggregate2D(n, r=1.0, len=1.0, output_text=True, plot=True)

    Parameters
    ----------
    n : number of particles in the aggregate
    r : particle radius
    len : step length (recommend keeping the same as 'r')
    output_text : whether to save a text file of the particle centers
        for external plotting or for other scripts
    plot : whether to plot the data directly

    Example
    -------
    aggregate2D(200, output_text=False)
    """

    # initialize arrays to hold center points, aggregate radius, dist btwn particles
    centers = np.zeros((n, 2))  # x, y
    r_max = np.zeros((n, 1))  # aggregate radius

    space = np.zeros((n, 1))  # dist btwn current particle and others

    # initialize the starting position, starting gap and index
    r_current = 2.5*r
    gap_current = 2.5*r
    gap_idx = 0

    for i in range(1, n):
        r_release = r_current + 5.0*r  # release radius
        theta = 2*np.pi*np.random.random()
        pos_current = r_release*np.array([[np.cos(theta), np.sin(theta)]])  # release pt
        centers[i, :] = pos_current  # update particle center to current position

        space = cdist(centers[:i, ], pos_current)

        gap_current = np.amin(space[:i])  # calculate current min distance
        gap_idx = np.argmin(space[:i])  # current min distance index

        while gap_current > (2*r):
            theta = 2*np.pi*np.random.random()
            step = len*np.array([[np.cos(theta), np.sin(theta)]])
            pos_current += step
            r_current = cdist(pos_current, centers[[0, ]])

            # reset distance if too far outside of range (kill radius)
            if r_current > 2.0*r_release:
                pos_current = r_release*np.array([[np.cos(theta), np.sin(theta)]])
                continue

            space = cdist(centers[:i, ], pos_current)

            gap_current = np.amin(space[:i])  # calculate new min distance
            gap_idx = np.argmin(space[:i])  # new min distance index

            # is there intersection?
            if gap_current <= (2*r):
                print('Found intersection ', i)
                p2 = pos_current  # current step location
                angle = len/r*step  # step vector direction
                p1 = p2 - angle  # prior step location

                # end point adjustment to correct location without overlap
                Asq = np.dot((p2[0, ] - p1[0, ]),
                             (p2[0, ] - p1[0, ]))
                Bsq = np.dot((p2[0, ] - centers[gap_idx, ]),
                             (p2[0, ] - centers[gap_idx, ]))
                AB = np.dot((p2[0, ] - p1[0, ]),
                            (p2[0, ] - centers[gap_idx, ]))
                alpha = np.roots([Asq, 2*AB, (Bsq - (2*r)**2)])

                # two pontetial corrected center locations (from the 2 roots)
                pf1 = p2 + alpha[0]*angle
                pf2 = p2 + alpha[1]*angle

                # select location closest to the prior step location
                if cdist(pf1, p1) > cdist(pf2, p1):
                    pf = pf2
                else:
                    pf = pf1

                # update the various arrays
                centers[i] = pf
                r_max[i] = cdist(pf, centers[[0, ]])
                r_current = np.amax(r_max)
                break

    # save text file of centers for external plotting or for other scripts
    if output_text:
        np.savetxt('centers.txt', centers, delimiter=' ')

    # plot of the aggregate
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        patches = []

        # draw circle patch at each particle center, adjust radius if you want overlap
        for i in range(n):
            circle = Circle((centers[i]), radius=r)
            patches.append(circle)

        p = PatchCollection(patches, color='black')
        ax.add_collection(p)
        plt.axis('off')
        plt.axis('equal')
        plt.show()


# --------------------------------------------------------------------------------------
def bounded_aggregate(path='centers.txt', r=1.0):
    """
    Plots bounding circle around previously generated 2D aggregate.

    bounded_aggregate(path='/centers.txt', r=1.0)

    Parameters
    ----------
    path : location of text file with particle centers (generated from aggregate2D)
    r : particle radius
        need to keep the same as the radius used to generate the text file

    Example
    -------
    bounded_aggregate(path='centers.txt')
    """
    # load the text file containing the centers
    cntrs = np.loadtxt(path)
    n = int(np.size(cntrs)/2)  # number of particles

    # determine the convex hull for the data set, use to find bounding points
    hull = ConvexHull(cntrs)

    hull_idx = np.unique(np.ravel(hull.simplices[:]))
    hull_pts = cntrs[hull_idx]

    # find distances from hull points to center
    r_hull = cdist(hull_pts, np.zeros((1, 2)))

    # set first bounding circle point to max hull radius
    p1 = hull_pts[[np.argmax(r_hull)]]

    # set second bounding circle point to max dist between 1st and other hull pts
    dist_p1 = cdist(hull_pts, p1)
    p2 = hull_pts[[np.argmax(dist_p1)]]

    # set third bounding circle point to max dist between midpoint and other hull pts
    dist_p_avg = cdist(hull_pts, (p1 + p2)/2)
    p3 = hull_pts[[np.argmax(dist_p_avg)]]

    # reorient point arrays for convenience
    p1, p2, p3 = p1.ravel(), p2.ravel(), p3.ravel()

    # find center of bounding circle
    A = np.array([[2*(p1[0] - p2[0]), 2*(p1[1] - p2[1])],
                  [2*(p1[0] - p3[0]), 2*(p1[1] - p3[1])]])
    b = np.array([[(p1[0]**2 - p2[0]**2) + (p1[1]**2 - p2[1]**2)],
                  [(p1[0]**2 - p3[0]**2) + (p1[1]**2 - p3[1]**2)]])

    # center of bounding circle
    x = np.linalg.solve(A, b)

    # bounding radius
    r_bounding = np.sqrt((p1[0] - x[0])**2 + (p1[1] - x[1])**2)

    # plot the aggregate and the bounding circle
    fig = plt.figure()
    ax = fig.add_subplot(111)
    patches = []

    # plot boundary circle patch
    bound = Circle((x[0], x[1]), radius=r_bounding + r, fill=False, ec='red')
    ax.add_patch(bound)

    # plot the particles
    for i in range(n):
        circle = Circle((cntrs[i]), radius=r)
        patches.append(circle)

    p = PatchCollection(patches, color='black')
    ax.add_collection(p)
    plt.axis('off')
    plt.axis('equal')
    plt.show()


# --------------------------------------------------------------------------------------
def aggregate_conn(path='centers.txt', r=1.0):
    """
    Plots connectivity diagram for previously generated 2D aggregate.

    aggregate_conn(path='centers.txt', r=1.0)

    Parameters
    ----------
    path : location of text file with particle centers (generated from aggregate2D)
    r : particle radius
        need to keep the same as the radius used to generate the text file

    Example
    -------
    aggregate_conn(path='/centers.txt')
    """
    # load the text file containing the centers
    cntrs = np.loadtxt(path)

    # calculate the distance between all particle centers
    dist_matrix = np.round(cdist(cntrs, cntrs), 1)

    # determine neighbor indices
    ind = np.argwhere(dist_matrix == 2)

    # plot the connectivity diagram for aggregate
    mpl.rcParams['lines.solid_capstyle'] = 'round'

    fig = plt.figure()
    ax = fig.add_subplot(111)
    lines = []

    for i in range(int(np.size(ind)/2)):
        seg = np.array([cntrs[ind[i, 0]],
                        cntrs[ind[i, 1]]])
        lines.append(seg)

    p = LineCollection(lines, linestyle='solid',
                       linewidth=0.2, color='black')
    ax.add_collection(p)
    plt.axis('off')
    plt.axis('equal')
    plt.show()
