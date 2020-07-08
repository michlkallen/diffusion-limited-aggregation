# Diffusion-Limited Aggregation

This set of files generates off-lattice diffusion-limited aggregates in 2D and 3D. The files are written in Python and use `numpy`, `scipy`, and `matplotlib` to generate the centers and images.

## Why Off-Lattice?
Gridded simulations of diffusion-limited aggregation create artifacts in the image (since the "particle" can only move orthogonally). Off-lattice simulations avoid this by allowing the "particles" to move anywhere in the region.

The scripts use patches instead of pixels, so the images can be saved in a vector format if desired for better scalability. The algorithm itself borrows heavily from the work by Braga and Ribeiro to correct the particle location after contact.

## `dla.py`
Generates a 2D diffusion-limited aggregate. Typical use would be:
```python
>>> from dla import *
>>> aggregate2D(1000)  # generate 1000 particle aggregate
>>> bounded_aggregate()  # plot aggregate + bounding circle for newly created aggregate
>>> aggregate_conn()  # plot connectivity diagram for newly created aggregate
```

### Functions in `dla.py`
`aggregate2D(n, r=1.0, len=1.0, output_text=True, plot=True)`: base function - use to generate a 2D aggregate with options for plotting and exporting the particle centers as a text file (used in other functions)

![500 particle DLA](https://github.com/michlkallen/diffusion_limited_aggregation/blob/master/images/dla_500_particles.png)

`bounded_aggregate(path='centers.txt', r=1.0)`: plots the aggregate and its bounding circle (uses text file generated from `aggregate2D`)

![Bounded 500 particle DLA](https://github.com/michlkallen/diffusion_limited_aggregation/blob/master/images/dla_500_bounded.png)

`aggregate_conn(path='centers.txt', r=1.0)`: plots connectivity diagram for a given aggregate (uses text file generated from `aggregate2D`)

![500 particle connectivity](https://github.com/michlkallen/diffusion_limited_aggregation/blob/master/images/dla_500_connect.png)

## `dla3d.py`
Generates a 3D diffusion-limited aggregate. Similar procedure to above.

### Functions in `dla3d.py`
`aggregate3D(n, r=1.0, len=1.0)`: outputs a text file of particle centers. Could view by generating OpenSCAD script, for example.
