# diffusion-limited-aggregation
This set of files generates diffusion-limited aggregates in 2D and 3D. The files are written in Python and use numpy, scipy, and matplotlib to generate the centers and images.

The main aggregate scripts are for off-lattice DLA simulations to avoid the artifacts present in a gridded simulation. The scripts also use patches instead of pixels, so the images can be saved in a vector format if desired for better scalability. The algorithm itself borrows heavily from the work by Braga and Ribeiro to correct the particle location after contact.

DLA.py generates a 2D diffusion-limited aggregate.

DLA3d.py generates a 3D diffusion-limited aggregate (no plotting).

DLAboundary.py calculates the bounding circle for a given aggregate and plots both.

DLAConnect.py connects the particle centers with lines instead of plotting the particles themselves.
