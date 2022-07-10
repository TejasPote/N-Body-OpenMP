# N-Body-OpenMP
The given C program simulates the dynamics trajectory of N bodies in a box moving under mutual attraction. We consider 1000 spherical bodies each of mass unity enclosed in a cuboid shaped container of dimensions 100X200X400. We are provided with the initial 3D coordinates of the bodies. The force acting on a pair of bodies is calculated as the product of their masses divided by the square of the Euclidean distance between them. We then update the positions of the particles by calculating the half-step and full step velocities. We carry out this procedure for a total of 720000 iteration steps. This serial implementation is then parallelized using OpenMP, thereby obtaining a significant speedup in execution time. 
