<Please submit this file with your solution.>

CSCI 520, Assignment 1

Christin Carter
================

<Description of what you have accomplished>
1. Shear Forces
Each mass point is connected 20 diagonal neighbors. Hook's law is calculated for the force that each mass point
has on each of those neighbors
   
2. Bend Forces
Each mass is connected its 6 neighbors that are two mass points away. Hook's law is also calculated for the force
that the mass point applies to these neighbors

3. Structural Forces
Each mass is connected its 6 direct neighbors in the x, y, and z directions. Hook's law is also calculated for the force
that the mass point applies to these neighbors

4. Collision Forces
Each mass point's position is examined to see if they go outside the boundaries of the box and applies a collision
force if so.

5. Force Field Forces
Each mass point's position is examined to see if it's inside a force field voxel, and applies a trilinear interpolation
of the 8 force field vectors making up the corners of the voxel.

6. Damping Forces
For every spring force (shear, bend, structural) and collision force, a damping force has been applied.

<Modified files>
physics.cpp
jello.cpp (For enabling screenshots)
createWorld.cpp (For testing purposes)