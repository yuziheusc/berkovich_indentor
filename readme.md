this is a lammps fix which does the nano indentation using berkovich indentor. 
This code is based on the original lammps fix and the work of B Wang.
the berkovich indentor consists of three rigid planes. 
The force is calculated using the distance from the plane. 
The functional form of the force is F(r) = - K (r - R)^2. 

** Format

fix id group myindent K bkvch v_x v_y v_z sharp dummy units box

dummy is not used for now. 
v_x,v_y,v_z is the ueser defined computations which provide the trajectory of berkovich tip.
sharp is the tangent of the angle between indentor edge and the z-axis. Large sharp means the indentor is sharp.

** Example

fix 2 all myindent 10.0 bkvch 57 50 v_z 0.4 0.0 units box


Y. He
