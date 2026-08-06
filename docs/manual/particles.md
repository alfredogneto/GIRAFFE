# Particles

*Starts a command block for creation of particles.*

### Syntax:

Particles N

Name ID data

-   N: number of particles

-   Name: current particle name

-   ID: current particle identification number

-   data: current particle data (depends on particle resources and requirements)

### Example:

Particles 2

Sphere 1 Mat 1 CS 1 Radius 2.5 Node 1

Sphere 2 Mat 1 CS 1 Radius 1.5 Node 1

### Additional information:

Each particle is defined by a specific keyword followed by the particle identification number (must be an ascending sequence starting from number one) and additional data. Each particle available and its input data is explained next.

## Sphere

*Creates a spherical particle.*

### Syntax:

Sphere PID Mat MID CS CSID Radius RV Node NID

-   PID: current particle identification number

-   MID: material identification number

-   CSID: coordinate system identification number

-   RV: sphere radius value

-   NID: identification number of the node defining the particle

### Example:

Sphere 1 Mat 1 CS 1 Radius 2.5 Node 1

### Additional information:

This is a spherical particle centered at a given nodal position. Three displacement and three rotation DOFs are defined for the node. In order to provide Giraffe the necessary data to evaluate the sphere mass and its moment of inertia, one has to choose a material identification number and a radius value for the particle. Furthermore, a coordinate system has to be chosen, as the reference orientation of the sphere. This CS is used for rendering plot of the sphere. The particle follows a rigid body kinematics behavior.

## Polyhedron

*Creates a particle with polyhedral shape.*

### Syntax:

Polyhedron PID Mat MID CS CSID CADData CID Node NID

-   PID: current particle identification number

-   MID: material identification number

-   CSID: coordinate system identification number

-   CID: CAD data identification number

-   NID: identification number of the node defining the particle

### Example:

Polyhedron 1 Mat 1 CS 1 CADData 1 Node 1

### Additional information:

This is a rigid polyhedral particle, with the geometrical shape defined by the CADData input and the material (homogeneous) defined by the Mat input. The origin in the CADData supplied will correspond to the particle nodal position. A convenient choice for the CAD origin is the center of mass of the particle. The orientation of the CAD coordinate system is interpreted as aligned with the directions defined by the coordinate system CS. Therefore, the CS definition permits to change the particle alignment.

Three displacement and three rotation DOFs are defined for the node. The particle follows a rigid body kinematics behavior.

##  VEMPolyhedron

*Creates a flexible particle with polyhedral shape. Flexibility is solved employing the Virtual Element Method.*

### Syntax:

VEMPolyhedron PID Mat MID CS CSID CADData CID SuperNode SNID StabFlag SFID StabBetaStiff SBS StabBetaMassLoad SBM LumpedMass LMF RayleighDamping Alpha RAV Beta RBV

-   PID: current particle identification number

-   MID: material identification number

-   CSID: coordinate system identification number

-   CID: CAD data identification number

-   SNID: identification number of the super node defining the particle

-   SFID: Virtual Element stabilization flag. Use 1 for a Finite Element – based stabilization

-   SBS: Stiffness stabilization constant. Use a number between 0 and 1

-   SBM: Mass stabilization constant. Use a number between 0 and 1

-   LMF: Flag to activate/inactivate the lumped mass matrix (1 or 0, respectively)

-   RAV: Rayleigh damping stiffness multiplier constant

-   RBV: Rayleigh damping mass multiplied constant

### Example:

VEMPolyhedron 1 Mat 1 CS 1 CADData 1 SuperNode 1 StabFlag 1 StabBetaStiff 0.25 StabBetaMassLoad 0.00 LumpedMass 0 RayleighDamping Alpha 0 Beta 0

### Additional information:

This is a flexible polyhedral particle, with the reference (stress-free) geometrical shape defined by the CADData input and the material (homogeneous) defined by the Mat input. The origin in the CADData supplied will correspond to the particle super node position. The orientation of the CAD coordinate system is interpreted as aligned with the directions defined by the coordinate system CS. Therefore, the CS definition permits to change the particle alignment.

Each vertex of the particle has three DOFs (displacements). A local virtual element mesh with a single element per particle is generated to handle the particle flexibility. When using a stabilization based on the finite element method, one needs to supply together with each CAD file a mesh file containing the connectivity for a internal tetrahedron mesh.
