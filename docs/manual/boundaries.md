# Boundaries

*Starts a command block for creation of boundaries (within a particle simulation context).*

### Syntax:

Boundaries N

Name ID data

-   N: number of boundaries

-   Name: current boundary name

-   ID: current boundary identification number

-   data: current boundary data (depends on boundary resources and requirements)

### Example:

Boundaries 1

STLBoundary 1 Mat 1 CS 1 CADData 1 Node 1

### Additional information:

Each boundary is defined by a specific keyword followed by its identification number (must be an ascending sequence starting from number one) and additional data. Each boundary type available and its input data is explained next.

## STLBoundary

*Creates a boundary based on a STL file geometry.*

### Syntax:

STLBoundary BID Mat MID CS CSID CADData CID Node NID

-   BID: current boundary identification number

-   MID: material identification number

-   CSID: coordinate system identification number

-   CID: CAD data identification number

-   NID: identification number of the node defining the boundary

### Example:

STLBoundary 1 Mat 1 CS 1 CADData 1 Node 1

### Additional information:

This is a rigid boundary, with the geometrical shape defined by a STL CAD file, which identification number is provided. The origin in the CADData supplied will correspond to the boundary nodal position. The orientation of the CAD coordinate system is interpreted as aligned with the directions defined by the coordinate system CS. Therefore, the CS definition permits to change the boundary alignment. The material definition is employed only to establish a proper interface property between the boundary and the particles.

Three displacement and three rotation DOFs are defined for the node. The has no contributions to inertial forces. Rather, it provides only physical limitations to the motion of particles.
