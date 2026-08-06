# Nodes

*Starts a command block for creation of nodes to be used to compound a finite element mesh or particles positions.*

### Syntax:

Nodes N

Node ID X Y Z

-   N: number of nodes

-   ID: current node identification number

-   X: current node X coordinate (on a global coordinate system)

-   Y: current node Y coordinate (on a global coordinate system)

-   Z: current node Z coordinate (on a global coordinate system)

### Example:

Nodes 3

Node 1 1.0 0.0 3.0

Node 2 0.0 2.5 -5.1

Node 3 0.0 0.0 -10.0

### Additional information:

Each node is defined by the keyword Node followed by the node identification number (must be an ascending sequence starting from number one), coordinates X, Y and Z.

After reading the input file, Giraffe checks all nodes, elements and particles connectivity. Based on this check, it evaluates how many degrees of freedom (DOFs) and which nature of DOFs have to be assigned for each node.
