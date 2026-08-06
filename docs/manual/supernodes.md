# SuperNodes

*Starts a command block for creation of super nodes to be used to create flexible particles reference positions.*

### Syntax:

SuperNodes N

SuperNode ID X Y Z

-   N: number of super nodes

-   ID: current super node identification number

-   X: current super node X coordinate (on a global coordinate system)

-   Y: current super node Y coordinate (on a global coordinate system)

-   Z: current super node Z coordinate (on a global coordinate system)

### Example:

SuperNodes 3

SuperNode 1 1.0 0.0 3.0

SuperNode 2 0.0 2.5 -5.1

SuperNode 3 0.0 0.0 -10.0

### Additional information:

Each super node is defined by the keyword SuperNode followed by the super node identification number (must be an ascending sequence starting from number one), coordinates X, Y and Z.

After reading the input file, Giraffe checks all super nodes and particles. Based on this check, it creates the necessary degrees of freedom (DOFs) for the analysis.
