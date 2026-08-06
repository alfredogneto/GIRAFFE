# Coordinate Systems

*Starts a command block for creation of coordinate systems.*

### Syntax:

CoordinateSystems N

CS ID E1 E1XV E1YV E1ZV E3 E3XV E3YV E3ZV

-   N: number of coordinate systems

-   ID: current coordinate system identification number

-   E1XV, E1YV and E1ZV: components of direction E1

-   E3XV, E3YV and E3ZV: components of direction E3

### Example:

CoordinateSystems 2

CS 1 E1 1 0 0 E3 0 1 0

CS 2 E1 0 1 0 E3 0 0 1

### Additional information:

Each coordinate system is defined by the keyword CS followed by an identification number (must be an ascending sequence starting from number one). The coordinate systems are Cartesian and are defined by three unit-vectors, named E1, E2 and E3. Only E1 and E3 components have to be defined. The orientation E2 is internally calculated using the cross product.
