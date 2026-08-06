# Pipe Sections

*Starts a command block for creation of pipe sections.*

### Syntax:

PipeSections N

PS ID EA EAV EI EIV GJ GJV GA GAV Rho RV CDt CDTV CDn CDNV CAt CATV CAn CANV De DEV Di DIV

-   N: number of pipe sections

-   ID: current pipe section identification number

-   EAV : equivalent axial stiffness product value

-   EIV: equivalent bending stiffness product value

-   GJV: equivalent torsional stiffness product value

-   GAV: equivalent shearing stiffness product value

-   RV: equivalent mass per unit reference length value

-   CDTV: drag coefficient in tangential direction

-   CDNV: drag coefficient in normal direction

-   CATV: added mass coefficient in tangential direction

-   CANV: added mass coefficient in normal direction

-   DEV: external diameter value

-   DIV: internal diameter value

### Example:

PipeSections 2

PS 1 EA 400000000 EI 15000 GJ 2250000 GA 200000000 Rho 120 CDt 0.0 CDn 1.0 CAt 0 CAn 1.0 De 0.25 Di 0.2

PS 2 EA 6080489749 EI 110821607.7 GJ 85247390.5 GA 2338649904 Rho 237 CDt 0.0 CDn 1.0 CAt 0 CAn 1.0 De 0.4 Di 0.2

### Additional information:

Each pipe cross section centroid must lie at the pipe axis, defined by the nodes of the mesh. More details about theoretical details of the pipe section attributes can be found in \[2\], together with applications for offshore risers simulations.

This cross-section is to be used with Pipe\_1 element type. Alternatively, the user may also use it with Truss\_1 element type. In this case, only axial stiffness and hydrodynamic coefficients are considered in the model.
