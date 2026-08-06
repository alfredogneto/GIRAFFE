# Constraints

*Starts a command block for creation of constraints.*

### Syntax:

Constraints N

Name ID data

-   N: number of constraints

-   Name: current constraint name

-   ID: current constraint identification number

-   data: current constraint data (depends on constraint resources and requirements)

### Example:

Constraints 2

NodalConstraint 1 NodeSet 1

UX BoolTable 1

UY BoolTable 1

UZ BoolTable 1

ROTX BoolTable 1 1

ROTY BoolTable 1 1 0 1

ROTZ BoolTable 1

NodalConstraint 2 NodeSet 2

UX BoolTable 0

UY BoolTable 0 0

UZ BoolTable 0 0

ROTX BoolTable 0 0

ROTY BoolTable 0 0

ROTZ BoolTable 0 0

### Additional information:

Each constraint is defined by a specific keyword followed by the constraint identification number (must be an ascending sequence starting from number one) and additional data. Each constraint available and its input data is explained next.

**Note:** Constraints imposition follow the increasing sequence defined in the input file. In case the user establishes conflicting input data (such as defining node sets with common nodes and applying for each one a distinct boundary condition), the **fixed** condition for a DOF prevails.

## NodalConstraint

*Creates a nodal constraint.*

### Syntax:

NodalConstraint ID NodeSet NSID

constraint data

-   ID: current constraint identification number

-   NSID: node set identification number

-   constraint data: specific keywords UX, UY, UZ, ROTX, ROTY and ROTZ, each followed by a BoolTable keyword to define nodal constraint solution sequence data

### Example:

Constraints 2

NodalConstraint 1 NodeSet 1

UX BoolTable 1

UY BoolTable 1

UZ BoolTable 1

ROTX BoolTable 1 1

ROTY BoolTable 1 1 0 1

ROTZ BoolTable 1

NodalConstraint 2 NodeSet 2

UX BoolTable 1

UY BoolTable 0 1

ROTZ BoolTable 1 0

### Additional information:

A nodal constraint is employed to define fixed DOFs in the model. When establishing a mesh, all DOFs are free. To establish Dirichlet boundary conditions in model regions, it is necessary to establish node sets and, on that regions, apply desirable constraints. For that, the user may separately operate on distinct DOFs of each node, within a chosen node set. These are:

-   UX: displacement in global X direction

-   UY: displacement in global Y direction

-   UZ: displacement in global Z direction

-   ROTX: rotation in global X direction

-   ROTY: rotation in global Y direction

-   ROTZ: rotation in global Z direction

The BoolTable resource is used to establish the behavior of each DOF along solution evolution (see BoolTable), that is, if the constraint is turned on/off. This permits to establish scenarios of alternating constraints along solution evolution.

DOFs with constraints turned on are fixed. In case of definition of displacements associated with such nodes, these will be prescribed along solution evolution, following the established table of displacement along time, for each DOF (see Displacements). Otherwise, Giraffe will simply consider a zero-displacement value.
