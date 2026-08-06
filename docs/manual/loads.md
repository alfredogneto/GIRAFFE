# Loads

*Starts a command block for creation of loads.*

### Syntax:

Loads N

Name ID data

-   N: number of loads

-   Name: current load name

-   ID: current load identification number

-   data: current load data (depends on load resources and requirements)

### Example:

Loads 1

NodalLoad 1 NodeSet 2 CS 1 NTimes 2

//Time FX FY FZ MX MY MZ

0 0 0 0 0 0 0

1 1000 0 0 0 0 0

### Additional information:

Each load is defined by a specific keyword followed by the load identification number (must be an ascending sequence starting from number one) and additional data. Each load available and its input data is explained next.

## NodalLoad

*Creates a nodal load.*

### Syntax:

NodalLoad ID NodeSet NSID CS CSID NTimes N

//Time FX FY FZ MX MY MZ

table data

or

NodalLoad ID NodeSet NSID CS CSID MathCode

//MathCodeData

mathcode data

-   ID: current load identification number

-   NSID: node set identification number

-   CSID: coordinate system identification number

-   N: number of lines to be input in table data

-   table data: table of nodal loads data following the rule:

column 1: time

columns 2-4: components of force vector

columns 5-7: components of the moment vector

-   mathcode data: pieces of code to evaluate a time-variable expression (6 pieces of code are to be entered as the time-series of:

force components (FX, FY, FZ)

moment components (MX, MY, MZ)

### Example:

NodalLoad 1 NodeSet 1 CS 1 NTimes 2

//Time FX FY FZ MX MY MZ

0 0 0 0 0 0 0

1 1000 0 0 0 0 0

NodalLoad 2 NodeSet 2 CS 1 MathCode

//MathCode example - provides expressions to evaluate time-series

//The variable "t" has to be used to represent the time quantity

//Each code piece has to provide to Giraffe a mathematical expression to evaluate a given time series

//FX

Begin sin(2\*pi\*t) End

//FY

Begin exp(t^2) End

//FZ

Begin 0 End

//MX

Begin 0 End

//MY

Begin 0 End

//MZ

Begin 0 End

### Additional information:

A nodal load employs a table or a code expression to define a time-varying force and a time-varying moment. Both are applied to the nodes of the node set. Force and moment components are to be defined using any defined coordinate system CSID. The values of forces and moments components are divided between the nodes of the node set and each node receives the same amount of force/moment, such that the total magnitude of force and moment fulfills the user input.

**Note 1**: Division of forces/moments on nodes of the node set is done checking if the DOFs are active /non-active in each node. Thus, for Shell\_1 element, for example, only mid-side nodes receive moment loads. Corners have only translational DOFs active. This checking is automatically done by Giraffe.

**Note 2** : Giraffe performs a linear interpolation between data provided as a tabular time-series. In case of a need of data outside the defined time-range in a table, Giraffe considers: (i) the initially defined value as constant for time-values prior to it and (ii) the lastly defined value as constant for time-values after it.

**Note 3**: More details of MathCode may be seen in section: “MathCode”.

## NodalFollowerLoad

*Creates a nodal follower load.*

### Syntax:

NodalFollowerLoad ID NodeSet NSID CS CSID NTimes N

//Time FX FY FZ MX MY MZ

table data

or

NodalFollowerLoad ID NodeSet NSID CS CSID MathCode

//MathCodeData

mathcode data

-   ID: current load identification number

-   NSID: node set identification number

-   CSID: coordinate system identification number

-   N: number of lines to be input in table data

-   table data: table of nodal loads data following the rule:

column 1: time

columns 2-4: components of force vector

columns 5-7: components of the moment vector

-   mathcode data: pieces of code to evaluate a time-variable expression (6 pieces of code are to be entered as the time-series of:

force components (FX, FY, FZ)

moment components (MX, MY, MZ)

### Example:

NodalFollowerLoad 1 NodeSet 1 CS 1 NTimes 2

//Time FX FY FZ MX MY MZ

0 0 0 0 0 0 0

1 1000 0 0 0 0 0

NodalFollowerLoad 2 NodeSet 2 CS 1 MathCode

//MathCode example - provides expressions to evaluate time-series

//The variable "t" has to be used to represent the time quantity

//Each code piece has to provide to Giraffe a mathematical expression to evaluate a given time series

//FX

Begin sin(2\*pi\*t) End

//FY

Begin exp(t^2) End

//FZ

Begin 0 End

//MX

Begin 0 End

//MY

Begin 0 End

//MZ

Begin 0 End

### Additional information:

A nodal follower load employs a table or a code expression to define a time-varying force and a time-varying moment. Both are applied to the nodes of the node set. Force and moment components are to be defined using any defined coordinate system CSID. The values of forces and moments components are divided between the nodes of the node set and each node receives the same amount of force/moment, such that the total magnitude of force and moment fulfills the user input.

**Note 1**: Division of forces/moments on nodes of the node set is done checking if the DOFs are active /inactive in each node. Thus, for Shell\_1 element, for example, only mid-side nodes receive moment loads. Corners have only translational DOFs active. This checking is automatically done by Giraffe.

The components of force/moment are kept in a local coordinate system that follows the rotations of each node. Then, the follower load updates according to the movement experienced by the node.

**Note 2**: Giraffe performs a linear interpolation between data provided as a tabular time-series. In case of a need of data outside the defined time-range in a table, Giraffe considers: (i) the initially defined value as constant for time-values prior to it and (ii) the lastly defined value as constant for time-values after it.

**Note 3**: More details of MathCode may be seen in section: “MathCode”.

## PipeLoad

*Creates a pipe load (to be used together with Pipe\_1 elements).*

### Syntax:

PipeLoad ID ElementSet ESID NTimes N

//Time P0I P0E RhoI RhoE

table data

or

PipeLoad ID ElementSet ESID MathCode

//MathCodeData

mathcode data

-   ID: current load identification number

-   ESID: element set identification number

-   N: number of lines to be input in table data

-   table data: table of pipe loads data following the rule:

column 1: time

column 2: internal pressure of the pipe

column 3: external pressure of the pipe

column 4: internal fluid specific mass (currently not used by Giraffe)

column 5: external fluid specific mass (currently not used by Giraffe)

-   mathcode data: pieces of code to evaluate a time-variable expression (4 pieces of code are to be entered as the time-series of:

(1) internal pressure of the pipe;

(2) external pressure of the pipe;

(3) internal fluid specific mass (currently not used by Giraffe);

(4) external fluid specific mass (currently not used by Giraffe).

### Example:

PipeLoad 1 ElementSet 1 NTimes 2

//Time P0I P0E RhoI RhoE

2 0 0 0 0

3 1200000 0 0 0

PipeLoad 2 ElementSet 2 MathCode

//MathCode example - provides expressions to evaluate time-series

//The variable "t" has to be used to represent the time quantity

//Each code piece has to provide to Giraffe a mathematical expression to evaluate a given time series

// P0I

Begin sin(2\*pi\*t) End

// P0E

Begin exp(t^2) End

// RhoI

Begin 0 End

// RhoE

Begin 0 End

### Additional information:

A pipe load employs a table or a code expression to define time-varying internal/external pressures in a pipe. Both are applied to each element of the element set. All elements of the element set must be of the type Pipe\_1.

**Note 1**: Giraffe performs a linear interpolation between data provided as a tabular time-series. In case of a need of data outside the defined time-range in a table, Giraffe considers: (i) the initially defined value as constant for time-values prior to it and (ii) the lastly defined value as constant for time-values after it.

**Note 2**: More details of MathCode may be seen in section: “MathCode”.

## ShellLoad

*Creates a shell load (to be used together with Shell\_1 elements).*

### Syntax:

ShellLoad ID ElementSet ESID AreaUpdate AUB NTimes N

//Time P

table data

or

ShellLoad ID ElementSet ESID AreaUpdate AUB MathCode

//MathCodeData

mathcode data

-   ID: current load identification number

-   ESID: element set identification number

-   AUB: Boolean variable to define if the area is to be updated for recalculation of the resultant due to pressure integration (1) or not (0)

-   N: number of lines to be input in table data

-   table data: table of shell loads data following the rule:

column 1: time

column 2: pressure applied on shell element surface

-   mathcode data: pieces of code to evaluate a time-variable expression (1 piece of code is to be entered as the time-series of the pressure applied on shell element surface.

### Example:

ShellLoad 1 ElementSet 1 AreaUpdate 1 NTimes 2

//Time P

0 0.0

2 +10.0

ShellLoad 2 ElementSet 2 AreaUpdate 1 MathCode

//MathCode example - provides expressions to evaluate time-series

//The variable "t" has to be used to represent the time quantity

//Each code piece has to provide to Giraffe a mathematical expression to evaluate a given time series

// Pressure

Begin t\*100\*sin(t) End

### Additional information:

A shell load employs a table or a code expression to define a time-varying pressure in a shell element. It is applied to each element of the element set. All elements of the element set must be of the type Shell\_1. A positive value of pressure acts opposite to the external normal direction of the element surface, which is associated with the sequence of nodes numbering in the element. The AreaUpdate keyword establishes if the resultant of pressure integrated along the element area should consider or not the area update, according to element deformations.

**Note 1**: Giraffe performs a linear interpolation between data provided as a tabular time-series. In case of a need of data outside the defined time-range in a table, Giraffe considers: (i) the initially defined value as constant for time-values prior to it and (ii) the lastly defined value as constant for time-values after it.

**Note 2**: More details of MathCode may be seen in section: “MathCode”.
