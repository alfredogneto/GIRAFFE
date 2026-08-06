# Displacements

*Starts a command block for creation of displacements (to be prescribed).*

### Syntax:

Displacements N

Name ID data

-   N: number of displacements

-   Name: current displacement name

-   ID: current displacement identification number

-   data: current displacement data (depends on displacement resources and requirements)

### Example:

Displacements 1

NodalDisplacement 1 NodeSet 1 CS 1 NTimes 3

//Time UX UY UZ ROTX ROTY ROTZ

0 0 0 0 0 0 0

2 0 1.25 0 0 0 0

3 0 0 0 0 0 6.28

### Additional information:

Each displacement is defined by a specific keyword followed by the displacement identification number (must be an ascending sequence starting from number one) and additional data. Each displacement available and its input data is explained next.

**Note:** displacement prescription in Giraffe has to be carefully carried out by the user. In case of simultaneous prescription of displacements/rotations to a given node, the last-defined value only will be considered, following the sequence of displacements input data.

## NodalDisplacement

*Creates a nodal displacement (to be prescribed).*

### Syntax:

NodalDisplacement ID NodeSet NSID CS CSID

NTimes N

//Time UX UY UZ ROTX ROTY ROTZ

table data

NodalDisplacement ID NodeSet NSID CS CSID AngularParameters AP BoolTable BTC

NTimes N

//Time UX UY UZ ROTX ROTY ROTZ

table data

or

NodalDisplacement ID NodeSet NSID CS CSID MathCode

//MathCodeData

mathcode data

NodalDisplacement ID NodeSet NSID CS CSID AngularParameters AP BoolTable BTC

MathCode

//MathCodeData

mathcode data

-   ID: current displacement identification number

-   NSID: node set identification number

-   CSID: coordinate system identification number

-   AP: keyword to specify the input for rotation parameters. Currently options are “EulerVector” or “EulerAngles”. The keyword AngularParameters is optional, and its default specification is “EulerAngles”

-   BTC: bool table data for the NodalDisplacement (see BoolTable)

-   N: number of lines to be input in table data

-   table data: table of nodal displacement data following the rule:

column 1: time

columns 2-4: components of displacement vector

columns 5-7: components of the rotation vector (Euler rotation vector or Euler angles, according to the specification of AP)

-   mathcode data: pieces of code to evaluate a time-variable expression (6 pieces of code are to be entered as the time-series of:

displacement components (UX, UY, UZ)

rotation components (ROTX, ROTY, ROTZ) (Euler rotation vector or Euler angles, according to the specification of AP)

### Examples:

NodalDisplacement 1 NodeSet 1 CS 1 AngularParameters EulerVector BoolTable 1 0

NTimes 3

//Time UX UY UZ ROTX ROTY ROTZ

0 0 0 0 0 0 0

2 0 1.25 0 0 0 0

3 0 0 0 0 0 6.28

NodalDisplacement 2 NodeSet 1 CS 1

MathCode

//MathCode example - provides expressions to evaluate time-series

//The variable "t" has to be used to represent the time quantity

//Each code piece has to provide to Giraffe a mathematical expression to evaluate a given time series

//UX

Begin sin(2\*pi\*t) End

//UY

Begin exp(t^2) End

//UZ

Begin 0 End

//ROTX

Begin 0 End

//ROTY

Begin 0 End

//ROTZ

Begin 0 End

### Additional information:

A nodal displacement employs a table or a code expression to define a time-varying displacement and a time-varying rotation vector. Both are applied (prescribed) to each node of the node set. Displacement and rotation vector components must be defined using the desirable coordinate system CSID. The rotation input is made in **radians** by a Euler rotation vector if AngularParameters is specified as “EulerVector” or by a set of three Euler angles, if AngularParameters is specified as “EulerAngles”.

Displacement/rotation imposition is done in an incremental way. Thus, each time-step of solution will prescribe an increment of displacement/rotation, following the table of Nodal displacement entry or evaluating the piece of code provided in MathCode. Giraffe evaluates, for each evaluated time-step, the desired displacement/rotation values.

For the case of displacements, the difference between the desired value at the end/beginning of the time-step is prescribed, following the table or the code provided in MathCode. Thus, data is interpreted as incremental. For the case of rotations, Giraffe specifies the proper rotation increment such that the rotation tensor evaluated at the end of each time-step matches the calculated using the specified time-series data. This procedure is done employing the rotation multiplicative decomposition. The values for rotations (EulerVetor or EulerAngles) is interpreted in an absolute way.

The BoolTable (optional) keyword can be used to enable/disable the NodalDisplacement, according to a sequence of solution steps. This can be used to avoid undesirable displacement impositions, particularly useful when rotations take place.

The increment of displacement/rotation is prescribed within the time-step. The procedure continues until the end of the simulation. Note that the value of prescribed displacement/rotation in the table or evaluated by MathCode is not necessarily imposed during the simulation. It depends on how the user sets constraints on nodes. To prescribe displacements/rotations the user has to define both displacements and constraints inputs, in order to choose if each DOF is free or fixed (which may vary along solution sequence). If the DOF is free (default), the corresponding nodal displacement is ignored.

**Note 1**: When using a table time-series, Giraffe performs a linear interpolation between data provided to evaluate desired time-values for the series. In case of a need of data outside the defined time-range in a table, Giraffe considers: (i) the initially defined value as constant for time-values prior to it and (ii) the last defined value as constant for time-values after it.

**Note 2**: More details of MathCode may be seen in section: “MathCode”.

**Note 3**: The optional keywords “AngularParameters”, followed by “AP” and “BoolTable”, followed by “BTC” must appear in this sequence, and prior to table data (such as designation of “NTimes” or “MathCode” data.

## DisplacementField

*Creates a linear time-varying displacement field at nodes during a given solution step.*

### Syntax:

DisplacementField ID NNodes NN CS CSID SolutionStep SS

//Node UX UY UZ ROTX ROTY ROTZ

table data

-   ID: current displacement identification number

-   NN: number of nodes that will have an assigned displacement data

-   CSID: coordinate system identification number

-   SS: solution step number associated with the displacement field prescription

-   table data: table of nodal displacement data following the rule:

column 1: node

columns 2-4: components of displacement vector

columns 5-7: components of the rotation vector (see explanation below)

### Example:

DisplacementField 1 NNodes 5 CS 1 SolutionStep 2

//Node UX UY UZ ROTX ROTY ROTZ

1 0 1.0 0 0 0 0

2 0 0.5 0 0 0 0

3 0 0.0 0 0 0 0

6 0 -0.5 0 0 0 0

7 0 -1.0 0 0 0 0

### Additional information:

A nodal displacement field employs a table to define a field of generally distinct linear time-varying displacement/rotation for a set of nodes. All nodes assigned have independent displacement/rotation. Displacement and rotation vector components have to be defined using the desirable coordinate system CSID.

The rotation input is made in **radians** by a Euler rotation vector if AngularParameters is specified as “EulerVector” or by a set of three Euler angles, if AngularParameters is specified as “EulerAngles”. The default option is “EulerAngles”.

Displacement/rotation prescription is done in an incremental way. Thus, at each time-step of solution an increment of displacement/rotation is prescribed following a linear increment along time, such that along the solution step chosen the whole amount of displacement/rotation is prescribed.

Note that the value of prescribed displacement/rotation in the table is not necessarily imposed during the simulation. It depends on how the user sets constraints on nodes. To prescribe displacements/rotations the user has to define both displacements and constraints inputs, in order to choose if each DOF is free or fixed (which may vary along solution sequence). If the DOF is free (default), the corresponding nodal displacement is ignored.
