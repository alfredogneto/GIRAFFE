# Initial Conditions

*Starts a command block for creation of initial conditions to be used in a transient dynamic analysis.*

### Syntax:

InitialConditions N

InitialCondition ICID Node NID DU DUX DUY DUZ OMEGA OX OY OZ SolutionStep SSID

-   N: number of initial conditions

-   ICID: current initial condition identification number

-   NID: node identification number

-   DUX, DUY and DUZ: components of velocity vector (on a global coordinate system)

-   OX, OY and OZ: components of angular velocity vector (on a global coordinate system)

-   SSID: solution step identification number

### Example:

InitialConditions 1

InitialCondition 1 Node 1 DU 0.0 1.0 0.0 OMEGA 0.0 0.0 0.0 SolutionStep 1

### Additional information:

Each initial condition is defined followed by the initial condition identification number (must be an ascending sequence starting from number one) and additional data. All data is interpreted on global coordinate system. The solution step input is necessary to associate the initial condition to a given solution step (dynamic).

**Remark**: when inserting initial conditions to nodes involved in a special constraint, please note that some of them may be ignored, according to the kind of special constraint.
