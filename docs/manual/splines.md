# Splines

*Starts a command block for creation of spline curves based on* NodeSets*.*

### Syntax:

Splines N

Spline ID R RV NodeSet NSID

-   N: number of splines

-   ID: current spline identification number

-   RV: radius value to define the spline surface

-   NSID: NodeSet identification

### Example:

Splines 2

Spline 1 R 0.1 NodeSet 1

Spline 2 R 0.3 NodeSet 2

### Additional information:

The current implementation of spline curves adopts a quadratic spline description of the centerline (curve space). This description has *C*<sup>1</sup> continuity. The first and last nodes of the splines are interpolatory (coincident to the nodes). The knot vector is automatically defined with uniform distribution, except at the first and last nodes due to interpolation. The nodes of the NodeSet are taken as control points of the spline. Note that the nodes in the NodeSet should be defined sequentially as the nodes used to define a single structure composed of beam elements.
