# Arcs

*Starts a command block for creation of arcs to be used to compound geometric entities (e.g.: extruded or revolved surfaces).*

### Syntax:

Arcs N

Arc ID InitialPoint XIP YIP FinalPoint XFP YFP CenterPoint XCP YCP

-   N: number of arcs

-   ID: current arc identification number

-   XIP: X coordinate of the initial point of the arc

-   YIP: Y coordinate of the initial point of the arc

-   XFP: X coordinate of the final point of the arc

-   YFP: Y coordinate of the final point of the arc

-   XCP: X coordinate of the center point of the arc

-   YCP: Y coordinate of the center point of the arc

### Example:

Arcs 2

Arc 1 InitialPoint 0.0 -1.0 FinalPoint 0.0 1.0 CenterPoint -1.0e5 0.0

Arc 2 InitialPoint 0.024 -0.132 FinalPoint 0.024 0.132 CenterPoint -0.075 0.0

### Additional information:

Each arc is defined by the keyword Arc followed by the arc identification number (must be an ascending sequence starting from number one), its initial point, final point and center point.

The arc is supposed to lie on a local XY plane. The coordinate parameters indicating its initial point, end point and center point are understood on a local coordinate system (to be defined by the user). Figure 23 illustrates the arc:

<img src="../../assets/manual/image28.png" style="width:1.91743in;height:1.81328in" alt="" />

Figure 23 – Arc definition on a local coordinate system (figure from \[8\])

The arc parameterization **a**(*θ*) is given by:

|                                                   |        |
|---------------------------------------------------|--------|
| $\mathbf{a}(\theta)\mathbf{=}\begin{bmatrix} r\cos\theta + \text{XCP} \\ r\sin\theta + \text{YCP} \\ 0 \\ \end{bmatrix}$                                   | ( 43 ) |

where *r* is the arc radius (evaluated by Giraffe automatically). The arc center is given by local coordinates **c** = (*X**C**P*,*Y**C**P*). In Figure 23 one may observe two particular evaluations of parameterizations, for the initial point and for the end point of the arc, given by **i** and **f**. The location of these points in the local coordinate system compose the input data, such that **i****=**(XIP**,**YIP) and **f****=**(XFP**,**YFP).
