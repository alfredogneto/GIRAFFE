# Surfaces

*Starts a command block for creation of surfaces.*

### Syntax:

Surfaces N

Name ID data optional

-   N: number of surfaces

-   Name: current surface name

-   ID: current surface identification number

-   data: current surface data (depends on surface resources and requirements)

-   optional: optional keywords for degenerating a surface

### Example:

Surfaces 2

RigidTriangularSurface\_1 1 Points 1 2 3 PilotNode 1

FlexibleSECylinder\_1 2 A 0.1 B 0.1 N 3.0 CS 1 NormalExterior Nodes 1 2

### Additional information:

Each surface is defined by a specific keyword followed by the surface identification number (must be an ascending sequence starting from number one) and additional data. Each surface available and its input data is explained next.

The user may choose a fixed value for a given convective coordinate or, alternatively, establish a number of divisions to be performed by Giraffe along the valid range of a given convective coordinate, thus establishing a set of fixed values for such coordinate. If this is the choice, Giraffe performs divisions uniformly along that coordinate.

When employing degeneration, the original surface turns into a curve (or a set of curves) or a point (or a set of points). Examples of degenerations are given next:

//Degeneration into a single point with Coord1 = C1V and Coord2 = C2V

Degeneration Coord1 C1V Coord2 C2V

//Degeneration into a set of points with Coord1 = C1V and ND2 divisions along Coord2 range

Degeneration Coord1 C1V Div2 ND2

//Degeneration into a set of points with ND1 divisions along Coord1 range and Coord2 = C2V

Degeneration Div1 ND1 Coord2 C2V

//Degeneration into a set of points with ND1 divisions along Coord1 range and ND2 divisions along Coord2 range

Degeneration Div1 ND1 Div2 ND2

//Degeneration into a single curve with Coord1 = C1V

Degeneration Coord1 C1V

//Degeneration into a single curve with Coord2 = C2V

Degeneration Coord2 C2V

//Degeneration into a set of curves with ND1 divisions along Coord1

Degeneration Div1 ND1

//Degeneration into a set of curves with ND2 divisions along Coord2

Degeneration Div2 ND2


## RigidTriangularSurface\_1

*Creates a rigid triangular surface.*

### Syntax:

RigidTriangularSurface\_1 SID Points ID1 ID2 ID3 PilotNode PNID

-   SID: current surface identification number

-   ID1, ID2 and ID3: points identification numbers

-   PNID: pilot node identification number

### Example:

RigidTriangularSurface\_1 1 Points 1 2 3 PilotNode 1

### Additional information:

The vertices of a rigid triangular region are defined by A, B and C, positioned at **x**<sub>*A*</sub>, **x**<sub>*B*</sub> and **x**<sub>*C*</sub>, respectively. The surface is parameterized by:

<table><tbody><tr class="odd"><td><blockquote><p><span class="math inline">$\Gamma\left( \zeta,\ \theta \right) = \begin{bmatrix}
N_{A}(\zeta,\ \theta) &amp; N_{B}(\zeta,\ \theta) &amp; N_{C}(\zeta,\ \theta) \\
\end{bmatrix}\begin{bmatrix}
\mathbf{x}_{A}^{} \\
\mathbf{x}_{B}^{} \\
\mathbf{x}_{C}^{} \\
\end{bmatrix}$</span>,</p></blockquote></td><td>( 44 )</td></tr></tbody></table>

with $N_{A}\left( \zeta,\ \theta \right) = - \frac{1}{2}(\zeta + \theta)$, $N_{B}\left( \zeta,\ \theta \right) = \frac{1}{2}(1 + \zeta)$ and $N_{C}\left( \zeta,\ \theta \right) = \frac{1}{2}(1 + \theta)$. The parameters *ζ* and *θ* map the points inside the triangular region. An example is given in Figure 24.

<img src="../../assets/manual/emf_image29.png" style="max-width:100%" alt="" />

Figure 24 – Rigid triangular surface example

This surface is rigidly connected to a pilot node, which rules its movement.

**Note**: this surface is currently available for using together with the contact NSSS.

## RigidOscillatorySurface\_1

*Creates a rigid oscillatory surface.*

### Syntax:

RigidOscillatorySurface\_1 SID A1 A1V A2 A2V A12 A12V Lambda1 L1V Lambda2 L2V Phi1 P1V Phi2 P2V Waves1 W1N Waves2 W2N CS CSID PilotNode PNID

-   SID: current surface identification number

-   A1V: amplitude to sin in direction *ζ* of surface parameterization

-   A2V: amplitude to sin in direction *θ* of surface parameterization

-   A12V: amplitude to the product of sines in directions ζ and θ of surface parameterization

-   L1V: wave length along direction *ζ* of surface parameterization

-   L2V: wave length along direction *θ* of surface parameterization

-   P1V: phase along direction *ζ* of surface parameterization

-   P2V: phase along direction *θ* of surface parameterization

-   Waves1: number of waves along direction *ζ* of surface parameterization

-   Waves2: number of waves along direction *θ* of surface parameterization

-   CSID: coordinate system identification number

-   PNID: pilot node identification number

### Example:

RigidOscillatorySurface\_1 1 A1 1.0 A2 1.0 A12 0.0 Lambda1 1.0 Lambda2 2.0 Phi1 0.0 Phi2 0.0 Waves1 2.5 Waves2 1.5 CS 1 PilotNode 1

### Additional information:

A rigid oscillatory surface may be used to define wave patterns on a surface, possibly in two directions. One may define it as aligned with an arbitrary direction, and rigidly attached to a pilot node, which will rule its movement along the model evolution. Let one define a function in a local coordinate system (P, ζ, θ), where P is the pilot node position – origin of the system:

|                                                                                                |        |
|------------------------------------------------------------------------------------------------|--------|
| $\Gamma\left( \zeta,\ \theta \right) = \mathbf{x}_{P}^{} + \mathbf{Q}\begin{bmatrix} \zeta \\ \theta \\ \Psi\left( \zeta,\theta \right) \\ \end{bmatrix}$,                                                                                | ( 45 ) |

where **x**<sub>*P*</sub> is the pilot node position, **Q** is a rotation matrix that rules the alignment of the surface in space, and depends on its initial orientation and on the pilot node rotation experienced during the model evolution. Finally, *Ψ*(*ζ*, *θ*) is a function used to describe the local geometry of the surface, given by:

<table><tbody><tr class="odd"><td><blockquote><p><span class="math inline">$\Psi\left( \zeta,\ \theta \right) = A_{1}\sin\left( \frac{2\pi\zeta}{\lambda_{1}} + \phi_{1} \right) + A_{2}\sin\left( \frac{2\pi\theta}{\lambda_{2}} + \phi_{2} \right) + A_{12}\sin{\left( \frac{2\pi\zeta}{\lambda_{1}} + \phi_{1} \right)\sin\left( \frac{2\pi\theta}{\lambda_{2}} + \phi_{2} \right)}$</span>,</p></blockquote></td><td>( 46 )</td></tr></tbody></table>

where *A*<sub>1</sub>, *A*<sub>2</sub> and *A*<sub>12</sub> are amplitudes, *λ*<sub>1</sub> and *λ*<sub>2</sub> are wave-lengths and *ϕ*<sub>1</sub> and *ϕ*<sub>2</sub> are phases, all used to define the desired surface.

<img src="../../assets/manual/emf_image30.png" style="max-width:100%" alt="" />

Figure 25 – Oscillatory surface example

**Note**: this surface is currently available for using together with the contact NSSS.

## FlexibleSECylinder\_1

*Creates a flexible super elliptical surface.*

### Syntax:

//Surface with normal pointing outwards super elliptical cylinder

FlexibleSECylinder\_1 SID A AV B BV N NV CS CSID NormalExterior Nodes ID1 ID2

//Surface with normal pointing inwards super elliptical cylinder

FlexibleSECylinder\_1 SID A AV B BV N NV CS CSID NormalInterior Nodes ID1 ID2

//Defining surface with two distinct coordinate system alignments

FlexibleSECylinder\_1 SID A AV B BV N NV CSA CSAID CSB CSBID NormalInterior Nodes ID1 ID2

-   SID: current surface identification number

-   AV: semi-axis value lying in direction E1

-   BV: semi-axis value lying in direction E2

-   NV: super ellipse exponent value

-   CSID: coordinate system identification number (single alignment option)

-   CSAID and CSBID: coordinate systems identification numbers (two distinct alignments option)

-   ID1 and ID2: nodes identification numbers

NormalExterior or NormalInterior keywords are used to assign that the normal direction of the surface points outwards/inwards of super elliptical cylinder.

### Example:

//Surface with normal pointing outwards super elliptical cylinder

FlexibleSECylinder\_1 1 A 0.1 B 0.1 N 3.0 CS 1 NormalExterior Nodes 1 2

//Surface with normal pointing inwards super elliptical cylinder

FlexibleSECylinder\_1 2 A 0.1 B 0.1 N 3.0 CS 1 NormalInterior Nodes 1 2

//Defining surface with two distinct coordinate system alignments

FlexibleSECylinder\_1 3 A 0.1 B 0.1 N 3.0 CSA 1 CSB 2 NormalInterior Nodes 1 2

### Additional information:

A flexible surface with super-elliptical cross section can be defined with this command. A schematic visualization of a deformed surface can be seen in Figure 26. The cross section of the surface is given by the expression:

|                                                                                  |        |
|----------------------------------------------------------------------------------|--------|
| $\left\vert  \frac{x}{a} \right\vert ^{n} + \left\vert  \frac{y}{b} \right\vert ^{n} = 1$, | ( 47 ) |

where “a” and “b” are super ellipse semi-axis. Parameter “n” is the curve exponent. The surface extreme positions are guided by the movement of nodes, located at each extreme cross section centroid. This includes translation and rotation. In Figure 26 one may see an example where one extreme cross section is twisted with respect to the other extreme cross section. This can be used to represent the external surface of a beam element in an approximated way.

<img src="../../assets/manual/image31.png" style="width:2.41096in;height:1.89067in" alt="" />

Figure 26 – Cylinder with super-elliptical cross section

The alignment of the surface extreme cross-sections at reference configuration is defined by two coordinate systems located at both nodes that define the surface limits. In case of a single direction is defined for both nodes, use a single coordinate system input using CS keyword. In case of distinct alignments, use two coordinate system inputs by keywords CSA and CSB.

**Remark**: the definition of two distinct coordinate systems (e.g.: CSA 1 CSB 2) is particularly useful when the reference configuration of two connecting surfaces FlexibleSECylinder\_1 are not aligned. For example, see the extract of input code to Giraffe defined next, where two pairs of connected surfaces are defined. The pair on the left size uses two distinct alignments for FlexibleSECylinder\_1, such that at the connection of surfaces the same coordinate system identification number is assigned for both. The pair on the right size assumes a single coordinate system for each surface, creating a distinct connection with discontinuities between surfaces.

In cases of creating a surface set to represent a single contact patch composed by many surfaces, the first option is desirable, since no holes are created between surfaces, which may lead to contact lose or bad convergence/divergence in models.

<table><tbody><tr class="odd"><td><p>...</p><p>Node 12 -1 -1 0.7</p><p>Node 13 -1 +1 0.7</p><p>Node 14 +1 +1 0.7</p><p>Node 15 +2 -1 0.7</p><p>Node 16 +2 +1 0.7</p><p>Node 17 +4 +1 0.7</p><p>...</p><p>CS 4</p><p>CSYS 1 E1 1 0 0 E3 0 0 1</p><p>CSYS 2 E1 0 0 1 E3 0 1 0</p><p>CSYS 3 E1 0 0 1 E3 1 1 0</p><p>CSYS 4 E1 0 0 1 E3 1 0 0</p><p>...</p><p>FlexibleSECylinder_1 1 A 0.3 B 0.3 N 2.4 CSA 2 CSB 3 NormalExterior Nodes 12 13</p><p>FlexibleSECylinder_1 2 A 0.3 B 0.3 N 2.4 CSA 3 CSB 4 NormalExterior Nodes 13 14</p><p>FlexibleSECylinder_1 3 A 0.3 B 0.3 N 2.4 CS 2 NormalExterior Nodes 15 16</p><p>FlexibleSECylinder_1 4 A 0.3 B 0.3 N 2.4 CS 4 NormalExterior Nodes 16 17</p></td></tr></tbody></table>

<img src="../../assets/manual/image32.png" style="width:5.90556in;height:3.19404in" alt="" />

Figure 27 – Example of distinct/single coordinate systems to align two FlexibleSECylinder\_1 surfaces


## FlexibleTriangularSurface\_2

*Creates a flexible triangular surface.*

### Syntax:

FlexibleTriangularSurface\_2 SID Nodes ID1 ID2 ID3 ID4 ID5 ID6

-   SID: current surface identification number

-   ID1, ID2, ID3, ID4, ID5 and ID6: nodes identification numbers

### Example:

FlexibleTriangularSurface\_2 1 Nodes 1 2 3 4 5 6

### Additional information:

The nodes A, B and C are the vertices of a triangular surface and D, E and F are mid-points of the edges of a reference triangle. They are located on **x**<sub>*K*</sub> , where K assumes any of such node’s indexes. A triangular surface *Γ* may be parameterized by:

<table><tbody><tr class="odd"><td><blockquote><p><span class="math inline">$\Gamma\left( \zeta,\ \theta \right) = \begin{bmatrix}
\begin{matrix}
\begin{matrix}
N_{A}^{2}(\zeta,\ \theta) &amp; N_{B}^{2}(\zeta,\ \theta) \\
\end{matrix} &amp; \begin{matrix}
N_{C}^{2}(\zeta,\ \theta) &amp; N_{D}^{2}(\zeta,\ \theta) \\
\end{matrix} \\
\end{matrix} &amp; \begin{matrix}
N_{E}^{2}(\zeta,\ \theta) &amp; N_{F}^{2}(\zeta,\ \theta) \\
\end{matrix} \\
\end{bmatrix}\begin{bmatrix}
\begin{matrix}
\mathbf{x}_{A}^{} \\
\mathbf{x}_{B}^{} \\
\mathbf{x}_{C}^{} \\
\end{matrix} \\
\begin{matrix}
\mathbf{x}_{D}^{} \\
\begin{matrix}
\mathbf{x}_{E}^{} \\
\mathbf{x}_{F}^{} \\
\end{matrix} \\
\end{matrix} \\
\end{bmatrix}$</span>,</p></blockquote></td><td>( 48 )</td></tr></tbody></table>

where *N*<sub>*K*</sub><sup>2</sup>(*ζ*, *θ*) are shape functions in plane ζθ. This surface experiences deformation following the interpolation of nodes, located at points A, B, C, D, E and F. Figure 28 shows an example of such parameterization, deformed.

<img src="../../assets/manual/image33.png" style="width:2.30908in;height:1.43056in" alt="" />

Figure 28 – Flexible triangular surface example

This surface can be used to establish a surface candidate to contact covering a shell element. The node sequence pattern follows exact the same rule of Shell\_1 element (see Shell\_1).

## FlexibleArcExtrusion\_1

*Creates a surface generated by the extrusion of an arc along the path defined by nodes.*

### Syntax:

//Surface with external normal pointing to the center of arc

FlexibleArcExtrusion\_1 SID Arc AID CS CSID Nodes ID1 ID2 Concave

//Surface with external normal pointing opposite to the center of arc

FlexibleArcExtrusion\_1 SID Arc AID CS CSID Nodes ID1 ID2 Convex

-   SID: current surface identification number

-   AID: arc identification number

-   CSID: coordinate system identification number

-   ID1, ID2: nodes identification numbers

-   Concave: indicates that the surface external normal points to the center of arc

-   Convex: indicates that the surface external normal points in the direction opposite to the center of arc

### Example:

FlexibleArcExtrusion\_1 1 Arc 1 CS 1 Nodes 1 2 Concave

FlexibleArcExtrusion\_1 2 Arc 1 CS 1 Nodes 1 2 Convex

### Additional information:

This surface is based on the definition of an arc on a local coordinate system, as shown in Figure 23 and presented in \[8\]. The flexible extruded arc surface is designed to be attached to two nodes, defining the direction and the limits of extrusion.

<table><tbody><tr class="odd"><td><blockquote><p><span class="math inline"><em>Γ</em>(<em>ζ</em>, <em>θ</em>)<strong>=</strong><em>h</em><sub>1</sub>(<strong>Q</strong><sub>1</sub><strong>a</strong><strong>(</strong><em>θ</em><strong>)</strong><strong>+</strong><strong>x</strong><sub>1</sub>)<strong>+</strong><em>h</em><sub>2</sub>(<strong>Q</strong><sub>2</sub><strong>a</strong><strong>(</strong><em>θ</em><strong>)</strong><strong>+</strong><strong>x</strong><sub>2</sub>)</span>.</p></blockquote></td><td>( 49 )</td></tr></tbody></table>

where **Q**<sub>1</sub> and **Q**<sub>2</sub> are operators, which encompass: (i) transformation between local and global coordinate systems – from local arc definition to the desired global orientation and (ii) rotation experienced by nodes 1 and 2, respectively. Vectors **x**<sub>1</sub> and **x**<sub>2</sub> are the positions of nodes 1 and 2, respectively (considered as local origins of the local coordinate system employed to define the arc). Shape functions *h*<sub>1</sub> and *h*<sub>2</sub> are given by:

<table><tbody><tr class="odd"><td><blockquote><p><span class="math inline">$h_{1} = \frac{1}{2}\left( 1 - \zeta \right)$</span> and</p><p><span class="math inline">$h_{2} = \frac{1}{2}\left( 1 + \zeta \right)$</span>,</p></blockquote></td><td>( 50 )</td></tr></tbody></table>

which define the extrusion parameter *ζ*. The local curve parameterization is given by:

<table><tbody><tr class="odd"><td><blockquote><p><span class="math inline">$\mathbf{a}(\theta)\mathbf{=}\begin{bmatrix}
r\cos\theta + \text{XCP} \\
r\sin\theta + \text{YCP} \\
0 \\
\end{bmatrix}$</span><strong>,</strong></p></blockquote></td><td>( 51 )</td></tr></tbody></table>

as defined in Arcs section.

The extrusion direction is E3 of the local coordinate system CSID. Figure 29 illustrates the extruded surface. In this figure the arc center point **c** = (*X**C**P*,*Y**C**P*) = (*c*<sub>1</sub>, *c*<sub>2</sub>).

If the user’s choice for the external normal direction is “Convex”, it will be given by:

<table><tbody><tr class="odd"><td><blockquote><p><span class="math inline">$\mathbf{n}_{\text{ext}} = \frac{\Gamma,_{\theta}\mathbf{\times}\Gamma,_{\zeta}}{\left\| \Gamma,_{\theta}\mathbf{\times}\Gamma,_{\zeta} \right\|}$</span>.</p></blockquote></td><td>( 52 )</td></tr></tbody></table>

If the choice is “Concave”, it will be given by:

<table><tbody><tr class="odd"><td><blockquote><p><span class="math inline">$\mathbf{n}_{\text{ext}} = - \frac{\Gamma,_{\theta}\mathbf{\times}\Gamma,_{\zeta}}{\left\| \Gamma,_{\theta}\mathbf{\times}\Gamma,_{\zeta} \right\|}$</span>.</p></blockquote></td><td>( 53 )</td></tr></tbody></table>

<img src="../../assets/manual/image34.png" style="width:5.90556in;height:4.70417in" alt="" />

Figure 29 – Flexible extruded arc surface (figure from \[8\])

**Note**: this surface is currently available for using together with the contact SSSS.

## RigidArcRevolution\_1

*Creates a surface generated by the revolution of an arc about a local axis.*

### Syntax:

//Surface with external normal pointing to the center of arc

RigidArcRevolution\_1 SID Arc AID CS CSID Node NID Concave RevolutionAngle RANG FactorX XF FactorZ ZF

//Surface with external normal pointing opposite to the center of arc

RigidArcRevolution\_1 SID Arc AID CS CSID Node NID Convex RevolutionAngle RANG FactorX XF FactorZ ZF

-   SID: current surface identification number

-   AID: arc identification number

-   CSID: coordinate system identification number

-   NID: node identification number

-   Concave: indicates that the surface external normal points to the center of arc

-   Convex: indicates that the surface external normal points in the direction opposite to the center of arc

> Optional data:

-   RANG: value for the maximum revolution angle coordinate *ϕ* (default *ϕ* ranges from 0 to 2*π* – full revolution)

-   XF: ovalization coefficient for local x direction (default value is 1.0)

-   ZF: ovalization coefficient for local z direction (default value is 1.0)

### Example:

RigidArcRevolution\_1 1 Arc 1 CS 1 Node 1 Concave

### Additional information:

This surface is based on the definition of an arc on a local coordinate system, as shown in Figure 23 and presented in \[8\]. The rigid arc revolution surface is designed to be attached to a single node.

<table><tbody><tr class="odd"><td><blockquote><p><span class="math inline"><em>Γ</em>(<em>θ</em>, <em>ϕ</em>) <strong>=</strong> <strong>Q</strong><strong>a</strong>(<em>θ</em>,<em>ϕ</em>)<strong>+</strong><strong>x</strong></span>.</p></blockquote></td><td>( 54 )</td></tr></tbody></table>

where **Q** is an operator that transforms between local and global coordinate systems – from local arc definition to the desired global orientation. Vector **x** is the position of the node, considered as the origin of the local coordinate system employed to define the arc.

The local surface parameterization **a**(*θ*,*ϕ*) is given by:

<table><tbody><tr class="odd"><td><blockquote><p><span class="math inline">$\mathbf{a}\left( \theta,\phi \right)\mathbf{=}\begin{bmatrix}
\left( r\cos\theta + c_{1} \right)\left( x_{\text{factor}}\cos\phi \right) \\
r\sin\theta + c_{2} \\
 - \left( \left( r\cos\theta + c_{1} \right)\left( z_{\text{factor}}\sin\phi \right) \right) \\
\end{bmatrix}$</span><strong>,</strong></p></blockquote></td><td>( 55 )</td></tr></tbody></table>

such that the revolution axis is E2 of the local coordinate system CSID. The definition of the arc is made following the guidelines presented in in Arcs section. Revolution angular parameter is given by *ϕ* .Optional data *x*<sub>factor</sub> and *z*<sub>factor</sub> are coefficients to rule an ovalization pattern. If not input by the user, their values are considered as *x*<sub>factor</sub> = 1 and *z*<sub>factor</sub> = 1. Figure 30 illustrates the revolved surface. In this figure the arc center point **c** = (*X**C**P*,*Y**C**P*) = (*c*<sub>1</sub>, *c*<sub>2</sub>).

<img src="../../assets/manual/image35.png" style="width:3.63983in;height:3.05206in" alt="" />

Figure 30 – Rigid arc revolution surface (figure from \[8\])

If the user’s choice for the external normal direction is “Convex”, it will be given by:

<table><tbody><tr class="odd"><td><blockquote><p><span class="math inline">$\mathbf{n}_{\text{ext}} = \frac{\Gamma,_{\phi}\mathbf{\times}\Gamma,_{\theta}}{\left\| \Gamma,_{\phi}\mathbf{\times}\Gamma,_{\theta} \right\|}$</span>.</p></blockquote></td><td>( 56 )</td></tr></tbody></table>

If the choice is “Concave”, it will be given by:

<table><tbody><tr class="odd"><td><blockquote><p><span class="math inline">$\mathbf{n}_{\text{ext}} = - \frac{\Gamma,_{\phi}\mathbf{\times}\Gamma,_{\theta}}{\left\| \Gamma,_{\phi}\mathbf{\times}\Gamma,_{\theta} \right\|}$</span>.</p></blockquote></td><td>( 57 )</td></tr></tbody></table>

**Note**: this surface is currently available for using together with the contact SSSS.

##  RigidNURBS\_1

*Creates a rigid NURBS surface attached to a node and oriented according to a local coordinate system.*

### Syntax:

RigidNURBS\_1 SID CS CSID PilotNode NID CADData CDID

-   SID: current surface identification number

-   CSID: coordinate system identification number

-   NID: pilot node identification number

-   CDID: CADData identification number

### Example:

RigidNURBS\_1 1 CS 1 PilotNode 1 CADData 1

### Additional information:

RigidNURBS\_1 follows the rigid body-surface parameterization presented in \[16\] given by:

<table><tbody><tr class="odd"><td><blockquote><p><span class="math inline"><em>Γ</em>(<em>u</em>,<em>v</em>)<strong>=</strong><strong>Q</strong><sub><em>O</em></sub><strong>s</strong>(<em>u</em>,<em>v</em>)<strong>+</strong><strong>x</strong><sub><em>O</em></sub></span>,</p></blockquote></td><td>(3)</td></tr></tbody></table>

where s(*u*,*v*) is a locally-defined NURBS surface parameterization, as shown in equation (2), **Q**<sub>*O*</sub> is a rotation tensor to align the NURBS surface with the desired coordinate system defined by the CS keyword and **x**<sub>*O*</sub> to translate it to a desired location, which is the pilot node position. The expression *Γ*(*u*,*v*) depends on the model degrees of freedom, leading to the possibility of updating the rigid body surface on a transient dynamics model evolution. The pilot node translation and rotation will rule the rigid surface position/orientation.

An example of RigidNURBS\_1 is shown in Figure 31:

<img src="../../assets/manual/image36.png" style="width:2.047in;height:1.84485in" alt="" />

Figure 31 – Rigid NURBS surface (figure from \[16\])

**Note**: this surface is currently available for using together with the contact SSSS.
