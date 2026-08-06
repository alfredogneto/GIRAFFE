# Elements

*Starts a command block for creation of elements to be used to compound a finite element mesh.*

### Syntax:

Elements N

Name ID data

-   N: number of elements

-   Name: current element name

-   ID: current element identification number

-   data: current element data (depends on element resources and requirements)

### Example:

Elements 2

Beam\_1 1 Mat 1 Sec 1 CS 1 Nodes 1 2 3

Shell\_1 2 Mat 1 Sec 1 Nodes 1 2 3 4 5 6

### Additional information:

Each element is defined by a specific keyword followed by the element identification number (must be an ascending sequence starting from number one) and additional data. Each element available and its input data is explained next.

## Beam\_1

*Creates an initially straight beam finite element defined by three nodes.*

### Syntax:

Beam\_1 EID Mat MID Sec SID CS CSID Nodes ID1 ID2 ID3

-   EID: current element identification number

-   MID: material identification number

-   SID: element cross-section identification number

-   CSID: coordinate system identification number

-   ID1, ID2 and ID3: identification number of nodes defining the element

### Example:

Beam\_1 1 Mat 1 Sec 1 CS 1 Nodes 1 2 3

### Additional information:

This is a 3D beam element, with three equally spaced nodes. Three displacement and three rotation DOFs are defined for each node. Then, each element has eighteen DOFs. The Beam\_1 element uses two Gauss points for integration. The nodes used to create the Beam\_1 element have to be established by the keyword Nodes, followed by the node numbers corresponding to nodes 1, 2 and 3 in a local reference (Figure 8). The only environmental field loading that can be used with this element is the self-weight induced by the gravitational field, defined by BoolTable environment data.

<img src="../../assets/manual/image9.png" style="width:1.71274in;height:1.09189in" alt="" />

Figure 8 – Beam\_1 element local nodes numbering reference

This element internally assumes its local framework containing **e<sub>3</sub>** axis aligned with axial direction of the beam. Then, local directions **e<sub>1 </sub>**and **e<sub>2</sub>** are orthogonal to the beam direction. For establishing the cross section correct alignment one has to choose the element coordinate system such that **e<sub>3 </sub>**lies at the element axial direction and **e<sub>1</sub>** is aligned with the direction used to define the cross section (see Sections). More details about theoretical details of this beam formulation can be found in \[1\] and \[2\]. For post-processing Beam\_1 element results, one has the following sequence (to be chosen in Paraview<sup>TM</sup> post-processing):

Table 1 – Beam\_1 element results

| **Element result index** | **Meaning**                                         |
|--------------------------|-----------------------------------------------------|
| 0                        | Shear force in direction **e**<sub>1</sub>          |
| 1                        | Shear force in direction **e**<sub>2</sub>          |
| 2                        | Axial force (in direction **e**<sub>3</sub>)        |
| 3                        | Bending moment around direction **e**<sub>1</sub>   |
| 4                        | Bending moment around direction **e**<sub>2</sub>   |
| 5                        | Torsion moment (around direction **e**<sub>3</sub>) |

## Pipe\_1

*Creates an initially straight pipe finite element defined by three nodes.*

### Syntax:

Pipe\_1 EID PipeSec PSID CS CSID Nodes ID1 ID2 ID3

-   EID: current element identification number

-   PSID: pipe cross-section identification number

-   CSID: coordinate system identification number

-   ID1, ID2 and ID3: identification number of nodes defining the element

### Example:

Pipe\_1 1 PipeSec 1 CS 1 Nodes 1 2 3

### Additional information:

This is a 3D pipe element. The structural behavior is the same as Beam\_1 element. However, the input attributes are different. It is possible to make use of Pipe\_1 with environmental loading, such as weight, Morison sea current drag loading, internal and external pressure loading.

For post-processing Pipe\_1 element results, one has the following sequence (to be chosen in Paraview<sup>TM</sup> post-processing):

> Table 2 – Pipe\_1 element results

| **Element result index** | **Meaning**                                         |
|--------------------------|-----------------------------------------------------|
| 0                        | Shear force in direction **e**<sub>1</sub>          |
| 1                        | Shear force in direction **e**<sub>2</sub>          |
| 2                        | Axial force (in direction **e**<sub>3</sub>)        |
| 3                        | Bending moment around direction **e**<sub>1</sub>   |
| 4                        | Bending moment around direction **e**<sub>2</sub>   |
| 5                        | Torsion moment (around direction **e**<sub>3</sub>) |

## Shell\_1 

*Creates an initially planar shell element defined by six nodes.*

### Syntax:

Shell\_1 EID Mat MID Sec SID CS CSID Nodes ID1 ID2 ID3 ID4 ID5 ID6

-   EID: current element identification number

-   MID: material identification number

-   SID: shell section identification number

-   CSID: optional coordinate system identification number. Necessary for composite shell structures

-   ID1, ID2, ID3, ID4, ID5 and ID6: identification number of nodes defining the element

### Example:

Shell\_1 1 Mat 1 Sec 1 Nodes 1 2 3 4 5 6

### Example for composite shell structures:

Shell\_1 1 Mat 1 Sec 1 CS 1 Nodes 1 2 3 4 5 6

### Additional information:

This is a triangular shell element with six nodes. It uses three points to integrate along the element area. The sequence of nodes must be provided according to the numbering sequence shown in Figure 9. Note that the direction chosen to increase the number of nodes implicitly defines the normal of the shell element, according to the right-hand rule. The user can establish the external normal direction **n** by performing the cross product between the vectors **v**<sub>**1**</sub> **=** **(****P**<sub>**2**</sub>**−****P**<sub>**1**</sub>**)** and **v**<sub>**2**</sub>**=****(****P**<sub>**3**</sub>**−****P**<sub>**1**</sub>**)**, such that $\mathbf{n =}\frac{\mathbf{v}_{\mathbf{1}}\mathbf{\times}\mathbf{v}_{\mathbf{2}}}{\left\Vert  \mathbf{v}_{\mathbf{1}}\mathbf{\times}\mathbf{v}_{\mathbf{2}} \right\Vert }$.This external normal direction is used to define pressure loading on shell elements.

<img src="../../assets/manual/image10.png" style="width:1.78472in;height:1.38861in" alt="" />

Figure 9 – Shell\_1 element local nodes numbering reference

More details about theoretical details of this shell formulation can be found in \[3\].

To establish a local coordinate system, Giraffe uses the reference configuration of the shell element. The direction **e**<sub>3</sub><sup>*r*</sup> is the normal direction of the shell, at reference configuration. The local **e**<sub>1</sub><sup>*r*</sup> is defined by the global x direction projection on the shell reference plane. If this projection is null, then **e**<sub>1</sub><sup>*r*</sup> is defined by the global y direction projection on the shell plane. Finally, **e**<sub>2</sub><sup>*r*</sup> = **e**<sub>3</sub><sup>*r*</sup> × **e**<sub>1</sub><sup>*r*</sup>.

When dealing with composite shell structures it is necessary to set a local coordinate system. Each local coordinate system will define the principal material axes orientation and the stacking order of the laminas in the shell element.

For post-processing Shell\_1 element results, one has the following sequence (to be chosen in Paraview<sup>TM</sup> post-processing):

> Table 3 – Shell\_1 element results

| **Element result index** | **Meaning**                                                                                   |
|--------------------------|-----------------------------------------------------------------------------------------------|
| 0                        | Force in direction **e**<sub>1</sub> (cutting plane with normal direction **e**<sub>1</sub>)  |
| 1                        | Force in direction **e**<sub>2</sub> (cutting plane with normal direction **e**<sub>1</sub>)  |
| 2                        | Force in direction **e**<sub>3</sub> (cutting plane with normal direction **e**<sub>1</sub>)  |
| 3                        | Moment in direction **e**<sub>1</sub> (cutting plane with normal direction **e**<sub>1</sub>) |
| 4                        | Moment in direction **e**<sub>2</sub> (cutting plane with normal direction **e**<sub>1</sub>) |
| 5                        | Moment in direction **e**<sub>3</sub> (cutting plane with normal direction **e**<sub>1</sub>) |
| 6                        | Force in direction **e**<sub>1</sub> (cutting plane with normal direction **e**<sub>2</sub>)  |
| 7                        | Force in direction **e**<sub>2</sub> (cutting plane with normal direction **e**<sub>2</sub>)  |
| 8                        | Force in direction **e**<sub>3</sub> (cutting plane with normal direction **e**<sub>2</sub>)  |
| 9                        | Moment in direction **e**<sub>1</sub> (cutting plane with normal direction **e**<sub>2</sub>) |
| 10                       | Moment in direction **e**<sub>2</sub> (cutting plane with normal direction **e**<sub>2</sub>) |
| 11                       | Moment in direction **e**<sub>3</sub> (cutting plane with normal direction **e**<sub>2</sub>) |

## Mass\_1

*Creates a single-node lumped mass element.*

### Syntax:

Mass\_1 EID Mass MV Node NID

-   EID: current element identification number

-   MV: mass value

-   NID: identification number of the node defining the element

### Example:

Mass\_1 1 Mass 150.3 Node 1

### Additional information:

This is a lumped mass element. It can be used to model a portion of mass not included in the finite element model, but that may affect the system response due to gravitational field loads and/or inertial loads.

This element has no direct influence in system’s stiffness, but only in the mass matrix coefficients related to translational DOFs and external loads, corresponding to inertial and gravitational field.

> This element has no specific results for post-processing.

## SpringDashpot\_1

*Creates a two-node spring and dashpot element.*

### Syntax:

SpringDashpot\_1 EID Stiffness SV Damping DV Nodes ID1 ID2

-   EID: current element identification number

-   SV: stiffness value

-   DV: damping value

-   ID1 and ID2: identification number of the nodes defining the element

### Example:

SpringDashpot\_1 1 Stiffness 200.2 Damping 1.34 Nodes 1 2

### Additional information:

This is a spring and dashpot element. It can be used in applications were the stiffness/damping coefficients are known a priori.

<img src="../../assets/manual/fig10_springdashpot.png" style="max-width:100%" alt="" />

Figure 10 – SpringDashpot\_1 element local nodes numbering reference

It is a two-node element (see local nodes numbering in Figure 10), with a linear stiffness “k” – proportional to the relative displacement of the nodes – and a linear damping “c” – proportional to the relative velocity of the nodes, projected on the direction of the line that connects the two nodes. The element is geometrically nonlinear. Then, the direction affected by the stiffness/damping follows the current position of the element nodes. It can handle large rigid body rotations and translations.

> Table 4 – SpringDashpot\_1 element results

| **Element result index** | **Meaning**                                                                                                                                                 |
|--------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 0                        | Spring elongation (positive value: augmenting the length, negative value: diminishing the length)                                                           |
| 1                        | Elastic force from spring (positive value: tension, negative value: compression)                                                                            |
| 2                        | Damping force from dashpot (positive value: relative velocity increasing the damper length, negative value: relative velocity decreasing the damper length) |

## RigidBody\_1

*Creates a single-node rigid body element.*

### Syntax:

RigidBody\_1 EID RigidBodyData RBID CS CID Node NID

-   EID: current element identification number

-   RBID: rigid body data identification number

-   CID: coordinate system identification number

-   NID: identification number of the node defining the element

### Example:

RigidBody\_1 1 RigidBodyData 1 CS 1 Node 1

### Additional information:

This is a rigid body element with mass and inertia properties. It can be used in multibody systems to represent relatively stiff components.

It is a single-node element that accounts for mass and inertia properties from a 3D solid body (provided within RigidBodyData). Each rigid body element has a unique identification number that follows the keyword “RigidBody\_1”. Each element receives its properties and therefore it is necessary to indicate the identification number of the RigidBodyData.

Rigid bodies are oriented in Giraffe according to local Coordinate System (CS). The objective of this local CS is to orient the global axes from the CAD file in the Giraffe three-dimensional space. The CAD origin is placed at the Rigid Body node. This is illustrated in Figure 11.

<img src="../../assets/manual/image12.png" style="width:5.51181in;height:2.75283in" alt="" />

Figure 11 – Orientation of RigidBody\_1 elements

In Figure 11 (a) a generic CAD model is represented. It was modelled in the yz plane. Suppose the orientation shown in Figure 11 (b) is desired in Giraffe, then the local CS represented by E1, E2 and E3 has to be specified. For example, in Figure 11 (b) the axes were oriented according to the following CS:

<table><tbody><tr class="odd"><td><p>CS 1</p><p>CSYS 1 E1 1 0 0 E3 0 1 0</p></td></tr></tbody></table>

This means that the x axis of the CAD file was aligned with the X axis in Giraffe, and the z axis was aligned with the Y axis. Note that the main objective of this CS is to orient the geometry in Giraffe. It is a very important feature since it is directly related to the inertia properties and it is used for postprocessing purposes (to visualize geometry in ParaView<sup>TM</sup>).

The Node used to create RigidBody\_1 can be placed anywhere. Therefore, the CAD model has to be generated in a way that its origin coincides with the position of this node. It is very convenient to create the RigidBody\_1 node at the center of mass or at some place that will be constrained during the simulation, Figure 12 and Figure 13 clarify this question.

<img src="../../assets/manual/image13.png" style="width:5.51181in;height:2.81314in" alt="" />

Figure 12 – RigidBody\_1 positioning in Giraffe platform: Example 1

<img src="../../assets/manual/image14.png" style="width:5.51181in;height:2.79972in" alt="" />

Figure 13 – RigidBody\_1 positioning in Giraffe platform: Example 2

RigidBody\_1 is usually used in conjunction with special constraints, especially rigid node set. For example, one may be interested in monitoring kinematic quantities of a point located anywhere in the body. In order to do that, it is necessary to create a node at the position of interest and then to define a rigid node set from the RigidBody\_1 Node (pilot node) to the node of interest (slave node). Nodes at different locations of the body can be added to rigid node set. The kinematic quantities can be obtained via Monitors.

To obtain quantities such as Kinetic Energy (T), Linear Momentum (L) and Angular Momentum about the center of mass (HG) one can request element Monitors.

This element has no specific results for post-processing using Paraview<sup>TM</sup>.

## Truss\_1

*Creates a two-node truss element.*

### Syntax:

This element type has two possible syntax entries shown below:

Truss\_1 EID Mat MID Sec SID Nodes ID1 ID2

Truss\_1 EID PipeSec PSID Nodes ID1 ID2

-   EID: current element identification number

-   MID: material identification number

-   CSID: cross-section data identification number

-   PSID: pipe cross-section data identification number

-   ID1 and ID2: identification number of the nodes defining the element

### Example:

Truss\_1 1 Mat 1 Sec 1 Nodes 1 2

Truss\_1 2 PipeSec 1 Nodes 2 3

### Additional information:

This is a 3D truss element. There are two possible entries for establishing such element: by material and cross-section data or by pipe cross-section data. When establishing material data, this element can handle large strain. Material model may be chosen between:

1.  linear-elastic (Hooke)

2.  elastic-plastic with isotropic hardening

Alternatively, if pipe cross-section data is defined, the element employs only axial stiffness and adopts is as constant. Additionally, it may be used to define environmental loading such as weight and Morison sea current drag and added mass loading.

For post-processing Truss\_1 element results one has the following sequence (to be chosen in Paraview<sup>TM</sup> post-processing):

> Table 5 – Truss\_1 element results

| **Element result index** | **Meaning**             |
|--------------------------|-------------------------|
| 0                        | Element axial tension   |
| 1                        | Cross-section area      |
| 2                        | Element length          |
| 3                        | Plastic deformed length |
| 4                        | Kirchhoff stress        |

## TwoNodeConnector\_1

*Creates a flexible connection between two nodes with defined stiffness and damping matrices.*

### Syntax:

TwoNodeConnector\_1 EID CS CSID Nodes ID1 ID2

StiffnessData

//Stiffness matrix data

DampingData

//Damping matrix data

-   EID: current element identification number

-   CSID: coordinate system identification number

-   ID1, ID2: identification number of nodes defining the element

### Example:

TwoNodeConnector\_1 6 CS 1 Nodes 1 12

StiffnessData

10000 0 0 0 0 0

0 10000 0 0 0 0

0 0 10000 0 0 0

0 0 0 10000 0 0

0 0 0 0 10000 0

0 0 0 0 0 10000

DampingData

100 0 0 0 0 0

0 100 0 0 0 0

0 0 100 0 0 0

0 0 0 100 0 0

0 0 0 0 100 0

0 0 0 0 0 100

### Additional information:

This element aims at creating a flexible connection between two nodes A and B. The relative movement between nodes A and B generates elastic and damping forces, according to user-defined stiffness/damping matrices. Figure 14 shows the degrees of freedom considered for each node, where (1,2,3) refer to **u**<sub>*A*</sub> and (7,8,9) refer to **u**<sub>*B*</sub>, which are total displacements of nodes A and B, respectively. Analogously, one may find (4,5,6) referring to **α**<sub>*A*</sub> and (10,11,12) to **α**<sub>*B*</sub>, which are Rodrigues rotation vectors of nodes A and B, respectively.

<img src="../../assets/manual/image15.png" style="width:2.52778in;height:1.12014in" alt="" />

Figure 14 – TwoNodeConnector\_1 element degrees of freedom

Let **K** be a 6x6 stiffness matrix and **C** be a 6x6 damping matrix, defined by the user in element input data. Indexes 1-6 refer to relative displacements (1,2,3) and relative rotations (4,5,6) between nodes A and B, following the directions illustrated in Figure 14 and considering the local coordinate system directions (CS input). With that, the element contribution to the model weak form is divided in two terms, as follows:

|                                                         |       |
|---------------------------------------------------------|-------|
| $\delta W_{1} = \mathbf{K}\left( \begin{bmatrix} \mathbf{u}_{A} \\ \mathbf{\alpha}_{A} \\ \end{bmatrix} - \begin{bmatrix} \mathbf{u}_{B} \\ \mathbf{\alpha}_{B} \\ \end{bmatrix} \right) \cdot \left( \begin{bmatrix} {\delta\mathbf{u}}_{A} \\ {\delta\mathbf{\alpha}}_{A} \\ \end{bmatrix} - \begin{bmatrix} {\delta\mathbf{u}}_{B} \\ {\delta\mathbf{\alpha}}_{B} \\ \end{bmatrix} \right)$,                                | ( 1 ) |

stemming from the virtual work of internal (elastic) loads and

|                                                         |       |
|---------------------------------------------------------|-------|
| $\delta W_{2} = \mathbf{C}\left( \begin{bmatrix} {\dot{\mathbf{u}}}_{A} \\ {\dot{\mathbf{\alpha}}}_{A} \\ \end{bmatrix} - \begin{bmatrix} {\dot{\mathbf{u}}}_{B} \\ {\dot{\mathbf{\alpha}}}_{B} \\ \end{bmatrix} \right) \cdot \left( \begin{bmatrix} {\delta\mathbf{u}}_{A} \\ {\delta\mathbf{\alpha}}_{A} \\ \end{bmatrix} - \begin{bmatrix} {\delta\mathbf{u}}_{B} \\ {\delta\mathbf{\alpha}}_{B} \\ \end{bmatrix} \right)$,                                | ( 2 ) |

stemming from the virtual work of damping (viscous) loads. The quantities *δ**W*<sub>1</sub> and *δ**W*<sub>2</sub> may be re-written as:

|                                                                          |       |
|--------------------------------------------------------------------------|-------|
| $\delta W_{1} = \begin{bmatrix} \begin{matrix} {\delta\mathbf{u}}_{A}^{T} & {\delta\mathbf{\alpha}}_{A}^{T} \\ \end{matrix} & \begin{matrix} {\delta\mathbf{u}}_{B}^{T} & {\delta\mathbf{\alpha}}_{B}^{T} \\ \end{matrix} \\ \end{bmatrix}\begin{bmatrix} \mathbf{+ K} & \mathbf{-}\mathbf{K} \\ \mathbf{-}\mathbf{K} & \mathbf{+ K} \\ \end{bmatrix}\begin{bmatrix} \begin{matrix} \mathbf{u}_{A} \\ \mathbf{\alpha}_{A} \\ \end{matrix} \\ \begin{matrix} \mathbf{u}_{B} \\ \mathbf{\alpha}_{B} \\ \end{matrix} \\ \end{bmatrix}$,                                                          | ( 3 ) |

and

|                                                                          |       |
|--------------------------------------------------------------------------|-------|
| $\delta W_{2} = \begin{bmatrix} \begin{matrix} {\delta\mathbf{u}}_{A}^{T} & {\delta\mathbf{\alpha}}_{A}^{T} \\ \end{matrix} & \begin{matrix} {\delta\mathbf{u}}_{B}^{T} & {\delta\mathbf{\alpha}}_{B}^{T} \\ \end{matrix} \\ \end{bmatrix}\begin{bmatrix} \mathbf{+ C} & \mathbf{-}\mathbf{C} \\ \mathbf{-}\mathbf{C} & \mathbf{+ C} \\ \end{bmatrix}\begin{bmatrix} \begin{matrix} {\dot{\mathbf{u}}}_{A} \\ {\dot{\mathbf{\alpha}}}_{A} \\ \end{matrix} \\ \begin{matrix} {\dot{\mathbf{u}}}_{B} \\ {\dot{\mathbf{\alpha}}}_{B} \\ \end{matrix} \\ \end{bmatrix}$.                                                          | ( 4 ) |

Therefore, one may see this element as a contribution involving 12 degrees of freedom, connecting nodes A and B by a user-defined stiffness and damping which rule the relative movement between these nodes.

**Note**: No updates in directions of stiffness/damping coefficients are considered. They follow the directions defined by the local coordinate system CS, independently of displacements/rotations experienced by nodes A and B.
