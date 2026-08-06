# Special Constraints

*Starts a command block for creation of special constraints.*

### Syntax:

SpecialConstraints N

Name ID data

-   N: number of special constraints

-   Name: current special constraint name

-   ID: current special constraint identification number

-   data: current special constraint data (depends on special constraint resources and requirements)

### Example:

SpecialConstraints 1

SameDisplacement 1 Nodes 1 2 BoolTable 1

### Additional information:

Each special constraint is defined by a specific keyword followed by the special constraint identification number (must be an ascending sequence starting from number one) and additional data. Each special constraint available and its input data is explained next. Models of special constraints are presented in \[7\] and \[8\].

Differently from Constraints keyword, SpecialConstraints are established keeping original system DOFs and including additional unknowns to enforce desired constraints (Lagrange Multipliers). Thus, inclusion of special constraints in the model will increase the number of unknows.

## SameDisplacement

*Creates a “same displacement” special constraint.*

### Syntax:

SameDisplacement SCID Nodes ID1 ID2 BoolTable BDSC

-   SCID: current special constraint identification number

-   ID1 and ID2: nodes identification numbers

-   BDSC: BoolTable data for current special constraint

### Example:

SameDisplacement 1 Nodes 1 2 BoolTable 1

### Additional information:

This constraint is used to enforce that the two selected nodes will present the same displacements (but not necessarily the same rotation). It can be used to represent a spherical joint in a mechanism, for example. The selected nodes should be initially coincident for that purpose.

BoolTable keyword is optional. It permits creating a scenario in which the special constraint is turned on/off along solution steps (see BoolTable). If the user does not include BoolTable, Giraffe assumes that the special constraint will be considered for all solution steps as turned on.

**Remark**: if this special constraint is used in dynamic simulations, velocity initial conditions are to be set only to the first node. If one sets different velocity initial conditions for both nodes, the first node velocity initial conditions are considered for both nodes. The second node velocity initial conditions are ignored.

## SameRotation

*Creates a “same rotation” special constraint.*

### Syntax:

SameRotation SCID Nodes ID1 ID2 BoolTable BDSC

-   SCID: current special constraint identification number

-   ID1 and ID2: nodes identification numbers

-   BDSC: BoolTable data for current special constraint

### Example:

SameRotation 1 Nodes 1 2 BoolTable 1

### Additional information:

This constraint is used to enforce that the two selected nodes will present the same rotation (but not necessarily the same displacement).

BoolTable keyword is optional. It permits creating a scenario in which the special constraint is turned on/off along solution steps (see BoolTable). If the user does not include BoolTable, Giraffe assumes that the special constraint will be considered for all solution steps as turned on.

**Remark**: if this special constraint is used in dynamic simulations, angular velocity initial conditions are to be set only to the first node. If one sets different angular velocity initial conditions for both nodes, the first node angular velocity initial conditions are considered for both nodes. The second node angular initial velocity conditions are ignored.

## RigidNodeSet

*Creates a “rigid node set” special constraint.*

### Syntax:

RigidNodeSet SCID PilotNode PID NodeSet NSID BoolTable BDSC

-   SCID: current special constraint identification number

-   PID: pilot node identification number

-   NSID: node set identification number

-   BDSC: BoolTable data for current special constraint

### Example:

RigidNodeSet 1 PilotNode 1 NodeSet 1 BoolTable 1

### Additional information:

This constraint is used to establish a rigid region, which is formed by the connection of all the nodes in the selected node set. The pilot node will always present 6 degrees of freedom, which will rule the movement of all the nodes in the node set. The pilot node can be a node inside the node set or an independent node.

BoolTable keyword is optional. It permits creating a scenario in which the special constraint is turned on/off along solution steps (see BoolTable). If the user does not include BoolTable, Giraffe assumes that the special constraint will be considered for all solution steps as turned on.

This special constraint imposes that general rigid body movement may play a role for the set of nodes assigned. Then, one can handle large displacements and large rotations with no kinematic limitations.

**Remark**: if this special constraint is used in dynamic simulations, displacement and rotation initial conditions are to be set only for the pilot node. Giraffe automatically evaluates, using rigid body’s equations, the proper initial conditions for each node of the node set. If one sets arbitrary initial conditions for the nodes, which are not compatible to pilot node’s conditions, these are ignored.

## HingeJoint

*Creates a “hinge joint” special constraint.*

### Syntax:

HingeJoint SCID Nodes ID1 ID2 CS CSID LinearStiffness LSV LinearDamping LDV QuadraticDamping QDV BoolTable BDSC

-   SCID: current special constraint identification number

-   ID1 and ID2: nodes identification numbers

-   CSID: coordinate system identification number

-   LSV: linear stiffness coefficient value

-   LDV: linear damping coefficient value

-   QDV: quadratic damping coefficient value

-   BDSC: BoolTable data for current special constraint

### Example:

HingeJoint 1 Nodes 1 2 CS 1 LinearStiffness 0.0 LinearDamping 0.0 QuadraticDamping 0.0 BoolTable 1

### Additional information:

This constraint is used to represent a hinge joint. It enforces that the two selected nodes will present the same displacements, as in SameDisplacement constraint. Furthermore, the direction **e**<sub>**3**</sub> from the chosen coordinate system CS will represent the direction of free relative rotation between the chosen nodes (hinge axis). The direction **e**<sub>**3**<sub>**A**</sub></sub> is updated, since it is attached to the rotation of the node A (first node defined for the hinge joint). The directions of **e**<sub>**1**<sub>**B**</sub></sub> and **e**<sub>**2**<sub>**B**</sub></sub> are also updated, following the node B rotations (second node defined for the hinge joint).

At the beginning, we assume that both nodes coordinate systems lie at the same directions. During the simulation, the hinge joint ensures that **e**<sub>**3**<sub>**A**</sub></sub>**≡****e**<sub>**3**<sub>**B**</sub></sub>. For that, two constraints *r*<sub>1</sub> and *r*<sub>2</sub> for rotations are enforce by:

| *r*<sub>1</sub>**=****e**<sub>**3**<sub>**A**</sub></sub> ⋅ **e**<sub>**1**<sub>**B**</sub></sub>**=**0  | ( 36 ) |
|---------------------------------------------------------------------------------------------------------|--------|
|                                                                                                         |        |
| *r*<sub>2</sub>**=****e**<sub>**3**<sub>**A**</sub></sub> ⋅ **e**<sub>**2**<sub>**B**</sub></sub>**=**0  | ( 37 ) |
|                                                                                                         |        |

<img src="../../assets/manual/image23.png" style="width:1.67277in;height:1.31507in" alt="" /> (a) <img src="../../assets/manual/image24.png" style="width:1.61243in;height:1.28017in" alt="" /> (b)

Figure 20 – (a) Example of a hinge joint between two beams. (b) The coordinate systems of nodes A and B after nodes rotation.

Figure 20(a) shows an example of hinge joint, located at nodes A and B (coincident). The coordinate systems of both nodes is initially the same, and defined by the keyword CS. After some movement, the rotation at nodes A and B may differ, such that the systems (which follows nodes rotations) are as in Figure 20(b). Note that the constraints ( 36 ) and ( 37 ) are obeyed in such transformation.

The parameter LinearStiffness permits entering a stiffness coefficient, such that the hinge presents a torsion stiffness. The moment *M*<sub>spring</sub> generated by the torsion stiffness spring is given by:

|                                             |        |
|---------------------------------------------|--------|
| *M*<sub>spring</sub> = *K*<sub>*θ*</sub>*θ*  | ( 38 ) |

where *K*<sub>*θ*</sub> is the linear stiffness coefficient and *θ* is the accumulated angle of relative rotation around the hinge direction. The variable *θ* may represent finite rotations, involving many turns around the hinge axis.

Analogously, one can enter linear and quadratic damping coefficients by the keywords LinearDamping and QuadraticDamping. Then, the moment *M*<sub>damper</sub> can be evaluated by:

| *M*<sub>damper</sub> = *C*<sub>1*θ*</sub>(**ω**<sub>*A*</sub>−**ω**<sub>*B*</sub>) ⋅ **e**<sub>3*A*</sub> + *C*<sub>2*θ*</sub>∥(**ω**<sub>*A*</sub>−**ω**<sub>*B*</sub>)∥(**ω**<sub>*A*</sub>−**ω**<sub>*B*</sub>) ⋅ **e**<sub>3*A*</sub>  | ( 39 ) |
|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------|
|                                                                                                                                                                                                                                           |        |

where *C*<sub>1*θ*</sub> and *C*<sub>2*θ*</sub> are respectively linear and quadratic damping coefficients and **ω**<sub>*A*</sub> and **ω**<sub>*B*</sub> are the instantaneous angular velocities of nodes A and B, respectively.

**Remark**: If one sets different velocity initial conditions for both nodes, the first node velocity initial conditions are considered for both nodes. The second node velocity initial conditions are, then, ignored. Furthermore, angular velocity conditions may have independent value for both nodes. However, their compatibility is done by:

1.  Giraffe evaluates the direction of the hinge axis **e**<sub>3*A*</sub> and projects the angular velocity initial condition of both nodes on this direction;

2.  The components that lie in the direction of the hinge axis **e**<sub>3*A*</sub>, are independent in both nodes;

The components that lie in direction orthogonal to the hinge axis are set, taking the first node values and copying then to the second one’s, which has its values ignored.

BoolTable keyword is optional. It permits creating a scenario in which the special constraint is turned on/off along solution steps (see BoolTable). If the user does not include BoolTable, Giraffe assumes that the special constraint will be considered for all solution steps as turned on.

## UniversalJoint

*Creates a “universal joint” (Cardan) special constraint.*

### Syntax:

UniversalJoint SCID Nodes ID1 ID2 CSA CSAID CSB CSBID BoolTable BDSC

-   SCID: current special constraint identification number

-   ID1 and ID2: nodes identification numbers

-   CSAID and CSBID: coordinate systems identification numbers

-   BDSC: BoolTable data for current special constraint

### Example:

UniversalJoint 1 Nodes 1 2 CSA 1 CSB 2 BoolTable 1

### Additional information:

This constraint is used to represent a universal (cardan) joint. It enforces that the two selected nodes will present the same displacements, as in SameDisplacement constraint. The selected nodes should be coincident for that purpose. Additionally, two coordinate systems are defined, associated with the movement of the first and the second defined nodes (A and B, respectively). The directions **e**<sub>**3**<sub>**A**</sub></sub> and **e**<sub>**3**<sub>**B**</sub></sub> should be used to represent the directions of the axles that are connected using the joint. Both coordinate systems CSA and CSB directions are updated during simulation, according to the rotations experienced by nodes A and B. The constraint imposed to the represent the rotation transmission of cardan joint is:

| *r*<sub>3</sub>**=****e**<sub>**1**<sub>**A**</sub></sub> ⋅ **e**<sub>**2**<sub>**B**</sub></sub>**=**0  | ( 40 ) |
|---------------------------------------------------------------------------------------------------------|--------|
|                                                                                                         |        |

<img src="../../assets/manual/image25.png" style="width:2.10422in;height:1.19289in" alt="" />(a) <img src="../../assets/manual/image26.png" style="width:1.01319in;height:0.94167in" alt="" /> (b)

Figure 21 – (a) Example of a universal joint between two beams. (b) The coordinate systems of nodes A and B.

Figure 21(a) shows an example of cardan joint, located at nodes A and B (coincident). The coordinate systems of both nodes are different and **must obey** at the beginning of the simulation the constraint ( 40 ), otherwise an error message will be shown in Giraffe output window.

**Remark**: If one sets different velocity initial conditions for both nodes, the first node velocity initial conditions are considered for both nodes. The second node velocity initial conditions are ignored. Furthermore, angular velocity conditions may have independent value for both nodes. However, their compatibility is done by:

1.  Giraffe evaluates the directions **e**<sub>**3**<sub>**A**</sub></sub> and **e**<sub>**3**<sub>**B**</sub></sub>, and projects the angular velocity initial condition of both nodes to respective directions;

2.  The component of angular velocity of node A that lies in direction **e**<sub>**3**<sub>**A**</sub></sub> is also imposed to direction **e**<sub>**3**<sub>**B**</sub></sub>;

3.  The components of angular velocity in other directions are independent.

BoolTable keyword is optional. It permits creating a scenario in which the special constraint is turned on/off along solution steps (see BoolTable). If the user does not include BoolTable, Giraffe assumes that the special constraint will be considered for all solution steps as turned on.

## TranslationalJoint

*Creates a “translational joint” special constraint.*

### Syntax:

TranslationalJoint SCID Nodes ID1 ID2 RotationNode ID3 CS CSID BoolTable BDSC

-   SCID: current special constraint identification number

-   ID1 and ID2: nodes identification numbers (to be connected by translational joint)

-   ID3: rotation node identification number (to rule translational joint direction update)

-   CSID: coordinate system identification number

-   BDSC: BoolTable data for current special constraint

### Example:

TranslationalJoint 1 Nodes 1 2 RotationNode 3 CS 1 BoolTable 1

### Additional information:

This constraint is used to represent a translational joint between the first and second chosen nodes (nodes A and B, identified by ID1 and ID2). It enforces that the two selected nodes will present relative displacements only along a direction **e**<sub>**3**</sub>. The direction **e**<sub>**3**</sub> is taken from the chosen coordinate system CS and is updated along simulation according to the rotations experienced by the rotation node assigned (ID3). The constraints enforced are the following:

| *r*<sub>1</sub>**=****e**<sub>**1**</sub> ⋅ (**u**<sub>**A**</sub>**−****u**<sub>**B**</sub>)**=**0  | ( 41 ) |
|-----------------------------------------------------------------------------------------------------|--------|
|                                                                                                     |        |
| *r*<sub>2</sub>**=****e**<sub>**2**</sub> ⋅ (**u**<sub>**A**</sub>**−****u**<sub>**B**</sub>)**=**0  | ( 42 ) |

This type of joint is very useful for establishing suspension systems, as in the example:

<img src="../../assets/manual/fig22_suspension.png" style="max-width:100%" alt="" />

Figure 22 – Example of a suspension system

In this example a spring/dashpot element is defined between nodes A and B. Also, one has a mass element defined at node B, while node A is embedded in frame structure, as a part meshed using beam elements. In order to avoid undesirable rotations of the spring/dashpot one may establish a translational joint constraint between nodes A and B, thus considering gobal direction Y as the initial direction for **e**<sub>**3**</sub>. In this case node A should be used as RotationNode. With that, only relative displacements along current direction **e**<sub>**3**</sub> are permitted and the suspension system may behave as desirable. Nodes A and B may translate or rotate with the frame structure, as direction **e**<sub>**3**</sub> is updated, accordingly.

BoolTable keyword is optional. It permits creating a scenario in which the special constraint is turned on/off along solution steps (see BoolTable). If the user does not include BoolTable, Giraffe assumes that the special constraint will be considered for all solution steps as turned on.
