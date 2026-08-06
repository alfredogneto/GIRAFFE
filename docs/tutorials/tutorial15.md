# Tutorial 15: Lateral buckling of a pipeline laid on a flat rigid surface

!!! example "Input files"
    [`tutorial15.inp`](../assets/tutorials/inputs/tutorial15/tutorial15.inp)
    (also in the repository under `inputs/tutorial15/`)

Input file name: tutorial15.inp

*Developed by Marina Vendl Craveiro & Alfredo Gay Neto*

## Objective

This tutorial presents a geometrically-nonlinear simulation of a pipeline laid on a flat rigid surface with an isolated oscillatory imperfection. The aim is to simulate the laying process of the pipeline on the surface and evaluate the buckling of the structure when it is subjected to an internal pressure load.

<img src="../../assets/tutorials/image120.png" style="width:3.18898in;height:1.94488in" alt="" /> (a) <img src="../../assets/tutorials/image27.png" style="width:1.1875in;height:1.54599in" alt="" />(b)

Figure 15‑1 – (a) Model of the pipeline laid on a surface with imperfection (b) Cross section of the pipeline subjected to internal pressure field.

## Description

The pipeline AB (2000 m length) is meshed considering 2001 equally-spaced nodes and 1000 Pipe\_1 elements. The flat surfaces and the oscillatory imperfection are modeled using rigid oscillatory surfaces. The pipeline is initially positioned at the level of the top of the oscillatory imperfection. In order to establish contact with the rigid surface, the pipeline is released under the action of its self-weight. The interaction between the pipeline and its support is described by the NSSS contact model. Note that laying the pipeline on the oscillatory imperfection generates a vertical imperfection. Here, in addition to the vertical imperfection, a lateral imperfection is also imposed by means of a prescribed displacement at the central node of the pipeline. The laying process and the imposition of the lateral imperfection are done through static solution steps and lead to an initial imperfect pipeline configuration that is then subjected dynamically to the internal pressure load. The internal pressure is responsible for triggering the buckling of the pipeline.

## Input file step-by-step

### Creation of nodes

The finite element mesh contains 2001 nodes. Other 3 nodes are created to be pilot nodes of rigid oscillatory surfaces. So, the total of nodes is 2004.

<table><tbody><tr class="odd"><td><p>//Creation of nodes</p><p>Nodes 2004</p><p>Node 1 -1000 0.25 0</p><p>Node 2 -999 0.25 0</p><p>Node 3 -998 0.25 0</p><p>Node 4 -997 0.25 0</p><p>Node 5 -996 0.25 0</p><p>Node 6 -995 0.25 0</p><p>Node 7 -994 0.25 0</p><p>Node 8 -993 0.25 0</p><p>Node 9 -992 0.25 0</p><p>Node 10 -991 0.25 0</p><p>Node 11 -990 0.25 0</p><p>Node 12 -989 0.25 0</p><p>Node 13 -988 0.25 0</p><p>Node 14 -987 0.25 0</p><p>Node 15 -986 0.25 0</p><p>Node 16 -985 0.25 0</p><p>Node 17 -984 0.25 0</p><p>Node 18 -983 0.25 0</p><p>Node 19 -982 0.25 0</p><p>Node 20 -981 0.25 0</p><p>…</p><p>Node 1982 981 0.25 0</p><p>Node 1983 982 0.25 0</p><p>Node 1984 983 0.25 0</p><p>Node 1985 984 0.25 0</p><p>Node 1986 985 0.25 0</p><p>Node 1987 986 0.25 0</p><p>Node 1988 987 0.25 0</p><p>Node 1989 988 0.25 0</p><p>Node 1990 989 0.25 0</p><p>Node 1991 990 0.25 0</p><p>Node 1992 991 0.25 0</p><p>Node 1993 992 0.25 0</p><p>Node 1994 993 0.25 0</p><p>Node 1995 994 0.25 0</p><p>Node 1996 995 0.25 0</p><p>Node 1997 996 0.25 0</p><p>Node 1998 997 0.25 0</p><p>Node 1999 998 0.25 0</p><p>Node 2000 999 0.25 0</p><p>Node 2001 1000 0.25 0</p><p>Node 2002 0 0.125 0</p><p>Node 2003 -753.125 0 0</p><p>Node 2004 753.125 0 0</p></td></tr></tbody></table>

### Creation of node sets

Ten node sets are defined to establish loads, displacements and constraints. Node set 1 groups all the pipeline nodes and is used to create the NSSS contact between the pipeline and the surfaces. Node sets 2 and 4 are associated with the pipeline end nodes. Node set 3 contains the central node of the pipeline, which establishes contact with the top of the oscillatory imperfection and at which the lateral displacement is prescribed. Node sets 5, 6 and 7 contain pilot nodes of rigid oscillatory surfaces.

In order to avoid rigid body movements during the laying process of the pipeline, the release of the pipeline nodes on the rigid surfaces is done in more than one step, also requiring the grouping of the pipeline nodes into three different node sets. Node set 8 contains 101 central nodes and node sets 9 and 10 contain the other nodes adjacent to node set 8, except for the pipeline end nodes (which are included in node sets 2 and 4).

<table><tbody><tr class="odd"><td><p>//Creation of node sets</p><p>NodeSets 10</p><p>NodeSet 1 Nodes 2001 Sequence Initial 1 Increment 1</p><p>NodeSet 2 Nodes 1 List 1</p><p>NodeSet 3 Nodes 1 List 1001</p><p>NodeSet 4 Nodes 1 List 2001</p><p>NodeSet 5 Nodes 1 List 2002</p><p>NodeSet 6 Nodes 1 List 2003</p><p>NodeSet 7 Nodes 1 List 2004</p><p>NodeSet 8 Nodes 101 Sequence Initial 951 Increment 1</p><p>NodeSet 9 Nodes 949 Sequence Initial 2 Increment 1</p><p>NodeSet 10 Nodes 949 Sequence Initial 1052 Increment 1</p></td></tr></tbody></table>

### Creation of elements

A total of 1000 elements is created. All the elements have the same cross section (identification number 1) and the same coordinate system (identification number 1). Each Pipe\_1 element is defined by 3 nodes, which must be equally spaced and contained in a straight line.

<table><tbody><tr class="odd"><td><p>//Creation of elements</p><p>Elements 1000</p><p>Pipe_1 1 PipeSec 1 CS 1 Nodes 1 2 3</p><p>Pipe_1 2 PipeSec 1 CS 1 Nodes 3 4 5</p><p>Pipe_1 3 PipeSec 1 CS 1 Nodes 5 6 7</p><p>Pipe_1 4 PipeSec 1 CS 1 Nodes 7 8 9</p><p>Pipe_1 5 PipeSec 1 CS 1 Nodes 9 10 11</p><p>Pipe_1 6 PipeSec 1 CS 1 Nodes 11 12 13</p><p>Pipe_1 7 PipeSec 1 CS 1 Nodes 13 14 15</p><p>Pipe_1 8 PipeSec 1 CS 1 Nodes 15 16 17</p><p>Pipe_1 9 PipeSec 1 CS 1 Nodes 17 18 19</p><p>Pipe_1 10 PipeSec 1 CS 1 Nodes 19 20 21</p><p>Pipe_1 11 PipeSec 1 CS 1 Nodes 21 22 23</p><p>Pipe_1 12 PipeSec 1 CS 1 Nodes 23 24 25</p><p>Pipe_1 13 PipeSec 1 CS 1 Nodes 25 26 27</p><p>Pipe_1 14 PipeSec 1 CS 1 Nodes 27 28 29</p><p>Pipe_1 15 PipeSec 1 CS 1 Nodes 29 30 31</p><p>…</p><p>Pipe_1 986 PipeSec 1 CS 1 Nodes 1971 1972 1973</p><p>Pipe_1 987 PipeSec 1 CS 1 Nodes 1973 1974 1975</p><p>Pipe_1 988 PipeSec 1 CS 1 Nodes 1975 1976 1977</p><p>Pipe_1 989 PipeSec 1 CS 1 Nodes 1977 1978 1979</p><p>Pipe_1 990 PipeSec 1 CS 1 Nodes 1979 1980 1981</p><p>Pipe_1 991 PipeSec 1 CS 1 Nodes 1981 1982 1983</p><p>Pipe_1 992 PipeSec 1 CS 1 Nodes 1983 1984 1985</p><p>Pipe_1 993 PipeSec 1 CS 1 Nodes 1985 1986 1987</p><p>Pipe_1 994 PipeSec 1 CS 1 Nodes 1987 1988 1989</p><p>Pipe_1 995 PipeSec 1 CS 1 Nodes 1989 1990 1991</p><p>Pipe_1 996 PipeSec 1 CS 1 Nodes 1991 1992 1993</p><p>Pipe_1 997 PipeSec 1 CS 1 Nodes 1993 1994 1995</p><p>Pipe_1 998 PipeSec 1 CS 1 Nodes 1995 1996 1997</p><p>Pipe_1 999 PipeSec 1 CS 1 Nodes 1997 1998 1999</p><p>Pipe_1 1000 PipeSec 1 CS 1 Nodes 1999 2000 2001</p></td></tr></tbody></table>

### Creation of element sets

One element set is defined, which contains all the elements created previously. It is used to create the internal pressure load field.

<table><tbody><tr class="odd"><td><p>//Creation of element sets</p><p>ElementSets 1</p><p>ElementSet 1 Elements 1000 Sequence Initial 1 Increment 1</p></td></tr></tbody></table>

### Creation of surfaces

Both the flat surfaces and the oscillatory imperfection are modeled using rigid oscillatory surfaces. So, 3 surfaces are created: the first represents the central imperfection and the other 2 represent the adjacent flat surfaces. Note that rigid triangular surfaces could also be used to represent the flat surfaces.

To create the oscillatory imperfection it is necessary to provide information on amplitude, wavelength, phase and number of waves for each direction of the surface. Such parameters must be adjusted to create the desired geometry of the imperfection, besides establishing a smooth connection with the adjacent flat surfaces. Using the coordinate system with identification number 2, the oscillatory imperfection illustrated in Figure 15‑1 (a) can be created by using one wave in each direction of the surface. Amplitude A1 must be equal to half the desired height for the imperfection and Lambda1 must be equal to the desired wavelength for the imperfection. Phi1 must be -*π*/2. Once the imperfection height along direction E2 does not vary, amplitudes A2 and A12 must be zero and Lambda2 must be equal to the desired width for the surface. Note that it is also necessary to specify a pilot node to rule the surface movement along the model evolution. In this tutorial, the pilot node must be located in the center of the surface, at half the height of the imperfection.

The flat surfaces can be created using the same procedure. For a surface to be flat, all amplitudes must be zero.

<table><tbody><tr class="odd"><td><p>//Creation of surfaces</p><p>Surfaces 3</p><p>RigidOscillatorySurface_1 1 A1 -0.125 A2 0 A12 0 Lambda1 12.5 Lambda2 1000 Phi1 -1.57 Phi2 1.57 Waves1 1 Waves2 1 CS 2 PilotNode 2002</p><p>RigidOscillatorySurface_1 2 A1 0 A2 0 A12 0 Lambda1 1493.75 Lambda2 1000 Phi1 -1.57 Phi2 1.57 Waves1 1 Waves2 1 CS 2 PilotNode 2003</p><p>RigidOscillatorySurface_1 3 A1 0 A2 0 A12 0 Lambda1 1493.75 Lambda2 1000 Phi1 -1.57 Phi2 1.57 Waves1 1 Waves2 1 CS 2 PilotNode 2004</p></td></tr></tbody></table>

### Creation of surface sets

Two surface sets are created: the first groups the flat surfaces and the second contains the oscillatory imperfection. The surfaces sets are used to create the NSSS contact between the pipeline and the surfaces.

<table><tbody><tr class="odd"><td><p>//Creation of surface sets</p><p>SurfaceSets 2</p><p>SurfaceSet 1 Surfaces 1 List 1</p><p>SurfaceSet 2 Surfaces 2 Sequence Initial 2 Increment 1</p></td></tr></tbody></table>

### Creation of contacts

Two NSSS contacts are created to test the interaction between the pipeline and the rigid surfaces along the simulation. The NSSS contact with identification number 1 tests the interaction between the pipeline nodes (node set 1) and the oscillatory imperfection (surface set 1) while the NSSS contact with identification number 2 tests the interaction between the pipeline nodes (node set 1) and the flat surfaces (surface set 2). The NSSS contact model uses a node to surface description in which constraint enforcements are done by penalty method. It is necessary to input the coefficient of friction, the normal and tangential penalty coefficients, besides the normal and tangential damping. It is assumed here that the surfaces are rigid. So, large values of penalty coefficients are used (within the possibilities of numerical convergence) in order to minimize the vertical penetrations of the pipeline into the surfaces and to approximate the tangential constitutive relation to the Coulomb’s law. The coefficient of friction is 0.4 and no damping is considered.

Once the pipeline axis is at the level of the surfaces, the radius is zero. Neither the rigid body rotation that occurs in rolling scenarios nor the moment of the friction force are evaluated. A pinball radius of 2500 m ensures that whenever it exists, the contact is detected by the model. Besides this, only a single pointwise contact is allowed between the pipeline and the surfaces.

The BoolTable is not included so the contact constraint is considered turned on for all solution steps.

<table><tbody><tr class="odd"><td><p>//Creation of contacts</p><p>Contacts 2</p><p>NSSS 1 NodeSet 1 SurfaceSet 1 MU 0.4 EPN 1e9 CN 0 EPT 1e8 CT 0 Pinball 2500 Radius 0 MaxPointwiseInt 1</p><p>NSSS 2 NodeSet 1 SurfaceSet 2 MU 0.4 EPN 1e7 CN 0 EPT 1e6 CT 0 Pinball 2500 Radius 0 MaxPointwiseInt 1</p></td></tr></tbody></table>

### Creation of cross-sections

One only cross section is created, which is shown in the Figure 15‑1 (b).

<table><tbody><tr class="odd"><td><p>//Creation of pipeline sections</p><p>PipeSections 1</p><p>PS 1 EA 5.9847E9 EI 3.0182E8 GJ 2.4146E8 GA 2.3939E9 Rho 234.9 CDt 0.1 CDn 1.2 CAt 0 CAn 1 De 0.65 Di 0.62</p></td></tr></tbody></table>

### Creation of coordinate systems

Three coordinate systems are created. The first one is used to define the element properties. Note again that, since employing Pipe\_1 element, direction E3 must lie in the pipeline alignment direction. The second coordinate system is used to define the rigid oscillatory surfaces. Directions E1 and E2 must be the horizontal directions of the oscillatory surfaces while direction E3 must be related to their amplitudes. Finally, the third coordinate system is used to define the direction of application of the nodal loads and prescription of the nodal displacements.

<table><tbody><tr class="odd"><td><p>//Creation of coordinate systems</p><p>CoordinateSystems 3</p><p>CS 1 E1 0 1 0 E3 1 0 0</p><p>CS 2 E1 1 0 0 E3 0 1 0</p><p>CS 3 E1 1 0 0 E3 0 0 1</p></td></tr></tbody></table>

### Creation of solution steps

This tutorial has eight solution steps. The first seven solution steps are static and are used to simulate the laying process of the pipeline on the rigid surfaces and to impose a lateral imperfection at the central node of the pipeline. The laying process is done in more than one step to avoid rigid body movements that can occur due to the low initial stiffness of the structure. As already mentioned, the pipeline is initially positioned at the level of the top of the oscillatory imperfection and has all nodes fixed. The first solution step consists in prescribing a vertical displacement of 0.25 m (imperfection height) at node sets 9 and 10 in order to position the pipeline at the level of the flat rigid surfaces. Note that the pipeline self-weight is not applied yet. So, no contact is established between the pipeline and the surfaces. Such step establishes the end time equal to 1. The time step used is 0.04, that is, the prescribed displacement is divided into increments of 4%. In the second solution step, the self-weight is applied to the structure and the node set 8 is released in order to establish contact with the rigid surfaces. The solution step 2 establishes the end time equal to 2. The time step used is 0.004, that is, the load is divided into increments of 0.4%. The third solution step releases the node sets 9 and 10 under the action of the pipeline self-weight in order to establish contact with the rigid surfaces. The solution step 3 establishes the end time equal to 3. Again, the time step used is 0.004. The fourth solution step releases the node sets 2 and 4 in order to establish contact with the rigid surfaces. Such node sets are fixed again in the fifth solution step. The solution steps 4 and 5 establish the end times equal to 4 and 5, respectively. The time step used is 1, that is, the load is applied in a single step.

Once finished the laying process of the pipeline, a lateral displacement of 0.501 m is imposed at the central node in the sixth solution step. Such node is released in the seventh solution step and, due to the friction between the pipeline and the surface, the lateral imperfection of the pipeline is maintained. The solution steps 6 and 7 establish the end times equal to 6 and 7, respectively. The time step used is 1, that is, the load is applied in a single step.

The eighth solution step is dynamic and is used to apply an internal pressure of 20 MPa to the pipeline. The internal pressure is applied over 11240 s. Such time is approximately eleven times the largest natural period of the structure (a modal analysis can be performed previously to determine the natural period of the structure) and avoids the excitation of any natural vibration mode of the structure. The eighth step establishes the end time equal to 11247. The time step used is 1, that is, the internal pressure is initially divided into increments of approximately 0.0089%.

<table><tbody><tr class="odd"><td><p>//Creation of solution steps</p><p>SolutionSteps 8</p><p>Static 1</p><p>EndTime 1</p><p>TimeStep 0.04</p><p>MaxTimeStep 1</p><p>MinTimeStep 1.0E-10</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.1</p><p>Sample 10</p><p>Static 2</p><p>EndTime 2</p><p>TimeStep 0.004</p><p>MaxTimeStep 1</p><p>MinTimeStep 1.0E-10</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.1</p><p>Sample 10</p><p>Static 3</p><p>EndTime 3</p><p>TimeStep 0.004</p><p>MaxTimeStep 1</p><p>MinTimeStep 1.0E-10</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.1</p><p>Sample 10</p><p>Static 4</p><p>EndTime 4</p><p>TimeStep 1</p><p>MaxTimeStep 1</p><p>MinTimeStep 1.0E-10</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.1</p><p>Sample 10</p><p>Static 5</p><p>EndTime 5</p><p>TimeStep 1</p><p>MaxTimeStep 1</p><p>MinTimeStep 1.0E-10</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.1</p><p>Sample 10</p><p>Static 6</p><p>EndTime 6</p><p>TimeStep 1</p><p>MaxTimeStep 1</p><p>MinTimeStep 1.0E-10</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.1</p><p>Sample 10</p><p>Dynamic 7</p><p>EndTime 7</p><p>TimeStep 1</p><p>MaxTimeStep 1</p><p>MinTimeStep 1.0E-10</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.1</p><p>Sample 1</p><p>RayleighDamping Alpha 0 Beta 0 Update 0</p><p>NewmarkCoefficients Beta 0.3 Gamma 0.5</p><p>Dynamic 8</p><p>EndTime 11247</p><p>TimeStep 1</p><p>MaxTimeStep 50</p><p>MinTimeStep 1e-8</p><p>MaxIt 15</p><p>MinIt 3</p><p>ConvIncrease 10</p><p>IncFactor 1.1</p><p>Sample 1</p><p>RayleighDamping Alpha 0 Beta 0 Update 0</p><p>NewmarkCoefficients Beta 0.3 Gamma 0.5</p></td></tr></tbody></table>

### Creation of loads

An internal pressure of 20 MPa is applied to the element set 1 over 11240 s. Such load is applied as a ramp-load in the solution step 8.

<table><tbody><tr class="odd"><td><p>//Creation of loads</p><p>Loads 1</p><p>PipeLoad 1 ElementSet 1 NTimes 2</p><p>//Time P0I P0E RhoI RhoE</p><p>7 0 0 0 0</p><p>11247 20000000 0 0 0</p></td></tr></tbody></table>

### Creation of displacements

Three nodal displacements are prescribed along the simulation. The nodal displacements 1 and 2 consist in vertical displacements of 0.25 m prescribed at node sets 9 and 10, respectively. The idea is to place the pipeline closer to the flat surfaces in order to avoid rigid body movements. Such displacements are prescribed as ramp-displacements in the first solution step. The nodal displacement 3 consists in the creation of the lateral imperfection of 0.501 m at node set 3 and is prescribed as a ramp-displacement in the sixth solution step.

<table><tbody><tr class="odd"><td><p>//Creation of displacements</p><p>Displacements 3</p><p>NodalDisplacement 1 NodeSet 9 CS 3 NTimes 2</p><p>//Time UX UY UZ ROTX ROTY ROTZ</p><p>0 0 0 0.000 0 0 0</p><p>1 0 -0.250 0 0 0 0</p><p>NodalDisplacement 2 NodeSet 10 CS 3 NTimes 2</p><p>//Time UX UY UZ ROTX ROTY ROTZ</p><p>0 0 0 0.000 0 0 0</p><p>1 0 -0.250 0 0 0 0</p><p>NodalDisplacement 3 NodeSet 3 CS 3 NTimes 2</p><p>//Time UX UY UZ ROTX ROTY ROTZ</p><p>5 0 0 0.000 0 0 0</p><p>6 0 0 0.501 0 0 0</p></td></tr></tbody></table>

### Creation of constraints

Nine nodal constraints are created. The first three constraints refer to the pilot nodes and fix their DOFs UX, UY, UZ, ROTX, ROTY and ROTZ throughout the simulation. The fourth constraint refers to the node set 8 and fixes its DOFs UX, UY, UZ, ROTX, ROTY and ROTZ in the first solution step. Such DOFs are released in the second solution step, which allows to the central region of the pipeline to move under the action of its self-weight. The fifth and sixth constraints refer to the node sets 9 and 10, respectively, and fix their DOFs UX, UY, UZ, ROTX, ROTY and ROTZ in the first solution step (to prescribe the vertical displacement of 0.25 m) as well as in the second solution step (while the node set 8 is released). Such DOFs are released in the third solution step, which allows to the end regions of the pipeline to move under the action of their self-weights. The seventh and eight constraints refer to the node sets 2 and 4, respectively. Such node sets are maintained fixed throughout the simulation, except for the solution step 4, in which the end nodes are released to establish contact with the surfaces. Finally, the last constraint refers to the node set 3 whose DOFs UX, UY, UZ, ROTX, ROTY and ROTZ are fixed in the solution steps 1 and 6. In the solution step 6, a lateral displacement of 0.501 m is prescribed at node set 3, requiring its constraints to be turned on.

<table><tbody><tr class="odd"><td><p>//Creation of constraint</p><p>Constraints 9</p><p>NodalConstraint 1 NodeSet 5</p><p>UX BoolTable 1</p><p>UY BoolTable 1</p><p>UZ BoolTable 1</p><p>ROTX BoolTable 1</p><p>ROTY BoolTable 1</p><p>ROTZ BoolTable 1</p><p>NodalConstraint 2 NodeSet 6</p><p>UX BoolTable 1</p><p>UY BoolTable 1</p><p>UZ BoolTable 1</p><p>ROTX BoolTable 1</p><p>ROTY BoolTable 1</p><p>ROTZ BoolTable 1</p><p>NodalConstraint 3 NodeSet 7</p><p>UX BoolTable 1</p><p>UY BoolTable 1</p><p>UZ BoolTable 1</p><p>ROTX BoolTable 1</p><p>ROTY BoolTable 1</p><p>ROTZ BoolTable 1</p><p>NodalConstraint 4 NodeSet 8</p><p>UX BoolTable 1 0 0 0 0 0 0 0</p><p>UY BoolTable 1 0 0 0 0 0 0 0</p><p>UZ BoolTable 1 0 0 0 0 0 0 0</p><p>ROTX BoolTable 1 0 0 0 0 0 0 0</p><p>ROTY BoolTable 1 0 0 0 0 0 0 0</p><p>ROTZ BoolTable 1 0 0 0 0 0 0 0</p><p>NodalConstraint 5 NodeSet 9</p><p>UX BoolTable 1 1 0 0 0 0 0 0</p><p>UY BoolTable 1 1 0 0 0 0 0 0</p><p>UZ BoolTable 1 1 0 0 0 0 0 0</p><p>ROTX BoolTable 1 1 0 0 0 0 0 0</p><p>ROTY BoolTable 1 1 0 0 0 0 0 0</p><p>ROTZ BoolTable 1 1 0 0 0 0 0 0</p><p>NodalConstraint 6 NodeSet 10</p><p>UX BoolTable 1 1 0 0 0 0 0 0</p><p>UY BoolTable 1 1 0 0 0 0 0 0</p><p>UZ BoolTable 1 1 0 0 0 0 0 0</p><p>ROTX BoolTable 1 1 0 0 0 0 0 0</p><p>ROTY BoolTable 1 1 0 0 0 0 0 0</p><p>ROTZ BoolTable 1 1 0 0 0 0 0 0</p><p>NodalConstraint 7 NodeSet 2</p><p>UX BoolTable 1 1 1 0 1 1 1 1</p><p>UY BoolTable 1 1 1 0 1 1 1 1</p><p>UZ BoolTable 1 1 1 0 1 1 1 1</p><p>ROTX BoolTable 1 1 1 0 1 1 1 1</p><p>ROTY BoolTable 1 1 1 0 1 1 1 1</p><p>ROTZ BoolTable 1 1 1 0 1 1 1 1</p><p>NodalConstraint 8 NodeSet 4</p><p>UX BoolTable 1 1 1 0 1 1 1 1</p><p>UY BoolTable 1 1 1 0 1 1 1 1</p><p>UZ BoolTable 1 1 1 0 1 1 1 1</p><p>ROTX BoolTable 1 1 1 0 1 1 1 1</p><p>ROTY BoolTable 1 1 1 0 1 1 1 1</p><p>ROTZ BoolTable 1 1 1 0 1 1 1 1</p><p>NodalConstraint 9 NodeSet 3</p><p>UX BoolTable 1 0 0 0 0 1 0 0</p><p>UY BoolTable 1 0 0 0 0 1 0 0</p><p>UZ BoolTable 1 0 0 0 0 1 0 0</p><p>ROTX BoolTable 1 0 0 0 0 1 0 0</p><p>ROTY BoolTable 1 0 0 0 0 1 0 0</p><p>ROTZ BoolTable 1 0 0 0 0 1 0 0</p></td></tr></tbody></table>

### Creation of environment data

Environment data are created to introduce the gravity field on the model. In this example, the gravity field is linearly imposed along solution step 2 and is maintained throughout the simulation.

<table><tbody><tr class="odd"><td><p>//Creation of environment data</p><p>Environment</p><p>GravityData</p><p>G 0 -9.81 0 BoolTable 0 1 1 1 1 1 1 1</p></td></tr></tbody></table>

### Creation of solver options

It is requested Giraffe to apply three threads to solve this case. The linear system solution chosen is Direct.

<table><tbody><tr class="odd"><td><p>//Creation of solver options</p><p>SolverOptions</p><p>Processors 3 LinSys Direct</p></td></tr></tbody></table>

### Creation of monitors

It is requested Giraffe to monitor 3 nodes along the simulation: node 1 (point A), node 1001 (midspan) and node 2001 (point B). It is also requested Giraffe to monitor the two NSSS contacts created.

<table><tbody><tr class="odd"><td><p>//Creation of monitors</p><p>Monitor Sample 1</p><p>MonitorNodes 1 1001 2001</p><p>MonitorContacts 1 2</p></td></tr></tbody></table>

### Creation of post files

The options WriteMesh, WriteRigidContactSurfaces and WriteContactForces are chosen for post-processing using Paraview<sup>TM</sup>. The parameter MagFactor assumes a unit value so the deformed shape of Paraview<sup>TM</sup> frames reflects the real scale.

<table><tbody><tr class="odd"><td><p>//Creation of post files</p><p>PostFiles</p><p>MagFactor 1.0</p><p>WriteMesh 1</p><p>WriteRenderMesh 0</p><p>WriteRigidContactSurfaces 1</p><p>WriteFlexibleContactSurfaces 0</p><p>WriteForces 0</p><p>WriteConstraints 0</p><p>WriteSpecialConstraints 0</p><p>WriteContactForces 1</p><p>WriteRenderRigidBodies 0</p><p>WriteRenderParticles 0</p></td></tr></tbody></table>

### Creation of convergence criteria

Here, it is defined the convergence criteria. This is an optional task, which may be done only by advanced users. For this case, however, it is necessary to re-stablish such criteria according to here provided data in order to ensure convergence.

<table><tbody><tr class="odd"><td><p>//Creation of convergence criteria</p><p>ConvergenceCriteria</p><p>ForceTolerance 1e-5</p><p>MomentTolerance 1e-5</p><p>ForceMinimumReference 1e0</p><p>MomentMinimumReference 1e0</p><p>ConstraintMinimumReference 1e-7</p><p>DisplacementTolerance 1e-4</p><p>RotationTolerance 1e-4</p><p>LagrangeTolerance 1e-4</p><p>DisplacementMinimumReference 1e-4</p><p>RotationMinimumReference 1e-4</p><p>LagrangeMinimumReference 1e-6</p><p>DivergenceReference 1e+15</p></td></tr></tbody></table>

## Results and post-processing

### Paraview<sup>TM</sup> post-processing

To post-process results using Paraview<sup>TM</sup>, first we must open the files “whole\_solution\_contactsurfaces.pvd” and “whole\_solution\_mesh.pvd”, which are located at the folder “post” in “tutorial16”.

<img src="../../assets/tutorials/image121.png" style="width:5.27559in;height:2.2874in" alt="" />

Figure 15‑2 – Opening the files in Paraview<sup>TM</sup>.

Then, we must click the preview icon and the Apply button on the Pipeline browser.

<img src="../../assets/tutorials/image122.png" style="width:1.92913in;height:2.9252in" alt="" />

Figure 15‑3 – Pipeline Browser.

Figure 15‑4 shows how Paraview<sup>TM</sup> displays the results.

<img src="../../assets/tutorials/image123.png" style="width:4.69167in;height:2.51686in" alt="" />

Figure 15‑4 – Results using Paraview<sup>TM</sup>.

It is interesting to perform an animation to observe how the structure deforms during the laying process and during the application of the internal pressure. Note that the model is three-dimensional and a better visualization of the deformation of the structure requires multiple views. This can be done by dividing the screen into multiple views. First, you can divide the screen horizontally into two parts. For that, use the menu located in the upper-right corner of the current view (Figure 15‑5). The left-hand view can be used, for example, to visualize the perspective of the model (RenderView1). The right-hand view, in its turn, can also be divided vertically into two parts. For that, use the menu located in the upper-right corner of the right-hand view (Figure 15‑6). Choose, then, two new render views. The top view, for example, can be used to visualize the xy plane (RenderView2) and the bottom view can be used to visualize the xz plane (RenderView3).

<img src="../../assets/tutorials/image124.png" style="width:2.27778in;height:1.20794in" alt="" />

<img src="../../assets/tutorials/image125.png" style="width:4.06667in;height:2.18158in" alt="" />

Figure 15‑5 – Dividing the screen horizontally into two views.

<img src="../../assets/tutorials/image126.png" style="width:2.76667in;height:1.26667in" alt="" />

<img src="../../assets/tutorials/image127.png" style="width:4.06578in;height:2.1811in" alt="" />

Figure 15‑6 – Dividing the screen vertically into two views.

It is possible to adjust the direction and the zoom of each of the views for a better visualization. Texts can also be inserted by clicking, in sequence, Sources, Alphabetical and Text. Figure 15‑7 shows the final layout. The animation can be performed by clicking the Play button located in the VCR controls toolbar.

<img src="../../assets/tutorials/image128.png" style="width:5.9in;height:3.15833in" alt="" />

Figure 15‑7 – Animation layout.

### Monitors post-processing

Monitors are useful, for example, to plot the time series of node coordinates and displacements. It allows to evaluate the buckling process.

<img src="../../assets/tutorials/image129.png" style="width:5.81111in;height:3.89167in" alt="" />

Figure 15‑8 – Coordinate z of the central node of the pipeline versus time.
