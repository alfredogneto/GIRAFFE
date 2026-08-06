# Tutorial 12: Simplified bogie dynamics

!!! example "Input files"
    [`Bolster_tutorial.stl`](../assets/tutorials/inputs/tutorial12/Bolster_tutorial.stl) · [`Eixo_tutorial.stl`](../assets/tutorials/inputs/tutorial12/Eixo_tutorial.stl) · [`Sideleft_tutorial.stl`](../assets/tutorials/inputs/tutorial12/Sideleft_tutorial.stl) · [`Sideright_tutorial.stl`](../assets/tutorials/inputs/tutorial12/Sideright_tutorial.stl) · [`tutorial12.inp`](../assets/tutorials/inputs/tutorial12/tutorial12.inp)
    (also in the repository under `inputs/tutorial12/`)

Input file name: tutorial12.inp

*Developed by Gabriel Maki*

## Objective

This tutorial presents a dynamic simulation of a bogie. The aim is to show how to model rigid bodies, suspensions and identify the contact between two surfaces. Besides, the tutorial has the objective of study the bogie’s movement and the hunting phenomenon in railway systems.

## Description

The bogie is meshed considering 34 nodes and 24 elements. It is made up of five rigid bodies, the bolster, the left-side and right-side frames and the frontal and back wheelsets. All the bodies are connected by suspensions. A lumped mass is included in the model to represent half of the wagon mass. It is rigidly linked to the bolster. For the simulation, a lateral nodal load is applied on the bolster barycenter to induce the hunting oscillations. The simulation will be done in 4 solution steps.

## Input file step-by-step

### Creation of nodes

A total of 34 nodes is created to establish a finite element mesh. The nodes 1 to 13 and 28 to 34 were created to represent the rigid bodies, while the nodes 14 to 27 were created to establish connections between bodies.

<table><tbody><tr class="odd"><td><p>Nodes 34</p><p>//Bolster</p><p>Node 1 -1 0 0.5</p><p>//Sideframe - left</p><p>Node 2 -1.00219 1.2 0.47834</p><p>//Sideframe - right</p><p>Node 3 -1.00219 -1.2 0.47834</p><p>//Wheelset - front</p><p>Node 4 0.01347 0 0.4123</p><p>//Wheelset - back</p><p>Node 5 -1.91533 0 0.4123</p><p>//Nodes - Sideframe - left</p><p>Node 6 -0.08653 1.2 0.8925</p><p>Node 7 -1.81533 1.2 0.8925</p><p>//Nodes - Sideframe - right</p><p>Node 8 -0.08653 -1.2 0.8925</p><p>Node 9 -1.81533 -1.2 0.8925</p><p>//Nodes - Wheelset - front</p><p>Node 10 0.01347 1.0095 0.4123</p><p>Node 11 0.01347 -1.0095 0.4123</p><p>//Nodes - Wheelset - back</p><p>Node 12 -1.91533 1.0095 0.4123</p><p>Node 13 -1.91533 -1.0095 0.4123</p><p>//Nodes - Secondary suspension longitudinal direction</p><p>Node 14 -1.00093 1.0095 0.5</p><p>Node 15 -1.00093 -1.0095 0.5</p><p>//Nodes - primary suspension</p><p>Node 16 0.01347 1.2 0.8925</p><p>Node 17 0.01347 -1.2 0.8925</p><p>Node 18 -1.91533 1.2 0.8925</p><p>Node 19 -1.91533 -1.2 0.8925</p><p>Node 20 0.01347 1.0095 0.8925</p><p>Node 21 -1.91533 1.0095 0.8925</p><p>Node 22 0.01347 -1.0095 0.8925</p><p>Node 23 -1.91533 -1.0095 0.8925</p><p>//Nodes - secondary suspension</p><p>Node 24 -1.00093 1.2 0.5</p><p>Node 25 -1.00093 -1.2 0.5</p><p>Node 26 -1.00093 1.2 0.1952</p><p>Node 27 -1.00093 -1.2 0.1952</p><p>//Rails</p><p>Node 28 -2.5 0.9662 -0.114205</p><p>Node 29 2500 0.9662 -0.114205</p><p>Node 30 -2.5 -0.9662 -0.114205</p><p>Node 31 2500 -0.9662 -0.114205</p><p>//Wagon barycenter</p><p>Node 32 -1.00093 0 2.1154</p><p>//Nodes - Bolster</p><p>Node 33 -1 1.0095 0.5</p><p>Node 34 -1 -1.0095 0.5</p></td></tr></tbody></table>

### Creation of node sets

13 node sets are created to establish the constraints, loads and special constraints.

<table><tbody><tr class="odd"><td><p>NodeSets 13</p><p>//Wheelset - front</p><p>NodeSet 1 Nodes 2 List 10 11</p><p>//Wheelset - back</p><p>NodeSet 2 Nodes 2 List 12 13</p><p>//Sideframe - left</p><p>NodeSet 3 Nodes 3 List 6 7 26</p><p>//Sideframe - right</p><p>NodeSet 4 Nodes 3 List 8 9 27</p><p>//Bolster</p><p>NodeSet 5 Nodes 2 List 33 34</p><p>//BC - wheelset front</p><p>NodeSet 6 Nodes 1 List 4</p><p>//BC - wheelset back</p><p>NodeSet 7 Nodes 1 List 5</p><p>//Load - bolster</p><p>NodeSet 8 Nodes 1 List 1</p><p>//Rails</p><p>NodeSet 9 Nodes 1 List 28</p><p>NodeSet 10 Nodes 1 List 29</p><p>NodeSet 11 Nodes 1 List 30</p><p>NodeSet 12 Nodes 1 List 31</p><p>NodeSet 13 Nodes 1 List 32</p></td></tr></tbody></table>

### Creation of special constraints

26 special constraints are defined to represent the rigid node sets and the suspensions. To define suspensions, it is necessary to create a translational joint between nodes that link the bodies. This joint allows that two nodes to have relative displacements in a single direction, defined by the direction E3 of a local coordinate system. Besides, we must define a rotational connection between the bodies, so we use the Same Rotation constraint.

<table><tbody><tr class="odd"><td><p>SpecialConstraints 26</p><p>//Wheelset - front</p><p>RigidNodeSet 1 PilotNode 4 NodeSet 1</p><p>//Wheelset - back</p><p>RigidNodeSet 2 PilotNode 5 NodeSet 2</p><p>//Sideframe - left</p><p>RigidNodeSet 3 PilotNode 2 NodeSet 3</p><p>//Sideframe - right</p><p>RigidNodeSet 4 PilotNode 3 NodeSet 4</p><p>//Bolster</p><p>RigidNodeSet 5 PilotNode 1 NodeSet 5</p><p>//Primary Suspension - vertical</p><p>TranslationalJoint 6 Nodes 10 20 RotationNode 2 CS 3</p><p>TranslationalJoint 7 Nodes 11 22 RotationNode 3 CS 3</p><p>TranslationalJoint 8 Nodes 12 21 RotationNode 2 CS 3</p><p>TranslationalJoint 9 Nodes 13 23 RotationNode 3 CS 3</p><p>//Primary Suspension - lateral</p><p>TranslationalJoint 10 Nodes 20 16 RotationNode 2 CS 4</p><p>TranslationalJoint 11 Nodes 22 17 RotationNode 3 CS 4</p><p>TranslationalJoint 12 Nodes 21 18 RotationNode 2 CS 4</p><p>TranslationalJoint 13 Nodes 23 19 RotationNode 3 CS 4</p><p>//Primary Suspension - longitudinal</p><p>TranslationalJoint 14 Nodes 6 16 RotationNode 2 CS 5</p><p>TranslationalJoint 15 Nodes 8 17 RotationNode 3 CS 5</p><p>TranslationalJoint 16 Nodes 7 18 RotationNode 2 CS 5</p><p>TranslationalJoint 17 Nodes 9 19 RotationNode 3 CS 5</p><p>//Secondary Suspension - vertical</p><p>TranslationalJoint 18 Nodes 26 24 RotationNode 2 CS 3</p><p>TranslationalJoint 19 Nodes 27 25 RotationNode 3 CS 3</p><p>//Secondary Suspension - lateral</p><p>TranslationalJoint 20 Nodes 14 24 RotationNode 2 CS 4</p><p>TranslationalJoint 21 Nodes 15 25 RotationNode 3 CS 4</p><p>//Hinge Joints</p><p>SameRotation 22 Nodes 1 2</p><p>SameRotation 23 Nodes 1 3</p><p>//Wagon barycenter - rigid node set fixed to bolster</p><p>RigidNodeSet 24 PilotNode 1 NodeSet 13</p><p>//Secondary Suspension - longitudinal</p><p>TranslationalJoint 25 Nodes 34 15 RotationNode 3 CS 5</p><p>TranslationalJoint 26 Nodes 33 14 RotationNode 2 CS 5</p></td></tr></tbody></table>

### Creation of elements

24 elements are created to compose the mesh. Here we use three different types of element. The first type is the RigidBody\_1, an element that create a rigid body based on a single node, normally chosen as the barycenter of the body (but not necessarily). The second type is the SpringDashpot\_1 element, used to create a spring and damper assembly in parallel. The third type is the Mass\_1, which is used to create a lumped mass.

<table><tbody><tr class="odd"><td><p>Elements 24</p><p>RigidBody_1 1 RigidBodyData 1 CS 1 Node 1</p><p>RigidBody_1 2 RigidBodyData 2 CS 1 Node 2</p><p>RigidBody_1 3 RigidBodyData 3 CS 1 Node 3</p><p>RigidBody_1 4 RigidBodyData 4 CS 1 Node 4</p><p>RigidBody_1 5 RigidBodyData 4 CS 1 Node 5</p><p>//Pad - vertical</p><p>SpringDashpot_1 6 Stiffness 7000000 Damping 200000 Nodes 10 20</p><p>SpringDashpot_1 7 Stiffness 7000000 Damping 200000 Nodes 11 22</p><p>SpringDashpot_1 8 Stiffness 7000000 Damping 200000 Nodes 12 21</p><p>SpringDashpot_1 9 Stiffness 7000000 Damping 200000 Nodes 13 23</p><p>//Pad - lateral</p><p>SpringDashpot_1 10 Stiffness 150000000 Damping 20000 Nodes 20 16</p><p>SpringDashpot_1 11 Stiffness 150000000 Damping 20000 Nodes 22 17</p><p>SpringDashpot_1 12 Stiffness 150000000 Damping 20000 Nodes 21 18</p><p>SpringDashpot_1 13 Stiffness 150000000 Damping 20000 Nodes 23 19</p><p>//Pad - longitudinal</p><p>SpringDashpot_1 14 Stiffness 140000000 Damping 20000 Nodes 16 6</p><p>SpringDashpot_1 15 Stiffness 140000000 Damping 20000 Nodes 17 8</p><p>SpringDashpot_1 16 Stiffness 140000000 Damping 20000 Nodes 18 7</p><p>SpringDashpot_1 17 Stiffness 140000000 Damping 20000 Nodes 19 9</p><p>//Secondary - vertical</p><p>SpringDashpot_1 18 Stiffness 1700000 Damping 0 Nodes 26 24</p><p>SpringDashpot_1 19 Stiffness 1700000 Damping 0 Nodes 27 25</p><p>//Secondary - lateral</p><p>SpringDashpot_1 20 Stiffness 130000000 Damping 0 Nodes 14 24</p><p>SpringDashpot_1 21 Stiffness 130000000 Damping 0 Nodes 15 25</p><p>Mass_1 22 Mass 55000 Node 32</p><p>//Secondary - longitudinal</p><p>SpringDashpot_1 23 Stiffness 170000000 Damping 200000 Nodes 34 15</p><p>SpringDashpot_1 24 Stiffness 170000000 Damping 200000 Nodes 33 14</p></td></tr></tbody></table>

### Creation of rigid body data

Here we insert the mass and inertia properties of each kind of rigid body. Besides, we can request Giraffe to use a CAD file in the simulation to represent the rigid body (for post-processing purposes). For that, one may use any CAD which is able to export in .stl ASCII file to develop the geometry.

<table><tbody><tr class="odd"><td><p>RBData 1</p><p>Mass 800</p><p>J11 930 J22 22 J33 480 J12 0.0 J13 0.0 J23 0.0</p><p>Barycenter 0.0 0.0 0.0</p><p>CADFile Bolster_tutorial.stl</p><p>RBData 2</p><p>Mass 700</p><p>J11 60 J22 340 J33 300 J12 0.0 J13 0.0 J23 0.0</p><p>Barycenter 0.0 0.0 0.0</p><p>CADFile Sideleft_tutorial.stl</p><p>RBData 3</p><p>Mass 700</p><p>J11 60 J22 340 J33 300 J12 0.0 J13 0.0 J23 0.0</p><p>Barycenter 0.0 0.0 0.0</p><p>CADFile Sideright_tutorial.stl</p><p>RBData 4</p><p>Mass 1550</p><p>J11 1050 J22 120 J33 1050 J12 0.0 J13 0.0 J23 0.0</p><p>Barycenter 0.0 0.0 0.0</p><p>CADFile Eixo_tutorial.stl</p></td></tr></tbody></table>

### Creation of coordinate systems

Nine coordinate systems are defined considering all the constraints and surfaces that will be created after. Note that CS 6 and CS 7 have some different directions because they are used to define the rail cant (inclination). They are made so that the direction E2 stays parallel to rail cross sections. In this case, the cant used is 1:20.

<table><tbody><tr class="odd"><td><p>CoordinateSystems 9</p><p>//Global CS</p><p>CS 1 E1 1 0 0 E3 0 0 1</p><p>//Wheelset CS - E1 aligned with wheelset axle</p><p>CS 2 E1 0 1 0 E3 0 0 1</p><p>//Supension - vertical</p><p>CS 3 E1 1 0 0 E3 0 0 1</p><p>//Suspension - lateral</p><p>CS 4 E1 1 0 0 E3 0 1 0</p><p>//Suspension - longitudinal</p><p>CS 5 E1 0 1 0 E3 1 0 0</p><p>//Rails</p><p>CS 6 E1 0 0.998752339179 4.99376109278e-2 E3 1 0 0</p><p>CS 7 E1 0 -0.998752339179 4.99376109278e-2 E3 -1 0 0</p><p>//Wheels</p><p>CS 8 E1 0 0 -1 E3 1 0 0</p><p>CS 9 E1 0 0 -1 E3 -1 0 0</p></td></tr></tbody></table>

### Creation of arcs

Here we insert the circumference arcs for the creation of contact surfaces. Note that input data include: initial, final and center points of the arc in a plan. The coordinate system and the node established as the system origin are defined on the creation of each surface.

<table><tbody><tr class="odd"><td><p>//Creation of arcs</p><p>Arcs 14</p><p>//Wheel arcs</p><p>Arc 1 InitialPoint 0.425796370582452 0.882054430691546 FinalPoint 0.443208542157583 0.889979653581409 CenterPoint 0.421579962014195 0.914408507213130</p><p>Arc 2 InitialPoint 0.443208542157585 0.889979653581408 FinalPoint 0.446742276607354 0.897460714881545 CenterPoint 0.436246087964145 0.897843543549290</p><p>Arc 3 InitialPoint 0.446742276607329 0.897460714881546 FinalPoint 0.442401755983930 0.908021749519942 CenterPoint 0.432948005503375 0.897963834881460</p><p>Arc 4 InitialPoint 0.442401755983959 0.908021749519973 FinalPoint 0.434723623779321 0.912588822523442 CenterPoint 0.428246647673200 0.892962025722420</p><p>Arc 5 InitialPoint 0.426919157227811 0.918935831513070 FinalPoint 0.434723623779279 0.912588822523315 CenterPoint 0.439150899172595 0.926004532741030</p><p>Arc 6 InitialPoint 0.421791943301005 0.935976728599755 FinalPoint 0.426919157227897 0.918935831513120 CenterPoint 0.460460922049591 0.938319555863530</p><p>Arc 7 InitialPoint 0.421791943301296 0.935976728599773 FinalPoint 0.415483689711215 1.007853111080860 CenterPoint -0.910691689680170 0.855245892541840</p><p>Arc 8 InitialPoint 0.415483689711318 1.007853111080870 FinalPoint 0.413972014807964 1.019925467434800 CenterPoint -0.188853650451404 0.938310091375720</p><p>Arc 9 InitialPoint 0.413972014807492 1.019925467434730 FinalPoint 0.400805033000113 1.027311927799740 CenterPoint 0.404605274257736 1.018657322933800</p><p>//Rails arcs</p><p>Arc 10 InitialPoint 0.03685759800334810 0.08526439800008370</p><p>FinalPoint 0.02779624699522180 0.09820644298783040</p><p>CenterPoint 0.022574560946790100 0.084907320972669800</p><p>Arc 11 InitialPoint 0.02779624699875580 0.09820644299683120</p><p>FinalPoint 0.01777999999465240 0.10036312589318200</p><p>CenterPoint 0.016192499777804300 0.068652837381108900</p><p>Arc 12 InitialPoint 0.01777999999754500 0.10036312595096100</p><p>FinalPoint -0.01778000000000110 0.10036312600002300</p><p>CenterPoint -0.000000000491233737 -0.254792066634406000</p><p>Arc 13 InitialPoint -0.01778000000076080 0.10036312601519700</p><p>FinalPoint -0.02779624703571790 0.09820644309097000</p><p>CenterPoint -0.016192499912296400 0.068652838315372000</p><p>Arc 14 InitialPoint -0.02779624701583830 0.09820644304033850</p><p>FinalPoint -0.03685759803196000 0.08526439802029900</p><p>CenterPoint -0.022574560818925100 0.084907320913286000</p></td></tr></tbody></table>

### Creation of surfaces

Here we request Giraffe to create surfaces based on the arcs defined previously. There are resources in Giraffe that uses arc definitions to establish extruded/revolved surfaces. For the wheel, we must use the RigidArcRevolution\_1 (revolution surface). For the rail, we must use the FlexibleArcExtrusion\_1 (extrusion surface).

<table><tbody><tr class="odd"><td><p>Surfaces 46</p><p>//Left front wheel</p><p>RigidArcRevolution_1 1 Arc 1 CS 8 Node 4 Convex //interior</p><p>RigidArcRevolution_1 2 Arc 2 CS 8 Node 4 Convex</p><p>RigidArcRevolution_1 3 Arc 3 CS 8 Node 4 Convex</p><p>RigidArcRevolution_1 4 Arc 4 CS 8 Node 4 Convex</p><p>RigidArcRevolution_1 5 Arc 5 CS 8 Node 4 Concave</p><p>RigidArcRevolution_1 6 Arc 6 CS 8 Node 4 Concave</p><p>RigidArcRevolution_1 7 Arc 7 CS 8 Node 4 Convex</p><p>RigidArcRevolution_1 8 Arc 8 CS 8 Node 4 Convex</p><p>RigidArcRevolution_1 9 Arc 9 CS 8 Node 4 Convex //exterior</p><p>//Right front wheel</p><p>RigidArcRevolution_1 10 Arc 1 CS 9 Node 4 Convex //interior</p><p>RigidArcRevolution_1 11 Arc 2 CS 9 Node 4 Convex</p><p>RigidArcRevolution_1 12 Arc 3 CS 9 Node 4 Convex</p><p>RigidArcRevolution_1 13 Arc 4 CS 9 Node 4 Convex</p><p>RigidArcRevolution_1 14 Arc 5 CS 9 Node 4 Concave</p><p>RigidArcRevolution_1 15 Arc 6 CS 9 Node 4 Concave</p><p>RigidArcRevolution_1 16 Arc 7 CS 9 Node 4 Convex</p><p>RigidArcRevolution_1 17 Arc 8 CS 9 Node 4 Convex</p><p>RigidArcRevolution_1 18 Arc 9 CS 9 Node 4 Convex //exterior</p><p>//Left rail</p><p>FlexibleArcExtrusion_1 19 Arc 10 CS 6 Nodes 28 29 Convex</p><p>FlexibleArcExtrusion_1 20 Arc 11 CS 6 Nodes 28 29 Convex</p><p>FlexibleArcExtrusion_1 21 Arc 12 CS 6 Nodes 28 29 Convex</p><p>FlexibleArcExtrusion_1 22 Arc 13 CS 6 Nodes 28 29 Convex</p><p>FlexibleArcExtrusion_1 23 Arc 14 CS 6 Nodes 28 29 Convex</p><p>//Right rail</p><p>FlexibleArcExtrusion_1 24 Arc 10 CS 7 Nodes 31 30 Convex</p><p>FlexibleArcExtrusion_1 25 Arc 11 CS 7 Nodes 31 30 Convex</p><p>FlexibleArcExtrusion_1 26 Arc 12 CS 7 Nodes 31 30 Convex</p><p>FlexibleArcExtrusion_1 27 Arc 13 CS 7 Nodes 31 30 Convex</p><p>FlexibleArcExtrusion_1 28 Arc 14 CS 7 Nodes 31 30 Convex</p><p>//Left back wheel</p><p>RigidArcRevolution_1 29 Arc 1 CS 8 Node 5 Convex //interior</p><p>RigidArcRevolution_1 30 Arc 2 CS 8 Node 5 Convex</p><p>RigidArcRevolution_1 31 Arc 3 CS 8 Node 5 Convex</p><p>RigidArcRevolution_1 32 Arc 4 CS 8 Node 5 Convex</p><p>RigidArcRevolution_1 33 Arc 5 CS 8 Node 5 Concave</p><p>RigidArcRevolution_1 34 Arc 6 CS 8 Node 5 Concave</p><p>RigidArcRevolution_1 35 Arc 7 CS 8 Node 5 Convex</p><p>RigidArcRevolution_1 36 Arc 8 CS 8 Node 5 Convex</p><p>RigidArcRevolution_1 37 Arc 9 CS 8 Node 5 Convex //exterior</p><p>//Right back wheel</p><p>RigidArcRevolution_1 38 Arc 1 CS 9 Node 5 Convex //interior</p><p>RigidArcRevolution_1 39 Arc 2 CS 9 Node 5 Convex</p><p>RigidArcRevolution_1 40 Arc 3 CS 9 Node 5 Convex</p><p>RigidArcRevolution_1 41 Arc 4 CS 9 Node 5 Convex</p><p>RigidArcRevolution_1 42 Arc 5 CS 9 Node 5 Concave</p><p>RigidArcRevolution_1 43 Arc 6 CS 9 Node 5 Concave</p><p>RigidArcRevolution_1 44 Arc 7 CS 9 Node 5 Convex</p><p>RigidArcRevolution_1 45 Arc 8 CS 9 Node 5 Convex</p><p>RigidArcRevolution_1 46 Arc 9 CS 9 Node 5 Convex //exterior</p></td></tr></tbody></table>

### Creation of surface sets

After the creation of surfaces, we must define surface sets to be used in contact pairs establishment. On this tutorial, we have six different surface sets (four wheels and two rails).

<table><tbody><tr class="odd"><td><p>SurfaceSets 6</p><p>SurfaceSet 1 Surfaces 5 List 19 20 21 22 23</p><p>SurfaceSet 2 Surfaces 9 List 1 2 3 4 5 6 7 8 9</p><p>SurfaceSet 3 Surfaces 5 List 24 25 26 27 28</p><p>SurfaceSet 4 Surfaces 9 List 10 11 12 13 14 15 16 17 18</p><p>SurfaceSet 5 Surfaces 9 List 29 30 31 32 33 34 35 36 37</p><p>SurfaceSet 6 Surfaces 9 List 38 39 40 41 42 43 44 45 46</p></td></tr></tbody></table>

### Creation of contacts

Here we define contact pairs These are composed by combination of surface sets for which wheel/rail possibility. The total is four contact pairs.

<table><tbody><tr class="odd"><td><p>Contacts 4</p><p>SSSS 1 SurfaceSet1 1 SurfaceSet2 2 MU 0.3 EPN 5e8 CN 5e5 EPT 5e7 CT 0 Pinball 1.0e6 BoolTable 1</p><p>SSSS 2 SurfaceSet1 3 SurfaceSet2 4 MU 0.3 EPN 5e8 CN 5e5 EPT 5e7 CT 0 Pinball 1.0e6 BoolTable 1</p><p>SSSS 3 SurfaceSet1 1 SurfaceSet2 5 MU 0.3 EPN 5e8 CN 5e5 EPT 5e7 CT 0 Pinball 1.0e6 BoolTable 1</p><p>SSSS 4 SurfaceSet1 3 SurfaceSet2 6 MU 0.3 EPN 5e8 CN 5e5 EPT 5e7 CT 0 Pinball 1.0e6 BoolTable 1</p></td></tr></tbody></table>

### Creation of solution steps

This tutorial has four solution steps. All the steps are dynamics. The first, ends at 0.5 s and the time step is 0.1 s because it is the beginning of simulation where the bogie’s movement is a free fall (prior to contact interaction with the rails). From the second step on the contact begins to occur and the time step is reduced to 0.002 s. During the second step, the bogie accomodates on the rail. At the beginning of the third step, an initial velocity is imposed on the bogie and the simulation finishes at the fourth step after 150 s of movement.

<table><tbody><tr class="odd"><td><p>SolutionSteps 4</p><p>Dynamic 1</p><p>EndTime 0.5</p><p>TimeStep 0.1</p><p>MaxTimeStep 0.1</p><p>MinTimeStep 1e-7</p><p>MaxIt 15</p><p>MinIt 3</p><p>ConvIncrease 2</p><p>IncFactor 1.4</p><p>Sample 100</p><p>RayleighDamping Alpha 0 Beta 0 Update 0</p><p>NewmarkCoefficients Beta 0.3 Gamma 0.5</p><p>Dynamic 2</p><p>EndTime 6</p><p>TimeStep 0.002</p><p>MaxTimeStep 0.02</p><p>MinTimeStep 1e-7</p><p>MaxIt 15</p><p>MinIt 3</p><p>ConvIncrease 2</p><p>IncFactor 1.4</p><p>Sample 100</p><p>RayleighDamping Alpha 0 Beta 0 Update 0</p><p>NewmarkCoefficients Beta 0.3 Gamma 0.5</p><p>Dynamic 3</p><p>EndTime 16.2</p><p>TimeStep 0.002</p><p>MaxTimeStep 0.002</p><p>MinTimeStep 1e-9</p><p>MaxIt 15</p><p>MinIt 3</p><p>ConvIncrease 2</p><p>IncFactor 1.4</p><p>Sample 100</p><p>RayleighDamping Alpha 0 Beta 0 Update 0</p><p>NewmarkCoefficients Beta 0.3 Gamma 0.5</p><p>Dynamic 4</p><p>EndTime 150</p><p>TimeStep 0.002</p><p>MaxTimeStep 0.002</p><p>MinTimeStep 1e-9</p><p>MaxIt 15</p><p>MinIt 3</p><p>ConvIncrease 2</p><p>IncFactor 1.4</p><p>Sample 10</p><p>RayleighDamping Alpha 0 Beta 0 Update 0</p><p>NewmarkCoefficients Beta 0.3 Gamma 0.5</p></td></tr></tbody></table>

### Creation of initial conditions

At the beginning of the third solution step, the bogie starts to move. We define initial conditions of translational velocity for all the nodes of the model, such as angular velocities to the wheelsets.

<table><tbody><tr class="odd"><td><p>//Creation of initial conditions</p><p>InitialConditions 30</p><p>InitialCondition 1 Node 1 DU 1.6 0 0 OMEGA 0 0 0 SolutionStep 3</p><p>InitialCondition 2 Node 2 DU 1.6 0 0 OMEGA 0 0 0 SolutionStep 3</p><p>InitialCondition 3 Node 3 DU 1.6 0 0 OMEGA 0 0 0 SolutionStep 3</p><p>…</p></td></tr></tbody></table>

### Creation of constraints

Four nodal constraints are created to fix all the degrees of freedom of rail nodes.

<table><tbody><tr class="odd"><td><p>//Constraints</p><p>Constraints 4</p><p>NodalConstraint 1 NodeSet 9</p><p>UX BoolTable 1</p><p>UY BoolTable 1</p><p>UZ BoolTable 1</p><p>ROTX BoolTable 1</p><p>ROTY BoolTable 1</p><p>ROTZ BoolTable 1</p><p>NodalConstraint 2 NodeSet 10</p><p>UX BoolTable 1</p><p>UY BoolTable 1</p><p>UZ BoolTable 1</p><p>ROTX BoolTable 1</p><p>ROTY BoolTable 1</p><p>ROTZ BoolTable 1</p><p>NodalConstraint 3 NodeSet 11</p><p>UX BoolTable 1</p><p>UY BoolTable 1</p><p>UZ BoolTable 1</p><p>ROTX BoolTable 1</p><p>ROTY BoolTable 1</p><p>ROTZ BoolTable 1</p><p>NodalConstraint 4 NodeSet 12</p><p>UX BoolTable 1</p><p>UY BoolTable 1</p><p>UZ BoolTable 1</p><p>ROTX BoolTable 1</p><p>ROTY BoolTable 1</p><p>ROTZ BoolTable 1</p></td></tr></tbody></table>

### Creation of loads

A single nodal load is created. This is done establishing a time varying force in direction Y. This lateral force is applied to induce an oscillatory movement of the bogie. The force is applied on the bolster barycenter.

<table><tbody><tr class="odd"><td><p>//Creation of loads</p><p>Loads 1</p><p>NodalLoad 1 NodeSet 8 CS 1 NTimes 3</p><p>//Time FX FY FZ MX MY MZ</p><p>Time 10.0 0 0 0 0 0 0</p><p>Time 20.0 0 20000 0 0 0 0</p><p>Time 30.0 0 0 0 0 0 0</p></td></tr></tbody></table>

### Creation of environment data

An environment data is created to introduce the gravity field on the model.

<table><tbody><tr class="odd"><td><p>Environment</p><p>GravityData</p><p>G 0 0 -9.81 BoolTable 1</p></td></tr></tbody></table>

### Creation of solver options

Here we request Giraffe to apply six threads to solve this case. The linear system solution chosen is Direct.

<table><tbody><tr class="odd"><td><p>SolverOptions</p><p>Processors 6 LinSys Direct</p></td></tr></tbody></table>

### Creation of monitors

We request Giraffe to monitor three nodes (Bolster and wheelsets barycenter) and the four contacts.

<table><tbody><tr class="odd"><td><p>Monitor Sample 10</p><p>MonitorNodes 1 4 5</p><p>MonitorContacts 1 2 3 4</p></td></tr></tbody></table>

### Creation of post files

We choose the MagFactor parameter to be a unit value so that the Paraview<sup>TM</sup> frames reflect the real scale. Besides, the WriteRigidContactSurfaces, WriteFlexibleContactSurface, WriteSpecialConstraints, WriteContactForces and WriteRenderRigidBodies are chosen for post-processing. If the user wants to reduce Giraffe’s runtime, she/he may not choose the WriteRenderRigidBodies because it has a high computational cost.

<table><tbody><tr class="odd"><td><p>PostFiles</p><p>MagFactor 1</p><p>WriteMesh 0</p><p>WriteRenderMesh 0</p><p>WriteRigidContactSurfaces 1</p><p>WriteFlexibleContactSurfaces 1</p><p>WriteForces 0</p><p>WriteConstraints 0</p><p>WriteSpecialConstraints 1</p><p>WriteContactForces 1</p><p>WriteRenderRigidBodies 1</p><p>WriteRenderParticles 0</p></td></tr></tbody></table>

## Results and post-processing

To post-process data using Paraview<sup>TM</sup>, first we must open the files “whole\_solution\_contactforces.pvd”, “whole\_solution\_contactsurfaces.pvd” and “whole\_solution\_rb\_particles.pvd”.

<img src="../../assets/tutorials/image93.png" style="width:4.24943in;height:2.41404in" alt="" />

Figure 12‑1 - Opening the files of tutorial 12 in Paraview<sup>TM</sup>.

<img src="../../assets/tutorials/image94.png" style="width:2.70235in;height:1.06288in" alt="" />

Figure 12‑2 - Showing the Pipeline Browser.

For a more convenient view with no distorting dimensions, one may adopt an isometric perspective. For that, go to Properties and check “Camera Parallel Projection in Orientation Axes.

An interesting thing to observe during the simulation is the contact point. To get a better view of that point, the user may select the “whole\_solution\_contactsurfaces.pvd” on the Pipeline Browser. Then, select “Surface with Edges” like Figure 12‑3.

<img src="../../assets/tutorials/image95.png" style="width:1.65625in;height:1.09375in" alt="" />

Figure 12‑3 - Changing the view to surface with edges.

The Paraview<sup>TM</sup> will show the model like the Figure 12‑4.

<img src="../../assets/tutorials/image96.png" style="width:2.84932in;height:1.68366in" alt="" />

Figure 12‑4 - Visualization of results of tutorial 12 in ParaviewTM.

It may be helpful to increase the size point. For that, select “whole\_solution\_contactforces.pvd” on Pipeline Browser and change the point size at the Properties.

<img src="../../assets/tutorials/image97.png" style="width:2.14229in;height:1.83324in" alt="" />

Figure 12‑5 - Changing the point size on tutorial 12.

Next, we can create a Glyph to see the symbol of the contact forces. For that, go to the Pipeline Browser, click in “whole\_solution\_contactforces.pvd” and click in Apply. Then, click the Glyph.

After creating the Glyph, we have to select what kind of force we want to see (normal or friction force). To do this, go to the Glyph Properties and change the Vectors in Active Attributes to normal or friction.

After that, we must adjust the scale of vectors. To do this, go to Scaling, at the Properties, and set the scale mode to vector, change the scale factor to a small number like 0.000001 and click Apply.

<img src="../../assets/tutorials/image98.png" style="width:2.58872in;height:2.15207in" alt="" />

Figure 12‑6 - Adjusting the Glyph of tutorial 12.

Finally, we can show an animation of the simulation results. If the user wants to know the values of displacements of some nodes or contact forces in a certain instant, we can gather such information in the “monitors” folder, since they are requested prior to simulation solving.

<img src="../../assets/tutorials/image99.png" style="width:3.62443in;height:2.2056in" alt="" />

Figure 12‑7 - Bogie simulation.

Figure 12‑8 shows the lateral displacement of bolster barycenter over time. Note that, after 30 s, when the lateral force stops being applied, the bolster presents an oscillatory movement, which characterizes the hunting phenomenon.

<img src="../../assets/tutorials/emf_image100.png" style="max-width:100%" alt="" />

Figure 12‑8 - Lateral displacement (UY) of bolster barycenter.
