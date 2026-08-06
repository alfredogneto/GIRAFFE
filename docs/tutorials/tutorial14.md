# Tutorial 14: Contact between cantilever beams

!!! example "Input files"
    [`tutorial14.inp`](../assets/tutorials/inputs/tutorial14/tutorial14.inp)
    (also in the repository under `inputs/tutorial14/`)

Input file name: tutorial14.inp

*Developed by Debora Naomi Higa & Alfredo Gay Neto*

## Objective

This tutorial presents a simulation of contact between cantilever beams. The aim is to simulate the process of contact between beams with super-elliptical cross section subject to a prescribed axial rotation - applying nodal displacements and boundary conditions and solving a static geometrically-nonlinear analysis.

<img src="../../assets/tutorials/image113.png" style="width:5.90556in;height:2.11667in" alt="" />

Figure 14‑1 – Beams rendering view.

## Description 

Beams AB and CD (both with 3.0 m length) are meshed using 31 nodes and 15 elements each, using the element Beam\_1. Beam AB is aligned to X direction and fixed at point A. Beam CD is aligned to Y direction.

In a first step, a vertical displacement is imposed at point C in order to establish contact between beams. In the sequence, a rotation is prescribed about Y axis at point C, thus creating an oscillating contact force pattern.

## Input file step-by-step

### Creation of nodes

62 nodes are created to establish a finite element mesh (31 nodes for each beam).

```text
//Creation of nodes
Nodes 62
//Number X Y Z
Node 1 0 2.5 0
Node 2 0.1 2.5 0
Node 3 0.2 2.5 0
Node 4 0.3 2.5 0
Node 5 0.4 2.5 0
...
Node 58 2.5 2.6 0.3
Node 59 2.5 2.7 0.3
Node 60 2.5 2.8 0.3
Node 61 2.5 2.9 0.3
Node 62 2.5 3 0.3
```

### Creation of node sets

Two node sets are created to establish loads and constraints.

```text
//Creation of node sets
NodeSets 2
//Point A
NodeSet 1 Nodes 1 List 1
//Point B
NodeSet 2 Nodes 1 List 32
```

### Creation of materials

Two material models are defined, which identification numbers are 1 and 2.

```text
//Creation of materials
Materials 2
//Material 1
Hooke 1 E 2e9 Nu 0.30 Rho 8000
//Material 2
Hooke 2 E 2e9 Nu 0.30 Rho 8000
```

### Creation of cross-sections

Two cross-sections are defined, both super ellipses. The semi-axis A and B are defined lying in direction E1 and E2, respectively. Both super ellipse exponent value N are set by 3. The parameter AMeshFDM establishes the discretization employed by Giraffe to evaluate the St. Venant torsion constant (by the finite-difference method). This is done as a pre-processing. The identification numbers of the cross-sections are 1 and 2.

```text
//Creation of sections
Sections 2
//Sesction 1
SuperEllipse 1 A 0.2 B 0.06 N 3 AMeshFDM 200
//Section 2
SuperEllipse 2 A 0.06 B 0.1 N 3 AMeshFDM 200
```

### Creation of coordinate systems

Three coordinate systems are input. These are used to define the element properties and the nodal displacements. The identification number of the coordinate systems are 1, 2 and 3. Note that, since employing Beam\_1 element (shown in sequence), the direction E3 must lie in the beam alignment direction

```text
//Creation of coordinate systems
CoordinateSystems 3
CS 1 E1 0 1 0 E3 1 0 0
CS 2 E1 1 0 0 E3 0 1 0
CS 3 E1 1 0 0 E3 0 0 1
```

### Creation of elements

Thirty elements are defined (fifteen for each beam). The elements are assigned to different material properties (identification number 1 or 2), different cross section (identification number 1 or 2) and different coordinate system (identification number 1 or 2). Each Beam\_1 element is defined by three nodes, which **must** be equally spaced and contained in a straight line.

**Note**: In case of creation of not-aligned beam elements (for example, for representing an approximation of a curved beam or a frame), each beam element must be associated with a defined coordinate system, which direction E3 is aligned with the element direction. This is done in order to permit creation of arbitrary alignment between successive cross sections of a beam-like structure, composed possibly by many cross-sections or alignments of these.

```text
//Creation of elements
Elements 30
Beam\_1 1 Mat 1 Sec 1 CS 1 Nodes 1 2 3
Beam\_1 2 Mat 1 Sec 1 CS 1 Nodes 3 4 5
Beam\_1 3 Mat 1 Sec 1 CS 1 Nodes 5 6 7
Beam\_1 4 Mat 1 Sec 1 CS 1 Nodes 7 8 9
Beam\_1 5 Mat 1 Sec 1 CS 1 Nodes 9 10 11
Beam\_1 6 Mat 1 Sec 1 CS 1 Nodes 11 12 13
Beam\_1 7 Mat 1 Sec 1 CS 1 Nodes 13 14 15
Beam\_1 8 Mat 1 Sec 1 CS 1 Nodes 15 16 17
Beam\_1 9 Mat 1 Sec 1 CS 1 Nodes 17 18 19
Beam\_1 10 Mat 1 Sec 1 CS 1 Nodes 19 20 21
Beam\_1 11 Mat 1 Sec 1 CS 1 Nodes 21 22 23
Beam\_1 12 Mat 1 Sec 1 CS 1 Nodes 23 24 25
Beam\_1 13 Mat 1 Sec 1 CS 1 Nodes 25 26 27
Beam\_1 14 Mat 1 Sec 1 CS 1 Nodes 27 28 29
Beam\_1 15 Mat 1 Sec 1 CS 1 Nodes 29 30 31
Beam\_1 16 Mat 2 Sec 2 CS 2 Nodes 32 33 34
Beam\_1 17 Mat 2 Sec 2 CS 2 Nodes 34 35 36
Beam\_1 18 Mat 2 Sec 2 CS 2 Nodes 36 37 38
Beam\_1 19 Mat 2 Sec 2 CS 2 Nodes 38 39 40
Beam\_1 20 Mat 2 Sec 2 CS 2 Nodes 40 41 42
Beam\_1 21 Mat 2 Sec 2 CS 2 Nodes 42 43 44
Beam\_1 22 Mat 2 Sec 2 CS 2 Nodes 44 45 46
Beam\_1 23 Mat 2 Sec 2 CS 2 Nodes 46 47 48
Beam\_1 24 Mat 2 Sec 2 CS 2 Nodes 48 49 50
Beam\_1 25 Mat 2 Sec 2 CS 2 Nodes 50 51 52
Beam\_1 26 Mat 2 Sec 2 CS 2 Nodes 52 53 54
Beam\_1 27 Mat 2 Sec 2 CS 2 Nodes 54 55 56
Beam\_1 28 Mat 2 Sec 2 CS 2 Nodes 56 57 58
Beam\_1 29 Mat 2 Sec 2 CS 2 Nodes 58 59 60
Beam\_1 30 Mat 2 Sec 2 CS 2 Nodes 60 61 62
```

### Creation of contact surfaces

Thirty contact surfaces are defined (fifteen for each beam), all flexible super elliptical cylinder surfaces. The cylinder surfaces semi-axes (A and B) and exponent (N) are set in a similar way from the creation of cross-sections. Coordinate systems (CS) chosen are same employed during the creation of elements (identification number 1 or 2). The normal direction of the surfaces points outwards of the super elliptical cylinder. Each cylinder surface is defined by two nodes.

```text
//Creation of contact surfaces
Surfaces 30
FlexibleSECylinder\_1 1 A 0.2 B 0.06 N 3 CS 1 NormalExterior Nodes 1 3
FlexibleSECylinder\_1 2 A 0.2 B 0.06 N 3 CS 1 NormalExterior Nodes 3 5
FlexibleSECylinder\_1 3 A 0.2 B 0.06 N 3 CS 1 NormalExterior Nodes 5 7
...
```

### Creation of surface sets

Surface sets are created to set the potential contact bodies. Two surface sets are defined (one for each beam) input in a sequence. Both surfaces sets have 15 surfaces, which are 1 to 15 (SurfaceSet 1) and 16 to 30 (SurfaceSet 2).

```text
//Creation of surface sets
SurfaceSets 2
SurfaceSet 1 Surfaces 15 Sequence Initial 1 Increment 1
SurfaceSet 2 Surfaces 15 Sequence Initial 16 Increment 1
```

### Creation of displacements

Two nodal displacements are created (identification numbers 1 and 2), both with three time steps (NTimes) and coordinate system identification 3. The first prescribes no displacement to node set 1. The second prescribes displacements and rotations to node set 2. From time 0 to time 1 a displacement is linearly prescribed in Z direction (-0.5). From time 1 to time 2 a linearly varying displacement in Y direction of two units is prescribed. The displacement in Z direction (-0.5) is kept. From time 1 to time 2 a large rotation in Y axis is also linearly prescribed (25 rad).

```text
//Creation of displacements
Displacements 2
NodalDisplacement 1 NodeSet 1 CS 3 NTimes 3
//Time UX UY UZ ROTX ROTY ROTZ
0 0 0 0 0 0 0
1 0 0 0 0 0 0
2 0 0 0 0 0 0
NodalDisplacement 2 NodeSet 2 CS 3 NTimes 3
//Time UX UY UZ ROTX ROTY ROTZ
0 0 0 0 0 0 0
1 0 0 -0.5 0 0 0
2 0 2.0 -0.5 0 25 0
```

### Creation of constraints

Two nodal constraints are created. Both are done fixing the following DOFs for the node sets 1 and 2: UX, UY, UZ, ROTX, ROTY and ROTZ.

```text
//Creation of constraints
Constraints 2
NodalConstraint 1 NodeSet 1
UX BoolTable 1
UY BoolTable 1
UZ BoolTable 1
ROTX BoolTable 1
ROTY BoolTable 1
ROTZ BoolTable 1
NodalConstraint 2 NodeSet 2
UX BoolTable 1
UY BoolTable 1
UZ BoolTable 1
ROTX BoolTable 1
ROTY BoolTable 1
ROTZ BoolTable 1
```

### Creation of contacts

One contact constraint for the interaction between two surface sets (SSSS) is defined. The surface sets are identified by numbers 1 and 2. The coefficient of friction value (MU) is set by 0.3.

Constraints enforcements are done by the Penalty Method. Thus, it is necessary for the user to input penalty parameters data. Usually these may be calibrated based on physical information related to the desired scenario, basing on equivalent local stiffness, leading to allowable penetration on each contact zone. In this case, penalty coefficient to enforce normal and tangential contact constraints (EPN and EPT) are respectively set to 1e7 and 1e6.

Normal and tangential damping parameter coefficients (CN and CT) are useful for dissipation of energy during impact simulations, avoiding high frequency oscillations on contact forces. In this case, there is no damping on contact since the coefficients are set to 0.

Pinball radius value is a rough search geometrical parameter used by Giraffe to establish probable and not probable contact interactions. The pinball radius value is set to 1.

The maximum number of contact pointwise interactions between surfaces (MaxPointwiseInt) is set to 1. Note that on non-convexity scenarios, surfaces have the possibility of seeking for more than one pointwise contact solution.

Bool table data for current contact constraint (BoolTable) is set to 1. It is an optional keyword and it permits to create a scenario in which the contact constraint is turned on/off along solution steps. If BoolTable is not included, Giraffe assumes that the contact constraint will be considered “on” for all solution steps.

```text
//Creation of contacts
Contacts 1
SSSS 1 SurfaceSet1 1 SurfaceSet2 2 MU 0.3 EPN 1e7 CN 0 EPT 1e6 CT 0 Pinball 1.0 MaxPointwiseInt 1 BoolTable 1
```

### Creation of solution steps

This tutorial has two solution steps (both are statics). The first step is created establishing the end time equal to “1”. The time-step is set to “0.05”, which means that the load is divided into increments of 5% along the first solution step. The maximum time-step is set by “0.05” and the minimum time-step is set by “1e-4” which permits Giraffe to automatically decrease the time-step in case of difficulties of convergence along the simulation. The second step is created establishing the end time equal to “2”, time-step “0.005”, maximum time-step “0.1” and minimum time-step “1e-6”. The sampling for saving post-processing files is set to “1” which ensures saving post-processing files for each converged configuration.

```text
//Creation of the solution steps
SolutionSteps 2
Static 1
EndTime 1
TimeStep 0.05
MaxTimeStep 0.05
MinTimeStep 1e-4
MaxIt 15
MinIt 3
ConvIncrease 2
IncFactor 1.4
Sample 1
Static 2
EndTime 2
TimeStep 0.005
MaxTimeStep 0.1
MinTimeStep 1e-6
MaxIt 15
MinIt 3
ConvIncrease 2
IncFactor 1.4
Sample 1
```

### Creation of post-files

The options WriteMesh, WriteRenderMesh, WriteFlexibleContactSurfaces and WriteContactForces are chosen for post-processing using Paraview<sup>TM</sup>. The parameter MagFactor assumes a unit value so the deformed shape of Paraview<sup>TM</sup> frames reflects the real scale.

//Creation of post files

PostFiles

MagFactor 1

WriteMesh 1

WriteRenderMesh 1

WriteRigidContactSurfaces 0

WriteFlexibleContactSurfaces 1

WriteForces 0

WriteConstraints 0

WriteSpecialConstraints 0

WriteContactForces 1

WriteRenderRigidBodies 0

WriteRenderParticles 0

### Convergence criteria

Here, we define some convergence and divergence tolerances. For detailed explanation on the meaning of each of these quantities, refer to Giraffe users’ manual. Altering default convergence criteria is recommended only for advanced users.

//Convergence criteria

ConvergenceCriteria

ForceTolerance 1e-4

MomentTolerance 1e-4

ForceMinimumReference 1e-3

MomentMinimumReference 1e-3

ConstraintMinimumReference 1e-7

DisplacementTolerance 1e-4

RotationTolerance 1e-4

LagrangeTolerance 1e-4

DisplacementMinimumReference 1e-5

RotationMinimumReference 1e-5

LagrangeMinimumReference 1e-6

DivergenceReference 1e+15

### Creation of monitors

It is requested Giraffe to monitor the contact 1 (SSSS) created.

<table><tbody><tr class="odd"><td><p>//Creation of monitors</p><p>Monitor Sample 1</p><p>MonitorContacts 1</p></td></tr></tbody></table>

## Results and post-processing

Result files are automatically saved during the solving process on the folder “/tutorial15/post/”.

### Paraview<sup>TM</sup> post-processing

To view the post-processing in Paraview<sup>TM</sup>, first we must open the files “whole\_solution\_contact\_forces.pvd”, “whole\_solution\_contact\_forces.pvd”, “whole\_solution\_mesh.pvd” and “whole\_solution\_rendermesh.pvd”, which are located at the folder “post” in “tutorial15”.

<img src="../../assets/tutorials/image114.png" style="width:5.89653in;height:3.27361in" alt="" />

Figure 14‑2 – Opening the files in Paraview<sup>TM</sup>

Then, we must click the preview icon and the Apply button on the Pipeline browser.

<img src="../../assets/tutorials/image115.png" style="width:2.08598in;height:2.966in" alt="" />

Figure 14‑3 – Pipeline browser in Paraview<sup>TM</sup>

As a particular feature, forces output may be used in Paraview<sup>TM</sup> together with Glyphs in order to establish special visualizations. For that, select “whole\_solution\_contactforces.pvd” in the Pipeline Browser. Then, click in the Glyph feature button and in the Apply button to create a Glyph.

<img src="../../assets/tutorials/image116.png" style="width:2.52989in;height:4.82626in" alt="" />

Figure 14‑4 – Using Glyphs in Paraview<sup>TM</sup>

After creating the Glyph, select it on Pipeline Browser to modify its default options as depicted next. Glyph Type is set to “Arrow”, Vectors is set to “normal” and Glyph Mode is set to “All Points”. Scale Mode is set to “vector” and scale factor is set to 0.0001. This value is arbitrary and may vary from model to model and controls the size of the created Glyph feature. It should be chosen based on the desired size of the Glyph feature on the visualization of the model.

<img src="../../assets/tutorials/image117.png" style="width:2.69388in;height:7.44241in" alt="" />

Figure 14‑5 – Setting Glyph attributes in Paraview<sup>TM</sup>

Visualization will look like the shown in Figure 14‑6:

<img src="../../assets/tutorials/image118.jpeg" style="width:5.90556in;height:3.19792in" alt="" />

Figure 14‑6 – Using Glyphs in Paraview<sup>TM</sup>

The user may animate results to see the Glyph feature changing its size, according to the force magnitude along time.

Giraffe generates a monitor text file for the requested contact containing the time series of the forces. For example we can plot the time series of the normal forces.

<img src="../../assets/tutorials/image119.png" style="width:4.984in;height:3.10788in" alt="" />

Figure 14‑7 – Normal force vs time.

## What to try now?

-   change the magnitude and the direction of the force

-   change the beams cross section

-   change material properties of one beam

-   change material properties of both beams

-   create different displacement compositions

-   change the boundary conditions (fixed and free nodes)
