# Tutorial 16: Suspension of a catenary mooring line

!!! example "Input files"
    [`tutorial16.inp`](../assets/tutorials/inputs/tutorial16/tutorial16.inp)
    (also in the repository under `inputs/tutorial16/`)

Input file name: tutorial16.inp

*Developed by Guilherme Rocha Martins*

## Objective

This tutorial presents a geometrically-nonlinear simulation of a catenary polymeric and homogeneous mooring line. The objective is to show how to suspend a catenary mooring line (or a catenary riser).

## Description 

The mooring line (835.35 m length) is meshed using 269 nodes and 134 elements, using the element Pipe\_1. There is an additional node, the pilot node attached to the seabed (to rule its behavior). The anchor (A) is fixed and has all your degrees of freedom blocked during all the analysis. The fairlead (B) is established along the simulation, while moving the mooring line tip to the desired location, at the floating unit. The mooring line is initially lying on the seabed, aligned in direction x. We consider here environment data, such as gravity and sea current loads. Next, a force of 1,000,000 N is applied at the point B in direction x and, lastly, a displacement of 186 m in direction z is prescribed to the point B to locate it close to surface level, attached to the floating unit. In this example, this point is located at 14 m depth.

## Input file step-by-step

### Creation of nodes

The mooring line and the pilot node of the seabed add up to 270 nodes. Node 1 is the anchor, node 269 is the fairlead and node 270 will be set, after, as the pilot node of the seabed. We made the choice here of considering the mooring line initially located at level z=0 and aligned in direction x.

Note that the discretization is more refined in the future expected touchdown zone of the catenary (between nodes 36 and 171).

```text
//Creation of nodes
Nodes 270
//Number X Y Z
//Anchor mooring line
Node 1 0.000000 0 0
Node 2 4.176750 0 0
Node 3 8.353500 0 0
Node 4 12.53025 0 0
Node 5 16.70700 0 0
Node 6 20.88375 0 0
Node 7 25.06050 0 0
Node 8 29.23725 0 0
Node 9 33.41400 0 0
Node 10 37.59075 0 0
...
Node 34 137.83275 0 0
Node 35 142.00950 0 0
Node 36 144.09788 0 0
Node 37 146.18625 0 0
Node 38 148.27463 0 0
Node 39 150.36300 0 0
Node 40 152.45138 0 0
...
Node 260 797.75925 0 0
Node 261 801.93600 0 0
Node 262 806.11275 0 0
Node 263 810.28950 0 0
Node 264 814.46625 0 0
Node 265 818.64300 0 0
Node 266 822.81975 0 0
Node 267 826.99650 0 0
Node 268 831.17325 0 0
//Fairlead node
Node 269 835.35000 0 0
//Pilot node
Node 270 0 0 0
```

### Creation of elements

A total of 134 elements is defined. All the elements are assigned to the same pipe cross-section (identification number 1) and same coordinate system (identification number 1). Each Pipe\_1 element is defined by three nodes, which must be equally spaced and contained in a straight line.

```text
//Creation of elements
Elements 134
Pipe\_1 1 PipeSec 1 CS 1 Nodes 1 2 3
Pipe\_1 2 PipeSec 1 CS 1 Nodes 3 4 5
Pipe\_1 3 PipeSec 1 CS 1 Nodes 5 6 7
Pipe\_1 4 PipeSec 1 CS 1 Nodes 7 8 9
Pipe\_1 5 PipeSec 1 CS 1 Nodes 9 10 11
Pipe\_1 6 PipeSec 1 CS 1 Nodes 11 12 13
Pipe\_1 7 PipeSec 1 CS 1 Nodes 13 14 15
Pipe\_1 8 PipeSec 1 CS 1 Nodes 15 16 17
Pipe\_1 9 PipeSec 1 CS 1 Nodes 17 18 19
Pipe\_1 10 PipeSec 1 CS 1 Nodes 19 20 21
...
Pipe\_1 125 PipeSec 1 CS 1 Nodes 249 250 251
Pipe\_1 126 PipeSec 1 CS 1 Nodes 251 252 253
Pipe\_1 127 PipeSec 1 CS 1 Nodes 253 254 255
Pipe\_1 128 PipeSec 1 CS 1 Nodes 255 256 257
Pipe\_1 129 PipeSec 1 CS 1 Nodes 257 258 259
Pipe\_1 130 PipeSec 1 CS 1 Nodes 259 260 261
Pipe\_1 131 PipeSec 1 CS 1 Nodes 261 262 263
Pipe\_1 132 PipeSec 1 CS 1 Nodes 263 264 265
Pipe\_1 133 PipeSec 1 CS 1 Nodes 265 266 267
Pipe\_1 134 PipeSec 1 CS 1 Nodes 267 268 269
```

### Creation of pipe sections

One only cross-section is defined. There is no need to define the material properties, because the pipe section already encompasses such data, together with cross section geometric information.

```text
//Creation of pipe sections
PipeSections 1
PS 1 EA 7.536000e+08 EI 2.763621e+05 GJ 1.905945e+05 GA 2.598621e+08 Rho 1.133500e+02 CDt 4.000e-01 CDn 2.000e+00 CAt 2.500e-01 CAn 8.000e-01 De 7.660e-02 Di 0.000e+00
```

### Creation of coordinate systems

Two coordinate systems are created. The first one is used to define the element properties. Note that, since employing Pipe\_1 element, the direction E3 must lie in the pipe alignment direction. The second coordinate system is the global one, used to define the direction of the application of the loads.

```text
//Creation of coordinate systems
CoordinateSystems 2
//Pipe\_1 Element CS (1)
CS 1 E1 0 1 0 E3 1 0 0
//Global aligned CS (2)
CS 2 E1 1 0 0 E3 0 0 1
```

### Creation of solution steps

This tutorial has four solution steps (all statics). The choices of all end-times of solution steps are arbitrary, since there is no direct meaning of physical time in statics.

The first solution step is established to impose the self-weight loading on the mooring line and its contact with the seabed. The end-time is equal to “1” and a time-step is set to “0.1”, which means that the load is divided into increments of 10%. IncFactor 2.0 is considered, which permits Giraffe to increase the time-step automatically, in case of easy convergence rate.

The second solution step is established to provide stress-stiffness to the line, by applying an axial force at point B (the node that will become the fair lead, at the end of the simulation). This step is defined with the end-time equal to “2” with time-step equal to “0.1”, which means that the load is also divided into increments of 10%. IncFactor 2.0 is considered, which permits Giraffe to increase the time-step automatically, in case of easy convergence rate.

The third load step is established to prescribe a displacement at point B, from the seabed level to the desired level of the fair lead. For that, the end-time choice is equal to “3” and the time-step is considered “0.01”, with IncFactor 1.5.

The fourth load step is established to consider the sea current load. The end-time is equal to “4”and the time-step is considered “0.01”, with IncFactor 1.5.

```text
//Creation of the solution steps
SolutionSteps 4
//Step 1
Static 1
EndTime 1
TimeStep 0.1
MaxTimeStep 1.0
MinTimeStep 0.0001
MaxIt 20
MinIt 3
ConvIncrease 1
IncFactor 2
Sample 10
//Step 2
Static 2
EndTime 2
TimeStep 0.1
MaxTimeStep 1.0
MinTimeStep 0.0001
MaxIt 30
MinIt 3
ConvIncrease 4
IncFactor 1.5
Sample 10
//Step 3
Static 3
EndTime 3
TimeStep 0.01
MaxTimeStep 0.1
MinTimeStep 0.0001
MaxIt 20
MinIt 3
ConvIncrease 4
IncFactor 1.5
Sample 10
//Step 4
Static 4
EndTime 4
TimeStep 0.01
MaxTimeStep 0.1
MinTimeStep 0.00001
MaxIt 20
MinIt 3
ConvIncrease 4
IncFactor 1.5
Sample 10
```

### Creation of points

Four points are created to establish the seabed, through two triangular surfaces.

```text
//Creation of points
Points 4
Point 1 -100 -200 0.01
Point 2 1100 -200 0.01
Point 3 1100 +200 0.01
Point 4 -100 +200 0.01
```

### Creation of surfaces

Two triangular surfaces are created to compose a rectangular surface. Both surfaces are considered attached to the same pilot node (node 270). The pilot node movement rules the movement of the surface. Both translational and rotational movements are transferred from the pilot node to all surface points, as a rigid region.

```text
//Creation of surfaces
Surfaces 2
RigidTriangularSurface\_1 1 Points 1 2 3 PilotNode 270
RigidTriangularSurface\_1 2 Points 1 3 4 PilotNode 270
```

### Creation of node sets

Four node sets are created to establish loads and constraints. The first one represents the anchor and the second represents the fairlead. The third one contains all the mooring line nodes. Lastly, the fourth node set is the pilot node, on the seabed.

```text
//Creation of node sets
NodeSets 4
//Tip - Point A
NodeSet 1 Nodes 1 List 1
//Tip - Point B
NodeSet 2 Nodes 1 List 269
//Mooring line nodes (all - except the pilot node)
NodeSet 3 Nodes 269 Sequence Initial 1 Increment 1
//Tip - Pilot node (seabed)
NodeSet 4 Nodes 1 List 270
```

### Creation of the surface sets

Surface sets are created to establish contacts. Only one surface set is created, to represent the seabed.

```text
//Creation of surface sets
SurfaceSets 1
//Seabed (ground)
SurfaceSet 1 Surfaces 2 List 1 2
```

### Creation of contacts

Only one contact region is defined: between the mooring line (node set identification number 3) and the seabed (surface set identification number 1). The contact is between a node set and a surface set, indicated in Giraffe as “NSSS”. A friction coefficient of 0.2 is considered.

```text
//Creation of contacts
Contacts 1
NSSS 1 NodeSet 3 SurfaceSet 1 MU 2.000000e-01 EPN 1.000000e+06 CN 1.000000e+02 EPT 1.000000e+05 CT 0.000000e+00 Pinball 1.762983e+03 Radius 0 MaxPointwiseInt 1 BoolTable 1
```

### Creation of constraints

Three nodal constraints are created. The BoolTable feature is here employed, to consider each constraint turned “on” or “off”, according to desire in each sequential solution step.

The first nodal constraint is applied to the anchor node (node set 1). All translations are fixed, as well as the rotation around direction x. Other rotations are free. This condition is kept along the whole simulation.

The second nodal constraint is applied to the fairlead node (node set 2). It is equal to the first nodal constraint during the first solution step. Next, then the translation in direction x is released in the steps two and three (when the axial forces and vertical translation are imposed). In the last step (when the sea current is imposed), all degrees of freedom are fixed – since it is considered fixed at a floating unit.

The last nodal constraint concerns the seabed (pilot node). All degrees of freedom are fixed along the whole the simulation.

```text
//Creation of constraints
Constraints 3
NodalConstraint 1 NodeSet 1
UX BoolTable 1
UY BoolTable 1
UZ BoolTable 1
ROTX BoolTable 1
NodalConstraint 2 NodeSet 2
UX BoolTable 1 0 0 1
UY BoolTable 1 1 1 1
UZ BoolTable 1 1 1 1
ROTX BoolTable 1 1 1 1
ROTY BoolTable 0 0 0 1
ROTZ BoolTable 0 0 0 1
NodalConstraint 3 NodeSet 4
UX BoolTable 1
UY BoolTable 1
UZ BoolTable 1
ROTX BoolTable 1
ROTY BoolTable 1
ROTZ BoolTable 1
```

### Creation of displacements

A single prescribed nodal displacement is created. It is a vertical translation at the fairlead (node set identification number 2). It is a time-series data. Giraffe automatically performs linear interpolation between the given data, if needed.

```text
//Creation of displacements
Displacements 1
NodalDisplacement 1 NodeSet 2 CS 2 NTimes 2
//Time UX UY UZ RX RY RZ
2 0 0 0 0 0 0
3 0 0 186 0 0 0
```

### Creation of loads

A single nodal load in direction x (global) is created. This is done establishing a time-varying force on the fairlead (node set identification number 2). Giraffe automatically performs linear interpolation for undefined time values following input data.

```text
//Creation of loads
Loads 1
NodalLoad 1 NodeSet 2 CS 2 NTimes 3
//Time FX FY FZ MX MY MZ
1 0 0 0 0 0 0
2 1000000 0 0 0 0 0
3 1000000 0 0 0 0 0
```

### Creation of environment data

Environment data here includes gravity and ocean data. Both are considered according to the BoolTable input to turn “on” or “off” the effects, as desired and described here in solution steps section. Gravity data is considered during all the simulation and the sea current load is considered only in the last solution step. This can be noted in the BoolTable. We consider here the surface position at 200 m level (z). The sea current is a 1.5 m/s speed aligned in 90 deg azimuth direction.

//Creation of environment

Environment

//Gravity data

GravityData

G 0 0 -9.81 BoolTable 1

OceanData

RhoFluid 1024 SurfacePosition 0 0 200

SeaCurrent N 2 BoolTable 0 0 0 1

Depth 0 Speed 1.500 Angle 90

Depth 100 Speed 1.500 Angle 90

Depth 200 Speed 1.500 Angle 90

### Creation of post-files

We choose WriteMesh, WriteRenderMesh and WriteRigidContactSurfaces options for post-processing using Paraview<sup>TM</sup>. The MagFactor parameters can be used to generate plots with exaggerated displacement patterns, to help to visualize in case of small displacements. This is not the case of this example. Then, we will assume a unit value, then, requesting that the deformed shape of Paraview<sup>TM</sup> frames will reflect the real scale. WriteForces and WriteContactForces are also considered for post-processing.

//Creation of post files

PostFiles

MagFactor 1

WriteMesh 1

WriteRenderMesh 1

WriteRigidContactSurfaces 0

WriteFlexibleContactSurfaces 1

WriteForces 1

WriteConstraints 0

WriteSpecialConstraints 0

WriteContactForces 1

WriteRenderRigidBodies 0

WriteRenderParticles 0

### Creation of monitors

We request Giraffe to monitor 6 nodes during the simulation (numbers 2, 36, 88, 90, 170, 268) and 6 elements (numbers 1, 18, 44, 45, 85, 134). Note that we entered a sampling of 1 for saving data for each converged time-step. In case we desire smaller files, we less data, we could set larger sample rates.

```text
//Creation of monitors
Monitor Sample 1
MonitorNodes 10 1 2 36 88 90 170 268
MonitorElements 6 1 18 44 45 85 134
```

### Creation of convergence criteria

Some convergence and divergences tolerances are defined.

//Creation of convergence criteria

ConvergenceCriteria

ForceTolerance 1.0e-4

MomentTolerance 1.0e-4

ForceMinimumReference 1.0e-3

MomentMinimumReference 1.0e-3

ConstraintMinimumReference 1.0e-7

DisplacementTolerance 1.0e-4

RotationTolerance 1.0e-4

LagrangeTolerance 1.0e-4

DisplacementMinimumReference 1.0e-6

RotationMinimumReference 1.0e-6

LagrangeMinimumReference 1.0e-6

DivergenceReference 1.0e+15

### Creation of solver options

Here we request Giraffe to apply four cores to solve this case. The linear system solution chosen is Direct.

//Creation of solver options

SolverOptions

Processors 4 LinSys Direct

Note that the Giraffe input file commands not necessarily has to follow the here proposed sequence. The important is that all the commands here presented are included.

## Results and post-processing

Result files are automatically saved during the solving process on the folder “/tutorial17/post/” and “/tutorial17/monitors/”

To post-process results in Paraview<sup>TM</sup>, first we have to open the files “whole\_solution\_mesh.pvd”, “whole\_solution\_rendermesh.pvd”, “whole\_solution\_forces.pvd” “whole\_solution\_contactsurfaces.pvd”, which are located at the folder “post” in “tutorial17”.

<img src="../../assets/tutorials/image130.png" style="width:5.90556in;height:4.12361in" alt="" />

Figure 16‑1 - Opening the files of tutorial 17 in Paraview<sup>TM</sup>

Then we can click the preview icon and the “Apply button” on the Pipeline browser. Then, click in the Glyph feature button to create a Glyph, after doing that, we are ready to post-process the results.

<img src="../../assets/tutorials/image131.png" style="width:3.01667in;height:2.77582in" alt="" />

Figure 16‑2 - Using Glyphs in Paraview<sup>TM</sup>

### Paraview<sup>TM</sup> post-processing

As we have done, in order to avoid dimensions distortions, we must check the “Camera Parallel Projection” in “Orientation Axes”.

<img src="../../assets/tutorials/image132.png" style="width:5.825in;height:3.98983in" alt="" />

Figure 16‑3 - Checking the "Camera Parallel Projection" in Paraview<sup>TM</sup>

Then we able to create glyphs. It is interesting to visualize the contact forces between seabed and mooring line and the force applied to the fairlead node. To do that, go to the pipeline browser and click in “whole\_solution\_contactforces.pvd” than click in the glyph button. Do the same for the forces, click in “whole\_solution\_forces.pvd” than click in the glyph button.

<img src="../../assets/tutorials/image133.png" style="width:3.02381in;height:3.30329in" alt="" />

Figure 16‑4 - Creating the contact force’s glyph

Note that the size of each glyph may be set as proportional to the magnitude of force. If this is the case, the contact forces in the more refined region of the mesh (nearby touchdown point) are smaller than in the coarser region. This reflects the contact modeling, which considers nodal detection and pointwise (nodal) action.

<img src="../../assets/tutorials/image134.png" style="width:2.24405in;height:3.7268in" alt="" />

Figure 16‑5 - Creating the force's glyph

After creating the glyphs, we have to adjust the scale of the vectors. To do this, first click in the glyph of interest preview icon and. Than in “Orientation” set the orientation array to “normal” and in “Scale” set the scale array also to “normal”. Finally, set the vector scale mode to “scale by magnitude” and change the scale factor to a small number like 0.005 and click Apply. For the force’s glyph the process is analogous, but we choose “force” in the orientation array and scale array, and we can set a smallest scale factor, like 0.0002, as indicated in Figure 16‑7.

<img src="../../assets/tutorials/image135.png" style="width:5.90556in;height:3.28542in" alt="" />

Figure 16‑6 - Changing the properties of contact force's glyph

<img src="../../assets/tutorials/image136.png" style="width:5.90556in;height:3.27361in" alt="" />

Figure 16‑7 - Changing the properties of force's glyph

For a better visualization, we can chance the background color. To do that, go to “Background”, in view. Here, we apply a gradient, observing that the color from below it’s best to be darker, to simulate the seabed.

<img src="../../assets/tutorials/image137.png" style="width:5.90556in;height:3.32222in" alt="" />

Figure 16‑8 - Changing the background color

For a better view, we can also enhance the line width and chance the color. To do that, click in “whole\_solution\_mesh.pvd” icon and go to properties. The color can be changing in “coloring” and the line width in “Styling”, as showing in Figure 16‑9.

<img src="../../assets/tutorials/image138.png" style="width:5.90556in;height:3.28571in" alt="" />

Figure 16‑9 - Changing the mesh properties

Finally, we can show an animation of the simulation results. If the user wants to know the values of forces and moments of some elements in a certain instant, we can gather such information in the “monitors” folder, since they are requested prior to simulation solving.

<img src="../../assets/tutorials/image139.png" style="width:5.90556in;height:3.29167in" alt="" />

Figure 16‑10 - Configuration established before apply sea current

### Other post-processing

It is interesting to evaluate the bending moment and the shear force along the mooring line length, especially close to the touchdown zone. First, we have to open the last solution file of the third solution step (before applying the sea current) and copy it to a spreadsheet. The same analysis can be done after sea current is considered, here we opted to analyze before, to avoid the three-dimensional analysis. So, we can now go to the folder “post”, then to “solution\_3” and, finally, open the last converged solution.

<img src="../../assets/tutorials/image140.png" style="width:5.90556in;height:3.36736in" alt="" />

Figure 16‑11 - Opening the solution report

To plot the graphics, attempt to the coordinate system and remember that the elements results are given in the local coordinate system (of the element).

Said that, it is our interest here to evaluate the shear force in local direction *e*<sub>2</sub> and the bending moment around the local axis *e*<sub>1</sub>.

<img src="../../assets/tutorials/chart_t16_f12.png" style="max-width:100%" alt="" />

Figure 16‑12 - Bending moment along the mooring line length

<img src="../../assets/tutorials/chart_t16_f13.png" style="max-width:100%" alt="" />

Figure 16‑13 - Shear force along the mooring line length
