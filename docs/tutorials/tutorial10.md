# Tutorial 10: Transient dynamics of a beam

!!! example "Input files"
    [`tutorial10.inp`](../assets/tutorials/inputs/tutorial10/tutorial10.inp)
    (also in the repository under `inputs/tutorial10/`)

Input file name: tutorial10.inp

*Developed by Alfredo Gay Neto & Gabriel Maki*

<img src="../../assets/tutorials/image74.png" style="width:2.71875in;height:1.15774in" alt="" /> (a) <img src="../../assets/tutorials/image56.png" style="width:2.53371in;height:1.09375in" alt="" /> (b)

Figure 10‑1 - (a) – Beam model (b) – Beam mesh

## Objective

This tutorial presents a geometrically nonlinear simulation of a cantilever beam. The objective is show how to work with prescribed displacements and rotations in Giraffe.

## Description

The beam AB (1 m of length) is meshed using 11 nodes and 5 elements, using the element Beam\_1. There is a clamp boundary condition at point A. Besides, there is a single solution step, which is dynamic and every 2.5 seconds, a displacement or rotation is prescribed in point A.

## Input file step-by-step

### Creation of nodes

A total of 11 nodes is created to establish a finite element mesh.

<table><tbody><tr class="odd"><td><p>//Creation of nodes</p><p>Nodes 11</p><p>//Number X Y Z</p><p>Node 1 0 0 0.0</p><p>Node 2 0 0 0.1</p><p>Node 3 0 0 0.2</p><p>Node 4 0 0 0.3</p><p>Node 5 0 0 0.4</p><p>Node 6 0 0 0.5</p><p>Node 7 0 0 0.6</p><p>Node 8 0 0 0.7</p><p>Node 9 0 0 0.8</p><p>Node 10 0 0 0.9</p><p>Node 11 0 0 1.0</p></td></tr></tbody></table>

### Creation of node sets

Two node sets are created to establish displacements and constraints.

<table><tbody><tr class="odd"><td><p>//Creation of node sets</p><p>NodeSets 2</p><p>//Point A</p><p>NodeSet 1 Nodes 1 List 1</p><p>//Point B</p><p>NodeSet 2 Nodes 1 List 11</p></td></tr></tbody></table>

### Creation of elements

A total of 5 elements is defined. All the elements are assigned to the same material properties, same cross section and same coordinate system.

<table><tbody><tr class="odd"><td><p>//Creation of elements</p><p>Elements 5</p><p>Beam_1 1 Mat 1 Sec 1 CS 1 Nodes 1 2 3</p><p>Beam_1 2 Mat 1 Sec 1 CS 1 Nodes 3 4 5</p><p>Beam_1 3 Mat 1 Sec 1 CS 1 Nodes 5 6 7</p><p>Beam_1 4 Mat 1 Sec 1 CS 1 Nodes 7 8 9</p><p>Beam_1 5 Mat 1 Sec 1 CS 1 Nodes 9 10 11</p></td></tr></tbody></table>

### Creation of materials

One only material model is defined, which identification number is 1.

<table><tbody><tr class="odd"><td><p>//Creation of materials</p><p>Materials 1</p><p>Hooke 1 E 1e7 Nu 0.3 Rho 2000</p></td></tr></tbody></table>

### Creation of cross-sections

One only cross-section is defined, that is a rectangle. Both B and H are set to 0.1 m, to represent a square cross-section. The identification number of the cross-section is 1.

<table><tbody><tr class="odd"><td><p>//Creation of sections</p><p>Sections 1</p><p>Rectangle 1 B 0.1 H 0.1</p></td></tr></tbody></table>

### Creation of coordinate systems

One only coordinate system is created. It is important to mention that, since employing Beam\_1 element, the direction E3 must lie in the beam alignment direction. The identification number of the coordinate system is 1.

<table><tbody><tr class="odd"><td><p>//Creation of coordinate systems</p><p>CoordinateSystems 1</p><p>CS 1 E1 1 0 0 E3 0 0 1</p></td></tr></tbody></table>

### Creation of solution steps

This tutorial has a single solution step, which is dynamic. This is created establishing the end time equal to “10”.

<table><tbody><tr class="odd"><td><p>//Creation of the solution steps</p><p>SolutionSteps 1</p><p>//First step: displacement prescribed at point A</p><p>Dynamic 1</p><p>EndTime 10.0</p><p>TimeStep 0.005</p><p>MaxTimeStep 0.005</p><p>MinTimeStep 0.001</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.0</p><p>Sample 10</p><p>RayleighDamping Alpha 0 Beta 0 Update 0</p><p>NewmarkCoefficients Beta 0.3 Gamma 0.5</p></td></tr></tbody></table>

### Creation of displacements

A nodal displacement is created at point A. The rotation input is made in radians. The displacement/rotation imposition is done in an incremental way, like the loads. It is important to mention that**, if the degree of freedom is free, the nodal displacement related to this degree is ignored by Giraffe** (there is no way to prescribe a given displacement value to a “free” degree of freedom of the system).

<table><tbody><tr class="odd"><td><p>//Creation of displacements</p><p>Displacements 1</p><p>NodalDisplacement 1 NodeSet 1 CS 1 NTimes 5</p><p>//Time UX UY UZ RX RY RZ</p><p>0 0 0 0 0 0 0</p><p>2.5 0 0.5 +0.25 0 0 0</p><p>5.0 0 0 0 0 0 3.14</p><p>7.5 0 -0.5 -0.25 0 0 0</p><p>10 0 0 0 0 0 0</p></td></tr></tbody></table>

### Creation of constraints

A clamp boundary condition is created. This is done by fixing all the degrees of freedom during the entire simulation. It is important to remember that, in BoolTable, there is a single “1”, but it means that the constraint is applied for all the steps.

<table><tbody><tr class="odd"><td><p>//Creation of constraints</p><p>Constraints 1</p><p>NodalConstraint 1 NodeSet 1</p><p>UX BoolTable 1</p><p>UY BoolTable 1</p><p>UZ BoolTable 1</p><p>ROTX BoolTable 1</p><p>ROTY BoolTable 1</p><p>ROTZ BoolTable 1</p></td></tr></tbody></table>

### Creation of solver options

Here we request Giraffe to apply a single thread to solve this case. The linear system solution chosen is Direct.

<table><tbody><tr class="odd"><td><p>//Creation of solver options</p><p>SolverOptions</p><p>Processors 1 LinSys Direct</p></td></tr></tbody></table>

### Creation of monitors

We request Giraffe to monitor 2 nodes during the simulation: numbers 1 and 11.

<table><tbody><tr class="odd"><td><p>//Creation of monitors</p><p>Monitor Sample 1</p><p>MonitorNodes 1 11</p></td></tr></tbody></table>

### Creation of post-files

The options used here are the same as the tutorial 01, with an exception that is the “WriteForces”. On this tutorial, there is no load, so we do not have to create a file for forces.

<table><tbody><tr class="odd"><td><p>//Creation of post files</p><p>PostFiles</p><p>MagFactor 1</p><p>WriteMesh 1</p><p>WriteRenderMesh 1</p><p>WriteRigidContactSurfaces 0</p><p>WriteFlexibleContactSurfaces 0</p><p>WriteForces 0</p><p>WriteConstraints 0</p><p>WriteSpecialConstraints 0</p><p>WriteContactForces 0</p><p>WriteRenderRigidBodies 0</p><p>WriteRenderParticles 0</p></td></tr></tbody></table>

## Results and post-processing

To view the post-processing in Paraview<sup>TM</sup>, first we must open the files “whole\_solution\_mesh.pvd” and “whole\_solution\_rendermesh.pvd”, which are located at the folder “post” in “tutorial09”.

<img src="../../assets/tutorials/image75.png" style="width:5.90556in;height:3.36806in" alt="" />

Figure 10‑2 - Opening the post-files of tutorial 9 in Paraview<sup>TM</sup>.

Then, we must click the preview icon and the Apply button on the Pipeline browser.

<img src="../../assets/tutorials/image76.png" style="width:3.16757in;height:1.11318in" alt="" />

Figure 10‑3 - Showing the Pipeline Browser.

For a better view, we must use the isometric projection. To do this, go to Properties and check “Camera Parallel Projection in Orientation Axes.

Unlike other tutorials, we have no chance to create Glyphs, because there is no load in this model. So, we can run the animation. If the user wants to know the values of displacements of some nodes or the reaction forces on the clamp in a certain instant, we can acquire such information in the “monitors” folder.

In Figure 10‑4 we show a comparison between the vertical displacement at the beam tip (node 11) varying Rayleigh stiffness proportional damping, the oscillation clearly decreases when increasing damping.

<img src="../../assets/tutorials/image77.png" style="width:6.10958in;height:3.99434in" alt="" />

Figure 10‑4 – Vertical displacement UY vs. time
