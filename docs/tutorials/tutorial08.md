# Tutorial 8: Evaluation of system modal analysis during a static simulation

!!! example "Input files"
    [`tutorial08.inp`](../assets/tutorials/inputs/tutorial08/tutorial08.inp)
    (also in the repository under `inputs/tutorial08/`)

Input file name: tutorial08.inp

*Developed by Alfredo Gay Neto & Gabriel Maki*

<img src="../../assets/tutorials/image64.png" style="width:2.5625in;height:1.07098in" alt="" />(a) <img src="../../assets/tutorials/image56.png" style="width:2.39583in;height:1.03423in" alt="" /> (b)

Figure 8‑1 - (a) – Beam model (b) Beam mesh render view

## Objective

This tutorial presents a geometrically-nonlinear simulation of a cantilever beam. The aim is to give an example on how to work with concomitant solutions.

## Description

The beam AB (1 m of length) is meshed using 11 nodes and 5 elements, using the element Beam\_1. There is a clamp boundary condition at point A. Two loads are applied sequentially. First, a force F is applied at point B. then, a moment M is applied at the same point. The simulation is done using three static solution steps and, during the first and third solution steps, a concomitant modal analysis is done to evaluate natural frequencies for each deformed (and stressed) configuration.

## Input file step-by-step

### Creation of nodes

A total of 11 nodes is created to establish a finite element mesh.

<table><tbody><tr class="odd"><td><p>//Creation of nodes</p><p>Nodes 11</p><p>//Number X Y Z</p><p>Node 1 0 0 0.0</p><p>Node 2 0 0 0.1</p><p>Node 3 0 0 0.2</p><p>Node 4 0 0 0.3</p><p>Node 5 0 0 0.4</p><p>Node 6 0 0 0.5</p><p>Node 7 0 0 0.6</p><p>Node 8 0 0 0.7</p><p>Node 9 0 0 0.8</p><p>Node 10 0 0 0.9</p><p>Node 11 0 0 1.0</p></td></tr></tbody></table>

### Creation of node sets

Two node sets are created to establish loads and constraints.

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

This tutorial has three solution steps. All the steps are statics and have a time step of 0.1 which means that the loads are applied or reduced in fractions of 10%.

<table><tbody><tr class="odd"><td><p>/Creation of the solution steps</p><p>SolutionSteps 3</p><p>//First step: FX applied at the tip of the beam</p><p>Static 1</p><p>EndTime 1</p><p>TimeStep 0.1</p><p>MaxTimeStep 0.1</p><p>MinTimeStep 0.01</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.0</p><p>Sample 2</p><p>//Second step: FX kept and MZ applied at the tip of the beam</p><p>Static 2</p><p>EndTime 2</p><p>TimeStep 0.1</p><p>MaxTimeStep 0.1</p><p>MinTimeStep 0.01</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.0</p><p>Sample 2</p><p>//Third step: all loads are removed (no loads at the end)</p><p>Static 3</p><p>EndTime 3</p><p>TimeStep 0.1</p><p>MaxTimeStep 0.1</p><p>MinTimeStep 0.01</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.0</p><p>Sample 2</p></td></tr></tbody></table>

### Creation of concomitant solutions

A concomitant solution is created to run a modal analysis during the first and third solution step. In this case, we request Giraffe to evaluate the smallest ten natural frequencies related eigenvalues. The BoolTable keyword is employed to choose only the first and the third solution steps, to have the concomitant solution evaluated with.

<table><tbody><tr class="odd"><td><p>//Creation of concomitant solution</p><p>ConcomitantSolution</p><p>Sample 3 BoolTable 1 0 1</p><p>Modal NumberModes 10 Tolerance 1E-6</p></td></tr></tbody></table>

### Creation of loads

A single nodal load is created. Note that the load application time is associated with the start and end time of the solution steps.

<table><tbody><tr class="odd"><td><p>//Creation of loads</p><p>Loads 1</p><p>NodalLoad 1 NodeSet 2 CS 1 NTimes 4</p><p>//Time FX FY FZ MX MY MZ</p><p>0 0 0 0 0 0 0</p><p>1 1000 0 0 0 0 0</p><p>2 1000 0 0 0 0 250</p><p>3 0 0 0 0 0 0</p></td></tr></tbody></table>

### Creation of constraints

A clamp boundary condition is created. This is done by fixing all the degrees of freedom during the entire simulation. This is done in node set 1.

<table><tbody><tr class="odd"><td><p>//Creation of constraints</p><p>Constraints 1</p><p>NodalConstraint 1 NodeSet 1</p><p>UX BoolTable 1</p><p>UY BoolTable 1</p><p>UZ BoolTable 1</p><p>ROTX BoolTable 1</p><p>ROTY BoolTable 1</p><p>ROTZ BoolTable 1</p></td></tr></tbody></table>

### Creation of solver options

Here we request Giraffe to apply a single thread to solve this case. The linear system solution chosen is Direct.

<table><tbody><tr class="odd"><td><p>//Creation of solver options</p><p>SolverOptions</p><p>Processors 1 LinSys Direct</p></td></tr></tbody></table>

### Creation of monitors

We request Giraffe to monitor 2 nodes during the simulation: numbers 1 and 11.

<table><tbody><tr class="odd"><td><p>//Creation of monitors</p><p>Monitor Sample 1</p><p>MonitorNodes 1 11</p></td></tr></tbody></table>

### Creation of post-files

The options used here are the same as the tutorial 01, with an exception that is the “WriteForces”. This was turned off to avoid saving many solution files, unnecessary for present objectives.

<table><tbody><tr class="odd"><td><p>//Creation of post files</p><p>PostFiles</p><p>MagFactor 1</p><p>WriteMesh 1</p><p>WriteRenderMesh 1</p><p>WriteRigidContactSurfaces 0</p><p>WriteFlexibleContactSurfaces 0</p><p>WriteForces 0</p><p>WriteConstraints 0</p><p>WriteSpecialConstraints 0</p><p>WriteContactForces 0</p><p>WriteRenderRigidBodies 0</p><p>WriteRenderParticles 0</p></td></tr></tbody></table>

## Results and post-processing

To post-process results in Paraview<sup>TM</sup>, first we must open the files “whole\_solution\_mesh.pvd” and “whole\_solution\_rendermesh.pvd”, which are located at the folder “post” in “tutorial08”.

<img src="../../assets/tutorials/image65.png" style="width:4.97761in;height:2.86458in" alt="" />

Figure 8‑2 - Opening the files of tutorial 7 in Paraview<sup>TM</sup>.

Then, we have to click the preview icon and the Apply button on the Pipeline browser.

<img src="../../assets/tutorials/image66.png" style="width:2.52083in;height:0.99734in" alt="" />

Figure 8‑3 - Showing the Pipeline Browser.

The user may note that there is no file for forces, so we do not create Glyphs. But, unlike other tutorials, there is a folder called “concomitant\_solution”. Inside this folder, we find a text file with the same name. This file contains all the eigenvalues associated with the vibration modes of the structure, evaluated along time. That is, the eigenvalues are calculated while the structure is deforming.

<img src="../../assets/tutorials/emf_image67.png" style="max-width:100%" alt="" />

Figure 8‑4 – Smallest eigenvalue vs. time.
