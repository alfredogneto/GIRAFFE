# Tutorial 9: Evaluation of a system modal analysis during buckling of a cantilever beam

!!! example "Input files"
    [`tutorial09.inp`](../assets/tutorials/inputs/tutorial09/tutorial09.inp)
    (also in the repository under `inputs/tutorial09/`)

Input file name: tutorial09.inp

*Developed by Alfredo Gay Neto & Gabriel Maki*

<img src="../../assets/tutorials/image68.png" style="width:2.55311in;height:1.17708in" alt="" />(a) <img src="../../assets/tutorials/image56.png" style="width:2.39583in;height:1.03423in" alt="" />(b)

Figure 9‑1 - (a) – Beam subjected to buckling (b) Beam mesh render view

## Objective

This tutorial presents a geometrically-nonlinear simulation of a cantilever beam. The aim is to perform a modal analysis during the nonlinear buckling simulation of the structure. For that, it is necessary to introduce an imperfection on the beam and the trigger load, as shown in Figure 9‑1(a). Besides, the resource of concomitant solution is used on this tutorial too.

## Description

The beam AB (1 m of length) is meshed using 51 nodes and 25 elements, using the element Beam\_1. There is a clamp boundary condition at point A. Two nodal loads are applied sequentially at point B. First, a lateral load is applied. After, a compressive load is applied. The lateral force is necessary because it is not possible to conduct a nonlinear buckling analysis without an initial imperfection on the structure.

## Input file step-by-step

### Creation of nodes

A total of 51 nodes is created to establish a finite element mesh.

<table><tbody><tr class="odd"><td><p>//Creation of nodes</p><p>Nodes 51</p><p>//Number X Y Z</p><p>Node 1 0 0 0.0</p><p>Node 2 0 0 0.02</p><p>Node 3 0 0 0.04</p><p>Node 4 0 0 0.06</p><p>Node 5 0 0 0.08</p><p>…</p><p>Node 45 0 0 0.88</p><p>Node 46 0 0 0.9</p><p>Node 47 0 0 0.92</p><p>Node 48 0 0 0.94</p><p>Node 49 0 0 0.96</p><p>Node 50 0 0 0.98</p><p>Node 51 0 0 1.0</p></td></tr></tbody></table>

### Creation of node sets

Two node sets are created to establish loads and constraints.

<table><tbody><tr class="odd"><td><p>//Creation of node sets</p><p>NodeSets 2</p><p>//Point A</p><p>NodeSet 1 Nodes 1 List 1</p><p>//Point B</p><p>NodeSet 2 Nodes 1 List 51</p></td></tr></tbody></table>

### Creation of elements

A total of 25 elements is defined. All the elements are assigned to the same material properties, same cross section and same coordinate system.

<table><tbody><tr class="odd"><td><p>//Creation of elements</p><p>Elements 25</p><p>Beam_1 1 Mat 1 Sec 1 CS 1 Nodes 1 2 3</p><p>Beam_1 2 Mat 1 Sec 1 CS 1 Nodes 3 4 5</p><p>Beam_1 3 Mat 1 Sec 1 CS 1 Nodes 5 6 7</p><p>Beam_1 4 Mat 1 Sec 1 CS 1 Nodes 7 8 9</p><p>…</p><p>Beam_1 22 Mat 1 Sec 1 CS 1 Nodes 43 44 45</p><p>Beam_1 23 Mat 1 Sec 1 CS 1 Nodes 45 46 47</p><p>Beam_1 24 Mat 1 Sec 1 CS 1 Nodes 47 48 49</p><p>Beam_1 25 Mat 1 Sec 1 CS 1 Nodes 49 50 51</p></td></tr></tbody></table>

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

We create only two solution steps, both statics. The first step is associated with the initial imperfection that we must introduce for conduct a nonlinear buckling analysis. The second step is associated with the compression force that leads to buckling.

<table><tbody><tr class="odd"><td><p>//Creation of the solution steps</p><p>SolutionSteps 2</p><p>//First step: Imperfection force at the tip of the beam</p><p>Static 1</p><p>EndTime 1</p><p>TimeStep 0.1</p><p>MaxTimeStep 0.1</p><p>MinTimeStep 0.01</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.0</p><p>Sample 2</p><p>//Second step: Compression force applied at the tip of the beam</p><p>Static 2</p><p>EndTime 2</p><p>TimeStep 0.1</p><p>MaxTimeStep 0.1</p><p>MinTimeStep 0.01</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.0</p><p>Sample 2</p></td></tr></tbody></table>

### Creation of loads

Two loads are created. The first is the force applied at the tip of the beam to introduce an imperfection and induce buckling. The second load is the compression force. It is important to mention that the compression force used is a bit larger than the Euler’s critical load. That one is calculated by the following expression:

$$F_{\text{cr}} = \ \frac{\pi^{2}\text{EI}}{{(2l)}^{2}}$$

Where F<sub>cr</sub> is the Euler’s critical load, E is the Young’s Module, *I* is the cross section’s moment of inertia and *l* is the length of the beam. Using that expression, we obtain the value of 205.62 N.

<table><tbody><tr class="odd"><td><p>//Creation of loads</p><p>Loads 2</p><p>NodalLoad 1 NodeSet 2 CS 1 NTimes 2</p><p>//Time FX FY FZ MX MY MZ</p><p>0 0 0 0 0 0 0</p><p>1 0 1 0 0 0 0</p><p>NodalLoad 2 NodeSet 2 CS 1 NTimes 2</p><p>//Time FX FY FZ MX MY MZ</p><p>1 0 0 0 0 0 0</p><p>2 0 0 -206 0 0 0</p></td></tr></tbody></table>

### Creation of concomitant solution

A concomitant solution is created to run a modal analysis during the second solution step. In this case, we request Giraffe to evaluate the 10 smallest eigenvalues, associated with natural frequencies of the system.

<table><tbody><tr class="odd"><td><p>//Creation of concomitant solution</p><p>ConcomitantSolution</p><p>Sample 3 BoolTable 0 1</p><p>Modal NumberModes 10 Tolerance 1E-6</p></td></tr></tbody></table>

### Creation of constraints

A clamp boundary condition is created. This is done by fixing all the degrees of freedom during the entire simulation.

<table><tbody><tr class="odd"><td><p>//Creation of constraints</p><p>Constraints 1</p><p>NodalConstraint 1 NodeSet 1</p><p>UX BoolTable 1</p><p>UY BoolTable 1</p><p>UZ BoolTable 1</p><p>ROTX BoolTable 1</p><p>ROTY BoolTable 1</p><p>ROTZ BoolTable 1</p></td></tr></tbody></table>

### Creation of solver options

Here we request Giraffe to apply a single thread to solve this case. The linear system solution chosen is Direct.

<table><tbody><tr class="odd"><td><p>//Creation of solver options</p><p>SolverOptions</p><p>Processors 1 LinSys Direct</p></td></tr></tbody></table>

### Creation of monitors

We request Giraffe to monitor 2 nodes during the simulation: numbers 1 and 11.

<table><tbody><tr class="odd"><td><p>//Creation of monitors</p><p>Monitor Sample 1</p><p>MonitorNodes 1 51</p></td></tr></tbody></table>

### Creation of post-files

The options used here are the same as the tutorial 01, with an exception that is the “WriteForces”.

<table><tbody><tr class="odd"><td><p>//Creation of post files</p><p>PostFiles</p><p>MagFactor 1</p><p>WriteMesh 1</p><p>WriteRenderMesh 1</p><p>WriteRigidContactSurfaces 0</p><p>WriteFlexibleContactSurfaces 0</p><p>WriteForces 0</p><p>WriteConstraints 0</p><p>WriteSpecialConstraints 0</p><p>WriteContactForces 0</p><p>WriteRenderRigidBodies 0</p><p>WriteRenderParticles 0</p></td></tr></tbody></table>

## Results and post-processing

To post-process results in Paraview<sup>TM</sup>, first we have to open the files “whole\_solution\_mesh.pvd” and “whole\_solution\_rendermesh.pvd”, which are located at the folder “post” in “tutorial09”.

<img src="../../assets/tutorials/image69.png" style="width:5.90556in;height:3.16875in" alt="" />

Figure 9‑2 - Opening the files of tutorial 09 in Paraview<sup>TM</sup>.

Then, we can click the preview icon and the Apply button on the Pipeline browser. After doing that, we are ready to post-process the results.

<img src="../../assets/tutorials/image70.png" style="width:4.40625in;height:1.95908in" alt="" />

Figure 9‑3 - Deformed shape of the structure (render view).

Like mentioned in the previous tutorial, there is a folder called “concomitant\_solution”. Inside this folder, we find a text file with the same name. This file contains all the eigenvalues assigned to the vibration modes requested, evaluated over time. One may plot the time-evolution of the smallest eigenvalue real part, leading to Figure 9‑4.

<img src="../../assets/tutorials/emf_image71.png" style="max-width:100%" alt="" />

Figure 9‑4 - Smallest Eigenvalue vs time.

Note that the smallest eigenvalue approaches zero at the end of simulation, which means that the structure loses stiffness and thus experiences buckling. Besides, there is another way to check it. If we do not introduce an imperfection force, the same eigenvalue analysis may be done. In this case, the eigenvalue does cross the abscissa axis. The instant at which there is the abscissa crossing is related to the Euler’s critical load.

<img src="../../assets/tutorials/emf_image72.png" style="max-width:100%" alt="" />

Figure 9‑5 - Smallest Eigenvalue vs Time (without imperfection force).

Giraffe generates a monitor text file for each node requested containing the values of displacement and force in each global direction. Therefore, it is possible to obtain the vertical displacement of point B along the simulation evolution. Besides, we can vary the imperfection load in the input file, generating a collection of plots. Then, we can make an equilibrium configuration plot, shown in Figure 9‑6.

<img src="../../assets/tutorials/emf_image73.png" style="max-width:100%" alt="" />

Figure 9‑6 - Vertical Displacement UY vs Lateral Force FY.
