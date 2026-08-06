# Tutorial 7: Using several solution steps in a cantilever beam model

!!! example "Input files"
    [`tutorial07.inp`](../assets/tutorials/inputs/tutorial07/tutorial07.inp)
    (also in the repository under `inputs/tutorial07/`)

Input file name: tutorial07.inp

*Developed by Alfredo Gay Neto & Gabriel Maki*

<img src="../../assets/tutorials/image55.png" style="width:2.78163in;height:1.17708in" alt="" />(a) <img src="../../assets/tutorials/image56.png" style="width:2.58333in;height:1.11517in" alt="" /> (b)

Figure 7‑1 - (a) – Cantilever beam model (b) Beam mesh render view

## Objective

This tutorial presents a geometrically nonlinear simulation of a cantilever beam. The objective is to introduce the gravitational field in Giraffe and show how to perform a static analysis with many solution steps.

## Description

The beam AB (1 m of length) is meshed using 11 nodes and 5 elements, using the element Beam\_1. There is a clamp boundary condition at point A. Three loads are applied sequentially. First, a force F is applied at point B. Then, a moment M is applied at point B. The third load is the gravitational field, inducing the weight (distributed along the whole structure volume). Figure 7‑1 illustrates the problem. Three loads are applied during three steps. The force and the moment are applied during step 2 to 4 and the gravitational field is applied during step 3 to 5. On the sixth step, all the loads are removed.

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

This tutorial has six solution steps. The first solution step is created establishing the end time equal to “1”. The second solution step is created establishing the end time equal to “2” and the others follow this pattern. Each solution step has a time step of 0.1, that is, the loads are applied in increments of 10%.

<table><tbody><tr class="odd"><td><p>//Creation of the solution steps</p><p>SolutionSteps 6</p><p>//First step: FX applied at the tip of the beam</p><p>Static 1</p><p>EndTime 1</p><p>TimeStep 0.1</p><p>MaxTimeStep 0.1</p><p>MinTimeStep 0.01</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.0</p><p>Sample 2</p><p>//Second step: FX kept and MZ applied at the tip of the beam</p><p>Static 2</p><p>EndTime 2</p><p>TimeStep 0.1</p><p>MaxTimeStep 0.1</p><p>MinTimeStep 0.01</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.0</p><p>Sample 2</p><p>//Third step: FX and MZ kept and gravitational field is turned on</p><p>Static 3</p><p>EndTime 3</p><p>TimeStep 0.1</p><p>MaxTimeStep 0.1</p><p>MinTimeStep 0.01</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.0</p><p>Sample 2</p><p>//Foutrh step: All kept (FX, MZ and gravity)</p><p>Static 4</p><p>EndTime 4</p><p>TimeStep 0.1</p><p>MaxTimeStep 0.1</p><p>MinTimeStep 0.01</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.0</p><p>Sample 2</p><p>//Fifth step: Only gravity is kept and FX and MZ are removed</p><p>Static 5</p><p>EndTime 5</p><p>TimeStep 0.1</p><p>MaxTimeStep 0.1</p><p>MinTimeStep 0.01</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.0</p><p>Sample 2</p><p>//Sixth step: gravity is removed (no loads at the end)</p><p>Static 6</p><p>EndTime 6</p><p>TimeStep 0.1</p><p>MaxTimeStep 0.1</p><p>MinTimeStep 0.01</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.0</p><p>Sample 2</p></td></tr></tbody></table>

### Creation of loads

A single nodal load is created. Note that the load application time is associated with the start and end time of the solution steps in this example. This is not strictly necessary. The user may create time-varying loads, according to arbitrary choices.

<table><tbody><tr class="odd"><td><p>//Creation of loads</p><p>Loads 1</p><p>NodalLoad 1 NodeSet 2 CS 1 NTimes 5</p><p>//Time FX FY FZ MX MY MZ</p><p>0 0 0 0 0 0 0</p><p>1 1000 0 0 0 0 0</p><p>2 1000 0 0 0 0 250</p><p>4 1000 0 0 0 0 250</p><p>5 0 0 0 0 0 0</p></td></tr></tbody></table>

### Creation of constraints

A clamp boundary condition is created. This is done by fixing all the degrees of freedom during the entire simulation. It is important to remember that, in BoolTable, there is a single “1”, but it means that the constraint is applied for all the steps.

<table><tbody><tr class="odd"><td><p>//Creation of constraints</p><p>Constraints 1</p><p>NodalConstraint 1 NodeSet 1</p><p>UX BoolTable 1</p><p>UY BoolTable 1</p><p>UZ BoolTable 1</p><p>ROTX BoolTable 1</p><p>ROTY BoolTable 1</p><p>ROTZ BoolTable 1</p></td></tr></tbody></table>

### Creation of Environment Data

An environment data is created to introduce the gravity field on the model. Note that it is not kept during the entire simulation. The BoolTable has the objective of establishing in which solution steps such load is going to be considered. In our example, the gravitational field is linearly imposed along solution step 3, but is linearly turned off along solution step 6.

<table><tbody><tr class="odd"><td><p>//Creation of Environment data</p><p>Environment</p><p>//Gravity data</p><p>GravityData</p><p>G 9.81 0 0 BoolTable 0 0 1 1 1 0</p></td></tr></tbody></table>

### Creation of solver options

Here we request Giraffe to apply a single thread to solve this case. The linear system solution chosen is Direct.

<table><tbody><tr class="odd"><td><p>//Creation of solver options</p><p>SolverOptions</p><p>Processors 1 LinSys Direct</p></td></tr></tbody></table>

### Creation of monitors

We request Giraffe to monitor 2 nodes during the simulation: numbers 1 and 11.

<table><tbody><tr class="odd"><td><p>//Creation of monitors</p><p>Monitor Sample 1</p><p>MonitorNodes 1 11</p></td></tr></tbody></table>

### Creation of post-files

The options used here are the same as the tutorial 01, where the WriteMesh, WriteRenderMesh and WriteForces are chosen for post-processing using Paraview<sup>TM</sup> and the MagFactor parameters assume a unit value so the deformed shape of Paraview<sup>TM</sup> frames will reflect the real scale.

<table><tbody><tr class="odd"><td><p>//Cration of post files</p><p>PostFiles</p><p>MagFactor 1</p><p>WriteMesh 1</p><p>WriteRenderMesh 1</p><p>WriteRigidContactSurfaces 0</p><p>WriteFlexibleContactSurfaces 0</p><p>WriteForces 1</p><p>WriteConstraints 0</p><p>WriteSpecialConstraints 0</p><p>WriteContactForces 0</p><p>WriteRenderRigidBodies 0</p><p>WriteRenderParticles 0</p></td></tr></tbody></table>

## Results and post-processing

To post-process data using Paraview<sup>TM</sup>, first we must open the files “whole\_solution\_forces.pvd”, “whole\_solution\_mesh.pvd” and “whole\_solution\_rendermesh.pvd”, which are located at the folder “post” in “tutorial07”.

<img src="../../assets/tutorials/image57.png" style="width:5.48712in;height:3.13844in" alt="" />

Figure 7‑2 - Opening the files in Paraview<sup>TM</sup>.

<img src="../../assets/tutorials/image58.png" style="width:1.98936in;height:3.14184in" alt="" />

Figure 7‑3 - Showing the Pipeline Browser.

For a more convenient view with no distorting dimensions, one may adopt an isometric perspective. For that, go to Properties and check “Camera Parallel Projection in Orientation Axes.

Next, we can create a Glyph to see the symbol of the load being applied. For that, go to the Pipeline Browser and click in “whole\_solution\_forces.pvd”. Then, click the Glyph twice.

<img src="../../assets/tutorials/image59.png" style="width:2.5797in;height:1.86242in" alt="" />

Figure 7‑4 - Creating the two Glyphs.

After creating the Glyphs, we have to adjust the scale of vectors. To do this, go to the Pipeline Browser and, first, click the preview icon of Glyph. Then, in Scaling, set the scale mode to vector, change the scale factor to a small number like 0.001 and click Apply. After doing this for both Glyphs, in Glyph 2, we must change, in Active Attributes, the Vectors to Moment and click Apply.

<img src="../../assets/tutorials/image60.png" style="width:2.1875in;height:3.22245in" alt="" />

Figure 7‑5 - Changing the Properties of Glyph 2.

Finally, we can show an animation of the simulation results. If the user wants to know the values of displacements of some nodes or the reaction forces on the clamp in a certain instant, we can gather such information in the “monitors” folder, since they are requested prior to simulation solving.

<img src="../../assets/tutorials/image61.png" style="width:3.40625in;height:2.84158in" alt="" />

Figure 7‑6 - Deformed shape at the end of third solution step.

In Figure 7‑6 it is possible to see that the background color has changed. To do this, click the tool indicated in Figure 7‑7 and choose the color.

<img src="../../assets/tutorials/image62.png" style="width:2.44426in;height:2.05208in" alt="" />

Figure 7‑7 - Changing the background color.

<img src="../../assets/tutorials/emf_image63.png" style="max-width:100%" alt="" />

Figure 7‑8 – Absolute Moment Reaction \|MY\| vs Time step

An interesting result can be seen in Figure 7‑8. It shows the importance of geometric nonlinearity. Using the reactive moment in direction Y as example, if we calculate this reaction considering the undisturbed structure, it results in 1098 Nm. If we consider the nonlinearity and calculate numerically, it results in 438,4 Nm as shown in Figure 7‑8.
