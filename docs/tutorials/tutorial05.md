# Tutorial 5: Static and Dynamic analysis of a cantilever plate

!!! example "Input files"
    [`tutorial05.inp`](../assets/tutorials/inputs/tutorial05/tutorial05.inp)
    (also in the repository under `inputs/tutorial05/`)

Input file name: tutorial05.inp

*Developed by Alfredo Gay Neto & Gabriel Maki*

<img src="../../assets/tutorials/image32.png" style="width:3in;height:1.65275in" alt="" /> (a) <img src="../../assets/tutorials/image33.png" style="width:1.32292in;height:1.65505in" alt="" /> (b)

Figure 5‑1 - (a) Model of a cantilever plate (b) Mesh used for modeling

## Objective

This tutorial presents a geometrically-nonlinear simulation of a cantilever plate, subjected to a uniformly distributed load along a line and clamped in the opposite side, as shown in Figure 5‑1(a). The tutorial has the purpose of showing how to solve an analysis that has two steps, one static and one dynamic.

## Description

The plate is 0.395 m long and 0.195 m wide. It is meshed considering 861 nodes and 400 elements, using the element Shell\_1. There is a clamp boundary condition on one of the edges of the plate. The section has a thickness of 2 mm. Besides, the simulation is solved in two steps. In the first, the load is applied on an edge of the plate. This deforms the whole structure, inducing bending. In the second step, which is dynamic, the nodal load is removed. It is possible to observe how the plate oscillates after dynamically removing the load.

The main difference between static and dynamic solution steps is that dynamics considers not only the internal and external loads on the structure, but also inertial and damping loads. This may be useful for determining a time-series response of a given structure to dynamic phenomena, such as time-varying forces, impacts and other possibilities.

## Input file step-by-step

### Creation of nodes

A total of 861 nodes are created to establish a finite element mesh.

<table><tbody><tr class="odd"><td><p>//Creation of nodes</p><p>Nodes 861</p><p>Node 1 0 0 0</p><p>Node 2 0.0195 0 0</p><p>Node 3 0.039000001 0 0</p><p>Node 4 0.058499999 0 0</p><p>Node 5 0.078000002 0 0</p><p>Node 6 0.097499996 0 0</p><p>Node 7 0.116999999 0 0</p><p>Node 8 0.136500001 0 0</p><p>Node 9 0.156000003 0 0</p><p>Node 10 0.175500005 0 0</p><p>Node 11 0.194999993 0 0</p><p>...</p><p>Node 850 0.136500001 0.385125011 0</p><p>Node 851 0.126749992 0.385125011 0</p><p>Node 852 0.126749992 0.395000011 0</p><p>Node 853 0.156000003 0.385125011 0</p><p>Node 854 0.14625001 0.385125011 0</p><p>Node 855 0.14625001 0.395000011 0</p><p>Node 856 0.175500005 0.385125011 0</p><p>Node 857 0.165749997 0.385125011 0</p><p>Node 858 0.165749997 0.395000011 0</p><p>Node 859 0.194999993 0.385125011 0</p><p>Node 860 0.185249999 0.385125011 0</p><p>Node 861 0.185249999 0.395000011 0</p></td></tr></tbody></table>

### Creation of node sets

Two node sets are defined to establish the constraints and the load on the structure. The node sets are represented in Figure 5‑2 by red lines. The nodes of node set 1 belong to the clamped line.

<img src="../../assets/tutorials/image34.png" style="width:4.19792in;height:2.48302in" alt="" />

Figure 5‑2 - Representation of node sets on the plate.

<table><tbody><tr class="odd"><td><p>//Creation of node sets</p><p>NodeSets 2</p><p>//Left tip</p><p>NodeSet 1 Nodes 21 List 1 2 3 4 5 6 7 8 9 10 11 232 237 241 245 249 253 257 261 265 269</p><p>//Right tip</p><p>NodeSet 2 Nodes 21 List 221 222 223 224 225 226 227 228 229 230 231 833 837 840 843 846 849 852 855 858 861</p></td></tr></tbody></table>

### Creation of elements

A total of 400 elements is created. All the elements have the same material properties and the same shell section.

<table><tbody><tr class="odd"><td><p>//Creation of elements</p><p>//Elements - in this case only Shell_1 elements are created</p><p>//Important: the 6 nodes of each shell element MUST lie in a plane at reference configuration (input file)</p><p>//After deformation, shells may deform to out-of-plane configurations.</p><p>//Each element has attributes to material (Mat) and section (Sec).</p><p>Elements 400</p><p>Shell_1 1 Mat 1 Sec 1 Nodes 1 2 13 232 233 234</p><p>Shell_1 2 Mat 1 Sec 1 Nodes 13 12 1 235 236 234</p><p>Shell_1 3 Mat 1 Sec 1 Nodes 2 3 14 237 238 239</p><p>Shell_1 4 Mat 1 Sec 1 Nodes 14 13 2 240 233 239</p><p>Shell_1 5 Mat 1 Sec 1 Nodes 3 4 15 241 242 243</p><p>…</p><p>Shell_1 395 Mat 1 Sec 1 Nodes 217 218 229 824 853 854</p><p>Shell_1 396 Mat 1 Sec 1 Nodes 229 228 217 855 850 854</p><p>Shell_1 397 Mat 1 Sec 1 Nodes 218 219 230 827 856 857</p><p>Shell_1 398 Mat 1 Sec 1 Nodes 230 229 218 858 853 857</p><p>Shell_1 399 Mat 1 Sec 1 Nodes 219 220 231 830 859 860</p><p>Shell_1 400 Mat 1 Sec 1 Nodes 231 230 219 861 856 860</p></td></tr></tbody></table>

### Creation of Materials

One only material model is defined, which identification number is 1.

<table><tbody><tr class="odd"><td><p>//Creation of Materials</p><p>Materials 1</p><p>Hooke 1 E 200e9 Nu 0.3 Rho 8000</p></td></tr></tbody></table>

### Creation of shell sections

One shell section is defined, homogeneous and with 2 mm thickness. Its identification number is 1.

<table><tbody><tr class="odd"><td><p>//Creation of ShellSections</p><p>ShellSections 1</p><p>Homogeneous 1 Thickness 0.002</p></td></tr></tbody></table>

### Creation of Coordinate Systems

One only coordinate system is created, which is used to define the direction of the nodal load.

<table><tbody><tr class="odd"><td><p>//Creation of Coordinate Systems</p><p>CoordinateSystems 1</p><p>CS 1 E1 1 0 0 E3 0 0 1</p></td></tr></tbody></table>

### Creation of solution steps

This tutorial has two solutions steps. The first is static and is created establishing the end time equal to “1”. The time step is set to “0.1”, which means that the load is divided into increments of 10% along the first solution step. The second solution step is dynamic and is created establishing the end time equal to “1.4”. Thus, it encompasses simulation from time 1 to 1.4 (from the end of previous solution step until the end time defined for it). In a dynamic analysis, the time step is the physical time, unlike in the static analysis, where time is a tracking scalar parameter, serving only as a guide to define varying loads along the simulation evolution.

<table><tbody><tr class="odd"><td><p>SolutionSteps 2</p><p>//First step: statics - from time 0 - 1</p><p>Static 1</p><p>EndTime 1.0</p><p>TimeStep 0.1</p><p>MaxTimeStep 0.1</p><p>MinTimeStep 0.01</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.0</p><p>Sample 1</p><p>//Second step: dynamics - from time 1 - 1.4</p><p>Dynamic 2</p><p>EndTime 1.4</p><p>TimeStep 0.005</p><p>MaxTimeStep 0.005</p><p>MinTimeStep 0.001</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 5</p><p>IncFactor 1.2</p><p>Sample 1</p><p>RayleighDamping Alpha 0 Beta 0 Update 0</p><p>NewmarkCoefficients Beta 0.3 Gamma 0.5</p></td></tr></tbody></table>

Note that here we set MinTimeStep as 0.001. This permits Giraffe to automatically decrease the time-step in case of difficulties of convergence along the simulation.

The last two lines of the definition of a dynamic solution step refers to damping and time-integration scheme setting. Rayleigh damping permits the user to set Alpha and Beta (multipliers of mass and stiffness matrices of a Rayleigh damping model). Also, it is possible to Update or not such matrices at the beginning of each time-step. In case of the option Update 0, the reference configuration matrices are considered for damping evaluation along the whole simulation. In present example no damping is considered, since both Alpha and Beta parameters are set zero.

Lastly, it is important to mention that Giraffe employs an implicit time-integration scheme of equations of motion. Newmark method is used, in a special version of it, when dealing with rotations (see e.g.: \[2\]). The values of Beta and Gamma are coefficients to rule the time-integration. The values Beta 0.3 and Gamma 0.5 are recommended to avoid numerical damping.

### Creation of loads

A single nodal load is defined. This is done establishing a time varying force in direction Z. Note that, after 1 s, the nodal load is removed very quickly (in 0.0001 s). Giraffe employs linear interpolation between the input tabled values of loads along time.

<table><tbody><tr class="odd"><td><p>//Creation of loads</p><p>Loads 1</p><p>NodalLoad 1 NodeSet 2 CS 1 NTimes 3</p><p>//Time FX FY FZ MX MY MZ</p><p>0.00 0 0 0 0 0 0</p><p>1.00 0 0 1000 0 0 0</p><p>1.0001 0 0 0 0 0 0</p></td></tr></tbody></table>

### Creation of constraints

A single nodal constraint is created. This is done fixing all the degrees of freedom, representing a clamp boundary condition at node set 1.

<table><tbody><tr class="odd"><td><p>//Creation of constraints</p><p>Constraints 1</p><p>//Clamping the left tip</p><p>NodalConstraint 1 NodeSet 1</p><p>UX BoolTable 1</p><p>UY BoolTable 1</p><p>UZ BoolTable 1</p><p>ROTX BoolTable 1</p><p>ROTY BoolTable 1</p><p>ROTZ BoolTable 1</p></td></tr></tbody></table>

### Creation of solver options

Here it is possible to choose the number of processors used (parallel processing). Besides, we must define the solver type for linear equation systems. In this case, four processors are employed and the direct type of linear system solver is our choice.

<table><tbody><tr class="odd"><td><p>//Creation of solver options</p><p>SolverOptions</p><p>Processors 4 LinSys Direct</p></td></tr></tbody></table>

### Creation of monitors

We request Giraffe to monitor the two node sets during the simulation. This is useful for further post-processing.

<table><tbody><tr class="odd"><td><p>//Creation of monitors</p><p>Monitor Sample 1</p><p>MonitorNodeSets 1 2</p></td></tr></tbody></table>

### Creation of post-files

The options used here are the same as the tutorial 01, where the WriteMesh, WriteRenderMesh and WriteForces are chosen for post-processing using Paraview<sup>TM</sup> and the MagFactor parameters assume a unit value so the deformed shape of Paraview<sup>TM</sup> frames will reflect the real scale.

<table><tbody><tr class="odd"><td><p>//Creation of post files</p><p>PostFiles</p><p>MagFactor 1</p><p>WriteMesh 1</p><p>WriteRenderMesh 1</p><p>WriteRigidContactSurfaces 0</p><p>WriteFlexibleContactSurfaces 0</p><p>WriteForces 1</p><p>WriteConstraints 0</p><p>WriteSpecialConstraints 0</p><p>WriteContactForces 0</p><p>WriteRenderRigidBodies 0</p><p>WriteRenderParticles 0</p></td></tr></tbody></table>

## Results and post-processing

To view the post-processing in Paraview<sup>TM</sup>, first we must open the files “whole\_solution\_forces.pvd” and “whole\_solution\_mesh.pvd”, which are located at the folder “post” in “tutorial05”.

<img src="../../assets/tutorials/image35.png" style="width:3.94776in;height:2.27238in" alt="" />

Figure 5‑3 - Opening the files of tutorial 05 in Paraview<sup>TM</sup>.

Then, we must click on the preview icon and the Apply button on the Pipeline browser.

<img src="../../assets/tutorials/image36.png" style="width:1.49254in;height:2.33582in" alt="" />

Figure 5‑4 - Where to click in Pipeline Browser.

Figure 5‑5 shows how Paraview<sup>TM</sup> displays the results.

<img src="../../assets/tutorials/image37.png" style="width:3.29851in;height:2.18569in" alt="" />

Figure 5‑5 - Results using Paraview<sup>TM</sup>.

It is possible to change the view by using the tools of “Set view”. Figure 5‑6 depicts the icons that may be used for that.

<img src="../../assets/tutorials/image38.png" style="width:5.02083in;height:2.83396in" alt="" />

Figure 5‑6 - Changing the view in Paraview<sup>TM</sup>.

In addition, to changing the view, it is possible to establish an isometric perspective, instead of a conic one (default in Paraview).For that, we must click on any of Pipeline browser items, go to Properties and check “Camera Parallel Projection in Orientation Axes” option (Figure 5‑7).

<img src="../../assets/tutorials/image39.png" style="width:2.4375in;height:3.29916in" alt="" />

Figure 5‑7 - Establishing isometric perspective view.

Next, we can create a Glyph for visualization of the load being applied. For that, go to the Pipeline Browser and click in “whole\_solution\_forces.pvd”. Then, click the Glyph. The Figure 5‑8 shows where to click.

<img src="../../assets/tutorials/image22.png" style="width:2.27083in;height:2.875in" alt="" />

Figure 5‑8 - Creating a Glyph.

After creating the Glyph, we have to adjust the scale of vectors. For that, go to the Pipeline Browser and, first, click the preview icon of Glyph. Then, in Scaling, set the scale mode to vector and change the scale factor to a small number like 0.001. Finally, click Apply.

<img src="../../assets/tutorials/image40.png" style="width:2in;height:3.79683in" alt="" />

Figure 5‑9 - Adjusting the Scaling in Glyph.

Finally, we can perform an animation to observe how the structure deforms along time. Note that this post-processing may be done individually for each solution step, just opening in Paraview, instead of the whole solution results, the desired step results.

When needing to plot time series of displacements, forces or other kinematical quantities, monitors may be helpful.

<img src="../../assets/tutorials/image41.png" style="width:4.30208in;height:2.85777in" alt="" />

Figure 5‑10 - Deformed shape of structure after first solution step.

<img src="../../assets/tutorials/image42.png" style="width:4.30208in;height:2.85878in" alt="" />

Figure 5‑11 – Deformed structure at the end of the simulation.

Next plots refer to node set monitors post-processing. One may find kinematic and load quantities in such monitors.

<img src="../../assets/tutorials/emf_image43.png" style="max-width:100%" alt="" />

Figure 5‑12 - Vertical Displacement (UZ) of the plate tip vs. time.

<img src="../../assets/tutorials/emf_image44.png" style="max-width:100%" alt="" />

Figure 5‑13 - Vertical Reaction Force FZ at the clamp vs. time.

Figure 5‑12 shows the oscillatory movement that the structure presents. It is only possible to observe such behavior because the second solution step considers dynamic effects. Furthermore, Figure 5‑13 shows the reaction on the plate root vs. time, which is also a useful result.

If we change the second solution step to static, we will see that the plate returns to its initial configuration with no oscillations. To do so, we have to change the second solution step type to “static”. Furthermore, it is convenient to decrease the time-step of the simulation to provide Giraffe the possibility of creating a smoother path between the bended plate configuration (the end of the first solution step) and the expected end configuration of the second solution step, corresponding to the plate stress-free configuration. After making that changes and post-processing the model in Paraview<sup>TM</sup>, note that the initial and final configurations are the exactly the same.
