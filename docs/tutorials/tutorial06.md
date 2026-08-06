# Tutorial 6: Static and Dynamic analysis of a cantilever plate with a rigid edge

!!! example "Input files"
    [`tutorial06.inp`](../assets/tutorials/inputs/tutorial06/tutorial06.inp)
    (also in the repository under `inputs/tutorial06/`)

Input file name: tutorial06.inp

*Developed by Alfredo Gay Neto & Gabriel Maki*

<img src="../../assets/tutorials/image45.png" style="width:2.92435in;height:1.5in" alt="" /> (a) <img src="../../assets/tutorials/image46.png" style="width:1.3125in;height:2.3851in" alt="" /> (b)

Figure 6‑1 - (a) Model of cantilever plate rigidly linked to a node (b) Mesh used for modeling with the rigid node set represented

## Objective

The objective of this tutorial is to introduce a new resource for modeling called special constraint. In this case, the structure is exactly the same plate simulated in last tutorial. However, now we establish a rigid region along the edge of the plate by means of a rigid node set, a special constraint available in Giraffe.

## Description

The plate size is 0.395 m long and 0.195 m wide. It is meshed considering 862 nodes (one of which is called pilot node to use the special constraint) and 400 elements (Shell\_1). There is a clamp boundary condition on one of the edges of the plate. The section has a thickness of 2 mm. As in previous example, the simulation runs in two steps. During the first, a nodal load is applied on the pilot node, which deforms the structure. During the second step, which is dynamic, the load is kept for 0.25 s, then it starts to reduce until zero along more 0.25 s.

## Input file step-by-step

### Creation of nodes

A total of 862 nodes is created to establish a finite element mesh. Like mentioned before, there is one node that is called pilot node and it is used to define the special constraint.

<table><tbody><tr class="odd"><td><p>//Creation of nodes</p><p>Nodes 862</p><p>Node 1 0 0 0</p><p>Node 2 0.0195 0 0</p><p>Node 3 0.039000001 0 0</p><p>Node 4 0.058499999 0 0</p><p>Node 5 0.078000002 0 0</p><p>Node 6 0.097499996 0 0</p><p>Node 7 0.116999999 0 0</p><p>Node 8 0.136500001 0 0</p><p>…</p><p>Node 855 0.14625001 0.395000011 0</p><p>Node 856 0.175500005 0.385125011 0</p><p>Node 857 0.165749997 0.385125011 0</p><p>Node 858 0.165749997 0.395000011 0</p><p>Node 859 0.194999993 0.385125011 0</p><p>Node 860 0.185249999 0.385125011 0</p><p>Node 861 0.185249999 0.395000011 0</p><p>//Node 862 - will be the pilot node of the rigid node set</p><p>Node 862 0.1 0.6 0.0</p></td></tr></tbody></table>

### Creation of node sets

Three node sets are created to establish the constraints, special constraints and the load applied in the pilot node. The node sets 1 and 2 are the same as the previously tutorial and the node set 3 contains only the pilot node.

<table><tbody><tr class="odd"><td><p>//Creation of node sets</p><p>NodeSets 3</p><p>//Left tip</p><p>NodeSet 1 Nodes 21 List 1 2 3 4 5 6 7 8 9 10 11 232 237 241 245 249 253 257 261 265 269</p><p>//Right tip</p><p>NodeSet 2 Nodes 21 List 221 222 223 224 225 226 227 228 229 230 231 833 837 840 843 846 849 852 855 858 861</p><p>//Pilot node</p><p>NodeSet 3 Nodes 1 List 862</p></td></tr></tbody></table>

### Creation of elements

400 elements are created to compose the mesh. All the elements have the same material properties and the same shell section.

<table><tbody><tr class="odd"><td><p>//Creation of elements</p><p>//Elements - in this case only Shell_1 elements are created</p><p>//Important: the 6 nodes of each shell element MUST lie in a plane at reference configuration (input file)</p><p>//After deformation, shells may deform to out-of-plane configurations.</p><p>//Each element has attributes to material (Mat) and section (Sec).</p><p>Elements 400</p><p>Shell_1 1 Mat 1 Sec 1 Nodes 1 2 13 232 233 234</p><p>Shell_1 2 Mat 1 Sec 1 Nodes 13 12 1 235 236 234</p><p>Shell_1 3 Mat 1 Sec 1 Nodes 2 3 14 237 238 239</p><p>Shell_1 4 Mat 1 Sec 1 Nodes 14 13 2 240 233 239</p><p>Shell_1 5 Mat 1 Sec 1 Nodes 3 4 15 241 242 243</p><p>…</p><p>Shell_1 395 Mat 1 Sec 1 Nodes 217 218 229 824 853 854</p><p>Shell_1 396 Mat 1 Sec 1 Nodes 229 228 217 855 850 854</p><p>Shell_1 397 Mat 1 Sec 1 Nodes 218 219 230 827 856 857</p><p>Shell_1 398 Mat 1 Sec 1 Nodes 230 229 218 858 853 857</p><p>Shell_1 399 Mat 1 Sec 1 Nodes 219 220 231 830 859 860</p><p>Shell_1 400 Mat 1 Sec 1 Nodes 231 230 219 861 856 860</p></td></tr></tbody></table>

### Creation of materials

One only material model is defined, which identification number is 1.

<table><tbody><tr class="odd"><td><p>//Creation of Materials</p><p>Materials 1</p><p>Hooke 1 E 200e9 Nu 0.3 Rho 8000</p></td></tr></tbody></table>

### Creation of shell sections

One shell section is defined, which has a thickness of 2 mm and is homogeneous. Its identification number is 1.

<table><tbody><tr class="odd"><td><p>//Creation of ShellSections</p><p>ShellSections 1</p><p>Homogeneous 1 Thickness 0.002</p></td></tr></tbody></table>

### Creation of coordinate systems

One only coordinate system is created, which is used to define the direction of nodal load.

<table><tbody><tr class="odd"><td><p>//Creation of Coordinate Systems</p><p>CoordinateSystems 1</p><p>CS 1 E1 1 0 0 E3 0 0 1</p></td></tr></tbody></table>

### Creation of solution steps

This tutorial has two solutions steps. The first is static and is created establishing the end time equal to “1”. The time step is “0.05”, which means that the load is divided into increments of 5%. The second step is dynamic and is created establishing the end time equal to “2.5”.

<table><tbody><tr class="odd"><td><p>//Creation of the solution steps</p><p>SolutionSteps 2</p><p>//First step: statics - from time 0 - 1</p><p>Static 1</p><p>EndTime 1.0</p><p>TimeStep 0.05</p><p>MaxTimeStep 0.05</p><p>MinTimeStep 0.01</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.0</p><p>Sample 1</p><p>//Second step: dynamics - from time 1 - 2.5</p><p>Dynamic 2</p><p>EndTime 2.5</p><p>TimeStep 0.01</p><p>MaxTimeStep 0.01</p><p>MinTimeStep 0.0001</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 5</p><p>IncFactor 1.2</p><p>Sample 1</p><p>RayleighDamping Alpha 0 Beta 0 Update 0</p><p>NewmarkCoefficients Beta 0.3 Gamma 0.5</p></td></tr></tbody></table>

### Creation of loads

A single nodal load is defined. This is done establishing a time varying moment in direction X and Y. Note that, after 1.25 s, the moments start to reduce until be totally removed at 1.5 s. After that time, even if the user does not input data in the time-series table, Giraffe considers the last established value kept.

<table><tbody><tr class="odd"><td><p>//Creation of loads</p><p>Loads 1</p><p>NodalLoad 1 NodeSet 3 CS 1 NTimes 4</p><p>//Time FX FY FZ MX MY MZ</p><p>0.00 0 0 0 0 0 0</p><p>1.00 0 0 0 100 300 0</p><p>1.25 0 0 0 100 300 0</p><p>1.50 0 0 0 0 0 0</p></td></tr></tbody></table>

### Creation of constraints

A single nodal constraint is created. This is done fixing all the degrees of freedom, representing a clamp boundary condition on node set 1.

<table><tbody><tr class="odd"><td><p>//Creation of constraints</p><p>Constraints 1</p><p>//Clamping the left tip</p><p>NodalConstraint 1 NodeSet 1</p><p>UX BoolTable 1</p><p>UY BoolTable 1</p><p>UZ BoolTable 1</p><p>ROTX BoolTable 1</p><p>ROTY BoolTable 1</p><p>ROTZ BoolTable 1</p></td></tr></tbody></table>

### Creation of special constraints

The special constraint used in this tutorial is the Rigid Node Set, which defines a rigid region on a specific node set, in other words, the nodes contained in that node set behaves kinematically as a group of material points within a rigid body. Besides, that node set is rigidly linked to the pilot node as if there were a rigid bar between them. Therefore, the loads applied on the pilot node are transferred to the structure. It is important to mention that, when there is a single “1” on the BoolTable and there are more than one solution steps, it means that the special constraint is activated during all the remaining solution steps.

<table><tbody><tr class="odd"><td><p>//Definition of a rigid region (rigid node set) on node set 2. This region is ruled by the DOFs of the pilot node</p><p>SpecialConstraints 1</p><p>RigidNodeSet 1 PilotNode 862 NodeSet 2 BoolTable 1</p></td></tr></tbody></table>

### Creation of solver options

Here it is possible to choose the number of processors used (parallel processing). Besides, we must define the solver type for linear equation systems. In this case, four processors are employed and the direct type of linear system solver is our choice.

<table><tbody><tr class="odd"><td><p>//Creation of solver options</p><p>SolverOptions</p><p>Processors 4 LinSys Direct</p></td></tr></tbody></table>

### Creation of monitors

We request Giraffe to monitor two node sets during the simulation.

<table><tbody><tr class="odd"><td><p>//Creation of monitors</p><p>Monitor Sample 1</p><p>MonitorNodeSets 1 2</p></td></tr></tbody></table>

### Creation of post-files

The options used here are the same as the tutorial 01, where the WriteMesh, WriteRenderMesh and WriteForces are chosen for post-processing using Paraview<sup>TM</sup> and the MagFactor parameters assume a unit value so the deformed shape of Paraview<sup>TM</sup> frames will reflect the real scale. It is interesting to choose the “WriteSpecialConstraints” because it allows the user to view the Special Constraints in Paraview<sup>TM</sup>.

<table><tbody><tr class="odd"><td><p>//Creation of post files</p><p>PostFiles</p><p>MagFactor 1</p><p>WriteMesh 1</p><p>WriteRenderMesh 1</p><p>WriteRigidContactSurfaces 0</p><p>WriteFlexibleContactSurfaces 0</p><p>WriteForces 1</p><p>WriteConstraints 0</p><p>WriteSpecialConstraints 1</p><p>WriteContactForces 0</p><p>WriteRenderRigidBodies 0</p><p>WriteRenderParticles 0</p></td></tr></tbody></table>

### Creation of convergence criteria

Here, we can define some convergence tolerances. This is an optional task (for advanced users only).

<table><tbody><tr class="odd"><td><p>ConvergenceCriteria</p><p>ForceTolerance 1e-4</p><p>MomentTolerance 1e-4</p><p>ForceMinimumReference 1e-4</p><p>MomentMinimumReference 1e-4</p><p>ConstraintMinimumReference 1e-7</p><p>DisplacementTolerance 1e-4</p><p>RotationTolerance 1e-4</p><p>LagrangeTolerance 1e-4</p><p>DisplacementMinimumReference 1e-6</p><p>RotationMinimumReference 1e-6</p><p>LagrangeMinimumReference 1e-6</p><p>DivergenceReference 1e+15</p></td></tr></tbody></table>

## Results and post-processing

To view the post-processing in Paraview<sup>TM</sup>, first we must open the files “whole\_solution\_forces.pvd” and “whole\_solution\_mesh.pvd”, which are located at the folder “post” in “tutorial06”.

<img src="../../assets/tutorials/image47.png" style="width:5.90556in;height:3.38264in" alt="" />

Figure 6‑2 - Opening results of tutorial 6 in Paraview<sup>TM</sup>.

Then, we must click the preview icon and the Apply button on the Pipeline browser.

For a better view, we must use the isometric perspective view. To do this, go to Properties and check “Camera Parallel Projection in Orientation Axes.

<img src="../../assets/tutorials/image48.png" style="width:2.42705in;height:3.22917in" alt="" />

Figure 6‑3 - Enabling isometric perspective view in Paraview<sup>TM</sup>.

Next, we must create a Glyph for seeing the load being applied. For that, go to the Pipeline Browser and click in “whole\_solution\_forces.pvd”. Then, click the Glyph.

After creating the Glyph, we must adjust the scale of vectors. To do this, go to the Pipeline Browser and, first, click the preview icon of Glyph. Then, in Scaling, set the scale mode to vector and change the scale factor to a small number like 0.0005. Besides, we need to change, in Active Attributes, the Vectors to Moment. Finally, click Apply.

<img src="../../assets/tutorials/image49.png" style="width:2.29721in;height:1.89314in" alt="" />

Figure 6‑4 - Changing Properties of Glyph.

If the user wants to know the values of displacements of some nodes, node sets or the reaction forces on the clamp at a given instant, monitors resource may be employed.

<img src="../../assets/tutorials/image50.png" style="width:4.16667in;height:4.65833in" alt="" />

Figure 6‑5 - Deformed shape of structure after the first solution step.

Note that the moment is not applied on the plate. It is applied on the pilot node. Besides, it is possible to see that the rigid edge does not deform because of the special constraint imposed which is represented by the white lines in Figure 6‑5.

Obs: If the user wants to change the color of the arrow, the one can do it by clicking in Edit, on Coloring at Properties.

<img src="../../assets/tutorials/image51.png" style="width:3.27083in;height:2.72917in" alt="" />

Figure 6‑6 - Changing the color of arrows.

Besides, if the user wants to see the Special Constraint on a better way, it is possible to increase the point size. To do this, we have to go to the Properties of “whole\_solution\_symbols.pvd” and change the Point Size.

<img src="../../assets/tutorials/image52.png" style="width:3.375in;height:1.22917in" alt="" />

Figure 6‑7 - Changing the Point Size.

By monitoring the vertical displacement of the structure tip (node set 2), the user can note that it does not move along the first 0.25 s of dynamic solution step because there is no change in the load. When the moments decrease their values, the structure starts to oscillate like observed in the last tutorial. It is possible to see such results in Figure 6‑8 and Figure 6‑9.

<img src="../../assets/tutorials/emf_image53.png" style="max-width:100%" alt="" />

Figure 6‑8 – Vertical Displacement UZ vs Time.

<img src="../../assets/tutorials/emf_image54.png" style="max-width:100%" alt="" />

Figure 6‑9 - Moment Reaction MX vs time.
