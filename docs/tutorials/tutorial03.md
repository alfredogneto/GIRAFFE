# Tutorial 3: Buckling of a pipe triggered by a compression follower load

!!! example "Input files"
    [`tutorial03.inp`](../assets/tutorials/inputs/tutorial03/tutorial03.inp)
    (also in the repository under `inputs/tutorial03/`)

Input file name: tutorial03.inp

*Developed by Alfredo Gay Neto & Gabriel Maki*

<img src="../../assets/tutorials/image17.png" style="width:2.76875in;height:1.26875in" alt="" /> (a) <img src="../../assets/tutorials/image18.png" style="width:0.96875in;height:1.42708in" alt="" />(b)

Figure 3‑1 – (a) Model of a pipe experiencing buckling (b) Cross section of the pipe

##  Objective

This tutorial presents a geometrically-nonlinear simulation of a pipe. The aim is to evaluate the buckling of the structure, when subjected to a compression load. For that, it is necessary to introduce an imperfection on the pipe and the compression (trigger) load, as shown in Figure 3‑1(a).

## Description

The pipe AB (100 m length) is meshed considering 101 equally-spaced nodes and 50 elements, (Pipe\_1). Two nodal loads are applied sequentially. First, at the middle-spam of the pipe, a lateral load is applied. After, at point B, a compressive nodal load is applied. The force in the middle-spam is necessary because it is not possible to perform a nonlinear buckling analysis without an initial imperfection on the structure. This imperfection induces certain buckling mode. It could be considered also a geometric imperfection directly in the model mesh. In our case, however, the lateral load introduces some curvature to the pipe that leads to subsequent buckling, when the structure experiences the compression load. The cross section of the pipe is shown in Figure 3‑1(b) and the pipe is initially straight and aligned in direction x.

## Input File Step by Step

### Creation of nodes 

A total of 101 nodes is created to establish a finite element mesh.

<table><tbody><tr class="odd"><td><p>//Creation of nodes</p><p>Nodes 101</p><p>Node 1 0 0 0</p><p>Node 2 1 0 0</p><p>Node 3 2 0 0</p><p>Node 4 3 0 0</p><p>Node 5 4 0 0</p><p>Node 6 5 0 0</p><p>Node 7 6 0 0</p><p>Node 8 7 0 0</p><p>Node 9 8 0 0</p><p>Node 10 9 0 0</p><p>Node 11 10 0 0</p><p>Node 12 11 0 0</p><p>Node 13 12 0 0</p><p>Node 14 13 0 0</p><p>Node 15 14 0 0</p><p>…</p><p>Node 90 89 0 0</p><p>Node 91 90 0 0</p><p>Node 92 91 0 0</p><p>Node 93 92 0 0</p><p>Node 94 93 0 0</p><p>Node 95 94 0 0</p><p>Node 96 95 0 0</p><p>Node 97 96 0 0</p><p>Node 98 97 0 0</p><p>Node 99 98 0 0</p><p>Node 100 99 0 0</p><p>Node 101 100 0 0</p></td></tr></tbody></table>

### Creation of node sets

Three node sets are defined to establish the constraints and loads. Node sets 1 and 2 are associate with nodes located at the ends of the pipe. Node set 3 contains the node located at the middle of the pipe, where the imperfection load takes place.

<table><tbody><tr class="odd"><td><p>//Creation of node sets</p><p>NodeSets 3</p><p>//First node</p><p>NodeSet 1 Nodes 1 List 1</p><p>//Last node</p><p>NodeSet 2 Nodes 1 List 101</p><p>//Middle node</p><p>NodeSet 3 Nodes 1 List 51</p></td></tr></tbody></table>

### Creation of elements

A total of 50 elements is created. All the elements have the same cross section (identification number 1) and the same coordinate system (identification number 1). Each Pipe\_1 element is defined by three nodes, which must be equally spaced and contained in a straight line.

<table><tbody><tr class="odd"><td><p>//Creation of elements</p><p>Elements 50</p><p>Pipe_1 1 PipeSec 1 CS 1 Nodes 1 2 3</p><p>Pipe_1 2 PipeSec 1 CS 1 Nodes 3 4 5</p><p>Pipe_1 3 PipeSec 1 CS 1 Nodes 5 6 7</p><p>Pipe_1 4 PipeSec 1 CS 1 Nodes 7 8 9</p><p>Pipe_1 5 PipeSec 1 CS 1 Nodes 9 10 11</p><p>Pipe_1 6 PipeSec 1 CS 1 Nodes 11 12 13</p><p>Pipe_1 7 PipeSec 1 CS 1 Nodes 13 14 15</p><p>Pipe_1 8 PipeSec 1 CS 1 Nodes 15 16 17</p><p>Pipe_1 9 PipeSec 1 CS 1 Nodes 17 18 19</p><p>Pipe_1 10 PipeSec 1 CS 1 Nodes 19 20 21</p><p>Pipe_1 11 PipeSec 1 CS 1 Nodes 21 22 23</p><p>Pipe_1 12 PipeSec 1 CS 1 Nodes 23 24 25</p><p>Pipe_1 13 PipeSec 1 CS 1 Nodes 25 26 27</p><p>Pipe_1 14 PipeSec 1 CS 1 Nodes 27 28 29</p><p>Pipe_1 15 PipeSec 1 CS 1 Nodes 29 30 31</p><p>…</p><p>Pipe_1 40 PipeSec 1 CS 1 Nodes 79 80 81</p><p>Pipe_1 41 PipeSec 1 CS 1 Nodes 81 82 83</p><p>Pipe_1 42 PipeSec 1 CS 1 Nodes 83 84 85</p><p>Pipe_1 43 PipeSec 1 CS 1 Nodes 85 86 87</p><p>Pipe_1 44 PipeSec 1 CS 1 Nodes 87 88 89</p><p>Pipe_1 45 PipeSec 1 CS 1 Nodes 89 90 91</p><p>Pipe_1 46 PipeSec 1 CS 1 Nodes 91 92 93</p><p>Pipe_1 47 PipeSec 1 CS 1 Nodes 93 94 95</p><p>Pipe_1 48 PipeSec 1 CS 1 Nodes 95 96 97</p><p>Pipe_1 49 PipeSec 1 CS 1 Nodes 97 98 99</p><p>Pipe_1 50 PipeSec 1 CS 1 Nodes 99 100 101</p></td></tr></tbody></table>

### Creation of the cross section

One only pipe cross section is created, which is shown in the Figure 3‑1(b). It is important to mention that, unlike the tutorials previously presented, there is no need to define the material properties. Instead, here we use equivalent stiffness properties of the pipe to rule the structural behavior. These are to be directly input by the user (axial, bending, torsion and shear stiffness). This element may also be employed for simulations involving drag forces, when using it to represent offshore pipelines, risers or similar structures.

<table><tbody><tr class="odd"><td><p>//Creation of pipe sections</p><p>PipeSections 1</p><p>PS 1 EA 5984734005.0 EI 301817616.8 GJ 301817616.8 GA 2992367003 Rho 239.3893602 CDt 0 CDn 0 CAt 0 CAn 0 De 0.65 Di 0.62</p></td></tr></tbody></table>

### Creation of coordinate systems

Two coordinate systems are created. The first one is used to define the element properties. Note that, since employing Pipe\_1 element, direction E3 must lie in the pipe alignment direction (exactly as previously discussed for Beam\_1 element). The second coordinate system is used to define the direction of application of the nodal loads. It is optional but may provide a convenient way to create loads. Alternatively, the user could also employ the first coordinate system to create the loads.

```text
//Creation of coordinate systems
CoordinateSystems 2
CS 1 E1 0 1 0 E3 1 0 0
CS 2 E1 1 0 0 E3 0 0 1
```

### Creation of solution steps

This tutorial has two solution steps (both are statics). The first step is created establishing the end time equal to “1”. The time step is “1”, which means that the load (imperfection force) is applied in a single step. The second step is created establishing the end time equal to “2”. The time step used is “0.005”, that is, the compression force is divided into increments of 0.5% (variable according to an automatic time-stepping scheme based on difficulties for achieving convergence along simulation evolution).

<table><tbody><tr class="odd"><td><p>//Creation of solution steps</p><p>SolutionSteps 2</p><p>//Step 1 - Static - introduction of the perturbation load. Time: 0-1</p><p>Static 1 EndTime 1 TimeStep 1.0 MaxTimeStep 1.0 MinTimeStep 0.01 MaxIt 15 MinIt 3 ConvIncrease 2 IncFactor 1.0 Sample 1</p><p>//Step 2 - Static - introduction of the follower load. Time: 1-2</p><p>Static 2 EndTime 2 TimeStep 0.005 MaxTimeStep 1 MinTimeStep 0.002 MaxIt 15 MinIt 3 ConvIncrease 2 IncFactor 1.0 Sample 1</p></td></tr></tbody></table>

### Creation of loads

One nodal load and one nodal follower load are defined. Both are established by a force time series. The first is defined in Y direction and the second in X direction (considering CS 2). It is important to remember that the first load is the imperfection force, used to make feasible the evaluation of buckling. The second load is the compression force, which causes the buckling phenomenon studied. Besides, the follower load has a particular characteristic: it updates direction according to the rotations experienced by the node. So, it is possible to state that the compression force is always aligned with the pipe tip, in current configuration. Here the user may also replace the triggering nodal follower load by a nodal load to see the differences in results.

<table><tbody><tr class="odd"><td><p>//Creation of loads</p><p>Loads 2</p><p>//Imperfection load in the middle of the pipe</p><p>NodalLoad 1 NodeSet 3 CS 2 NTimes 2</p><p>//Time FX FY FZ MX MY MZ</p><p>0 0 0 0 0 0 0</p><p>1 0 1000 0 0 0 0</p><p>//Follower load in the last node of the pipe to compress it</p><p>NodalFollowerLoad 2 NodeSet 2 CS 2 NTimes 2</p><p>//Time FX FY FZ MX MY MZ</p><p>1 0 0 0 0 0 0</p><p>2 -362288.4648 0 0 0 0 0</p></td></tr></tbody></table>

### Creation of constraints

Two nodal constraints are created. The first, at the point A of the Figure 3‑1, is done by fixing DOFs: UX, UY, UZ, ROTX. The second constraint, at the point B, is done by fixing DOFs: UY and UZ. Note that Figure 3‑1 represents a bi-dimensional modeling of the problem, but the modeling in GIRAFFE is 3D, so the degree of freedom corresponding to a rotation in X direction must be fixed to restrict such rigid body movement, but it is not represented in Figure 3‑1. Note that the degrees of freedom that were not fixed do not need to appear on the file.

<table><tbody><tr class="odd"><td><p>//Creation of constraints</p><p>Constraints 2</p><p>NodalConstraint 1 NodeSet 1</p><p>UX BoolTable 1</p><p>UY BoolTable 1</p><p>UZ BoolTable 1</p><p>ROTX BoolTable 1</p><p>NodalConstraint 2 NodeSet 2</p><p>UY BoolTable 1</p><p>UZ BoolTable 1</p></td></tr></tbody></table>

### Creation of solver options

Here it is possible to choose the number of processors (threads) to be used (parallel processing). Besides, we must define the solver type for systems of linear equations. In this case, four processors are employed and the direct type of linear system solver is our choice.

<table><tbody><tr class="odd"><td><p>//Creation of solver options</p><p>SolverOptions</p><p>Processors 4 LinSys Direct</p></td></tr></tbody></table>

### Creation of monitors

We request Giraffe to monitor 3 nodes during the simulation: numbers 1 (point A), 51 (middle-spam) and 101 (point B).

<table><tbody><tr class="odd"><td><p>Monitor Sample 1</p><p>MonitorNodes 1 101 51</p></td></tr></tbody></table>

### Creation of post files

The options used here are the same as the tutorial 01, where the WriteMesh, WriteRenderMesh and WriteForces are chosen for post-processing using Paraview<sup>TM</sup> and the MagFactor parameters assume a unit value so the deformed shape of Paraview<sup>TM</sup> frames will reflect the real scale.

<table><tbody><tr class="odd"><td><p>PostFiles</p><p>MagFactor 1</p><p>WriteMesh 1</p><p>WriteRenderMesh 1</p><p>WriteRigidContactSurfaces 0</p><p>WriteFlexibleContactSurfaces 0</p><p>WriteForces 1</p><p>WriteConstraints 0</p><p>WriteSpecialConstraints 0</p><p>WriteContactForces 0</p><p>WriteRenderRigidBodies 0</p><p>WriteRenderParticles 0</p></td></tr></tbody></table>

### Creations of convergence criteria

Here, we must define some convergence and divergence tolerances. For detailed explanation on the meaning of each of these quantities, refer to Giraffe users’ manual.

**Note**: altering default convergence criteria is recommended only for advanced users, since one may create incorrect simulation results by doing so.

<table><tbody><tr class="odd"><td><p>ConvergenceCriteria</p><p>ForceTolerance 1.000000e-004</p><p>MomentTolerance 1.000000e-004</p><p>ForceMinimumReference 1.000000e-001</p><p>MomentMinimumReference 1.000000e-001</p><p>ConstraintMinimumReference 1.000000e-007</p><p>DisplacementTolerance 1.000000e-004</p><p>RotationTolerance 1.000000e-004</p><p>LagrangeTolerance 1.000000e-004</p><p>DisplacementMinimumReference 1.000000e-006</p><p>RotationMinimumReference 1.000000e-006</p><p>LagrangeMinimumReference 1.000000e-006</p><p>DivergenceReference 1.000000e+015</p></td></tr></tbody></table>

## Results and Post-Processing

To post-process results using Paraview<sup>TM</sup>, first, we must open the files “whole\_solution\_forces.pvd” and “whole\_solution\_mesh.pvd”, which are located at the folder “post” in “tutorial03”.

<img src="../../assets/tutorials/image19.png" style="width:5.90278in;height:3.34306in" alt="" />

Figure 3‑2 – Opening the files in Paraview<sup>TM</sup>.

Then, we must click the preview icon and the Apply button on the Pipeline browser. Note that the user has also the possibility of opening individual solution steps results, according to necessity.

<img src="../../assets/tutorials/image20.png" style="width:1.70903in;height:2.62708in" alt="" />

Figure 3‑3 – Pipeline Browser.

Figure 3‑4 shows how Paraview<sup>TM</sup> displays the results.

<img src="../../assets/tutorials/image21.png" style="width:5.90278in;height:2.70139in" alt="" />

Figure 3‑4 – Results using Paraview<sup>TM</sup>

Next, we can create a Glyph for observing the load being applied. For that, go to the Pipeline Browser and click in “whole\_solution\_forces.pvd”. Then, click the Glyph, as shown in Figure 3‑5.

<img src="../../assets/tutorials/image22.png" style="width:2.26875in;height:2.87292in" alt="" />

Figure 3‑5 – Creation of a Glyph.

After creating the Glyph, we must adjust the scale of vectors. For that, go to the Pipeline Browser and, first, click the preview icon of Glyph. Then, in Scaling, set the scale mode to vector and change the scale factor to a small number, e.g.: 0.0001. Finally, click Apply.

<img src="../../assets/tutorials/image23.png" style="width:1.52222in;height:2.84306in" alt="" />

Figure 3‑6 – Adjusting the Glyph.

Finally, we can make an animation and view the structure of the pipe experiencing buckling. If the user wants to know the values of displacement of selected nodes or forces in a certain instant, it is possible to obtain such information in the “monitors” folder. In this case, we requested monitors for three nodes (points where the forces were applied) but the user can request another node. Giraffe generates a text file for each node requested containing the values of displacement and force/moment in each global direction. It is possible to transfer the content of the file to a spreadsheet, getting easier to plot the values.

<img src="../../assets/tutorials/image24.png" style="width:5.90556in;height:1.35069in" alt="" />

Figure 3‑7 – Deformed shape of pipe at the end of simulation.

<img src="../../assets/tutorials/emf_image25.png" style="max-width:100%" alt="" />

Figure 3‑8 – Lateral displacement (UY) of the pipe middle spam vs. time
