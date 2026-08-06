# Tutorial 11: Two degrees of freedom lumped system

!!! example "Input files"
    [`tutorial11.inp`](../assets/tutorials/inputs/tutorial11/tutorial11.inp)
    (also in the repository under `inputs/tutorial11/`)

Input file name: tutorial11.inp

*Developed by Gabriel Maki and Alfredo Gay Neto*

## Objective

This tutorial presents a damped mass-spring system with two degrees of freedom. It aims at showing how to create a spring/dashpot element and how to apply a time-varying load. Besides, this tutorial proposes a simple example of a tuned mass damper.

<img src="../../assets/tutorials/image78.png" style="width:2.74792in;height:0.88195in" alt="" /> (a) <img src="../../assets/tutorials/image79.png" style="width:2.51242in;height:0.88194in" alt="" />(b)

Figure 11‑1 - (a) Geometric description of the system (b) 2D rendering for visualization purposes

## Description

The system consists of two lumped masses M<sub>1</sub> and M<sub>2</sub>, where M<sub>1</sub> = 22515.818 kg and M<sub>2</sub> = 500 kg, linked by springs and dashpots. Figure 11‑1 shows the model scheme and the degrees of freedom representation. In Giraffe Finite Element environment, the system is represented by three nodes and four elements. The first node is fixed, while the others are free in a single direction (global x). Two elements are created to represent the spring/dashpot and another two are created for the point masses. U<sub>1</sub>(t) is defined as the node 2 x-displacement, U<sub>2</sub>(t) is defined as the node 3 x-displacement and F(t) is the time varying (harmonic) load, applied at node 2 in direction x.

## Input file step-by-step

### Creation of nodes

Three nodes are created to represent the system modeled.

<table><tbody><tr class="odd"><td><p>//Creation of nodes</p><p>Nodes 3</p><p>//Number X Y Z</p><p>Node 1 0 0 0</p><p>Node 2 0.5 0 0</p><p>Node 3 1.0 0 0</p></td></tr></tbody></table>

### Creation of node sets

Four node sets are created to establish loads and constraints.

<table><tbody><tr class="odd"><td><p>NodeSets 4</p><p>NodeSet 1 Nodes 1 List 1</p><p>NodeSet 2 Nodes 1 List 2</p><p>NodeSet 3 Nodes 1 List 3</p><p>NodeSet 4 Nodes 2 List 2 3</p></td></tr></tbody></table>

### Creation of elements

Four elements are defined, using two different types of elements. The first is SpringDashpot\_1 element. It is used to define a spring/dashpot set. The second element is Mass\_1, used to establish a lumped mass.

<table><tbody><tr class="odd"><td><p>Elements 4</p><p>SpringDashpot_1 1 Stiffness 500000 Damping 5000 Nodes 1 2</p><p>SpringDashpot_1 2 Stiffness 11103.305 Damping 0 Nodes 2 3</p><p>Mass_1 3 Mass 22515.818 Node 2</p><p>Mass_1 4 Mass 500 Node 3</p></td></tr></tbody></table>

### Creation of coordinate systems

One coordinate system is created, using the same directions of the global coordinate system. It is used to define the nodal load.

<table><tbody><tr class="odd"><td><p>CoordinateSystems 1</p><p>CS 1 E1 1 0 0 E3 0 0 1</p></td></tr></tbody></table>

### Creation of solution steps

This tutorial has a single dynamic solution step. The time of simulation is defined as 70 seconds because it should be enough to establish a steady state oscillation pattern.

<table><tbody><tr class="odd"><td><p>SolutionSteps 1</p><p>Dynamic 1</p><p>EndTime 70</p><p>TimeStep 0.01</p><p>MaxTimeStep 0.01</p><p>MinTimeStep 0.00001</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 4</p><p>IncFactor 1.0</p><p>Sample 10</p><p>RayleighDamping Alpha 0 Beta 0 Update 0</p><p>NewmarkCoefficients Beta 0.3 Gamma 0.5</p></td></tr></tbody></table>

### Creation of loads

The external load F(t) applied on the first mass (node 1) is created establishing a time series.

$$F\left( t \right) = \ \rho_{0}\sin\left( \overline{\omega}t + \varphi \right);$$

Where *ρ*<sub>0</sub> = 25, 000 N is the load’s amplitude, $\overline{\omega} = 4.712\ rad/s$ is the load’s angular frequency and *φ* = 0 *r**a**d* is the load’s initial phase.

The time step used for time series is 0.05 s. Note that Giraffe performs automatically a linear interpolation between input data in a table form, as shown in the input file extract, next. Thus, if one needs a higher resolution in a time-series, there is a need for establishing an input file with more data times. Note also that an arbitrary time-series may be input. Our choice, for this tutorial, was a harmonic time-series.

<table><tbody><tr class="odd"><td><p>Loads 1</p><p>NodalLoad 1 NodeSet 2 CS 1 NTimes 1401</p><p>//Time FX FY FZ MX MY MZ</p><p>Time 0 0 0 0 0 0 0</p><p>Time 0.05 5836.134096 0 0 0 0 0</p><p>Time 0.1 11349.76249 0 0 0 0 0</p><p>Time 0.15 16236.20121 0 0 0 0 0</p><p>Time 0.2 20225.42486 0 0 0 0 0</p><p>…</p><p>Time 69.85 16236.20121 0 0 0 0 0</p><p>Time 69.9 11349.76249 0 0 0 0 0</p><p>Time 69.95 5836.134097 0 0 0 0 0</p><p>Time 70 1.53355E-07 0 0 0 0 0</p></td></tr></tbody></table>

### Creation of Solver Options

Here it is possible to choose the number of processors to be used. Besides, we must define the solver type for systems of linear equations. In this case, only one processor is employed, and the direct type of linear system solver is chosen.

<table><tbody><tr class="odd"><td><p>SolverOptions</p><p>Processors 1 LinSys Direct</p></td></tr></tbody></table>

### Creation of Monitors

Here we request Giraffe to monitor all the nodes during the simulation, for post-processing purposes.

<table><tbody><tr class="odd"><td><p>Monitor Sample 10</p><p>MonitorNodes 1 2 3</p></td></tr></tbody></table>

### Creation of post files

The options used here are the same as the tutorial 01, where the WriteMesh, WriteRenderMesh and WriteForces are chosen for post-processing using Paraview<sup>TM</sup> and the MagFactor parameters assume a unit value so the deformed shape of Paraview<sup>TM</sup> frames will reflect the real geometric scale.

<table><tbody><tr class="odd"><td><p>PostFiles</p><p>MagFactor 1</p><p>WriteMesh 1</p><p>WriteRenderMesh 1</p><p>WriteRigidContactSurfaces 0</p><p>WriteFlexibleContactSurfaces 0</p><p>WriteForces 1</p><p>WriteConstraints 0</p><p>WriteSpecialConstraints 0</p><p>WriteContactForces 0</p><p>WriteRenderRigidBodies 0</p><p>WriteRenderParticles 0</p></td></tr></tbody></table>

## Results and post-processing

To post-process results using Paraview<sup>TM</sup>, first we must open the files “whole\_solution\_forces.pvd”, “whole\_solution\_mesh.pvd” and “whole\_solution\_rendermesh.pvd”, which are located at the folder “post” in “tutorial11”.

<img src="../../assets/tutorials/image80.png" style="width:5.74943in;height:3.28104in" alt="" />

Figure 11‑2 - Opening files of tutorial 11 for post-processing

Then, we must click the preview icon and the Apply button on the Pipeline browser.

<img src="../../assets/tutorials/image81.png" style="width:1.77541in;height:2.40558in" alt="" />

Figure 11‑3 - Pipeline browser for tutorial 11.

We may use an isometric perspective and set a two-dimensional view to obtain a better view.

<img src="../../assets/tutorials/image82.png" style="width:1.99943in;height:2.18623in" alt="" />

Figure 11‑4 - Setting an isometric perspective.

<img src="../../assets/tutorials/image83.png" style="width:1.60417in;height:1.21875in" alt="" />

Figure 11‑5 - Defining the view direction, for a two-dimensional visualization.

After doing that, the Figure 11‑6 shows how would be the view.

<img src="../../assets/tutorials/image84.png" style="width:4.99943in;height:3.28573in" alt="" />

Figure 11‑6 - Results of tutorial 11 using Paraview<sup>TM</sup>.

Next, we can create a Glyph for observing the load being applied. For that, go to the Pipeline Browser and click in “whole\_solution\_forces.pvd”. Then, click the Glyph, as shown in Figure 11‑7.

<img src="../../assets/tutorials/image85.png" style="width:2.25in;height:1.84375in" alt="" />

Figure 11‑7 - Creation of a Glyph on tutorial 11.

After creating the Glyph, we must adjust the scale of vectors. For that, go to the Pipeline Browser and, first, click the preview icon of Glyph. Then, in Scaling, set the scale mode to vector and change the factor to a small number, e.g.: 0.00002.

<img src="../../assets/tutorials/image86.png" style="width:2.31924in;height:1.35473in" alt="" />

Figure 11‑8 - Adjusting the Glyph on tutorial 11.

Obs: If the user wants to change the arrow’s color, the one can do it by clicking in Edit, on Coloring at Glyph’s Properties.

<img src="../../assets/tutorials/image87.png" style="width:2.88542in;height:1.94792in" alt="" />

Figure 11‑9 - Changing the arrow’s color on tutorial 11.

Finally, we can perform an animation to observe the system’s movement. When needing to plot time series of displacement, forces or other kinematical quantities, monitors may be helpful.

<img src="../../assets/tutorials/image88.png" style="width:5.90556in;height:1.15694in" alt="" />

Figure 11‑10 - System’s oscillation.

<img src="../../assets/tutorials/emf_image89.png" style="max-width:100%" alt="" />

Figure 11‑11 - U<sub>1</sub> (M<sub>1</sub> Displacement).

<img src="../../assets/tutorials/emf_image90.png" style="max-width:100%" alt="" />

Figure 11‑12 - U<sub>2</sub> (M<sub>2</sub> Displacement).

One interesting thing to observe is the order of magnitude of displacements. M<sub>2</sub> has displacements much larger than M<sub>1</sub>. Besides, it is possible to view that U<sub>1</sub> experiences greater damping, compared to U<sub>2</sub>.

If we model only M<sub>1</sub> like shown in Figure 11‑13, the system will resonate because the load’s frequency is equal to the natural frequency. Figure 11‑14 presents the U<sub>1</sub> behavior in this case.

<img src="../../assets/tutorials/image91.png" style="width:2.87443in;height:1.48893in" alt="" />

Figure 11‑13 - Geometric description of the system of one degree of freedom.

<img src="../../assets/tutorials/emf_image92.png" style="max-width:100%" alt="" />

Figure 11‑14 - Displacement of M<sub>1</sub> – System of one degree of freedom.

In Figure 11‑14, it is possible to see that the displacement increases until establish a steady state (due to damping), which is characteristic of damped 1 degree-of-freedom resonant system. By including the M<sub>2</sub>, M<sub>1</sub> stops vibrating on resonance and the displacements decreases over time, like Figure 11‑11. Due to this effect that M<sub>2</sub> has on the system, we call it a tuned mass damper or harmonic absorber. It is an auxiliary device used to reduce the amplitude of vibration of a structure to an acceptable level.
