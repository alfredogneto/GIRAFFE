# Tutorial 13: Baseball bat hitting the ball

!!! example "Input files"
    [`tutorial13.inp`](../assets/tutorials/inputs/tutorial13/tutorial13.inp)
    (also in the repository under `inputs/tutorial13/`)

Input file name: tutorial13.inp

*Developed by Lucas da Silva*

<img src="../../assets/tutorials/image101.png" style="width:5.90556in;height:4.82778in" alt="" />

Figure 13‑1 – Baseball bat and ball visualization

## Objective

This tutorial presents a geometrically-nonlinear simulation of a baseball bat hitting a rigid ball, the two objects have given initial velocities, the bat has initial rotation velocity and the ball has initial linear velocity and spin. The tutorial has the purpose of showing how to set and solve an impact simulation.

## Description

The bat (0.84 m of length) is meshed using 41 nodes and 20 elements, using the element Beam\_1. The cross-sections of the elements are circular with varying radius to represent the geometry shown in Figure 13‑1. The base node of the bat is constrained, only the rotation about the *y* axis is permitted, our intention is that the motion takes place in the xz plane. The ball is meshed with a single node used in the creation of a rigid spherical particle with given radius.

## Input file step-by-step

### Creation of nodes

A total of 42 nodes is created to establish a finite element mesh.

<table><tbody><tr class="odd"><td><p>//Creation of nodes</p><p>Nodes 42</p><p>//Number X Y Z</p><p>//Bat nodes</p><p>Node 1 0 0 0.000</p><p>Node 2 0 0 0.0021</p><p>Node 3 0 0 0.0042</p><p>Node 4 0 0 0.0063</p><p>Node 5 0 0 0.0084</p><p>…</p><p>Node 38 0 0 0.777</p><p>Node 39 0 0 0.798</p><p>Node 40 0 0 0.819</p><p>Node 41 0 0 0.840</p><p>//Ball node</p><p>Node 42 -0.08 0.000 0.71</p></td></tr></tbody></table>

### Creation of node sets

Two node sets are created: the first one to establish the base constraints, the second one to establish contact.

<table><tbody><tr class="odd"><td><p>//Creation of node sets</p><p>NodeSets 2</p><p>//Base of the bat</p><p>NodeSet 1 Nodes 1 List 1</p><p>//Ball node</p><p>NodeSet 2 Nodes 1 List 42</p></td></tr></tbody></table>

### Creation of elements

A total of 20 elements is defined for the bat. All the elements are assigned to the same material properties and same coordinate system, the cross-section is different for each part of the bat.

<table><tbody><tr class="odd"><td><p>//Creation of elements</p><p>Elements 20</p><p>//”Knob”</p><p>Beam_1 1 Mat 1 Sec 1 CS 1 Nodes 1 2 3</p><p>//”Handle”</p><p>Beam_1 2 Mat 1 Sec 2 CS 1 Nodes 3 4 5</p><p>Beam_1 3 Mat 1 Sec 2 CS 1 Nodes 5 6 7</p><p>Beam_1 4 Mat 1 Sec 2 CS 1 Nodes 7 8 9</p><p>…</p><p>Beam_1 10 Mat 1 Sec 2 CS 1 Nodes 19 20 21</p><p>//”Transition”</p><p>Beam_1 11 Mat 1 Sec 3 CS 1 Nodes 21 22 23</p><p>Beam_1 12 Mat 1 Sec 4 CS 1 Nodes 23 24 25</p><p>Beam_1 13 Mat 1 Sec 5 CS 1 Nodes 25 26 27</p><p>Beam_1 14 Mat 1 Sec 6 CS 1 Nodes 27 28 29</p><p>Beam_1 15 Mat 1 Sec 7 CS 1 Nodes 29 30 31</p><p>Beam_1 16 Mat 1 Sec 8 CS 1 Nodes 31 32 33</p><p>//”Barrel”</p><p>Beam_1 17 Mat 1 Sec 9 CS 1 Nodes 33 34 35</p><p>…</p><p>Beam_1 20 Mat 1 Sec 9 CS 1 Nodes 39 40 41</p></td></tr></tbody></table>

### Creation of particles

The ball (node 42) is represented in the model with a rigid spherical particle of radius 0.037 (m). Giraffe evaluates the ball mass automatically, given the properties of the assigned material (2).

<table><tbody><tr class="odd"><td><p>//Creation of particles</p><p>Particles 1</p><p>Sphere 1 Mat 2 CS 1 Radius 0.037 Node 42</p></td></tr></tbody></table>

### Creation of materials

Two elastic materials are defined: the first one is assigned to the bat elements, it has properties of wood. The second material is assigned to the particle. Since the particle is modeled as a rigid body, its elastic properties are ignored.

<table><tbody><tr class="odd"><td><p>//Creation of materials</p><p>Materials 2</p><p>//Wood</p><p>Hooke 1 E 18.14e9 Nu 0.3 Rho 650</p><p>//Ball</p><p>Hooke 2 E 18.14e9 Nu 0.3 Rho 670</p></td></tr></tbody></table>

### Creation of cross-sections

A total of 9 cross sections is defined. All of them are SuperEllipse type sections (semi-axes A and B, exponent N) reduced to circular sections (A=B, N=2).

<table><tbody><tr class="odd"><td><p>//Creation of sections</p><p>Sections 9</p><p>SuperEllipse 1 A 0.025 B 0.025 N 2 AMeshFDM 100</p><p>SuperEllipse 2 A 0.015 B 0.015 N 2 AMeshFDM 100</p><p>SuperEllipse 3 A 0.018 B 0.018 N 2 AMeshFDM 100</p><p>SuperEllipse 4 A 0.021 B 0.021 N 2 AMeshFDM 100</p><p>SuperEllipse 5 A 0.024 B 0.024 N 2 AMeshFDM 100</p><p>SuperEllipse 6 A 0.027 B 0.027 N 2 AMeshFDM 100</p><p>SuperEllipse 7 A 0.030 B 0.030 N 2 AMeshFDM 100</p><p>SuperEllipse 8 A 0.0315 B 0.0315 N 2 AMeshFDM 100</p><p>SuperEllipse 9 A 0.033 B 0.033 N 2 AMeshFDM 100</p></td></tr></tbody></table>

### Creation of coordinate systems

A single coordinate system is created. It is important to mention that, since employing Beam\_1 element, the direction E3 must lie in the beam alignment direction. The identification number of the coordinate system is 1.

<table><tbody><tr class="odd"><td><p>//Creation of coordinate systems</p><p>CoordinateSystems 1</p><p>CS 1 E1 1 0 0 E3 0 0 1</p></td></tr></tbody></table>

### Creation of solution steps

There are two solution steps, both dynamics. The first step is associated with the impact between the bat and the ball, we force the time steps to be 0.035 ms or less, because the contact time is of about 0.7 ms, so the simulation gets nearly 20 time-steps to model the impact.

The second step takes place right after the first one, it is associated with the system evolution after the collision, the purpose of this step is to permit a greater time step after the impact, when we needed a really small time-step.

<table><tbody><tr class="odd"><td><p>//Creation of solution steps</p><p>/*</p><p>Two steps:</p><p>1st. Bat with given rotation velocity hitting a ball with given</p><p>initial linear velocity and spin</p><p>2nd. Evolution of the system post collision</p><p>*/</p><p>SolutionSteps 2</p><p>Dynamic 1</p><p>EndTime 0.001</p><p>TimeStep 0.000035</p><p>MaxTimeStep 0.000035</p><p>MinTimeStep 0.00000035</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 5</p><p>IncFactor 1.2</p><p>Sample 1</p><p>RayleighDamping Alpha 0 Beta 0 Update 0</p><p>NewmarkCoefficients Beta 0.3 Gamma 0.5</p><p>Dynamic 2</p><p>EndTime 0.025</p><p>TimeStep 0.000035</p><p>MaxTimeStep 0.00035</p><p>MinTimeStep 0.00000175</p><p>MaxIt 20</p><p>MinIt 3</p><p>ConvIncrease 5</p><p>IncFactor 1.2</p><p>Sample 1</p><p>RayleighDamping Alpha 0 Beta 0 Update 0</p><p>NewmarkCoefficients Beta 0.3 Gamma 0.5</p></td></tr></tbody></table>

### Creation of initial conditions

A total of 42 initial conditions is created, one for each node. The first initial condition is associated with the ball (node 42), it has an initial velocity of 40.0 (m/s) on the *x* direction and a rotation of 230 (rad/s) on the *y* direction.

The following 41 conditions are associated with the bat. To all nodes is assigned the same rotation of ω=-50 (rad/s) on the *y* direction, and a velocity given by:

*v* = *ω* ⋅ *r* ,

where *r* is the distance between the node and the base of the bat. With these velocities the bat has rigid body motion of rotation about the Oy axis.

<table><tbody><tr class="odd"><td><p>//Creation of initial conditions</p><p>InitialConditions 42</p><p>// Ball initial conditions</p><p>InitialCondition 1 Node 42 DU 40.0 0 0 OMEGA 0 230 0 SolutionStep 1</p><p>// Bat initial conditions</p><p>InitialCondition 2 Node 1 DU -0.00 0 0 OMEGA 0 -50 0 SolutionStep 1</p><p>InitialCondition 3 Node 2 DU -1.05 0 0 OMEGA 0 -50 0 SolutionStep 1</p><p>InitialCondition 4 Node 3 DU -2.10 0 0 OMEGA 0 -50 0 SolutionStep 1</p><p>InitialCondition 5 Node 4 DU -3.15 0 0 OMEGA 0 -50 0 SolutionStep 1</p><p>…</p><p>InitialCondition 39 Node 38 DU -38.85 0 0 OMEGA 0 -50 0 SolutionStep 1</p><p>InitialCondition 40 Node 39 DU -39.90 0 0 OMEGA 0 -50 0 SolutionStep 1</p><p>InitialCondition 41 Node 40 DU -40.95 0 0 OMEGA 0 -50 0 SolutionStep 1</p><p>InitialCondition 42 Node 41 DU -42.00 0 0 OMEGA 0 -50 0 SolutionStep 1</p></td></tr></tbody></table>

### Creation of constraints

The single constraint is applied to node set 1 (the base of the bat). It restricts all motion except for the rotation about the *y* axis.

<table><tbody><tr class="odd"><td><p>//Creation of constraints</p><p>Constraints 1</p><p>NodalConstraint 1 NodeSet 1</p><p>UX BoolTable 1</p><p>UY BoolTable 1</p><p>UZ BoolTable 1</p><p>ROTX BoolTable 1</p><p>ROTY BoolTable 0</p><p>ROTZ BoolTable 1</p></td></tr></tbody></table>

### Creation of contact surfaces

One contact surface is defined to enclose the barrel of the bat, the region the ball is supposed to hit. It has similar geometrical properties as the cross-section 9 assigned to the barrel part. The normal to the surface is set to be exterior, pointing to the region where the ball is supposed to be. Just like the barrel elements, the surface goes from node 33 to node 41.

<table><tbody><tr class="odd"><td><p>//Creation of contact surfaces</p><p>Surfaces 1</p><p>FlexibleSECylinder_1 1 A 0.033 B 0.033 N 2 CS 1 NormalExterior Nodes 33 41</p></td></tr></tbody></table>

### Creation of surface sets

One surface set is required to establish contact, it contains the only surface.

<table><tbody><tr class="odd"><td><p>//Creation of surface sets</p><p>SurfaceSets 1</p><p>SurfaceSet 1 Surfaces 1 List 1</p></td></tr></tbody></table>

### Creation of contact pairs

There are two interacting bodies, one of them is a spherical object. The Node Set to Surface Set (NSSS) contact is perfect for this problem because it considers that the nodes in the given node set are centers of spherical surfaces of given radius, therefore the only surface to be concerned is the barrel surface that was already defined and linked to the surface set.

The node set 2 (associated with the ball) and the surface set 1 are chosen. Then, some contact parameters are defined:

-   friction coefficient (MU)

-   normal penalty parameter (EPN)

-   normal dissipative coefficient (CN)

-   tangential penalty parameter (EPT)

-   tangential dissipative coefficient (CT)

The EPN parameter is strongly related with the collision time, the greater its value the smaller the period of contact. The CN parameter is also related with the collision time, but more importantly, with the coefficient of restitution (COR) of the collision, since it establishes a dissipation coefficient, proportional to the relative velocity of contacting bodies’ surfaces.

The values of these two parameters were set to establish a contact time about 0.0007 (s) and a COR about 0.33.

The value of EPT is the same as the EPN value with one order of magnitude less. The CT parameter is set equal to the CN parameter.

The pinball radius just needs to be large enough for contact detection to occur before the penetration of the sphere and the barrel surface.

The Radius parameter is the ball radius of 0.037 (m).

<table><tbody><tr class="odd"><td><p>//Creation of contact pairs</p><p>Contacts 1</p><p>NSSS 1 NodeSet 2 SurfaceSet 1 MU 0.3 EPN 3.2e6</p><p>CN 2e2 EPT 3.2e5 CT 2e2 Pinball 1 Radius 0.037</p><p>MaxPointwiseInt 1 BoolTable 1</p></td></tr></tbody></table>

### Creation of solver options

Here we request Giraffe to apply a single thread to solve this case. The linear system solution chosen is Direct.

<table><tbody><tr class="odd"><td><p>//Creation of solver options</p><p>SolverOptions</p><p>Processors 1 LinSys Direct</p></td></tr></tbody></table>

### Creation of monitors

We request Giraffe to monitor node 42 (the ball) and the contact pair during the simulation.

<table><tbody><tr class="odd"><td><p>//Creation of monitors</p><p>Monitor Sample 1</p><p>MonitorNodes 42</p><p>MonitorContacts 1</p></td></tr></tbody></table>

### Creation of post-files

We choose WriteRenderMesh, WriteFlexibleContactSurfaces, WriteContactForces and WriteRenderParticles options for post-processing using Paraview<sup>TM</sup>.

<table><tbody><tr class="odd"><td><p>//Creation of post files</p><p>PostFiles</p><p>MagFactor 1</p><p>WriteMesh 0</p><p>WriteRenderMesh 1</p><p>WriteRigidContactSurfaces 0</p><p>WriteFlexibleContactSurfaces 1</p><p>WriteForces 0</p><p>WriteConstraints 0</p><p>WriteSpecialConstraints 0</p><p>WriteContactForces 1</p><p>WriteRenderRigidBodies 0</p><p>WriteRenderParticles 1</p></td></tr></tbody></table>

## Results and post-processing

To visualize our model in Paraview<sup>TM</sup>, first we have to open the files “whole\_solution\_contactsurfaces.pvd”, “whole\_solution\_rb\_particles.pvd” and “whole\_solution\_rendermesh.pvd”, which are located at the folder “post” in “tutorial14”.

<img src="../../assets/tutorials/image102.png" style="width:5.90556in;height:3.15278in" alt="" />

Figure 13‑2 - Opening the files of tutorial 14 in Paraview<sup>TM</sup>.

Then, we can click the preview icon and the Apply button on the Pipeline browser. After doing that, we are ready to post-process the results.

<img src="../../assets/tutorials/image103.png" style="width:2.61495in;height:5.407in" alt="" />

Figure 13‑3 – Initial position of the system (render view).

It is possible to color the different elements of the model improving the visualization. To do that, we first select the object on the Pipeline browser (Figure 13‑4), then we go to Properties tab, click edit in the Coloring section and pick a color from the Basic colors or from the color panel or typing the RGB or HSV coordinates of the color (Figure 13‑5).

<img src="../../assets/tutorials/image104.png" style="width:2.92749in;height:3.87554in" alt="" />

Figure 13‑4 – Choosing object on Pipeline browser.

<img src="../../assets/tutorials/image105.png" style="width:5.90194in;height:3.16931in" alt="" />

**Figure** **13‑5 – Changing object color property.**

We can also change the Opacity of the object on the Styling section. After some changes, we get the following result.

<img src="../../assets/tutorials/image107.png" style="width:3.27129in;height:7.78234in" alt="" />

Figure 13‑6 – Choosing object on Pipeline browser.

Notice that the contact surface is present only on the top part of the bat.

We can now visualize the animation clicking on the play button on the Main Controls toolbar. In Figure 13‑7, we can see some frames of the animation, arrows indicating the velocity of the ball were added.

<img src="../../assets/tutorials/image107.png" style="width:1.15694in;height:2.75556in" alt="" /><img src="../../assets/tutorials/image108.png" style="width:0.95694in;height:2.75556in" alt="" /><img src="../../assets/tutorials/image109.png" style="width:1.20849in;height:2.75556in" alt="" /><img src="../../assets/tutorials/image110.png" style="width:2.18788in;height:2.75556in" alt="" /><img src="../../assets/tutorials/image111.png" style="width:2.63152in;height:2.75556in" alt="" /><img src="../../assets/tutorials/image112.png" style="width:3.13333in;height:2.75556in" alt="" />

Figure 13‑7 – Some selected frames of the simulation.

We could also plot a Glyph representing the contact force, but this is only for visualization purposes since we will inspect the forces by means of the contact monitor. The monitor files can be found in folder monitors in tutorial14.

Giraffe generates a monitor text file for the requested contact containing the time series of the normal force and friction force (Cartesian components and magnitude). Therefore, we can plot the time series of the forces and verify the impact time.

<img src="../../assets/tutorials/chart_t13_f8.png" style="max-width:100%" alt="" />

Figure 13‑8 – Normal Force vs Time.

<img src="../../assets/tutorials/chart_t13_f9.png" style="max-width:100%" alt="" />

Figure 13‑9 – Friction Force vs Time.

Giraffe also generates a file for the requested node containing the values of displacement and velocity in each global direction. Therefore, it is possible to plot the components of the velocity and obtain the exit velocity of the ball after the collision.

Figure 13‑10 shows the *x* component of the velocity of the ball, the dominant coordinate. The ball approaches the bat with speed 40 m/s and exits with a larger speed.

<img src="../../assets/tutorials/chart_t13_f10.png" style="max-width:100%" alt="" />

Figure 13‑10 – Velocity **x** component vs Time.

There is some exit velocity component in the *z* direction, gained because of the friction, which transforms the rotation energy from the spin into kinetic energy, as show in Figure 13‑11.

<img src="../../assets/tutorials/chart_t13_f11.png" style="max-width:100%" alt="" />

Figure 13‑11 – Velocity **z** component vs Time.
