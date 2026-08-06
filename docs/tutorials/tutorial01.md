# Tutorial 1: Cantilever beam statics

!!! example "Input files"
    [`tutorial01.inp`](../assets/tutorials/inputs/tutorial01/tutorial01.inp)
    (also in the repository under `inputs/tutorial01/`)

Input file name: tutorial01.inp

*Developed by Alfredo Gay Neto*

## Objective

This tutorial presents a geometrically-nonlinear simulation of a cantilever beam. It aims at showing how to apply nodal loads and nodal boundary conditions using Giraffe to solve a simple static geometrically-nonlinear analysis.

<img src="../../assets/tutorials/image2.png" style="width:2.28346in;height:1.71528in" alt="" /> (a)<img src="../../assets/tutorials/image3.png" style="width:2.78353in;height:1.83889in" alt="" />(b)

Figure 1‑1 – (a) geometric description (b) 3D rendering for visualization purposes

## Description 

The beam AB (1.0 m length) is meshed using 11 nodes and 5 elements, using the element Beam\_1. There is a clamp boundary condition at point A (origin of the global coordinate system). A force F is applied at point B. The cross section of the beam is square. The beam is initially straight and aligned in direction z.

## Input file step-by-step

### Creation of nodes

Eleven nodes are created to establish a finite element mesh.

//Creation of nodes

Nodes 11

//Number X Y Z

Node 1 0 0 0.0

Node 2 0 0 0.1

Node 3 0 0 0.2

Node 4 0 0 0.3

Node 5 0 0 0.4

Node 6 0 0 0.5

Node 7 0 0 0.6

Node 8 0 0 0.7

Node 9 0 0 0.8

Node 10 0 0 0.9

Node 11 0 0 1.0

### Creation of node sets

Node sets are created to establish loads and constraints. Node sets may be created by distinct syntaxes. In present case, the more convenient is to directly list the nodes included in each node set.

//Creation of node sets

NodeSets 2

//Point A

NodeSet 1 Nodes 1 List 1

//Point B

NodeSet 2 Nodes 1 List 11

### Creation of elements

Five elements are defined. All the elements are assigned to the same material properties (identification number 1), same cross section (identification number 1) and same coordinate system (identification number 1). Each Beam\_1 element is defined by three nodes, which **must** be equally spaced and contained in a straight line.

//Creation of elements

Elements 5

Beam\_1 1 Mat 1 Sec 1 CS 1 Nodes 1 2 3

Beam\_1 2 Mat 1 Sec 1 CS 1 Nodes 3 4 5

Beam\_1 3 Mat 1 Sec 1 CS 1 Nodes 5 6 7

Beam\_1 4 Mat 1 Sec 1 CS 1 Nodes 7 8 9

Beam\_1 5 Mat 1 Sec 1 CS 1 Nodes 9 10 11

### Creation of materials

One only material model is defined, which identification number is 1.

//Creation of materials

Materials 1

Hooke 1 E 1e7 Nu 0.3 Rho 2000

### Creation of cross-sections

One only cross-section is defined, a rectangle. Both B and H are set to 0.1 m (base and height), to represent a square cross-section. The identification number of the cross-section is 1.

//Creation of sections

Sections 1

Rectangle 1 B 0.1 H 0.1

### Creation of coordinate systems

One only coordinate system is input. This is used to define the element properties. Note that, since employing Beam\_1 element, the direction E3 must lie in the beam alignment direction. The identification number of the coordinate system is 1.

**Note**: In case of creation of not-aligned beam elements (for example, for representing an approximation of a curved beam or a frame), each beam element must be associated with a defined coordinate system, which direction E3 is aligned with the element direction. This is done in order to permit creation of arbitrary alignment between successive cross sections of a beam-like structure, composed possibly by many cross-sections or alignments of these.

//Creation of coordinate systems

CoordinateSystems 1

CS 1 E1 1 0 0 E3 0 0 1

### Creation of solution steps

This tutorial has a single solution step (static). This is created establishing the end time equal to “1”. This is an arbitrary tracking scalar parameter for statics. Time-steps are established to enforce total load to be divided in increments of 10%. Note that IncFactor 1.0 permits no increasing of time-stepping. The option Sample 2 ensures saving post-processing files for each 2 converged configurations. The larger the Sample set, the smaller number of files for post-processing that Giraffe will save along simulation evolution.

//Creation of the solution steps

SolutionSteps 1

Static 1

EndTime 1

TimeStep 0.1

MaxTimeStep 0.1

MinTimeStep 0.01

MaxIt 20

MinIt 3

ConvIncrease 4

IncFactor 1.0

Sample 2

### Creation of loads

A single nodal load is created. This is done establishing a time-varying force in direction X. Giraffe automatically performs linear interpolation for undefined time values following input data. Note that the user may choose local coordinate systems to define vector quantities of loads (forces/moments). In current case CS 1 was the choice.

//Creation of loads

Loads 1

NodalLoad 1 NodeSet 2 CS 1 NTimes 2

//Time FX FY FZ MX MY MZ

0 0 0 0 0 0 0

1 1000 0 0 0 0 0

### Creation of constraints

A single nodal constraint is created. This is done fixing DOFs: UX, UY, UZ, ROTX, ROTY and ROTZ. Note that BoolTable resource is employed. In this case the nodal constraint is established for the solution step 1 (the only step in this tutorial). So, just number “1” is entered, as a meaning of a Boolean variable “true” (to activate such constraint in the first solution step).

**Note**: In case of simulations with multiple solution steps the user may create a sequence of “1”and “0” to activate or inactivate constraints. Each “true” or “false” will refer to corresponding solution steps, in a solution sequence.

//Creation of constraints

Constraints 1

NodalConstraint 1 NodeSet 1

UX BoolTable 1

UY BoolTable 1

UZ BoolTable 1

ROTX BoolTable 1

ROTY BoolTable 1

ROTZ BoolTable 1

### Creation of solver options

Here we request Giraffe to apply a single processor thread to solve this case. If the user sets a larger number of processors, Giraffe will automatically use a shared memory parallel processing based on OpenMP<sup>TM</sup>. The linear system solution method chosen is Direct. Giraffe employs PARDISO<sup>TM</sup> and MKL<sup>TM</sup> Intel library to solve linear systems of equations along solution evolution.

//Creation of solver options

SolverOptions

Processors 1 LinSys Direct

### Creation of monitors

We request Giraffe to monitor 2 nodes during the simulation: numbers 1 and 11. This will create automatically a folder of monitors for post-processing data. Monitors are very useful for analyzing particular time-series of nodes, node-sets, contact regions or elements. When the user is interested in analyzing the whole mesh results at a given instant, the more convenient is not using monitors, but looking at solution post-processing files, instead. These are located in “post” folder and are organized according to the sequence of solution steps established for the simulation.

//Creation of monitors

Monitor Sample 1

MonitorNodes 1 11

### Creation of post-files

We choose WriteMesh, WriteRenderMesh and WriteForces options for post-processing using Paraview<sup>TM</sup>. The MagFactor parameters can be used to generate plots with exaggerated displacement patterns, to help to visualize in case of small displacements. This is not the case of this example. Then, we will assume a unit value, then, requesting that the deformed shape of Paraview<sup>TM</sup> frames will reflect the real scale. Use the MagFactor resource with care. It may create unrealistic plots and not-straightforward interpretations. If in doubt, the best is to set MagFactor 1.

//Creation of post files

PostFiles

MagFactor 1

WriteMesh 1

WriteRenderMesh 1

WriteRigidContactSurfaces 0

WriteFlexibleContactSurfaces 0

WriteForces 1

WriteConstraints 0

WriteSpecialConstraints 0

WriteContactForces 0

WriteRenderRigidBodies 0

WriteRenderParticles 0

Note that the Giraffe input file commands not necessarily has to follow the here proposed sequence. The important is that all the commands here presented are included in this tutorial.

## Results and post-processing

Result files are automatically saved during the solving process on the folder “/tutorial01/post/” and “/tutorial01/monitors/”

### Paraview<sup>TM</sup> post-processing

Opening in Paraview<sup>TM</sup> the “.pvd” files generated during simulation leads to visualization of the model results. The user may perform animations and extract high-quality images for reports.

<img src="../../assets/tutorials/image4.png" style="width:5.90556in;height:3.33264in" alt="" />

Figure 1‑2 – Results rendered using Paraview<sup>TM</sup>

As a particular feature, forces output may be used in Paraview<sup>TM</sup> together with Glyphs in order to establish special visualizations. For that, select “whole\_solution\_forces.pvd” in the Pipeline Browser. Then, click in the Glyph feature button to create a Glyph.

<img src="../../assets/tutorials/image5.png" style="width:2.79097in;height:2.0568in" alt="" />

Figure 1‑3 – Using Glyphs in Paraview<sup>TM</sup>

After creating the Glyph, select it on Pipeline Browser to modify its default options as depicted next. Glyph type is set to “Arrow”. Scale Mode is set to “vector” and scale factor is set to 0.001. This value is arbitrary and may vary from model to model and controls the size of the created Glyph feature. It should be chosen based on the desired size of the Glyph feature on the visualization of the model.

<img src="../../assets/tutorials/image6.png" style="width:2.84722in;height:3.15043in" alt="" />

Figure 1‑4 – Setting Glyph attributes in Paraview<sup>TM</sup>

Visualization will look like the shown in Figure 1‑5 (changing the color of the Glyph feature to red):

<img src="../../assets/tutorials/image7.png" style="width:5.90556in;height:3.33264in" alt="" />

Figure 1‑5 – Using Glyphs in Paraview<sup>TM</sup>

The user may animate results to see the Glyph feature changing its size, according to the force magnitude along time.

### Monitors post-processing

Monitors results are text files that saves time-series of specific results. In this example, they are associated with nodes. An example of post-processing is to plot the point B position along “time”, from 0-1.

<img src="../../assets/tutorials/emf_image8.png" style="max-width:100%" alt="" />

Figure 1‑6 – Monitoring the point B coordinates time-evolution

## What to try now?

-   change the magnitude and the direction of the force

-   create a load composed by moments, instead of forces

-   change the boundary conditions and apply loads in distinct nodes

-   change the time-step size to smaller or larger values to see the influence on simulation evolution

-   change the sampling rate for results saving both on monitor and on solution steps

-   change the beam cross-section and material properties

-   refine the mesh
