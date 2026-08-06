# Post Files

*Creates post files for post-processing results using Paraview<sup>TM</sup> post-processor interface.*

### Syntax:

PostFiles

MagFactor MFV

WriteMesh MF

WriteRenderMesh RMF

WriteRigidContactSurfaces RCSF

WriteFlexibleContactSurfaces FCSF

WriteForces FF

WriteConstraints CF

WriteSpecialConstraints SCF

WriteContactForces CFF

WriteRenderRigidBodies RBF

WriteRenderParticles PF

-   MFV: value of the magnification factor for displacements (for visualization purposes)

-   MF: Boolean flag to write (1) or not (0) the mesh file

-   RMF: Boolean flag to write (1) or not (0) the render mesh file

-   RCSF: Boolean flag to write (1) or not (0) the rigid contact surfaces file

-   FCSF: Boolean flag to write (1) or not (0) the flexible contact surfaces file

-   FF: Boolean flag to write (1) or not (0) the forces file

-   CF: Boolean flag to write (1) or not (0) the constraints file

-   SCF: Boolean flag to write (1) or not (0) the special constraints file

-   CFF: Boolean flag to write (1) or not (0) the contact forces file

-   RBF: Boolean flag to write (1) or not (0) the rigid bodies file

-   RBF: Boolean flag to write (1) or not (0) the particles file

### Example:

PostFiles

MagFactor 1.0

WriteMesh 1

WriteRenderMesh 1

WriteRigidContactSurfaces 0

WriteFlexibleContactSurfaces 0

WriteForces 0

WriteConstraints 0

WriteSpecialConstraints 0

WriteContactForces 0

WriteRenderRigidBodies 0

WriteRenderParticles 0

### Additional information:

PostFiles keyword activates saving of output files containing information for post-processing the simulation using PARAVIEW<sup>TM</sup>. The sampling for saving post files is established on solution steps definition.

The MagFactor keyword is a magnification factor that will be used to multiply all the displacements experienced in the model, in the visualization of deformed shape frames. If the user enters “1.0”, the deformed shape will show deformation patterns in real scale. A larger value than “1.0” can be used in case of simulations involving very small displacements/rotations, to help on visualization of results.

Some write control flags have to be set by the user. Each one may be turned on/off, by the values “1” or “0”, respectively. The choice of adequate save outputs permit to visualize more details of the model and are very useful for generating high-quality animations and good post-processing interpretations.

Each write control flag is described below:

-   WriteMesh: to write the base mesh information. E.g.: beams are represented by lines passing representing the axis. Particles are represented by points. Shells are represented my mid-surfaces.

-   WriteRenderMesh: to write the rendered mesh, every element as a 3D solid. E.g.: beams are represented by the chosen cross section extruded along the axis direction. Shells are represented including the thickness information.

-   WriteRigidContactSurfaces: to write the rigid contact surfaces. E.g.: RigidTriangularSurface\_1.

-   WriteFlexibleContactSurfaces: to write the flexible contact surfaces. E.g.: FlexibleTriangularSurface\_2.

-   WriteForces: to write data containing information of external applied forces on nodes. It covers NodalLoads and NodalFollowerLoads. It may be used to construct arrow glyphs in Paraview<sup>TM</sup> interface.

-   WriteConstraints: to write constraints symbols.

-   WriteSpecialConstraints: to write special constraints symbols.

-   WriteContactForces: to write data with contact forces locations and associated normal and friction values. It may be used to construct arrow glyphs in Paraview<sup>TM</sup> interface.

-   WriteRenderRigidBodies: to write data with rigid bodies rendering points.

-   WriteRenderParticles: to write particles data. External surfaces are represented.

Useful data is written by Giraffe when requesting WriteRenderMesh. A vector data associated with each cell, named ElementProperties contains information about:

-   Element type associated number (according to Table 7)

-   Associated material number

-   Associated section number (for beams, shells and trusses)

-   Associated coordinate system number

Table 7 – Element types and associated numbers

| **Element type** | **number** |
|------------------|------------|
| Beam\_1          | 1          |
| Pipe\_1          | 2          |
| Shell\_1         | 3          |
| Mass\_1          | 4          |
| SpringDashpot\_1 | 5          |
| RigidBody\_1     | 6          |
| Truss\_1         | 7          |

The objective of this data is to provide the user the possibility of creating selections for better post-processing in Paraview<sup>TM</sup> (e.g.: selecting only cells associated with beam elements, or with a given material number, etc.). See Appendix for more information.

**Note**: post files are always created in the model, even if the user does not request then. In such case, the only output will be the base mesh.
