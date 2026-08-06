# Rigid Body Data

*Starts a command block for creation of rigid body data, to be used together with RigidBody\_1 element type.*

### Syntax:

RigidBodyData N

RBData RBID

Mass MV

J11 J11V J22 J22V J33 J33V J12 J12V J13 J13V J23 J23V

Barycenter XGV YGV ZGV

//Optional keywords block below:

CADData CADN

-   N: number of rigid body data

-   RBID: current rigid body data identification number

-   MV: rigid body mass value

-   J11V, J22V, J33V, J12V, J13V and J23V: rigid body inertia tensor components values

-   XGV, YGV, and ZGV: coordinates of the barycenter position

-   CADN: number of the CAD data ID for post-processing purposes (see CADData keyword)

### Example:

RigidBodyData 1

RBData 1

Mass 3.63e-5

J11 0.00663 J22 0.00480 J33 0.00663 J12 0.0 J13 0.0 J23 0.0

Barycenter 0.0 0.0 0.0

CADData 1

### Additional information:

Each rigid body data is defined by the keyword RBData followed by an identification number (must be an ascending sequence starting from number one).

The mass value has to be provided according to the unit system adopted. In the example, length is defined in millimeters, force in Newton and time in seconds, so that mass has to be provided in tonnes. Users are allowed to choose any other consistent unit system for the whole model.

Inertia properties (J11, J22, J33, J12, J13, J23) must be provided with respect to barycentric axes, parallel to de CAD coordinate system. This is illustrated in Figure 19.

<img src="../../assets/manual/image22.png" style="width:3.14961in;height:3.29959in" alt="" />

Figure 19 – Inertia properties for Rigid Bodies

In Figure 19 we have the CAD coordinate system (P<sub>P</sub> xyz) and the body is already positioned in Giraffe (system OXYZ). Inertia properties should be provided with respect to axes x\*, y\* and z\* illustrated in the same picture, these axes are parallel to x, y and z, but have their origin located at the center of mass of the body (G). Inertia values have units of Mass x Length² and they are easily obtained from any 3D CAD software.

It is important to check how your CAD system computes inertia properties[^1]. Considering Figure 19, the input necessary to Giraffe could be computed using the following expressions:

| *J*<sub>11</sub> = ∫<sub>*V*</sub>*ρ*(*κ*<sub>2</sub><sup>2</sup>+*κ*<sub>3</sub><sup>2</sup>)dV  | ( 30 ) |
|--------------------------------------------------------------------------------------------------|--------|
|                                                                                                  |        |
| *J*<sub>22</sub> = ∫<sub>*V*</sub>*ρ*(*κ*<sub>3</sub><sup>2</sup>+*κ*<sub>1</sub><sup>2</sup>)dV  | ( 31 ) |
|                                                                                                  |        |
| *J*<sub>33</sub> = ∫<sub>*V*</sub>*ρ*(*κ*<sub>1</sub><sup>2</sup>+*κ*<sub>2</sub><sup>2</sup>)dV  | ( 32 ) |
|                                                                                                  |        |
| *J*<sub>12</sub> = ∫<sub>*V*</sub>*ρ*(*κ*<sub>1</sub>*κ*<sub>2</sub>)dV                           | ( 33 ) |
|                                                                                                  |        |
| *J*<sub>13</sub> = ∫<sub>*V*</sub>*ρ*(*κ*<sub>1</sub>*κ*<sub>3</sub>)dV                           | ( 34 ) |
|                                                                                                  |        |
| *J*<sub>23</sub> = ∫<sub>*V*</sub>*ρ*(*κ*<sub>2</sub>*κ*<sub>3</sub>)dV                           | ( 35 ) |

With *κ*<sub>1</sub>, *κ*<sub>2</sub> and *κ*<sub>3</sub> being the components of vector **κ** illustrated in Figure 19 and *ρ* is volumetric mass density function of the material considered for the body.

The barycenter position (G in Figure 19) has to be provided in the CAD coordinate system (O xyz). This information is also easily obtained from the 3D CAD software. The origin of the CAD must coincide with the RigidBody\_1 node position.

The last parameter (optional) to define the CADData is the graphic file for postprocessing purposes. It is not used for computing the system physics, which is provided via parameters (mass and inertia tensor entries).

[^1]: Some of the main CAD packages provide the inertia tensor as an output, which means that the products of inertia are shown as -J12, -J13 and -J23 (with negative signs). This is not the correct input for Giraffe. Instead of this, Giraffe expects the values directly obtained by the expressions above (without changing signs).
