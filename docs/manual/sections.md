# Sections

*Starts a command block for creation of cross-sections (here referenced as “sections”).*

### Syntax:

Sections N

Name ID data

-   N: number of sections

-   Name: current section name

-   ID: current section identification number

-   data: current section data (depends on section resources and requirements)

### Example:

Sections 1

Rectangle 1 B 1.0 H 2.5

### Additional information:

Each section is defined by a specific keyword followed by the section identification number (must be an ascending sequence starting from number one) and additional data. Each section available and its input data is explained next.

Sections are to be used as parameters for elements Beam\_1 and Truss\_1. Note that the local coordinate system used to define the cross section has origin on the cross-section intersection with the beam axis position, defined by the nodes of the mesh. The plane of the cross-section is parallel to the plane formed by **e<sub>1 </sub>**and **e<sub>2</sub>**, defined in the beam element coordinate system (see Beam\_1 input data).

## General

*Creates a general cross-section for usage with beam and truss elements.*

### Syntax:

General SID A AV I11 I11V I22 I22V I12 I12V JT JTV

-   SID: current section identification number

-   AV: cross-section area value

-   I11V: moment of inertia around E1 axis

-   I22V: moment of inertia around E2 axis

-   I12V: product of inertia related to axis E1 and E2

-   JTV: cross-section moment of torsion value

### Example:

General 1 A 0.1 I11 0.01 I22 0.01 I12 0.0 JT 0.02

### Additional information:

It is assumed that the element axis passes through cross-section centroids and shear centers. Thus, if this is not the case, use UserDefined cross-section, instead.

## Rectangle

*Creates a rectangular cross-section for usage with beam and truss elements.*

### Syntax:

Rectangle SID B BV H HV

-   SID: current section identification number

-   BV: base dimension value of the cross-section lying in direction E1

-   HV: height dimension value of the cross-section lying in direction E2

### Example:

Rectangle 1 B 1.0 H 2.5

### Additional information:

It is assumed that the element axis passes through cross-section centroids and shear centers. Thus, if this is not the case, use UserDefined cross-section, instead.

## SuperEllipse

*Creates a super elliptical cross-section for usage with beam and truss elements.*

### Syntax:

SuperEllipse SID A AV B BV N NV AMeshFDM MV

-   SID: current section identification number

-   AV: semi-axis value lying in direction E1

-   BV: semi-axis value lying in direction E2

-   NV: super ellipse exponent value

-   MV: number of divisions in a Finite Difference Method mesh discretization in the direction of radius A (performed to calculate Saint-Venant moment of torsion, prior to FEM simulation)

### Example:

SuperEllipse 1 A 1.0 B 2.0 N 3 AMeshFDM 200

### Additional information:

It is assumed that the element axis passes through cross-section centroids and shear centers. Thus, if this is not the case, use UserDefined cross-section, instead.

## Tube

*Creates a tubular cross-section for usage with beam and truss elements.*

### Syntax:

Tube SID De DEV Di DIV

-   SID: current section identification number

-   DEV: external diameter value

-   DIV: internal diameter value

### Example:

Tube 1 De 0.2 Di 0.1

### Additional information:

It is assumed that the element axis passes through cross-section centroids and shear centers. Thus, if this is not the case, use UserDefined cross-section, instead.

A null internal diameter can be used, and leads to a solid cylinder.

## UserDefined

*Creates a user defined cross-section for usage with beam and truss elements.*

### Syntax:

UserDefined SID

GA GAV

EA EAV

ES1 ES1V

ES2 ES2V

EI11 EI11V

EI22 EI22V

EI12 EI12V

GS1 GS1V

GS2 GS2V

GS1S GS1SV

GS2S GS2SV

GJT GJTV

J11 J11V

J22 J22V

J12 J12V

A AV

SC SC1V SC2V

BC BC1V BC2V

Rho RV

SD SDID

//Optional keywords block below:

AD ADID

AC AC1V AC2V

AeroLength ALV

-   SID: current section identification number

-   GAV: equivalent shear stiffness product

-   EAV: equivalent axial stiffness product

-   ES1V: equivalent *E**S*<sub>1</sub>

-   ES2V: equivalent *E**S*<sub>2</sub>

-   EI11V: equivalent bending stiffness *E**I*<sub>11</sub>

-   EI22V: equivalent bending stiffness *E**I*<sub>22</sub>

-   EI12V: equivalent bending stiffness *E**I*<sub>12</sub>

-   GS1V: equivalent *G**S*<sub>1</sub>

-   GS2V: equivalent *G**S*<sub>2</sub>

-   GS1SV: equivalent *G**S*<sub>1</sub><sup>*s*</sup>

-   GS2SV: equivalent *G**S*<sub>2</sub><sup>*s*</sup>

-   GJTV: equivalent torsion stiffness with respect to origin O (*G**J*<sub>*t*</sub>)

-   J11V: mass Moment of inertia per unit reference length *J*<sub>11</sub>

-   J22V: mass Moment of inertia per unit reference length *J*<sub>22</sub>

-   J12V: mass Product of inertia per unit reference length *J*<sub>12</sub>

-   AV: Cross-section area

-   SC1V and SC2V: shear center coordinates (s<sub>1</sub>, s<sub>2</sub>)

-   BC1V and BCV2: barycenter coordinates (g<sub>1</sub>, g<sub>2</sub>)

-   RV: mass <u>per unit reference length</u> $\mathbf{(}\overline{\rho})$

-   SDID: section details (defines the external contour of the cross section – to be used for rendering purposes)

Optional keywords block:

-   ADID: aerodynamic data identification number (defines aerodynamic curves to be used to evaluate environment wind forces)

-   AC1V and AC2V: aerodynamic center position (c<sub>1</sub>, c<sub>2</sub>)

-   ALV: aerodynamic reference length value (usually cross-section profile chord is employed)

### Example:

Sections 1

UserDefined 1

GA 656713326.769231

EA 1707454649.6

ES1 144398244.5

ES2 0.0

EI11 20026109.28137880

EI22 1200819.338227490

EI12 0.0

GS1 55537786.35

GS2 0.0

GS1S 27768893.1727725

GS2S 0.0

GJT 1373989.013

J11 0.787026095

J22 0.0471921999923403

J12 0.0

A 0.008537273248

SC 0.0 0.04228465

BC 0.0 0.08456930

Rho 67.10296773

SD 1

### Additional information:

For this cross-section, no assumptions are made with respect to the beam axis position. On 3D space, the beam element is defined by its nodes, normally, which are given in section Nodes. The properties to define the beam constitutive behavior will vary according to the chosen position of the axis. Accordingly, the results of internal loads have to be re-interpreted.

Many properties are required to use this cross section, but it permits a myriad of applications, such as composite beams, thin-walled cross-sections and complex-shape cross-sections.

<u>Important: when the UserDefined cross section is used, the material data assigned to the beam element is ignored.</u>

Figure 17 shows an example of UserDefined cross section. The beam axis intersection with the cross-section occurs at point O, origin of a local coordinate system that has to be used to define all quantities defined from now on. As in other Giraffe’s cross-sections for beams, directions **e**<sub>1</sub> and **e**<sub>2</sub> define the cross-section plane. A general material point on the cross-section may be described by coordinates (*x*<sub>1</sub>,*x*<sub>2</sub>) using the system (*O*,**e**<sub>1</sub>,**e**<sub>2</sub>). Let *ρ* be the material specific mass function *ρ* = *ρ̂*(*x*<sub>1</sub>,*x*<sub>2</sub>). The cross-section area domain is given by A.

<img src="../../assets/manual/emf_image19.png" style="max-width:100%" alt="" />

Figure 17 – UserDefined cross-section

The cross-section barycenter is located at material point G, such that

| **b** = (*G*−*O*) = *g*<sub>1</sub>**e**<sub>1</sub> + *g*<sub>2</sub>**e**<sub>2</sub>, | ( 8 ) |
|------------------------------------------------------------------------------------------|-------|
|                                                                                          |       |

with

|                                                                                                                                                                                 |       |
|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------|
| $g_{1} = \frac{\int_{A}^{}{\text{ρx}_{1}\text{dA}}}{\int_{A}^{}\text{ρdA}}$;$\ g_{2} = \frac{\int_{A}^{}{\text{ρx}_{2}\text{dA}}}{\int_{A}^{}\text{ρdA}}$. | ( 9 ) |

The shear center is located at material point S, such that

| **s** = (*S*−*O*) = *s*<sub>1</sub>**e**<sub>1</sub> + *s*<sub>2</sub>**e**<sub>2</sub>. | ( 10 ) |
|------------------------------------------------------------------------------------------|--------|
|                                                                                          |        |

For some applications involving fluid-structure interaction, it is also necessary to define the aerodynamic center, located at material point C, such that

| **c** = (*C*−*O*) = *c*<sub>1</sub>**e**<sub>1</sub> + *c*<sub>2</sub>**e**<sub>2</sub>. | ( 11 ) |
|------------------------------------------------------------------------------------------|--------|
|                                                                                          |        |

Next, one finds the convention used for all geometric properties evaluation for Giraffe. (see more details in \[6\]).

-   Cross-section area:

|                         |        |
|-------------------------|--------|
| *A* = ∫<sub>*A*</sub>dA  | ( 12 ) |

-   Moments of inertia with respect to area:

<table><tbody><tr><td><span class="math display"><em>I</em><sub>11</sub> = ∫<sub><em>A</em></sub><sup></sup><em>x</em><sub>2</sub><sup>2</sup>dA</span></td><td>( 13 )</td></tr><tr><td><span class="math display"><em>I</em><sub>22</sub> = ∫<sub><em>A</em></sub><sup></sup><em>x</em><sub>1</sub><sup>2</sup>dA</span></td><td>( 14 )</td></tr></tbody></table>

-   Product of inertia with respect to area:

<table><tbody><tr><td><span class="math display"><em>I</em><sub>12</sub> =  − ∫<sub><em>A</em></sub><sup></sup><em>x</em><sub>1</sub><em>x</em><sub>2</sub>dA</span></td><td>( 15 )</td></tr></tbody></table>

-   First-order moments (static moments) with respect to area:

| *S*<sub>1</sub> = ∫<sub>*A*</sub>*x*<sub>2</sub>dA     | ( 16 ) |
|-------------------------------------------------------|--------|
| *S*<sub>2</sub> =  − ∫<sub>*A*</sub>*x*<sub>1</sub>dA  | ( 17 ) |

Note: if the specific mass function *ρ* is constant on A, one may write:

| *S*<sub>1</sub> = *A**g*<sub>2</sub>     | ( 18 ) |
|-----------------------------------------|--------|
|                                         |        |
| *S*<sub>2</sub> =  − *A**g*<sub>1</sub>  | ( 19 ) |
|                                         |        |

-   First-order moments (static moments) with respect to the shear center position:

| *S*<sub>1</sub><sup>*s*</sup> = *S*<sub>1</sub> − As<sub>2</sub>   | ( 20 ) |
|-------------------------------------------------------------------|--------|
|                                                                   |        |
| *S*<sub>2</sub><sup>*s*</sup> = *S*<sub>2</sub> + As<sub>1</sub>   | ( 21 ) |
|                                                                   |        |

Note: if the specific mass function *ρ* is constant on A, one may write:

| *S*<sub>1</sub><sup>*s*</sup> = *A*(*g*<sub>2</sub> − *s*<sub>2</sub>)     | ( 22 ) |
|---------------------------------------------------------------------------|--------|
|                                                                           |        |
| *S*<sub>2</sub><sup>*s*</sup> =  − *A*(*g*<sub>1</sub> − *s*<sub>1</sub>)  | ( 23 ) |
|                                                                           |        |

One may also define properties that depend on mass distribution:

-   Mass per unit length

|                                                  |        |
|--------------------------------------------------|--------|
| $\overline{\rho} = \int_{A}^{}\text{ρdA}$  | ( 24 ) |

-   Moments of inertia with respect to mass (per unit reference length of beam, since they are integrated in the area and not in the volume):

| *J*<sub>11</sub> = ∫<sub>*A*</sub>*ρ**x*<sub>2</sub><sup>2</sup>dA  | ( 25 ) |
|--------------------------------------------------------------------|--------|
|                                                                    |        |
| *J*<sub>22</sub> = ∫<sub>*A*</sub>*ρ**x*<sub>1</sub><sup>2</sup>dA  | ( 26 ) |

-   Product of inertia with respect to mass (per unit reference length of beam, since it is integrated in the area and not in the volume):

|                                                                          |        |
|--------------------------------------------------------------------------|--------|
| *J*<sub>12</sub> =  − ∫<sub>*A*</sub>*ρ**x*<sub>1</sub>*x*<sub>2</sub>dA  | ( 27 ) |

The constitutive equation that Giraffe uses to evaluate the beam generalized stress (internal loads) **σ** relation with generalized strains **ϵ** is given by

|                                                                            |        |
|----------------------------------------------------------------------------|--------|
| $\mathbf{\sigma} = \begin{bmatrix} \text{GA} & 0 & 0 & 0 & 0 & G(S_{1}^{S} - S_{1}) \\ 0 & \text{GA} & 0 & 0 & 0 & G(S_{2}^{S} - S_{2}) \\ 0 & 0 & \text{EA} & ES_{1} & ES_{2} & 0 \\ 0 & 0 & ES_{1} & EI_{11} & EI_{12} & 0 \\ 0 & 0 & \text{ES}_{2} & EI_{12} & EI_{22} & 0 \\ G(S_{1}^{S} - S_{1}) & G(S_{2}^{S} - S_{2}) & 0 & 0 & 0 & GJ_{t} \\ \end{bmatrix}\mathbf{\epsilon}$.                                         | ( 28 ) |

All equivalent stiffness coefficients in this relation have to be input for using UserDefined cross-section. These values represent equivalent quantities already integrated on cross-section area. If handling a single material homogeneous cross-section, such stiffness coefficients may be evaluated using all geometric quantities previously defined, multiplied by material data (Young Modulus E and Shear Modulus G). If handling a composite material cross-section, equivalent stiffness coefficients evaluation may be non-straightforward and may be obtained by using a third-part software.

A coefficient needs a special attention: the torsion stiffness *G**J*<sub>*t*</sub>. When disregarding cross-section warping, the moment of torsion is given by *J*<sub>*t*</sub> = *I*<sub>11</sub> + *I*<sub>22</sub>. This gives exact results for circular or tubular cross sections, that do not experience warping under torsion. For all other shapes of cross section, such approximation gives larger values for *J*<sub>*t*</sub> than expected, when considering warping properly.

If one uses Saint-Venant torsion theory, it is possible to enhance evaluation of *J*<sub>*t*</sub> by defining a warping function. With that, one may define the shear-center and the moment of torsion, with respect to the shear-center, named *J*<sub>*t*</sub><sup>*S*</sup>. This value usually may be obtained by a CAD software or other third-part software. However, since the beam axis here considered is general, the moment of torsion for input in Giraffe is not *J*<sub>*t*</sub><sup>*S*</sup>, but may be obtained by considering a proper transport such that

|                                                                                                                     |        |
|---------------------------------------------------------------------------------------------------------------------|--------|
| *J*<sub>*t*</sub> = *J*<sub>*t*</sub><sup>*S*</sup> + *A*(*s*<sub>1</sub><sup>2</sup>+*s*<sub>2</sub><sup>2</sup>).  | ( 29 ) |

Table 6 shows the SI units for all quantities needed for input when using UserDefined cross-section.

Table 6 – SI units for all quantities needed for input in UserDefined cross-section

| **Coefficient**                                   | **Unit** | **Coefficient**                                                  | **Unit** | **Coefficient**              | **Unit** |
|---------------------------------------------------|----------|------------------------------------------------------------------|----------|------------------------------|----------|
| GA                                                | N        | *G**S*<sub>1</sub>, *G**S*<sub>2</sub>                           | N.m      | A                            | m²       |
| EA                                                | N        | *G**S*<sub>1</sub><sup>*s*</sup>, G*S*<sub>2</sub><sup>*s*</sup> | N.m      | s<sub>1</sub>, s<sub>2</sub> | m        |
| ES<sub>1</sub>, ES<sub>2</sub> | N.m | *G**J*<sub>*t*</sub> | N.m² | g<sub>1</sub>, g<sub>2</sub> | m |
| EI<sub>11</sub>, EI<sub>22</sub>, EI<sub>12</sub> | N.m² | *J*<sub>11</sub>, *J*<sub>22</sub>, *J*<sub>12</sub> | Kg.m | $\overline{\rho}$ | kg/m |
