# Materials

*Starts a command block for creation of materials.*

### Syntax:

Materials N

Name ID data

-   N: number of materials

-   Name: current material name

-   ID: current material identification number

-   data: current material data (depends on material resources and requirements)

### Example:

Materials 1

Hooke 1 E 210E9 Nu 0.3 Rho 7800

### Additional information:

Each material is defined by a specific keyword followed by the material identification number (must be an ascending sequence starting from number one) and additional data. Each material available and its input data is explained next.

## Hooke

*Creates a linear-elastic material (Hooke’s law)*

### Syntax:

Hooke MID E EV Nu NV Rho RV

-   MID: current material identification number

-   EV: Young’s Modulus value

-   NV: Poisson’s ratio value

-   RV: specific mass value

### Example:

Hooke 1 E 210E9 Nu 0.3 Rho 7800

### Additional information:

Even defining a Hooke material behavior, each element formulation makes use of different techniques to mount the constitutive equation. Different strain energy functions can be used. More details on the constitutive equation assumed for each element formulation can be found in the papers referenced in this manual, such as \[3\] and \[1\].

## ElasticPlasticIsoHardening

*Creates an elastic-plastic material with isotropic hardening rule.*

### Syntax:

ElasticPlasticIsoHardening MID E EV Nu NV Rho RV H HV YieldingStrength YSV

-   MID: current material identification number

-   EV: Young’s Modulus value

-   NV: Poisson’s ratio value

-   RV: specific mass value

-   HV: linear hardening slope value

-   YSV: yielding strength value

### Example:

ElasticPlasticIsoHardening 1 E 210000 Nu 0.3 Rho 8E-9 H 10000 YieldingStrength 250

### Additional information:

This material model is not available for all elements. Please, check availability for the element of interest prior to usage.

## Orthotropic

*Creates an orthotropic material*

### Syntax:

Orthotropic MID E1 E1V E2 E2V G12 G12V G23 G23V Nu12 N12V Rho RV

-   MID: current material identification number

-   E1V: Young’s Modulus value at direction principal direction 1

-   E2V: Young’s Modulus value at direction principal direction 2

-   G12V: Shear Modulus value associated with directions 1 and 2

-   G23V: Shear Modulus value associated with directions 2 and 3

-   N12V: Poisson’s ratio value associated with directions 1 and 2

-   RV: specific mass value

### Example:

Orthotropic 1 E1 41e9 E2 10.4e9 G12 4.3e9 G23 4.3e9 Nu12 0.28 Rho 1970

### Additional information:

This material model is not available for all elements. Please, check availability for the element of interest prior to usage. Orthotropic material constants orientation as defined in Figure 15.

<table><thead><tr class="header"><th><img src="../../assets/manual/image16.png" style="width:1.63724in;height:1.33766in" alt="" /></th><th><img src="../../assets/manual/image17.png" style="width:1.38312in;height:0.73562in" alt="" /></th></tr></thead><tbody><tr class="odd"><td><ol type="a"><li><p>Directions 1 and 2</p></li></ol></td><td><ol start="2" type="a"><li><p>Direction 3</p></li></ol></td></tr><tr class="even"><td>Figure 15 – Orthotropic material orientation</td><td></td></tr></tbody></table>
