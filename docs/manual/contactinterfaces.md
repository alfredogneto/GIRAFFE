# Contact Interfaces

*Starts a command block for creation of contact interfaces.*

### Syntax:

ContactInterfaces N

Name ID Materials MID1 MID2

data

-   N: number of contact interfaces

-   Name: current contact interface name

-   ID: current contact interface identification number

-   MID1 and MID2: material IDs to associate with the contact interface data

-   data: current contact interface data

### Example:

ContactInterfaces 1

Interface\_1 1 Materials 1 1

EPN1 1.0000E+8

N1 2

N2 -2

GNB 2.00E-04

Factor 0.3

ZetaN 8.0000E-01

MUS 0.3

MUD 0.3

EPT 1.00E+03

CT 0

### Additional information:

Each contact interface is defined by a specific keyword followed by its identification number (must be an ascending sequence starting from number one) and additional data. Each contact interface available and its input data is explained next.

## Interface\_1

*Creates a hybrid contact interface for a Barrier-based approach.*

### Syntax:

Interface\_1 IID Materials MID1 MID2

EPN1 EP1V

N1 N1V

N2 N2V

GNB GNBV

Factor FV

ZetaN ZV

MUS MUSV

MUD MUDV

EPT EPTV

CT CTV

Each coefficient is defined next, according to the explanation on the Additional information section.

-   EP1V: value of the coefficient *ϵ*<sub>1</sub> for the evaluation of the normal elastic contact force

-   N1V: value of the coefficient *n*<sub>1</sub> for the evaluation of the normal elastic contact force

-   N2V: value of the coefficient *n*<sub>2</sub> for the evaluation of the normal elastic contact force

-   GNBV: value of the coefficient ${\overline{g}}_{n}$ for the evaluation of the normal elastic contact force

-   FV: value of the factor  f to evaluate $\ {\overline{\overline{g}}}_{n} = f{\overline{g}}_{n}$

-   ZV: value of the normal damping coefficient *ζ*

-   MUSV: value of the static coefficient of friction

-   MUDV: value of the dynamic coefficient of friction

-   EPTV: value of the coefficient *ϵ*<sub>*t*</sub>

-   CTV: value of the coefficient *c*<sub>*t*</sub>

### Example:

Interface\_1 1 Materials 1 1

EPN1 1.0000E+8

N1 2

N2 -2

GNB 2.00E-04

Factor 0.3

ZetaN 8.0000E-01

MUS 0.3

MUD 0.3

EPT 1.00E+03

CT 0

### Additional information:

This interface law is provided to evaluate a contact force, which involves normal and tangential components.

The normal contribution depends on the normal gap evaluation *g*<sub>*n*</sub> and on its time derivative *ġ*<sub>*n*</sub>. This interface law is based on the description given in \[5\], in which a hybrid proposal is given to evaluate the normal elastic contact force *f*<sub>*n*</sub>:

|                                                                                                                                                                                                     |       |
|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------|
| $f_{n} = \left\{ \begin{matrix} 0\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \text{if}\ g_{n} \geq {\overline{g}}_{n} \\ \epsilon_{1}\left( {\overline{g}}_{n} - g_{n} \right)^{n_{1}}\text{\ \ \ \ \ \ \ \ \ \ }\text{if}\ {\overline{\overline{g}}}_{n} &lt; g_{n} &lt; {\overline{g}}_{n} \\ \epsilon_{2}g_{n}^{n_{2}} + c_{2}\text{\ \ \ \ \ \ }\text{if}\ g_{n} &lt; {\overline{\overline{g}}}_{n} \\ \end{matrix} \right.$                                                                                                                                                                          | ( 5 ) |

where the parameters to evaluate *f*<sub>*n*</sub> has to be given, with the exception of *ϵ*<sub>2</sub> and *n*<sub>2</sub> that are evaluated automatically. The reader is invited to have a look at \[5\], where the motivation for this expression is provided, such as examples of parameters for this interface law.

The normal damping force is also evaluated, as following:

|                                                                                                           |       |
|-----------------------------------------------------------------------------------------------------------|-------|
| $f_{d} = 2\zeta\sqrt{- \frac{df_{n}}{dg_{n}}\frac{m_{A}m_{B}}{m_{A} + m_{B}}}{\dot{g}}_{n}$. | ( 6 ) |

where *m*<sub>*A*</sub> and *m*<sub>*B*</sub> are the masses of contacting particles. As the composed contact normal force is given by the elastic + damping contributions, the result can be a compressive or a tractive force. When it gives a tractive contact force, Giraffe saturates the total normal contact force to zero, as here we assume no adhesion effects.

The tangential contact force is based on a classical Coulomb model, as described in \[5\]. To evaluate a trial friction force *f*<sub>*t*</sub><sup>try</sup> (to be tested against the Coulomb limit), we adopt a linear penalty model together with a linear tangential damping contribution, as:

|                                                                                                            |       |
|------------------------------------------------------------------------------------------------------------|-------|
| *f*<sub>*t*</sub><sup>try</sup> = *ϵ*<sub>*t*</sub>*g*<sub>*t*</sub> + *c*<sub>*t*</sub>*ġ*<sub>*t*</sub>, | ( 7 ) |

where *ϵ*<sub>*t*</sub> and *c*<sub>*t*</sub> are respectively the tangential penalty and the tangential damping coefficients.
