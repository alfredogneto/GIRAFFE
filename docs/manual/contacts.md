# Contacts

*Starts a command block for creation of contact constraints.*

### Syntax:

Contacts N

Name ID data

-   N: number of contacts

-   Name: current contact name

-   ID: current contact identification number

-   data: current contact data (depends on contact resources and requirements)

### Example:

Contacts 2

NSSS 1 NodeSet 1 SurfaceSet 1 MU 0.0 EPN 1e8 CN 0.0 EPT 1e7 CT 0.0 Pinball 1000 Radius 0.0 MaxPointwiseInt 1 BoolTable 1

SSSS 2 SurfaceSet1 1 SurfaceSet2 2 MU 0.0 EPN 1e8 CN 0.0 EPT 1e7 CT 0.0 Pinball 1000 BoolTable 1

### Additional information:

Each contact is defined by a specific keyword followed by the contact identification number (must be an ascending sequence starting from number one) and additional data. Each contact available and its input data is explained next.

Details about Giraffe contact formulations can be found in papers: \[9\] \[10\] \[11\] \[12\] \[13\] \[14\]. A complete explanation on most methods implemented in Giraffe may be found in \[15\].

## NSSS

*Creates a contact constraint for the interaction between a node set and a surface set (NSSS).*

### Syntax:

NSSS CID NodeSet NSID SurfaceSet SSID MU MUV EPN EPNV CN CNV EPT EPTV CT CTV Pinball PV Radius RV MaxPointwiseInt NP BoolTable BTC

-   CID: current contact constraint identification number

-   NSID: node set identification number

-   SSID: surface set identification number

-   MUV: coefficient of friction value

-   EPNV: penalty coefficient to enforce normal contact constraint (no penetration)

-   CNV: normal damping parameter coefficient

-   EPTV: penalty coefficient to enforce tangential contact constraint (sticking condition)

-   CTV: tangential damping parameter coefficient

-   PV: pinball radius value (contact rough searching)

-   RV: sphere radius value surrounding each node in the node set

-   NP: maximum number of contact pointwise interactions between each sphere and surface

-   BTC: bool table data for current contact constraint (see BoolTable)

### Example:

NSSS 1 NodeSet 1 SurfaceSet 1 MU 0.0 EPN 1e8 CN 0.0 EPT 1e7 CT 0.0 Pinball 1000 Radius 0.0 MaxPointwiseInt 1 BoolTable 1

### Additional information:

This contact formulation uses developments detailed in \[13\]. It is an enhanced master-slave, which considers a spherical surface around each node defined in the chosen node set (sphere). Each sphere interacts with surfaces in the chosen surface set, in case of contact occurrence.

Constraints enforcements are done by penalty method. Thus, it is necessary for the user to input penalty parameters data. Usually these may be calibrated based on physical information related to the desired scenario, basing on equivalent local stiffness, leading to allowable penetration on each contact zone.

Damping coefficients are useful for dissipation of energy during impact simulations, avoiding high frequency oscillations on contact forces.

Pinball radius is a rough search geometrical parameter that is used by Giraffe to establish probable and not probable contact interactions. The larger the penalty parameter, the heavier will be the model, since Giraffe will spend time for a larger number of possible contact interactions. However, small pinball radii lead to loosing contact detection. Thus, it is a compromise between accuracy and solution speed. In case the user is in doubt about this parameter, it is better to test it with high values and, afterwards, decrease it.

Usually only a single pointwise contact interaction is permitted between each sphere and each surface. However, some surfaces have the possibility of seeking for more than one pointwise contact solution (typically on non-convexity scenarios).

BoolTable keyword is optional. It permits creating a scenario in which the contact constraint is turned on/off along solution steps (see BoolTable). If the user does not include BoolTable, Giraffe assumes that the contact constraint will be considered for all solution steps as turned on.

## SSSS

*Creates a contact constraint for the interaction between two surface sets (SSSS).*

### Syntax:

//Syntax 1 – a single friction coefficient value

SSSS CID SurfaceSet1 SS1ID SurfaceSet2 SS2ID MU MUV EPN EPNV CN CNV EPT EPTV CT CTV Pinball PV BoolTable BTC

//Syntax 2 – static and dynamic friction coefficient values

SSSS CID SurfaceSet1 SS1ID SurfaceSet2 SS2ID MUS MUSV MUD MUDV EPN EPNV CN CNV EPT EPTV CT CTV Pinball PV BoolTable BTC

//Optional keywords:

WriteReport

-   CID: current contact constraint identification number

-   SS1ID: surface set 1 identification number

-   SS2ID: surface set 2 identification number

-   MUV: coefficient of friction value

-   MUSV: static coefficient of friction value

-   MUDV: dynamic coefficient of friction value

-   EPNV: penalty coefficient to enforce normal contact constraint (no penetration)

-   CNV: normal damping parameter coefficient

-   EPTV: penalty coefficient to enforce tangential contact constraint (sticking condition)

-   CTV: tangential damping parameter coefficient

-   PV: pinball radius value (contact rough searching)

-   BTC: bool table data for current contact constraint (see BoolTable)

-   WriteReport: keyword to instruct Giraffe to produce reports for each local contact problem solved (only use it for debugging purposes because it takes time for writing)

### Example:

SSSS 1 SurfaceSet1 1 SurfaceSet2 2 MU 0.0 EPN 1e8 CN 0.0 EPT 1e7 CT 0.0 Pinball 1000 BoolTable 1

### Additional information:

This contact formulation uses developments detailed in \[10\] and \[11\]. It is master-master contact formulation, which considers interaction between surfaces with no election of slave points.

Constraints enforcements are done by penalty method. Thus, it is necessary for the user to input penalty parameters data. Usually these may be calibrated based on physical information related to the desired scenario, basing on equivalent local stiffness, leading to allowable penetration on each contact zone.

Damping coefficients are useful for dissipation of energy during impact simulations, avoiding high frequency oscillations on contact forces.

Pinball radius is a rough search geometrical parameter that is used by Giraffe to establish probable and not probable contact interactions. The larger the penalty parameter, the heavier will be the model, since Giraffe will spend time for a larger number of possible contact interactions. However, small pinball radii lead to loosing contact detection. Thus, it is a compromise between accuracy and solution speed. In case the user is in doubt about this parameter, it is better to test it with high values and, afterwards, decrease it.

When the user performs a degeneration of surfaces involved in the SSSS contact, Giraffe automatically considers all degenerated cases for contact. With that the user may construct a set of curve/surface or point/surface contact pairs automatically within the same contact creation.

BoolTable keyword is optional. It permits creating a scenario in which the contact constraint is turned on/off along solution steps (see BoolTable). If the user does not include BoolTable, Giraffe assumes that the contact constraint will be considered for all solution steps as turned on.

## SPSP

*Creates a contact constraint for the interaction between two splines (SPSP).*

### Syntax:

//Syntax 1 – a single friction coefficient value

SPSP CID Spline1 SP1ID Spline2 SP2ID MU MUV EPN EPNV CN CNV EPT EPTV CT CTV Pinball PV BoolTable BTC

//Syntax 2 – static and dynamic friction coefficient values

SPSP CID Spline1 SP1ID Spline2 SP2ID MUS MUSV MUD MUDV EPN EPNV CN CNV EPT EPTV CT CTV Pinball PV BoolTable BTC

//Syntax 3 – nonlinear normal contact law

SPSP CID Spline1 SP1ID Spline2 SP2ID MU MUV EPN EPNV EPNN EPNNV CN CNV EPT EPTV CT CTV Pinball PV BoolTable BTC

//Optional keywords:

WriteReport

-   CID: current contact constraint identification number

-   SP1ID: spline 1 identification number

-   SP2ID: spline 2 identification number

-   MUV: coefficient of friction value

-   MUSV: static coefficient of friction value

-   MUDV: dynamic coefficient of friction value

-   EPNV: penalty coefficient to enforce normal contact constraint (no penetration)

-   EPNNV: penalty exponent to enforce nonlinear normal contact constraint (no penetration)

-   CNV: normal damping parameter coefficient

-   EPTV: penalty coefficient to enforce tangential contact constraint (sticking condition)

-   CTV: tangential damping parameter coefficient

-   PV: pinball radius value (contact rough searching)

-   BTC: bool table data for current contact constraint (see BoolTable)

-   WriteReport: keyword to instruct Giraffe to produce reports for each local contact problem solved (only use it for debugging purposes because it takes time for writing)

### Example:

SPSP 1 Spline 1 Spline 2 MU 0.0 EPN 1e8 CN 0.0 EPT 1e7 CT 0.0 Pinball 1000 BoolTable 1

### Additional information:

This contact formulation is based on the formulation presented in \[10\] and \[11\]. It is master-master contact formulation, which considers interaction between spline-based surfaces with no election of slave points.

Constraints enforcements are done by penalty method. Thus, it is necessary for the user to input penalty parameters data. Usually these may be calibrated based on physical information related to the desired scenario, basing on equivalent local stiffness, leading to allowable penetration on each contact zone.

Damping coefficients are useful for dissipation of energy during impact simulations, avoiding high frequency oscillations on contact forces.

Pinball radius is a rough search geometrical parameter that is used by Giraffe to establish probable and not probable contact interactions. The larger the pinball, the heavier will be the model, since Giraffe will spend time for a larger number of possible contact interactions. However, small pinball radii lead to loosing contact detection. Thus, it is a compromise between accuracy and solution speed. In case the user is in doubt about this parameter, it is better to test it with high values and, afterwards, decrease it.

BoolTable keyword is optional. It permits creating a scenario in which the contact constraint is turned on/off along solution steps (see BoolTable). If the user does not include BoolTable, Giraffe assumes that the contact constraint will be considered for all solution steps as turned on.
