# Environment

*Creates environment data.*

### Syntax:

Environment

//Optional block to define gravity data:

GravityData

G GXV GYV GZV BoolTable BDG

//Optional block to define ocean data:

OceanData

RhoFluid RV SurfacePosition XSV YSV ZSV

SeaCurrent N NSC BoolTable BDO

Depth DV Speed SP Angle AV

-   GXV, GYV and GZV: components of gravity field vector

-   BDG: bool table data for gravity loads (see BoolTable)

-   RV: specific mass of ocean water

-   XSV, YSV and ZSV: coordinates of an arbitrary point located on ocean surface

-   NSC: number of points employed to define the sea current velocity field

-   BDO: bool table data for sea current velocity field (see BoolTable)

To define each depth data for sea current:

-   DV: depth value

-   SP: water speed value

-   AV: water speed azimuth angle orientation value (in degrees)

### Example:

Environment

//Optional block to define gravity data:

GravityData

G 0 0 -9.81 BoolTable 1

//Optional block to define ocean data:

OceanData

RhoFluid 1024 SurfacePosition 0 0 1200

SeaCurrent N 5 BoolTable 0 0 1

Depth 0 Speed 1.3 Angle 30

Depth 100 Speed 1.2 Angle 10

Depth 250 Speed 0.5 Angle -30

Depth 500 Speed 0.3 Angle -20

Depth 1200 Speed 0.1 Angle 45

### Additional information:

The Environment keyword is used to define environment data, which is given in blocks: GravityData and OceanData. Each block may be defined in arbitrary sequence after Environment keyword. The presence of all blocks in not mandatory. For example, it is possible to define only the GravityData block if one wants only the effect of gravity in the model. For OceanData block, the sea current is defined using a table with N points, which has to be assigned. OceanData employs buoyancy effect when the element is inside the water (applicable for Beam\_1, Pipe\_1 and Truss\_1 elements)

Both GravityData and OceanData are considered in the simulation according to the BoolTable data (see BoolTable).
