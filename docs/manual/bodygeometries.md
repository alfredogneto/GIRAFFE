# Body Geometries

*Starts a command block for creation of body geometries (within a multibody simulation context).*

### Syntax:

BodyGeometries N

BodyGeometry ID Geometries NG List LG

//Input method 1:

BodyGeometry ID Geometries NG List LG

//Input method 2:

BodyGeometry ID Geometries NG Sequence Initial BIN Increment IN

-   N: number of body geometries

-   ID: current body geometry identification number

-   NG: number of individual geometries that define the body geometry

-   LG: list with identification numbers of individual geometries that compose the body geometry

-   BIN: initial identification number of the individual geometry that compose the body geometry

-   IN: increment for the identification numbers of the individual geometries that compose the body geometry

### Example:

BodyGeometries 2

BodyGeometry 1 Geometries 3 List 1 2 3

BodyGeometry 2 Geometries 5 Sequence Initial 4 Increment 1

### Additional information:

Each body geometry is defined by a set of individual geometric entities, using the keywords of the kinds of “Geometries”. Together, such entities compose the whole-body external boundary, where contact can take place. Geometries can also include degenerations into curves/points, to handle singularities.
