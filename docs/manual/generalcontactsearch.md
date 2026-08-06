# General Contact Search

*Creates the automatic contact search in Giraffe (available for particle-particle, particle-boundary, body-body and particle-body interactions)*

### Syntax:

GeneralContactSearch Method MNAME

//Optional keywords block below:

BVFactor BVV

Domain

XMINV XMAXV

YMINV YMAXV

ZMINV ZMAXV

BoolTable BTD

-   MNAME: name of the method employed on the automatic contact search. The available options are “AlltoAll”, “LinkedCells”, “Verlet” or “VerletLinkedCells”.

-   BVV: bounding volume inflation factor. Default value is 0.1

-   XMINV: minimum value for the coordinate x global considered for contact search

-   XMAXV: maximum value for the coordinate x global considered for contact search

-   YMINV: minimum value for the coordinate y global considered for contact search

-   YMAXV: maximum value for the coordinate y global considered for contact search

-   ZMINV: minimum value for the coordinate z global considered for contact search

-   ZMAXV: maximum value for the coordinate z global considered for contact search

-   BTD: bool table data to activate/inactivate the contact search as desired (see BoolTable section).

### Example:

GeneralContactSearch Method Verlet

BVFactor 0.1

Domain

-100 100

-100 100

-100 100

### Additional information:

The general contact search technique used by Giraffe to handle particle-particle and particle-boundary interactions is explained in \[5\]. The basic idea is the creation of inflated bounding volumes surrounding geometrical entities of particles and, as the simulation evolves, search for bounding volume overlaps. With that a list of probable contact candidates is created and investigated, accordingly. The current approaches for evaluating the bounding volume overlap search are:

-   AlltoAll: all the particles are checked against the others

-   Verlet: a list of proximity is constructed first to perform the bounding volume overlap search (known as Verlet list)

-   LinkedCells: a spatial division is made, and the particles are associated with cells. Based on that, one creates a more efficient searching scheme than the all-to-all search

-   VerletLinkedCells: the combination of both previous approaches
