# Surface Sets

*Starts a command block for creation of surface sets.*

### Syntax:

SurfaceSets N

//Input method 1:

SurfaceSet SSID Surfaces NS List S1 S2 …

//Input method 2:

SurfaceSet SSID Surfaces NS Sequence Initial NIS Increment IS

-   N: number of surface sets

-   SSID: current surface set identification number

-   NS: number of surfaces defined in the current surface set

-   S1, S2, …: list with NS surface identification numbers

-   NIS: initial surface identification number

-   IS: increment for the surface identification number

### Example:

SurfaceSets 2

//Input method 1:

SurfaceSet 1 Surfaces 3 List 12 27 21

//Input method 2:

SurfaceSet 2 Surfaces 4 Sequence Initial 3 Increment 2

### Additional information:

Each surface set is defined by the keyword SurfaceSet followed by an identification number (must be an ascending sequence starting from number one).

There are two input methods to define the surface sets. The user must choose one of the following options:

-   List: it indicates to Giraffe that a list of surfaces will be provided as input. For example, the SurfaceSet 1 has three surfaces which are listed after the keyword List;

-   Sequence: it indicates to Giraffe that a sequence of surfaces will be provided. In the example, SurfaceSet 2 has 4 surfaces. The sequence generated automatically will be 3, 5, 7, 9.
