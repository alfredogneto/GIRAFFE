# CAD Data

*Starts a command block for creation of Computer Aided Design (CAD) data.*

### Syntax:

CADData N

Name ID data

-   N: number of CAD data inputs

-   Name: current CAD data type name

-   ID: current CAD data identification number

-   data: current CAD data information (depends on CAD data resources and requirements)

### Example:

CADData 2

NURBSSurface 1 box.txt

STLSurface 2 body.stl

### Additional information:

Each CAD data is defined by a specific keyword followed by the CAD data identification number (must be an ascending sequence starting from number one) and additional data. Each CAD data available and its input data is explained next.

## STLSurface

*Creates a STL (stereolithography) CAD information for usage with other Giraffe resources.*

### Syntax:

STLSurface CID file

-   CID: current CAD data identification number

-   file: *STL* file name (the file must be located inside a folder named “CAD” located in the same directory of Giraffe input file)

### Example:

STLSurface 2 body.stl

### Additional information:

The stl file must be input using the ASCII syntax.

(see [<u>https://en.wikipedia.org/wiki/STL\_(file\_format)</u>](https://en.wikipedia.org/wiki/STL_(file_format)) for details)

**Remark**:

When using the STLSurface input to define a flexible particle (after referred by the VEMPolyhedron command), one has to provide a supplementary file, with the same name as the *STL* file, but with the extension “.msh”. This file has information on a triangular mesh of the surface. One has to define vertices numbering, triangular faces definitions and an internal tetrahedra mesh, for the volume. In this kind of application, we assume that the STL file represents a closed volume.

Example of a file “example.msh”, provided together with a file “example.stl”:

VERTICES

Vertex 0 0.005773502 0.01 -0.015115226

Vertex 1 0.005773502 -0.01 -0.015115226

Vertex 2 0.018683447 0.00 -0.003568220

…

Vertex 11 -0.005773502 -0.01 0.01511522

FACES

Face 0 0 5 4

Face 1 0 3 5

Face 2 0 4 2

…

Face 19 8 10 11

TETRAHEDRA

Tetrahedron 0 0 5 4 3

Tetrahedron 1 0 2 1 4

…

Tetrahedron 15 6 10 11 9

Note that for each vertex, face and tetrahedron a first number is provided as an identification number (starting at zero). For each vertex, one has to define its coordinates. For each face one has to define its connectivity (following the right-hand rule to define the exterior side of the body). For each tetrahedron one has to define its connectivity, following the right-hand rule to define a tetrahedron face, followed by the last vertex in the opposite direction of the external normal defined in the face numbering. An example is depicted in Figure 16.

> <img src="../../assets/manual/image18.png" style="width:2.78932in;height:3.07159in" alt="" />

Figure 16 – Numbering of some vertices to define: face “0” connecting the vertices 0, 5 and 4 and tetrahedron “0” connecting the vertex 0, 5, 4 and 3.

## NURBSSurface

*Creates a NURBS (non uniform rational basis spline) surface CAD information for usage with other Giraffe resources.*

### Syntax:

NURBSSurface CID file

-   CID: current CAD data identification number

-   file: NURBS file name (the file must be located inside a folder named “CAD” located in the same directory of Giraffe input file)

### Example:

STLSurface 2 body.stl

### Additional information:

The input file for a NURBS surface in Giraffe has to follow a specific syntax. It is a text file with information provided as explained next. More information on the theory of NURBS surfaces is well-presented in \[4\], but the basic idea is herein presented next.

First, we have to define a net of control points **P**<sub>*i*, *j*</sub>, where 0 ≤ *i* ≤ *n* and 0 ≤ *j* ≤ *m*. Thus, we consider a total of (*n*+1)(*m*+1) control points. Each one is associated with a weight, given by *w*<sub>*i*, *j*</sub>.

A polynomial basis is constructed to represent points on the surface. Each point is given by a mapping procedure, starting from a parametric plane with coordinates *u* and *v*. The basis is independent for *u* and *v*. The polynomial degrees are also independent and are *p* and *q*, respectively for *u* and *v*. The so-called knot-vectors are non-decreasing sequences of real numbers, organized as **U** and **V** by:

<table><tbody><tr class="odd"><td><blockquote><p><span class="math inline"><strong>U</strong> = [<em>u</em><sub>0</sub>,<em>u</em><sub>1</sub>,…,<em>u</em><sub><em>n</em> + <em>p</em> + 1</sub>]</span> ,</p><p><span class="math inline"><strong>V</strong> = [<em>v</em><sub>0</sub>,<em>v</em><sub>1</sub>,…,<em>v</em><sub><em>m</em> + <em>q</em> + 1</sub>]</span> .</p></blockquote></td><td>(1)</td></tr></tbody></table>

The knot inputs in **U** and **V** are used to establish the polynomial functions *N*<sub>*i*, *p*</sub>(*u*) and *N*<sub>*j*, *q*</sub>(*v*) by de Cox de Boor recursive formula. At the end, the NURBS surface is given by the parameterization:

<table><tbody><tr class="odd"><td><blockquote><p><span class="math inline">$\mathbf{s}\left( u,v \right) = \frac{\sum_{i = 0}^{n}{\sum_{j = 0}^{m}{N_{i,p}\left( u \right)N_{j,q}\left( v \right)w_{i,j}\mathbf{P}_{i,j}}}}{\sum_{i = 0}^{n}{\sum_{j = 0}^{m}{N_{i,p}\left( u \right)N_{j,q}\left( v \right)w_{i,j}}}}$</span></p></blockquote></td><td>(2)</td></tr></tbody></table>

The range for *u* and *v* is defined by the coordinates organized in each knot vector. The control points and weights have direct influence on the surface location in space. NURBS surfaces are not usually interpolatory on control points.

The Giraffe NURBS input file contains all the information needed to establish the parameterization (2). Therefore, one needs to input keywords followed by numeric input data. The keywords are defined as follows:

-   UDim: dimension *n* + 1 (for direction *u*);

-   VDim: dimension *m* + 1 (for direction *m*);

-   UOrder: polynomial degree for the basis along direction *u*;

-   VOrder: polynomial degree for the basis along direction *v*;

-   UKnotVector: **U** knot vector;

-   VKnotVector: **V** knot vector;

-   Weights: weights;

-   ControlPoints: control points.

Weights *w*<sub>*i*, *j*</sub> and control points **P**<sub>*i*, *j*</sub>**=**(*x*<sub>ij</sub>,*y*<sub>ij</sub>,*z*<sub>ij</sub>) are read assuming that they are in a sequence associated with (*u*, *v*) as follows:

(*v*<sub>0</sub>, *u*<sub>0</sub>), (*v*<sub>0</sub>, *u*<sub>1</sub>) … (*v*<sub>0</sub>, *u*<sub>*n* + 1</sub>)

(*v*<sub>1</sub>, *u*<sub>0</sub>), (*v*<sub>1</sub>, *u*<sub>1</sub>) … (*v*<sub>1</sub>, *u*<sub>*n* + 1</sub>)

…

(*v*<sub>*m* + 1</sub>, *u*<sub>0</sub>), (*v*<sub>*m* + 1</sub>, *u*<sub>1</sub>) … (*v*<sub>*m* + 1</sub>, *u*<sub>*n* + 1</sub>)

An example of a simple NURBS input file for Giraffe is:

<table><tbody><tr class="odd"><td><p>UDim</p><p>2</p><p>UOrder</p><p>1</p><p>UKnotVector</p><p>0.0000000000000000e+00</p><p>0.0000000000000000e+00</p><p>1.0000000000000000e+00</p><p>1.0000000000000000e+00</p><p>VDim</p><p>2</p><p>VOrder</p><p>1</p><p>VKnotVector</p><p>0.0000000000000000e+00</p><p>0.0000000000000000e+00</p><p>1.0000000000000000e+00</p><p>1.0000000000000000e+00</p><p>Weights</p><p>1.0000000000000000e+00</p><p>1.0000000000000000e+00</p><p>1.0000000000000000e+00</p><p>1.0000000000000000e+00</p><p>ControlPoints</p><p>-1.0000000000000000e+00 -1.0000000000000000e+00 0.0000000000000000e+00</p><p>1.0000000000000000e+00 -1.0000000000000000e+00 0.0000000000000000e+00</p><p>-1.0000000000000000e+00 1.0000000000000000e+00 0.0000000000000000e+00</p><p>1.0000000000000000e+00 1.0000000000000000e+00 0.0000000000000000e+00</p></td></tr></tbody></table>
