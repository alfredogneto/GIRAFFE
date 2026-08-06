# General information

## Installing Giraffe

Giraffe runs on Windows<sup>TM</sup> 64 bit (the original platform) and, since the *Mac-version* branch, also on macOS.

!!! info "macOS version"
    The macOS port was contributed by Dr. Ajay B Harish (Harish Research Lab, The University of Manchester). For questions about the Mac version, please [open an issue on GitHub](https://github.com/Harish-Research-Lab/GIRAFFE/issues).

=== "Windows (original)"

    Building from source requires [Git](https://git-scm.com/), [CMake](https://cmake.org/), [Intel oneMKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html) and [Microsoft Visual Studio](https://visualstudio.microsoft.com/):

    1.  Clone the repository.

    2.  Run the provided `install_dependencies.bat`. It creates a `/dependencies` folder, clones [exprtk](https://github.com/ArashPartow/exprtk.git), [Eigen](https://gitlab.com/libeigen/eigen.git) and [vcpkg](https://github.com/microsoft/vcpkg.git), installs [arpack-ng](https://github.com/opencollab/arpack-ng.git), and creates the environment variables Giraffe needs.

    3.  Manually add `%GIRAFFE_PATH%` and `%MKL_PATH%` to the **user** environment variable `Path` (Control Panel → System → Advanced System Settings → Environment Variables).

    4.  Run `build.bat`. It generates the Visual Studio project files and builds Giraffe (Release mode by default). The executable is produced in `./build/Release`.

    The original manual also described installing Giraffe from a prebuilt distribution folder ("/Giraffe 2.0/Giraffe Install/", DLL patches "vcredist\_x64" and "ww\_icl\_redist\_intel64\_2016.4.246"); that distribution channel predates this repository.

=== "macOS"

    Building uses [Homebrew](https://brew.sh) instead of vcpkg/oneMKL:

    1.  Install [Homebrew](https://brew.sh), CMake (`brew install cmake`) and Git (ships with the Xcode Command Line Tools).

    2.  Clone the repository (branch *Mac-version*).

    3.  Run `./install_dependencies.sh`. It installs `arpack`, `libomp` and `eigen` via Homebrew, clones `exprtk` into `/dependencies`, and sets `GIRAFFE_INSTALL` in your shell profile. Restart your terminal (or `source` your shell profile) afterwards.

    4.  Run `./build.sh`. The executable is produced at `./build/Giraffe`.

    On macOS, Giraffe links against OpenBLAS/LAPACK and arpack-ng (via Homebrew) instead of Intel oneMKL/PARDISO, and sparse linear systems are solved with Eigen's `SparseLU` instead of PARDISO — direct solves only; the PARDISO "Iterative" solver option has no equivalent on this path. The `PSY.exe` file-coupling feature (`CoupleByFile`) is also unavailable, as `PSY.exe` is a separate Windows-only tool.

**Note**: Giraffe usage is restricted. No user is permitted to supply third-part people with copies of Giraffe with no previous authorization of the developer.

## Running Giraffe

To perform a simulation, launch the Giraffe executable:

=== "Windows (original)"

    Double-click `Giraffe.exe` (in `./build/Release` when built from source), or run it from the Visual Studio environment.

=== "macOS"

    Run the executable from a terminal:

    ```bash
    ./build/Giraffe
    ```

The instruction: “Enter the name of input file” will be given in the screen. Type the desired file name. Then, wait until the simulation finishes. In Figure 1 the file name typed is “beam01”. Do not include the extension of the input file in this typing procedure. If you type “beam01.inp” Giraffe will not find the file. The input-file naming rules and the folders Giraffe searches are described in the [Input file](#input-file) section below.

<img src="../../assets/manual/image2.png" style="width:5.90556in;height:3.27708in" alt="" />

Figure 1 – Giraffe command window

## Input file

Prior to perform a simulation, Giraffe creates a model database with all needed setup. This is done by reading a user input file containing all the commands to construct necessary data for the model, such as nodes, elements, loads, constraints, options for solution, etc. After reading and verifying if input data is consistent, the model is solved.

Giraffe reads a single input text file[^1]. It must be located inside a folder with the same name of the input file. It is mandatory the usage of the file extension “\*.inp” for the input file. Files with different extensions or with no extensions will result in error messages when Giraffe tries to read them.

Figure 2 shows Giraffe directory and some input files folders. For example, the folder “beam01”. The input file named “beam01.inp” is located inside the “beam01” folder (Figure 3).

The folder with the input file can be located in three possible directories:

-   The directory of the Giraffe executable file,

-   The `/inputs` folder located in the installation directory of the Giraffe software, or

-   The public “/Documents/Giraffe/” folder.

When trying to read an input file, Giraffe seeks it in the sequence presented above. If none of the tries succeed, an error message is prompted to the user.

<img src="../../assets/manual/image3.png" style="width:5.90556in;height:2.06111in" alt="" />

Figure 2 – Giraffe executable file directory example

<img src="../../assets/manual/image4.png" style="width:5.90556in;height:2.06389in" alt="" />

Figure 3 – Input file folder example

This user’s manual presents and gives examples of each keyword to be used as a part of Giraffe input file. Giraffe is prepared to read the input file keywords independently of any pre-defined sequence. For example, one may first define the finite element nodes and, after, define elements. Alternatively, one may first define the elements and, after, the nodes. Therefore, the sequence of commands here provided is not mandatory for compounding the input file structure.

To interrupt Giraffe reading process of an input file, the user optionally is allowed to introduce the keyword EOF, indicating “end of file”. This causes Giraffe to read the input file only up to that position. If this keyword is not included in input file, Giraffe will read the whole input file contents.

Explanatory user-comments can be included in some parts of the input file. Comments syntax for Giraffe input file is analogous to C and C++ language, as follows. Single line comments can be included by starting a new line with “//…”. Comments can also span multiple lines, for this purpose one may use “/\*…\*/”. An example of a commented input file is shown below:

```text
//Comment 0
Nodes 5
//Comment 1
Node 1 0 2.5 0
Node 2 0.1 2.5 0
Node 3 0.2 2.5 0
//Comment 2
Node 4 0.3 2.5 0
Node 5 0.4 2.5 0
//Comment 3
CoordinateSystems 2
CS 1 E1 0 1 0 E3 1 0 0
CS 2 E1 1 0 0 E3 0 1 0
Materials 1
Hooke 1 E 2e9 Nu 0.30 Rho 8000
//Comment 4
Sections 2
SuperEllipse 1 A 0.1 B 0.06 N 2 AMeshFDM 100
SuperEllipse 2 A 0.06 B 0.1 N 2 AMeshFDM 100
/*
Large comment
with multiple lines
*/
Elements 2
//Comment 5
Beam_1 1 Mat 1 Sec 1 CS 1 Nodes 1 2 3
//Comment 6
Beam_1 2 Mat 1 Sec 1 CS 1 Nodes 3 4 5
```

Note that it is possible to make comments between first level keywords (e.g.: Nodes, Elements, Sections). It is also possible to make comments between second-level keywords (e.g.: Node, CS, Hooke, Beam\_1).

**Note**: the user cannot introduce comments between lower-level keywords, for example:

Node 1 X /\*not allowed\*/ 0.1 Y 2.5 Z 0 CONSTR 0

Next, an example of Giraffe input file is shown. The reader finds on it a basic structure to establish a simple finite element model of a cantilever beam initially aligned in global Z direction (tutorial01).

/\* Example of an input file for Giraffe \*/

//Creation of nodes

Nodes 11

//Number X Y Z

Node 1 0 0 0.0

Node 2 0 0 0.1

Node 3 0 0 0.2

Node 4 0 0 0.3

Node 5 0 0 0.4

Node 6 0 0 0.5

Node 7 0 0 0.6

Node 8 0 0 0.7

Node 9 0 0 0.8

Node 10 0 0 0.9

Node 11 0 0 1.0

//Creation of node sets

NodeSets 2

NodeSet 1 Nodes 1 List 1

NodeSet 2 Nodes 1 List 11

//Creation of elements

Elements 5

Beam\_1 1 Mat 1 Sec 1 CS 1 Nodes 1 2 3

Beam\_1 2 Mat 1 Sec 1 CS 1 Nodes 3 4 5

Beam\_1 3 Mat 1 Sec 1 CS 1 Nodes 5 6 7

Beam\_1 4 Mat 1 Sec 1 CS 1 Nodes 7 8 9

Beam\_1 5 Mat 1 Sec 1 CS 1 Nodes 9 10 11

//Creation of materials

Materials 1

Hooke 1 E 1e7 Nu 0.3 Rho 2000

//Creation of sections

Sections 1

Rectangle 1 B 0.1 H 0.1

//Creation of coordinate systems

CoordinateSystems 1

CS 1 E1 1 0 0 E3 0 0 1

//Creation of the solution steps

SolutionSteps 1

Static 1

EndTime 1

TimeStep 0.1

MaxTimeStep 0.1

MinTimeStep 0.01

MaxIt 20

MinIt 3

ConvIncrease 4

IncFactor 1.0

Sample 2

//Creation of loads

Loads 1

NodalLoad 1 NodeSet 2 CS 1 NTimes 2

//Time FX FY FZ MX MY MZ

0 0 0 0 0 0 0

1 1000 0 0 0 0 0

//Creation of constraints

Constraints 1

NodalConstraint 1 NodeSet 1

UX BoolTable 1

UY BoolTable 1

UZ BoolTable 1

ROTX BoolTable 1

ROTY BoolTable 1

ROTZ BoolTable 1

//Creation of solver options

SolverOptions

Processors 4 LinSys Direct

//Creation of monitors

Monitor Sample 1

MonitorNodes 1 11

//Creation of post files

PostFiles

MagFactor 1.0

WriteMesh 1

WriteRenderMesh 1

WriteRigidContactSurfaces 0

WriteFlexibleContactSurfaces 0

WriteForces 0

WriteConstraints 0

WriteSpecialConstraints 0

WriteContactForces 0

WriteRenderRigidBodies 0

WriteRenderParticles 0

## Output files

During solution process Giraffe automatically saves all requested output files. They are saved inside the folder where the input file was read.

<img src="../../assets/manual/image5.png" style="width:5.90556in;height:2.13681in" alt="" />

Figure 4 – Example of a simulation folder (from file “beam01.inp”)

Each output file type is described next.

### File “output.inp” 

It is a text file that simply reflects the read information from the input file and may be used to check if Giraffe read input file correctly in some cases.

### File “simulation\_report.txt” 

It is a text file that simply reflects the Giraffe screen output, showing convergence history of all time-steps of simulation.

### Post files

Folder “/post” is always created after solving a simulation. Inside it, Giraffe creates sub-folders with outputs for each solution step. Figure 5 shows an example of “/post” folder contents for the example “beam01.inp”. Note that in this case a single folder “/post/solution\_1” was created. This is because a single solution step was requested in the input file. In case the user requests more solution steps, more folders will be automatically created inside the “/post” folder. A concomitant solution also creates a sub-folder on “/post”.

Note that there are “.pvd” files inside the “/post” folder. These are to be used together with Paraview<sup>TM</sup> post processor. Paraview<sup>TM</sup> “.pvd” files creates links to other files located inside solution folders, which contain the simulation results established by PostFiles keyword. Depending on which results are requested, more or less “.pvd” files will be saved. For this example, only mesh and render mesh were requested, which lead for only two types of “.pvd” files. Paraview<sup>TM</sup> “.pvd” files are very useful for creating animations or high-quality images. These files links Paraview<sup>TM</sup> to read the whole time series of subsequent node positions, node and element results. Inside PARAVIEW<sup>TM</sup><sub>,</sub> when referring to element results meaning, one may look at each element results sequence list, contained in each element presented on this user’s manual.

The “whole\_solution\_mesh.pvd” and “whole\_solution\_rendermesh.pvd” files contains the same contents of the solutions “.pvd” files, but encompassing all solution steps. These are useful for a visualization in Paraview<sup>TM</sup> of the whole time-history of the simulation. For current example, since a single solution step was requested, “.pvd” files for solution 1 and whole “.pvd” files will be the same.

<img src="../../assets/manual/image6.png" style="width:5.90556in;height:2.13681in" alt="" />

Figure 5 – “/post” folder example

Also, inside the “/post/solution\_i” folder (for i = 1,…,n – solution steps), Giraffe writes configuration text files, with nodes and elements results, adopting the same sampling employed for saving Paraview<sup>TM</sup> post-processing files. An example of inside contents of “/post/solution\_1” is shown in Figure 6.

<img src="../../assets/manual/image7.png" style="width:5.90556in;height:3.52083in" alt="" />

Figure 6 – “/post/solution\_1” folder

Note that in “/post/solution\_1” there are also “.vtu” files, which are to be read by Paraview<sup>TM</sup>, being referred in already mentioned “.pvd” files.

### Monitor files

When requested, monitors are very useful for analyzing and creating time series with information about a given node, element, contact region or node set of interest. Monitors are created from the beginning of the simulation. They are updated until the simulation finishes, during the whole solution process. The Figure 7 shows an example of “monitors” folder.

<img src="../../assets/manual/image8.png" style="width:5.90556in;height:1.67986in" alt="" />

Figure 7 – Monitors folder example

## Tutorials

This document has no tutorials. A specific tutorials document containing examples of Giraffe input files is available.

[^1]: Some exceptions are treated locally in this manual, when extra inputs are to be provided by the user.
