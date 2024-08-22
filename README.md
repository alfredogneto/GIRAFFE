
# GIRAFFE: Generic Interface Readilly Accessible for Finite Elements
<img src="./images/Giraffe.png" width="300">

**Table of contents**
- [Introduction](#introduction)
- [Dependencies](#dependencies)
- [Building and compiling](#building-and-compiling)
- [Executing](#executing)
- [Documentation](#documentation)

#  Introduction
**Giraffe** is the acronym of “Generic Interface Readily Accessible for Finite Elements”. It is a platform coded using C++ language, with the objective of generating a base-interface to be used by researchers, to implement their own finite element formulations. It was already employed for modeling structural problems, with beam and shell elements, rigid bodies and recently it includes discrete element method capabilities.

The platform is coded in a way that permits embracing new elements, new contact formulations, new constraint equations, among other features. With that aim, the **Giraffe** code was started in 2014 by Prof. Alfredo Gay Neto, at University of Sao Paulo, Brazil.

**Giraffe** has started as a generalization of a previous-developed finite element code, named “FemCable”, which had the objective of simulating offshore structures: risers for oil exploitation. It had implementations of geometric nonlinear beam elements and classical node to surface contact formulation. Since a natural expansion required including new contact models, new structural elements and other resources, Giraffe was designed to have all the models included in “FemCable”. Furthermore, it was thought to embrace easy inclusion of new resources, using object orientation programming. Giraffe is under continuous development by Prof. Alfredo Gay Neto and co-workers.

**Giraffe** has resources to straightforwardly switch on/off boundary conditions, loads, joints, contacts, etc. This leads to the possibility of creating scenarios where load sequence is an issue. Furthermore, it provides numerical strategies to achieve solution of challenging nonlinear problems. Also, post-processing possibilities are convenient, with an organized set of post files, which is automatically generated for using [Paraview](https://www.paraview.org/). Giraffe has a proper input file format, documented in its [user's manual](/documentation/Giraffe%20User's%20Manual%20v2.0.124.pdf).

# Dependencies

To compile **Giraffe** one needs:
- [Intel oneMKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html), as part of the [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html). You need to install this library, which will be linked to Giraffe via [CMake](https://cmake.org/);
- [exprtk](https://github.com/ArashPartow/exprtk.git), as a library used for mathematical expressions interpretation in Giraffe input files;
- [Eigen](https://gitlab.com/libeigen/eigen.git), as a library used for mathematical operations, matrix organizations, etc.;
- [arpack-ng](https://github.com/opencollab/arpack-ng.git), as a library used to evaluate large-scale eigenvalue problems.

# Building and compiling

You can build the solution with [CMake](https://cmake.org/) and the auxiliary batch files. Thus, install the CMake before proceeding.
Furthermore, the provided batch files automatically download some prerequisites from public repositories. For that, you need to install [Git](https://git-scm.com/).

Do not forget to install - [Intel oneMKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html) prior to run the batch files to install Giraffe.

The provided batch files are developed to Windows system. They will install the [vcpkg](https://github.com/microsoft/vcpkg.git)package, as an installer for [arpack-ng](https://github.com/opencollab/arpack-ng.git) that automatically provides the necessary files to use it in Windows environment.

In summary, Giraffe installing requires the following steps:

  1. Clone [this repository](https://github.com/alfredogneto/GIRAFFE.git)

```cmd
git clone https://github.com/alfredogneto/GIRAFFE.git
```
  2. Install [Intel oneMKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html)
  3. Run the batch files in the sequence, always waiting for the conclusion of each one prior ro run the next one: `build1.bat`, `build2.bat`, `build3.bat` and `build4.bat`. This will generate and build (Debug mode by default) the solution without opening Visual Studio.

:warning: Developer: the VS debugger is changed by `<root folder>/inputs`. See the [main CMake file](./CMakeLists.txt). Thus, create a folder called "inputs" in the root directory and put your files there. This is done just to avoid unnecessary uploads to Git, but you are free to change the `.gitignore` file in your forked/personal repository, just keep in mind that when you run the Giraffe simulation, the size may be quite large.

# Executing 



# Some disclaimers...

  - IDE and OS supported: 
    
    Currently, **Giraffe** is built and tested only with MSVC ([Microsoft Visual Studio](https://visualstudio.microsoft.com/)) and Windows platform. 
    We plan to expand it to Linux and Mac.
