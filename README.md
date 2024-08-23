
# **GIRAFFE**: Generic Interface Readilly Accessible for Finite Elements
<img src="./images/Giraffe.png" width="300">

## Table of contents
- [Introduction](#introduction)
- [Dependencies](#dependencies)
- [Building and compiling](#building-and-compiling)
- [Executing](#executing)
- [Documentation](#documentation)
- [Disclaimer](#disclaimer)

#  Introduction
**Giraffe** is the acronym of “Generic Interface Readily Accessible for Finite Elements”. It is a platform coded in C++, with the objective of generating a base-interface to be used by researchers, to implement their own finite element formulations. It was already employed for modeling structural problems, with beam, shell and rigid body elements. It also has discrete element method capabilities, particularly handling polyhedral elements. **Giraffe** has resources to straightforwardly switch on/off boundary conditions, loads, joints, contacts, etc. This leads to the possibility of creating scenarios where load sequence is an issue. Furthermore, it provides numerical strategies to achieve solution of challenging nonlinear problems. Also, post-processing possibilities are convenient, with an organized set of post files, which is automatically generated for using [Paraview](https://www.paraview.org/). **Giraffe** has a proper input file format, documented in its [user's manual](/documentation/Giraffe%20User's%20Manual%20v2.0.124.pdf).

The platform is coded in a way that permits embracing new elements, new contact formulations, new constraint equations, among other features. With that aim, the **Giraffe** code was started in 2014 by [Prof. Alfredo Gay Neto](http://sites.poli.usp.br/p/alfredo.gay/), at the University of Sao Paulo, Brazil.

**Giraffe** was started as a generalization of a previous-developed finite element code, named “FemCable”, which had the objective of simulating offshore structures: risers for oil exploitation. It had implementations of geometric nonlinear beam elements and classical node-to-surface contact formulation. Since a natural expansion required including new contact models, new structural elements and other resources, **Giraffe** was designed from scratch to have all the models included in “FemCable” and to embrace easy inclusion of new resources, using object orientation programming. **Giraffe** is under continuous development by [Prof. Alfredo Gay Neto](http://sites.poli.usp.br/p/alfredo.gay/), co-workers and oriented students.

# Dependencies

To compile and execute **Giraffe** you need some additional resources:
- [Intel oneMKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html), as part of the [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html). You need to install this library, which will be linked to **Giraffe** via [CMake](https://cmake.org/);
- [exprtk](https://github.com/ArashPartow/exprtk.git), as a library used for mathematical expressions interpretation in Giraffe input files;
- [Eigen](https://gitlab.com/libeigen/eigen.git), as a library used for mathematical operations, matrix organizations, etc.;
- [arpack-ng](https://github.com/opencollab/arpack-ng.git), as a library used to evaluate large-scale eigenvalue problems.

# Building and compiling

You can build the solution with [CMake](https://cmake.org/) and the auxiliary batch files here provided.
Furthermore, the provided batch files automatically download some prerequisites from public repositories. For that, you need to install [Git](https://git-scm.com/).

Do not forget to install - [Intel oneMKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html) prior to run the batch files to install **Giraffe**.

The provided batch files are developed to MS Windows system. They will automatically install the [vcpkg](https://github.com/microsoft/vcpkg.git) package, as an installer for [arpack-ng](https://github.com/opencollab/arpack-ng.git) that automatically provides the necessary files to use it in MS Windows environment.

In summary, **Giraffe** installing requires the following steps:

  1. Install [CMake](https://cmake.org/) and [Git](https://git-scm.com/) software.
  2. Clone [GIRAFFE repository](https://github.com/alfredogneto/GIRAFFE.git). For that, you can use git command or download it manually using the web browser.
```cmd
git clone https://github.com/alfredogneto/GIRAFFE.git
```
  3. Install [Intel oneMKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html) library.
  4. Run the here provided batch files in the sequence, always waiting for the conclusion of each one prior ro run the next one. Next we describe what each file does:
  
  `build1.bat`: clones [exprtk](https://github.com/ArashPartow/exprtk.git), [Eigen](https://gitlab.com/libeigen/eigen.git) and [vcpkg](https://github.com/microsoft/vcpkg.git), such as installs the last one.

  `build2.bat`: using [vcpkg](https://github.com/microsoft/vcpkg.git), installs [arpack-ng](https://github.com/opencollab/arpack-ng.git).

  `build3.bat`: builds and compiles **Giraffe** for [Microsoft Visual Studio](https://visualstudio.microsoft.com/) (Debug mode by default).

  `build4.bat`: creates environment variables for **GIRAFFE** execution. All variables are created in the user's profile and are based in the folder where Giraffe lies in your computer.

  5. The last step is to manually add to the environment variable `Path` of your operational system the following paths: `%GIRAFFE_PATH%` and `%MKL_PATH%`. This can be done, for example, using the procedure described [here](https://www.java.com/en/download/help/path.html). PS: you can add these paths either to system variables or to variables for the user. Both work fine. Take care not to overwrite some of the available paths defined in your system/user, to avoid problems.

After finishing all the steps, you will find the [Microsoft Visual Studio](https://visualstudio.microsoft.com/) project created in the folder `/build`, such as another files generated automatically during installation.

# Executing 

**Giraffe** execution is based on reading input files, with a proper syntax. Examples of input files are found in the folder `/inputs`. 

To execute **Giraffe**, you just have to open its executable file either in [Microsoft Visual Studio](https://visualstudio.microsoft.com/) environment or directly in the folder created after compiling it (`/build/Debug`) or (`/build/Release`), depending on the option chosen (default is Debug). The only instruction **Giraffe** needs is the name of the input file.

The **Giraffe** input  file must be located inside a folder with the same name of the input file. It is mandatory the usage of the file extension `*.inp` for the input file. Files with different extensions or with no extensions will result in error messages when **Giraffe** tries to read them.
Example: the input file named `tutorial01.inp` is located inside a folder named `tutorial01`.

The folder with the input file can be located in the possible locations:

  1. In the directory of `Giraffe.exe` executable file;
  2. In the folder `/inputs` located in the installation directory of **Giraffe** software;
  3. In the public `/Documents/Giraffe/` folder.

When trying to read an input file, **Giraffe** seeks for it in the sequence here presented. If not succeed, an error message is prompted to the user.

# Documentation

A folder `/documentation` is provided, containing **Giraffe** [user's manual](/documentation/Giraffe%20User's%20Manual%20v2.0.124.pdf) and [tutorials](/documentation/Giraffe%20Tutorials%202.0%20v21.pdf).

# Disclaimer

   Currently, **Giraffe** is built and tested only with [Microsoft Visual Studio](https://visualstudio.microsoft.com/) and Windows platform.
 We plan to expand it to Linux and Mac in a near future.
   
