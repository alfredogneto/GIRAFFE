
# **GIRAFFE**: Generic Interface Readilly Accessible for Finite Elements
<img src="./images/Giraffe.png" width="300">

## Table of contents
- [Introduction](#introduction)
- [Building and compiling](#building-and-compiling)
- [Executing](#executing)
- [Documentation](#documentation)
- [Dependencies](#dependencies)
- [Disclaimer](#disclaimer)

#  Introduction
#  Introduction
**Giraffe** is the acronym of “Generic Interface Readily Accessible for Finite Elements”. It is a platform coded in C++, with the objective of generating a base interface to be used by researchers, to implement their finite element formulations. It was already employed for modeling structural problems, with beam, shell, and rigid body elements. It also has discrete element method capabilities, particularly handling polyhedral elements. **Giraffe** has resources to switch on/off boundary conditions, loads, joints, contacts, etc. straightforwardly. This leads to creating scenarios where load sequence is an issue. Furthermore, it provides numerical strategies to achieve solutions to challenging nonlinear problems. Also, post-processing possibilities are convenient, with an organized set of post files, which is automatically generated for using [Paraview](https://www.paraview.org/). **Giraffe** has a proper input file format, documented in its [user's manual](/documentation/Giraffe%20User's%20Manual%20v2.0.124.pdf).

The platform is coded in a way that permits embracing new elements, novel contact formulations, and proposition of constraint equations, among other features. With that aim, the **Giraffe** code was started in 2014 by [Prof. Alfredo Gay Neto](http://sites.poli.usp.br/p/alfredo.gay/), at the University of Sao Paulo, Brazil.

**Giraffe** was started as a generalization of a previous-developed finite element code, named “FemCable”, which had the objective of simulating offshore structures: risers for oil exploitation. It had implementations of geometric nonlinear beam elements and classical node-to-surface contact formulation. Since a natural expansion required including new contact models, new structural elements, and other resources, **Giraffe** was designed from scratch to have all the models included in “FemCable” and to embrace easy inclusion of new resources, using object orientation programming. **Giraffe** is under continuous development by [Prof. Alfredo Gay Neto](http://sites.poli.usp.br/p/alfredo.gay/), co-workers, and supervised students.

# Building and compiling

You can build the Giraffe project and compile it automatically with the auxiliary batch files here provided.
The provided batch files automatically download prerequisites from public repositories. For that, you need to install Git.
The following instructions can be followed to install **Giraffe** in MS Windows OS:

1. Install the following software:

	[Git](https://git-scm.com/)

	[CMake](https://cmake.org/)

	[Intel oneMKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html)

	[Microsoft Visual Studio](https://visualstudio.microsoft.com/)

  2. Clone [GIRAFFE repository](https://github.com/alfredogneto/GIRAFFE.git). For that, you can use the git command or download it manually using the web browser.
```cmd
git clone https://github.com/alfredogneto/GIRAFFE.git
```
  3. Run the provided `install_dependencies.bat`.
  
  `install_dependencies.bat`: creates a `/dependencies` folder, clones [exprtk](https://github.com/ArashPartow/exprtk.git), [Eigen](https://gitlab.com/libeigen/eigen.git) and [vcpkg](https://github.com/microsoft/vcpkg.git), such as installs the last one. Afer, it installs [arpack-ng](https://github.com/opencollab/arpack-ng.git). It creates environment variables for **GIRAFFE** execution. All variables are created in the user's profile and are based in the folder where Giraffe is located on your computer.

  4. After `install_dependencies.bat` finishes, manually add the following entries to the **user** environment variable `Path`: `%GIRAFFE_PATH%` and `%MKL_PATH%`.  This can be done, for example, using the procedure described [here](https://www.java.com/en/download/help/path.html). PS: If you want to add these paths to the **system** version of `Path`, you have to copy the variables into the system scope too.
  
  5. The last step is to run `build.bat`.

   `build.bat`: sets up the project files for [Microsoft Visual Studio](https://visualstudio.microsoft.com/) and builds **Giraffe** (Release mode by default). Both the project files and executable can be found in `./build` and `./build/Release` respectively.

# Executing 

**Giraffe** execution is based on reading input files, with a proper syntax. Examples of input files are found in the folder `/inputs`. 

To execute **Giraffe**, you simply have to open its executable file either in [Microsoft Visual Studio](https://visualstudio.microsoft.com/) environment or directly in the folder created after compiling it (`/build/Debug`) or (`/build/Release`), depending on the option chosen (default is Release). The only instruction **Giraffe** needs is the name of the input file.

The **Giraffe** input  file must be located inside a folder with the same name as the input file. It is mandatory to use the file extension `*.inp` for the input file. Files with different extensions or with no extensions will result in error messages when **Giraffe** tries to read them.
Example: the input file named `tutorial01.inp` is located inside a folder named `tutorial01`.

The folder with the input file can be located:

  1. In the directory of `Giraffe.exe` executable file;
  2. In the folder `/inputs` located in the installation directory of **Giraffe** software;
  3. In the public `/Documents/Giraffe/` folder.

When trying to read an input file, **Giraffe** seeks it in the sequence here presented. If not succeed, an error message is prompted to the user.

# Documentation

A folder `/documentation` is provided, containing **Giraffe** [user's manual](/documentation/Giraffe%20User's%20Manual%20v2.0.124.pdf) and [tutorials](/documentation/Giraffe%20Tutorials%202.0%20v21.pdf).

# Dependencies

To compile and execute **Giraffe** you need some additional resources:
- [Intel oneMKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html), as part of the [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html). You need to install this library, which will be linked to **Giraffe** via [CMake](https://cmake.org/);
- [exprtk](https://github.com/ArashPartow/exprtk.git), as a library used for mathematical expressions interpretation in Giraffe input files;
- [Eigen](https://gitlab.com/libeigen/eigen.git), as a library used for mathematical operations, matrix organizations, etc.;
- [arpack-ng](https://github.com/opencollab/arpack-ng.git), as a library used to evaluate large-scale eigenvalue problems.

# Disclaimer

   Currently, **Giraffe** is built and tested only with [Microsoft Visual Studio](https://visualstudio.microsoft.com/) and the MS Windows OS.
 We plan to expand it to Linux and Mac soon.
   
