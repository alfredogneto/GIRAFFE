setx GIRAFFE_INSTALL "%cd%\dependencies\vcpkg\packages"
setx GIRAFFE_PATH %%GIRAFFE_INSTALL%%\arpack-ng_x64-windows\bin\;%%GIRAFFE_INSTALL%%\lapack-reference_x64-windows\bin\;%%GIRAFFE_INSTALL%%\openblas_x64-windows\bin\;%%GIRAFFE_INSTALL%%\vcpkg-gfortran_x64-windows\bin;
setx MKL_PATH "C:\Program Files (x86)\Intel\oneAPI\mkl\2024.2\bin"

::Include in the system PATH:
::..\dependencies\vcpkg\packages\arpack-ng_x64-windows\bin\
::..\dependencies\vcpkg\packages\lapack-reference_x64-windows\bin\
::..\dependencies\vcpkg\packages\openblas_x64-windows\bin\
::..\dependencies\vcpkg\packages\vcpkg-gfortran_x64-windows\bin\
::C:\Program Files (x86)\Intel\oneAPI\mkl\2024.2\bin