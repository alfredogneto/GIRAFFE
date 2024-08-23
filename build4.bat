::Creates environment variables for Giraffe execution
::Giraffe installation path
setx GIRAFFE_INSTALL "%cd%"
::Giraffe dependencies path
setx GIRAFFE_DEPENDENCIES "%cd%/dependencies/vcpkg/packages"
::Giraffe path dependencies appended in a single variable
setx GIRAFFE_PATH %%GIRAFFE_DEPENDENCIES%%/arpack-ng_x64-windows/bin/;%%GIRAFFE_DEPENDENCIES%%/lapack-reference_x64-windows/bin/;%%GIRAFFE_DEPENDENCIES%%/openblas_x64-windows/bin/;%%GIRAFFE_DEPENDENCIES%%/vcpkg-gfortran_x64-windows/bin/;
::MKL and Intel compiler paths
setx MKL_PATH "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/bin/;C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin/;"
::The last step is to include in the system PATH the environment variables: %GIRAFFE_PATH% and %MKL_PATH%
::This is not done automatically by this script routine to avoid changing the OS variables.
