setx GIRAFFE_INSTALL "%cd%\dependencies\vcpkg\packages"
setx GIRAFFE_PATH %%GIRAFFE_INSTALL%%\arpack-ng_x64-windows\bin\;%%GIRAFFE_INSTALL%%\lapack-reference_x64-windows\bin\;%%GIRAFFE_INSTALL%%\openblas_x64-windows\bin\;%%GIRAFFE_INSTALL%%\vcpkg-gfortran_x64-windows\bin\;
setx MKL_PATH "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\bin\;C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\;"

::The last step is to include in the system PATH the environment variables: %GIRAFFE_PATH% and %MKL_PATH%
::This is not done automatically by this script routine.
