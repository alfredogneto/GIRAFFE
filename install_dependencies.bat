@echo off
::Creation of Giraffe installation path environment variable
setx GIRAFFE_INSTALL "%cd%"
set GIRAFFE_INSTALL="%cd%"
:: Prerequisites:
:: OpenMKL
:: ---------------------------
:: Dependencies
:: ---------------------------
::
:: exprtk
::
echo Searching for exprtk...
if not exist dependencies/exprtk (
echo exprtk not found, cloning exprtk...& ^
git clone https://github.com/ArashPartow/exprtk.git dependencies/exprtk
) else (
echo exprtk found.
)
::
:: eigen
::
echo.
echo Searching for eigen...
if not exist dependencies/eigen (
echo eigen not found, cloning eigen...& ^
git clone https://gitlab.com/libeigen/eigen.git dependencies/eigen
) else (
echo eigen found.
)
::
:: vcpkg
::
echo.
echo Searching for vcpkg...
IF NOT DEFINED vcpkg_root (
if not exist dependencies/vcpkg (
echo vcpkg not found, cloning vcpkg...& ^
git clone https://github.com/microsoft/vcpkg.git dependencies/vcpkg & ^
echo installing vcpkg...& ^
cd depedencies/vcpkg & bootstrap-vcpkg.bat & ^
echo vcpkg installed.& ^
echo.&^
echo installing arpack-ng...& ^
vcpkg install arpack-ng & ^
echo arpack-ng installed.& ^
cd ../..
)) else (
echo vcpkg found.& ^
mklink /j "dependencies/vcpkg" "%vcpkg_root%"&^
echo.&^
echo installing arpack-ng...& ^
%vcpkg_root%/vcpkg install arpack-ng & ^
echo arpack-ng installed.
)
:: 
:: Set GIRAFFE_PATH
::
echo.
echo Setting GIRAFFE_PATH...

setx GIRAFFE_PATH "%GIRAFFE_INSTALL%\dependencies\vcpkg\packages\arpack-ng_x64-windows\bin\;%GIRAFFE_INSTALL%\dependencies\vcpkg\packages\lapack-reference_x64-windows\bin\;%GIRAFFE_INSTALL%\dependencies\vcpkg\packages\openblas_x64-windows\bin\;%GIRAFFE_INSTALL%\dependencies\vcpkg\packages\vcpkg-gfortran_x64-windows\bin;"
set GIRAFFE_PATH="%GIRAFFE_INSTALL%\dependencies\vcpkg\packages\arpack-ng_x64-windows\bin\;%GIRAFFE_INSTALL%\dependencies\vcpkg\packages\lapack-reference_x64-windows\bin\;%GIRAFFE_INSTALL%\dependencies\vcpkg\packages\openblas_x64-windows\bin\;%GIRAFFE_INSTALL%\dependencies\vcpkg\packages\vcpkg-gfortran_x64-windows\bin;"

::set GIRAFFE_DEPENDENCIES="%cd%/packages"
::setx GIRAFFE_PATH "%GIRAFFE_DEPENDENCIES%/arpack-ng_x64-windows/bin/;%GIRAFFE_DEPENDENCIES%/lapack-reference_x64-windows/bin/;%GIRAFFE_DEPENDENCIES%/openblas_x64-windows/bin/;%GIRAFFE_DEPENDENCIES%/vcpkg-gfortran_x64-windows/bin/;"
::set GIRAFFE_PATH="%%GIRAFFE_DEPENDENCIES%%/arpack-ng_x64-windows/bin/;%%GIRAFFE_DEPENDENCIES%%/lapack-reference_x64-windows/bin/;%%GIRAFFE_DEPENDENCIES%%/openblas_x64-windows/bin/;%%GIRAFFE_DEPENDENCIES%%/vcpkg-gfortran_x64-windows/bin/;"
echo GIRAFFE_PATH set.
:: 
:: Set MKL_PATH
::
echo.
echo Setting MKL_PATH...
setx MKL_PATH "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/bin/;C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin/;"
set MKL_PATH="C:/Program Files (x86)/Intel/oneAPI/mkl/latest/bin/;C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin/;"
echo MKL_PATH set.
::
::
::
echo.
echo.
echo Step 1 of the setup is ready
echo Please update your PATH variable (user or system) to contain ^%%GIRAFFE_PATH%% and ^%%MKL_PATH%%.
echo After this, proceed with the CMake build by executing build.bat