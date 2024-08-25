@echo off
::Creation of Giraffe installation path environment variable
setx GIRAFFE_INSTALL "%cd%"
:: Prerequisites:
:: OpenMKL
:: ---------------------------
:: Dependencies
:: ---------------------------
if not exist dependencies (mkdir dependencies)
cd dependencies
::
:: exprtk
::
echo Searching for exprtk...
if not exist exprtk (
echo exprtk folder not found, cloning exprtk...& ^
git clone https://github.com/ArashPartow/exprtk.git
) else (
echo exprtk found.
)
::
:: eigen
::
echo.
echo Searching for eigen...
if not exist eigen (
echo eigen folder not found, cloning eigen...& ^
git clone https://gitlab.com/libeigen/eigen.git
) else (
echo eigen found.
)
::
:: vcpkg
::
echo.
echo Searching for vcpkg...
if not exist vcpkg (
echo vcpkg folder not found, cloning vcpkg...& ^
git clone https://github.com/microsoft/vcpkg.git & ^
echo installing vcpkg...& ^
cd vcpkg & bootstrap-vcpkg.bat & ^
echo vcpkg installed.& ^
echo.& ^
echo installing arpack-ng...& ^
vcpkg install arpack-ng & ^
echo arpack-ng installed.& ^
cd ../..
) else (
echo vcpkg found.& ^
cd vcpkg &^
echo.&^
echo installing arpack-ng...& ^
vcpkg install arpack-ng & ^
echo arpack-ng installed.
)
:: 
:: Set GIRAFFE_PATH
::
echo.
echo Setting GIRAFFE_PATH...
set GIRAFFE_DEPENDENCIES="%cd%/packages"& ^
setx GIRAFFE_PATH "%GIRAFFE_DEPENDENCIES%/arpack-ng_x64-windows/bin/;%GIRAFFE_DEPENDENCIES%/lapack-reference_x64-windows/bin/;%GIRAFFE_DEPENDENCIES%/openblas_x64-windows/bin/;%GIRAFFE_DEPENDENCIES%/vcpkg-gfortran_x64-windows/bin/;"
echo GIRAFFE_PATH set.
:: 
:: Set MKL_PATH
::
echo.
echo Setting MKL_PATH...
setx MKL_PATH "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/bin/;C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin/;"
echo MKL_PATH set.
::
::
::
echo.
echo.
echo Step 1 of the setup is ready
echo Please update your PATH variable (user or system) to contain ^%%GIRAFFE_PATH%% and ^%%MKL_PATH%%.
echo After this, proceed with the CMake build by executing build.bat