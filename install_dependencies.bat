@echo off
::Creation of Giraffe installation path environment variable
echo Creating GIRAFFE_INSTALL variable.
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
if exist dependencies/vcpkg goto begininstallarpack
if defined vcpkg_root goto beginvcpkgdefined
if not defined vcpkg_root goto begincheckpath

:beginvcpkgdefined
echo VCPKG_ROOT found.
mklink /j "dependencies/vcpkg" "%vcpkg_root%"
goto begininstallarpack

:begincheckpath
for %%G in ("%path:;=" "%") do @echo %%G >> temp.txt
findstr /c:"vcpkg" temp.txt > nul
if %errorlevel%==0 (
goto beginfoundinpath
) else (
goto beginclonevcpkg
)
:endcheckpath
del temp.txt
goto begininstallarpack

:beginfoundinpath
echo vcpkg found in path
for /f %%G in ('findstr vcpkg temp.txt') do set vcpkgroot=%%G & goto setjunction
:setjunction
mklink /j "dependencies/vcpkg" "%vcpkgroot%"
goto endcheckpath

:beginclonevcpkg
echo vcpkg not found, cloning vcpkg...
git clone https://github.com/microsoft/vcpkg.git dependencies/vcpkg
echo installing vcpkg...
cd dependencies/vcpkg
call bootstrap-vcpkg.bat
echo vcpkg installed.

cd ../..
goto endcheckpath

:begininstallarpack
echo.
echo installing arpack-ng...
cd dependencies/vcpkg
vcpkg install arpack-ng
cd ../..
echo arpack-ng installed.
:: 
:: Set GIRAFFE_PATH
::
echo.
echo Setting GIRAFFE_PATH...
setx GIRAFFE_PATH "%GIRAFFE_INSTALL%\dependencies\vcpkg\packages\arpack-ng_x64-windows\bin\;%GIRAFFE_INSTALL%\dependencies\vcpkg\packages\lapack-reference_x64-windows\bin\;%GIRAFFE_INSTALL%\dependencies\vcpkg\packages\openblas_x64-windows\bin\;%GIRAFFE_INSTALL%\dependencies\vcpkg\packages\vcpkg-gfortran_x64-windows\bin;"
set GIRAFFE_PATH="%GIRAFFE_INSTALL%\dependencies\vcpkg\packages\arpack-ng_x64-windows\bin\;%GIRAFFE_INSTALL%\dependencies\vcpkg\packages\lapack-reference_x64-windows\bin\;%GIRAFFE_INSTALL%\dependencies\vcpkg\packages\openblas_x64-windows\bin\;%GIRAFFE_INSTALL%\dependencies\vcpkg\packages\vcpkg-gfortran_x64-windows\bin;"
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







::del temp.txt
