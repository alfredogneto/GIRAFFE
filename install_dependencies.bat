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
:: Look if it is already adequately set up
if exist dependencies/vcpkg echo vcpkg found.& goto begininstallarpack
:: If not, see if there is VCPKG_ROOT definition
if defined vcpkg_root goto beginvcpkgdefined
:: If not, check the PATH variable for likely location of the executable
goto begincheckpath

:beginvcpkgdefined
echo VCPKG_ROOT found.
:: Create directory junction
mklink /j "dependencies/vcpkg" "%vcpkg_root%"
goto begininstallarpack

:begincheckpath
:: Write PATH to temporary file
for %%G in ("%path:;=" "%") do @echo %%G >> temp.txt
:: Look for the string vcpkg in it
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
:: Check, for the entries with the string "vcpkg", if they contain the executable
:: If they do, create the junction, if not, clone vcpkg and install it.
for /f %%G in ('findstr vcpkg temp.txt') do if exist %%G/vcpkg.exe (set vcpkgroot=%%G & goto setjunction)
goto beginclonevcpkg

:setjunction
echo vcpkg found in path.
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
setx MKL_PATH "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/bin/;C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin/;C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\redist\intel64_win\compiler"
set MKL_PATH="C:/Program Files (x86)/Intel/oneAPI/mkl/latest/bin/;C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin/;C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\redist\intel64_win\compiler"
echo MKL_PATH set.
::
:: Final instructions
::
echo.
echo.
echo Step 1 of the setup is ready
echo Please update your PATH variable (user or system) to contain ^%%GIRAFFE_PATH%% and ^%%MKL_PATH%%.
echo After this, proceed with the CMake build by executing build.bat
pause