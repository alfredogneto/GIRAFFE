::Pre requisits: 
::OpenMKL

:: Window title
title Generate and build Giraffe
::Dependencies
md dependencies & cd .\dependencies
::exprtk
git clone https://github.com/ArashPartow/exprtk.git
::eigen
git clone https://gitlab.com/libeigen/eigen.git
::vcpkg
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg && bootstrap-vcpkg.bat
