:: Create directory to build the solution and projects
md build & cd .\build
:: Generate project
cmake .. 
:: Build Giraffe (Debug)
cmake --build . --config Debug

::Useful references:
::https://github.com/opencollab/arpack-ng#windows-support
::https://www.dealii.org/developer/external-libs/arpack.html
::https://github.com/opencollab/arpack-ng.git
::https://github.com/microsoft/vcpkg
::https://github.com/microsoft/vcpkg?tab=readme-ov-file
::https://learn.microsoft.com/pt-br/vcpkg/get_started/get-started?pivots=shell-cmd
