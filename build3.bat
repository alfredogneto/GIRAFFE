:: Create directory to build the solution and projects
md build & cd .\build
:: Generate project
cmake .. 
:: Build Giraffe (Debug)
cmake --build . --config Debug

