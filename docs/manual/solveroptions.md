# Solver Options

*Sets solver options (for parallel processing).*

### Syntax:

SolverOptions

Processors NP LinSys ST

-   NP: number of processors (cores) to be used for processing the model

-   ST: solver type for systems of linear equations (“Direct” or “Iterative”).

### Example:

SolverOptions

Processors 4 LinSys Direct

### Additional information:

The SolverOptions keyword is used by Giraffe to set the parallel processing solver options. It rules how Giraffe will use OpenMP<sup>TM</sup> parallel processing routines, which can be really useful for speeding up model solution. Linear systems of equations are solved by PARDISO<sup>TM</sup> library routines.
