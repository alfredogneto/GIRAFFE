# Solution Steps

*Establishes a sequence of solution steps to be solved by Giraffe.*

### Syntax:

SolutionSteps N

Name ID data

-   N: number of solution steps

-   Name: current solution step name

-   ID: current solution step identification number

-   data: current solution step data (depends on solution step resources and requirements)

### Example:

SolutionSteps 2

Static 1

EndTime 2

TimeStep 0.1

MaxTimeStep 0.2

MinTimeStep 0.01

MaxIt 12

MinIt 3

ConvIncrease 2

IncFactor 1.2

Sample 1

Dynamic 2

EndTime 3

TimeStep 0.1

MaxTimeStep 0.2

MinTimeStep 0.01

MaxIt 12

MinIt 3

ConvIncrease 2

IncFactor 1.2

Sample 2

RayleighDamping Alpha 0 Beta 0 Update 0

NewmarkCoefficients Beta 0.3 Gamma 0.5

### Additional information:

Each solution step is defined by a specific keyword followed by the solution step identification number (must be an ascending sequence starting from number one) and additional data. Each solution step available and its input data is explained next.

Solution steps are used to create a sequence of solutions. There is a global “time” tracking parameter to rule all solution steps. The default start-time is zero. Then, each solution step has a definition of end-time. Note that the end-time of a given solution step must be larger than the end-time of the previous solution step, otherwise Giraffe will prompt an error message prior to solution start. Exceptions are modal analysis solution steps that have no end-time parameter as input. In this case, time is considered frozen during modal analysis.

Multiple solution steps may be created to establish starting/ending of constraint actions, contacts, special constraints, or even to split between statics and dynamics, according to the nonlinear model convenience.

The final converged model configuration at the end of a solution step is always taken as the start point for the next solution step. When modal analysis is performed, no changes exist on the model configuration for a next solution step.

Constraints, special constraints and contacts are considered according to the defined BoolTable in each particular creation of such resources. Definition of loads or displacements is made for each solution step following the time variable as a global tracking.

After simulation is finished, the user will find requested result files for each solution step in separate folders: “/post/solution\_i/”, where “i” is the solution identification number.

## Static

*Creates a solution step to solve a nonlinear static analysis.*

### Syntax:

Static SID

EndTime EV

TimeStep TS

MaxTimeStep MAX

MinTimeStep MIN

MaxIt MAXIT

MinIt MINIT

ConvIncrease CONV

IncFactor INCF

Sample SA

-   SID: current solution step identification number

-   EV: end time of current solution step

-   TS: time-step of current solution step

-   MAX: maximum time-step of current solution step

-   MIN: minimum time-step of current solution step

-   MAXIT: maximum number of iterations to be performed during Newton-Raphson routine, for each time-step within the current solution step

-   MINIT: minimum number of iterations, to indicate convergence easiness

-   CONV: number of sequential converged time-steps to indicate convergence easiness

-   INCF: time-step increasing factor

-   SA: sampling variable to rule post-processing files generation.

### Example:

Static 1

EndTime 2.5

TimeStep 0.10

MaxTimeStep 0.25

MinTimeStep 0.01

MaxIt 15

MinIt 3

ConvIncrease 4

IncFactor 1.5

Sample 1

### Additional information:

A static solution step is defined by a sequence of attributes, with the objective of creating an auto-adaptive scheme for nonlinear solution, capable of increasing or decreasing the time-step automatically. In the context of a static analysis, the time variable may be understood as a scalar tracking parameter that permits evaluation of a sequence of loads, constraints and boundary conditions. The time-period of a given static solution step goes from the previously converged time value until the end time of the current solution step, set by the user. Depending on the easiness or hardness of convergence, time-step may be updated automatically according to the parameters, explained below:

-   EndTime: defines the final instant for the current solution step.

-   TimeStep: defines the initial time step to be used for the evolution of the nonlinear model.

-   MaxTimeStep: defines the maximum time step to be used for the evolution of the nonlinear model.

-   MinTimeStep: defines the minimum time step to be used for the evolution of the nonlinear model. In case of convergence difficulties, Giraffe automatically performs bisections (i.e., decreases automatically the time-step). In case of many unsuccessful bisections, when the minimum time step is approached, the simulation stops with an error message.

-   MaxIt: defines the maximum number of iterations, which will be performed prior to assume that divergence occurred. There are two possibilities of divergence: by achieving the maximum number of iterations or by achieving a very high residual value, defined by ConvergenceCriteria keyword.

-   MinIt: defines the minimum number of iterations. Once a time step converges with less or equal to this number of iterations, the next time step will be increased by the IncFactor coefficient. Thus, it may be seen as an identifier of easy solution.

-   ConvIncrease: defines the sequential number of converged time steps that, once achieved, will also show that there is the possibility of increasing the time step value. Then, IncFactor coefficient is used to increase the time-step. Once applied, the converged partial solutions counting process re-starts.

-   IncFactor: defines the factor for increasing the time step, in case of easy convergence.

-   Sample: defines the sampling for saving post-processing files with partial converged solutions along time evolution. Entering the number “1” claims Giraffe to save all converged steps (many files can be generated).

## Dynamic

*Creates a solution step to solve a nonlinear dynamic analysis (transient dynamics).*

### Syntax:

Dynamic SID

EndTime EV

TimeStep TS

MaxTimeStep MAX

MinTimeStep MIN

MaxIt MAXIT

MinIt MINIT

ConvIncrease CONV

IncFactor INCF

Sample SA

RayleighDamping Alpha AD Beta BD Update UD

NewmarkCoefficients Beta BN Gamma GN

//Optional keyword to make null all the initial conditions stemming from previous load steps:

ZeroIC

-   SID: current solution step identification number

-   EV: end time of current solution step

-   TS: time-step of current solution step

-   MAX: maximum time-step of current solution step

-   MIN: minimum time-step of current solution step

-   MAXIT: maximum number of iterations to be performed during Newton-Raphson routine, for each time-step within the current solution step

-   MINIT: minimum number of iterations, to indicate convergence easiness

-   CONV: number of sequential converged time-steps to indicate convergence easiness

-   INCF: time-step increasing factor

-   SA: sampling variable to rule post-processing files generation.

-   AD: coefficient that multiplies mass matrix for Rayleigh damping evaluation

-   BD: coefficient that multiplies stiffness matrix for Rayleigh damping evaluation

-   UD: flag to update (1) or not (0) the Rayleigh damping matrix in each time step beginning

-   BN and GN: Newmark time-integrator β and γ parameters

-   ZeroIC: makes null all the initial conditions stemming from previous load steps

### Example:

Dynamic 1

EndTime 2.5

TimeStep 0.10

MaxTimeStep 0.25

MinTimeStep 0.01

MaxIt 15

MinIt 3

ConvIncrease 4

IncFactor 1.5

Sample 1

RayleighDamping Alpha 0.0 Beta 0.0 Update 0

NewmarkCoefficients Beta 0.3 Gamma 0.5

### Additional information:

A dynamic solution step is defined by a sequence of attributes, with the objective of creating an auto-adaptive scheme for nonlinear solution, capable of increasing or decreasing the time-step automatically. In the context of a dynamic analysis, the time variable is the physical time, differently from static solution steps. The time-period of a given dynamic solution step goes from the previously converged time value until the end time of the current solution step, set by the user. Depending on the easiness or hardness of convergence, time-step may be updated automatically according to the same parameters explained for “Static” type of solution step.

Dynamic simulations also include damping control. Rayleigh damping model is implemented, through usage of the following instruction example:

RayleighDamping Alpha 0 Beta 0 Update 0

The attributes are explained next:

-   Alpha: coefficient that multiplies mass matrix for compounding damping matrix

-   Beta: coefficient that multiplies stiffness matrix for compounding damping matrix

-   Update: flag, which can assume “1” or “0”. If updating is turned on then the damping matrix is updated in each time step beginning, with updated information about stiffness and mass matrices. If updating is turned off then the initial calculated damping (with initial stiffness and mass matrices) is kept during the whole solution step.

Newmark method is used to integrate equations along time. Two coefficients are defined in Newmark method. One can refer to \[17\] for more details. These coefficients are input through NewmarkCoefficients Beta 0.3 Gamma 0.5. These are the recommended values to be used for time-integration. The user can change such values in some particular simulations to induce numerical damping. For example, increasing Gamma from 0.5 to a value up to 0.6 and keeping Beta 0.3 usually introduces high-frequency numerical damping.

When performing a dynamic load step after previous load steps (static/dynamic), the previous velocity/acceleration state are considered as initial conditions. If one wants to make these conditions null, it is possible to include the optional keyword “ZeroIC”. With that, previous load steps final states are not propagated to the current dynamic load step. Initial conditions set by the keyword “InitialCondition” – and referring to such a load step, however, are included normally.

## Modal

*Creates a solution step to solve a modal analysis.*

### Syntax:

Modal SID

ExportMatrices EMF

NumberModes NM

Tolerance TV

ComputeEigenvectors CEF

NumberFrames NF

-   SID: current solution step identification number

-   EMF: a flag to export (1) or not (0) mass and stiffness matrices as text files (sparse matrix formats)

-   NM: number of modes required

-   TV: tolerance for ARPACK<sup>TM</sup> eigenvalues/eigenvectors extraction (use 0 for the Machine Precision based tolerance)

-   CEF: a flag to compute (1) or not (0) the model eigenvectors

-   NF: number of frames exported to animate each model eigenvector

### Example:

Modal 1

ExportMatrices 0

NumberModes 12

Tolerance 1E-6

ComputeEigenvectors 1

NumberFrames 6

### Additional information:

Modal analysis has no end-time information. Thus, during modal analysis time is considered frozen.

When performing modal analysis in a model with special constraints the results will have no meaning, due to Lagrange multipliers present in the model (still not treated for modal analysis in current Giraffe version).

Giraffe evaluates always the lowest magnitude eigenvalues of the system (possibly complex numbers). Depending on the requested results, the available files will be:

-   DOF\_table\_i.txt: a table containing the system connectivity. It contains, for each node and local DOF, the global DOF number.

-   eigenvalues\_solution\_i.txt: a list with the evaluated eigenvalues (real and imaginary parts). In case the eigenvalue is a real number, the natural frequency may be evaluated as the square root of it (rad/s).

-   m\_mass\_i.txt: mass matrix of the system

-   m\_stiffness\_i.txt: stiffness matrix of the system.

Post-processing of modal analysis is detailed in Appendix: Post-processing modal analysis using Paraview<sup>TM</sup>.
