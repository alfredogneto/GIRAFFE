# Concomitant Solution

*Establishes a modal concomitant solution to be solved repeated times within a given solution step (or solution steps).*

### Syntax:

ConcomitantSolution Sample SV BoolTable BDC

Modal NumberModes NM Tolerance TV

-   SV: sampling variable to rule concomitant solution call

-   BDC: bool table data for concomitant solution

-   NM: number of modes required

-   TV: tolerance for ARPACK<sup>TM</sup> eigenvalues/eigenvectors extraction (use 0 for the Machine Precision based tolerance)

### Example:

ConcomitantSolution Sample 5 BoolTable 1

Modal NumberModes 10 Tolerance 1E-6

### Additional information:

A concomitant solution may be created to ask Giraffe to solve extra solutions within a given solution step, or even along more than one solution step. For example, during a static or dynamic solution steps, the user may be interested in evaluating system modal analysis along time evolution. In this case a concomitant solution may be created. It does not influence in time-evolution, neither in solution steps sequence.

According to the sample variable a modal concomitant solution will be called. If SV is “1”, all converged time-steps will lead to a call of concomitant solution. Otherwise, larger integer SV will lead to less concomitant solution calls, always at each SV converged time-steps. Currently only concomitant modal analysis is available in Giraffe.

At the end of the simulation, the user will find as the result of concomitant solution a text file containing the time-series of evaluated eigenvalues along time. It is located inside the folder “/post/concomitant\_solution/”.
