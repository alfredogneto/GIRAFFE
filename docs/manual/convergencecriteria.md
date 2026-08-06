# Convergence Criteria

*Establishes convergence criteria.*

### Syntax:

ConvergenceCriteria

ForceTolerance FTV

MomentTolerance MTV

ForceMinimumReference FMRV

MomentMinimumReference MMRV

ConstraintMinimumReference CMRV

DisplacementTolerance DTV

RotationTolerance RTV

LagrangeTolerance LTV

DisplacementMinimumReference DMRV

RotationMinimumReference RMRV

LagrangeMinimumReference LMRV

DivergenceReference DRV

-   FTV: force tolerance value

-   MTV: moment tolerance value

-   FMRV: force minimum reference value

-   MMRV: moment minimum reference value

-   CMRV: constraint minimum referentece value

-   DTV: displacement tolerance value

-   RTV: rotation tolerance value

-   LTV: Lagrange multiplier tolerance value

-   DMRV: displacement minimum reference value

-   RMRV: rotation minimum reference value

-   LMRV: Lagrange multiplier minimum reference value

-   DRV: divergence reference value

### Example:

ConvergenceCriteria

ForceTolerance 1e-4

MomentTolerance 1e-4

ForceMinimumReference 1e-5

MomentMinimumReference 1e-5

ConstraintMinimumReference 1e-7

DisplacementTolerance 1e-4

RotationTolerance 1e-4

LagrangeTolerance 1e-4

DisplacementMinimumReference 1e-6

RotationMinimumReference 1e-6

LagrangeMinimumReference 1e-6

DivergenceReference 1e+15

### Additional information:

Since Giraffe was designed to solve nonlinear finite element models, one has to define convergence criteria, in order to guide the Newton-Raphson iterative method to stop, according to some rules. Next, we describe each individual convergence criterion that Giraffe applies. A solution is considered as “converged” if all the applied criteria are obeyed simultaneously. Stricter criteria will need more iterations to achieve convergence, but will have more precision.

Default convergence criteria usually works properly for general nonlinear simulations. In such cases, the user does not need to re-establish then. The usage of the ConvergenceCriteria command in Giraffe input file should be made with care, and is proper for advanced users.

-   ForceTolerance: a factor that multiplies the norm of the external forces vector, used to establish a criterion of maximum allowable error for the norm of the unbalanced forces vector. Default value is 0.01%.

-   MomentTolerance: a factor that multiplies the norm of the external moments vector, used to establish a criterion of maximum allowable error for the norm of the unbalanced moments vector. Default value is 0.01%.

-   ForceMinimumReference: a value of force, taken as very small, to avoid the establishment of a never achievable null convergence criterion. This would occur in cases for which no external forces are applied. In such situations, the criterion of maximum allowable error for the norm of the unbalance forces vector is evaluated by: ForceMinimumReference \* ForceTolerance. Default value is 1e-5.

-   MomentMinimumReference: a value of moment, taken as very small, to avoid the establishment of a never achievable null convergence criterion. This would occur in cases for which no external moments are applied. In such situations, the criterion of maximum allowable error for the norm of the unbalance moments vector is evaluated by: MomentMinimumReference \* MomentTolerance. Default value is 1e-5.

-   ConstraintMinimumReference: a small value, taken as the maximum residual for the constraints established by SpecialConstraints command. Default value is 1e-7.

-   DisplacementTolerance: a factor that multiplies the norm of the displacements vector (experienced during the current time-step – for dynamics or sub step – for statics). It is used to establish a criterion of maximum allowable norm for the iterative displacements increment (Newton Raphson). Default value is 0.01%.

-   RotationTolerance: a factor that multiplies the norm of the rotations vector (experienced during the current time-step – for dynamics or sub step – for statics). It is used to establish a criterion of maximum allowable norm for the iterative rotations increment (Newton Raphson). Default value is 0.01%.

-   LagrangeTolerance: a factor that multiplies the norm of the Lagrange multipliers vector (experienced during the current time-step – for dynamics or sub step – for statics). It applies only when the simulation has SpecialConstraints. It is used to establish a criterion of maximum allowable norm for the iterative Lagrange multipliers increment (Newton Raphson). Default value is 0.01%.

-   DisplacementMinimumReference: a value of displacement, taken as very small, to avoid the establishment of a never-achievable null criterion for the maximum allowable norm of the iterative displacements increment. This would occur for cases in which no displacements occur in the system. In such situations, the criterion for maximum allowable norm of the iterative displacements increment is evaluated by: DisplacementMinimumReference \* DisplacementTolerance. Default value is 1e-6.

-   RotationMinimumReference: a value of rotation, taken as very small, to avoid the establishment of a never-achievable null criterion for the maximum allowable norm of the iterative rotations increment. This would occur for cases in which no rotations occur in the system. In such situations, the criterion for maximum allowable norm of the iterative rotations increment is evaluated by: RotationMinimumReference \* RotationTolerance. Default value is 1e-6.

-   LagrangeMinimumReference: a value of Lagrange multiplier, taken as very small, to avoid the establishment of a never-achievable null criterion for the maximum allowable norm of the iterative Lagrange multipliers increment. This would occur for cases in which only null Lagrange multipliers occur in the system. In such situations, the criterion for maximum allowable norm of the iterative Lagrange multipliers increment is evaluated by: LagrangeMinimumReference \* LagrangeTolerance. It applies only when the simulation has SpecialConstraints. Default value is 1e-6.

-   DivergenceReference: defines a very large residual number, which once achieved, means “divergence”. Then, Giraffe automatically will perform a bisection (dividing the last load factor by two) in order to try to achieve convergence in next time-step or sub step. Default value is 1e+15.
