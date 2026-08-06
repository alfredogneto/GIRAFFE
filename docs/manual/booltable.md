# Bool Table

*Creates a Boolean table to rule the behavior of some loads, displacements, contacts and other resources along solution steps. This keyword is used as parameter for many Giraffe resources.*

### Syntax:

BoolTable B1 B2 B3 … BN

-   B1, B2, B3, …, BN: sequence of 1 or 0 values composing a Boolean table

### Example:

BoolTable 1 0 0 1

### Additional information:

Bootable command is used to provide to Giraffe a sequence of numbers 1 or 0. It is used to define if a given resource is active (1) or inactive (0) in a sequence of solution steps. It can be used with many Giraffe resources, such as environment data, constraints data, contact and special constraints.

As example, in a scenario of 3 sequential solution steps requested by the user, one may be interested in including environment data to define gravity field and also define some constraints, as depicted below:

Constraints 1

NodalConstraint 1 NodeSet 1

UX BoolTable 1 1 1

UY BoolTable 0 1 1

UZ BoolTable 1 0 0

Environment

GravityData

G 0 0 -9.81 BoolTable 0 1 1

Giraffe would interpret nodal constraints as:

-   UX constraint would be considered during the first, second and third solution steps;

-   UY constraint would be considered during the second and third solution steps, but would be disregarded during first solution step;

-   UZ constraint would be considered only during the first solution step and disregarded during second and third solution steps.

Gravity data would be interpreted to be not applied during the first solution step, but to be applied during the second solution step by a linear ramp function along time evolution. During the third solution step, gravity data would be kept.

Note that BoolTable can be defined with less data than the number of solution steps requested by the user. In this case, Giraffe uses the last provided value as the same for subsequent undefined data. For example, Giraffe interprets “BoolTable 1 1 1” data as the same as “BoolTable 1 1” and also “BoolTable 1”. Thus, if the user wants to include some resource using BoolTable from beginning and just keeping it defined along arbitrary solution sequence, it is enough to provide simply: “BoolTable 1”. Another example: “BoolTable 1 0 0” data is interpreted as the same as “BoolTable 1 0”.

Gravity and ocean data employs BoolTable, as provided in the example. Always the insertion of such loads within a given solution step is done by a linear ramp function along time evolution, both increasing or decreasing the load. For example, if one defines “BoolTable 0 1 0” for GravityData, during the first solution step gravity loads would not be applied, during the second solution step gravity data would be included in the model (by a linear ramp increasing function along time) and during the third solution step gravity data would be, again, switched off (by a linear ramp decreasing function along time).

When using BoolTable with another resources, such as contact or special constraints, the interpretation is straightforward, as one may see in following example:

SameDisplacement 1 Nodes 1 2 BoolTable 1 0 1

The SameDisplacement special constraint would be considered during the first and third solution steps, but not during the second solution step. Thus, BoolTable permits to switch on/off constraints, special constraints and contacts along solution, when changing between solution steps.
