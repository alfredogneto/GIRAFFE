# Element Sets

*Starts a command block for creation of element sets.*

### Syntax:

ElementSets N

//Input method 1:

ElementSet ESID Elements NE List E1 E2 …

//Input method 2:

ElementSet ESID Elements NE Sequence Initial EIN Increment IN

-   N: number of element sets

-   ESID: current element set identification number

-   NE: number of elements defined in the current element set

-   E1, E2, …, : list with NE element identification numbers

-   EIN: initial element identification number

-   IN: increment for the element identification number

### Example:

ElementSets 2

//Input method 1:

ElementSet 1 Elements 3 List 12 27 21

//Input method 2:

ElementSet 2 Elements 4 Sequence Initial 3 Increment 2

### Additional information:

Each element set is defined by the keyword ElementSet followed by an identification number (must be an ascending sequence starting from number one).

There are two input methods to define the element sets. The user must choose one of the following options:

-   List: it indicates to Giraffe that a list of elements will be provided as input. For example, the ElementSet 1 has three elements which are listed after the keyword List;

-   Sequence: it indicates to Giraffe that a sequence of elements will be provided. In the example, ElementSet 2 has 4 elements. The sequence generated automatically will be 3, 5, 7, 9.
