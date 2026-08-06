# Node Sets

*Starts a command block for creation of node sets.*

### Syntax:

NodeSets N

//Input method 1:

NodeSet NSID Nodes NN List N1 N2 …

//Input method 2:

NodeSet NSID Nodes NN Sequence Initial NIN Increment IN

-   N: number of node sets

-   NSID: current node set identification number

-   NN: number of nodes defined in the current node set

-   N1, N2, …: list with NN node identification numbers

-   NIN: initial node identification number

-   IN: increment for the node identification number

### Example:

NodeSets 2

//Input method 1:

NodeSet 1 Nodes 3 List 12 27 21

//Input method 2:

NodeSet 2 Nodes 4 Sequence Initial 3 Increment 2

### Additional information:

Each node set is defined by the keyword NodeSet followed by an identification number (must be an ascending sequence starting from number one).

There are two input methods to define the node sets. The user must choose one of the following options:

-   List: it indicates to Giraffe that a list of nodes will be provided as input. For example, the NodeSet 1 has three nodes which are listed after the keyword List;

-   Sequence: it indicates to Giraffe that a sequence of nodes will be provided. In the example, NodeSet 2 has 4 nodes. The sequence generated automatically will be 3, 5, 7, 9.
