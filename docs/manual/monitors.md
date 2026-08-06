# Monitors

*Creates monitors for post-processing results of nodes, elements, contacts and node sets.*

### Syntax:

Monitor Sample SV

//Optional keyword to monitor nodes:

MonitorNodes nodes IDs

//Optional keyword to monitor elements:

MonitorElements elements IDs

//Optional keyword to monitor contacts:

MonitorContacts contacts IDs

//Optional keyword to monitor node sets:

MonitorNodeSets node sets IDs

-   SV: sampling for saving data in monitor output files (use 1 to save all converged solutions and larger integer numbers for decreasing data size)

-   nodes IDs: list of node identification numbers (to be monitored)

-   elements IDs: list of element identification numbers (to be monitored)

-   contacts IDs: list of contact identification numbers (to be monitored)

-   node sets IDs: list of node set identification numbers (to be monitored)

### Example:

Monitor Sample SV

MonitorNodes 1 2

MonitorElements 1 2

MonitorContacts 1

MonitorNodeSets 1

### Additional information:

Monitors are extremely useful for analyzing time series. Each monitor has specific file formats, depending on the entity chosen. Giraffe saves monitor files during the simulation evolution with a sampling frequency ruled by the attribute input after Sample keyword. NodeSets monitor evaluates along time the following quantities:

-   the average position of the nodes in the node set;

-   the total force applied on all nodes in the node set;

-   the total moment applied on the nodes in the node set (the total moment is composed by the moments applied at each node and the transport moment of the force in each node to the average position of the nodes – taken as pole)
