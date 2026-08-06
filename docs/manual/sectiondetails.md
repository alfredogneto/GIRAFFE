# Section Details

*Starts a command block for creation of section details (cross-section details).*

### Syntax:

SectionDetails N

Name ID data

-   N: number of section details

-   Name: current section detail name

-   ID: current section detail identification number

-   data: current section detail data (depends on section detail resources and requirements)

### Example:

SectionDetails 1

SolidSection 1 AxisPosition 0.1 0.12 NPoints 4

Point 1 0.4 0

Point 2 0 0.2

Point 3 -0.4 0

Point 4 0 -0.2

### Additional information:

Each section detail is defined by a specific keyword followed by the section detail identification number (must be an ascending sequence starting from number one) and additional data. Each section detail available and its input data is explained next.

## SolidSection

*Creates a solid section details (for post-processing and visualization purposes).*

### Syntax:

SolidSection SDID AxisPosition XAV YAV NPoints NP

Point PID XPV YPV

-   SDID: current section detail identification number

-   XAV and YAV: coordinates of the axis position in element local coordinate system (according to the properties employed to establish the element cross-section)

-   NP: number of points employed to define the cross-section

-   PID: identification number of each point employed to define the cross-section external boundary

-   XPV and YPV: coordinates of each point employed to define the cross-section external boundary

### Example:

SolidSection 1 AxisPosition 0.1 0.12 NPoints 4

Point 1 0.4 0

Point 2 0 0.2

Point 3 -0.4 0

Point 4 0 -0.2

### Additional information:

The objective of the section details is to be used only for post-processing purposes. This data is not used by Giraffe mathematical model. When using beam elements with UserDefined cross-sections, the section details are then addressed by the PostFiles keyword.

The cross-section points are defined in a local coordinate system, associated with the beam element, as presented in UserDefined cross-section explanation. The cross-section local plane is given by directions E1 and E2 assigned to the beam element (see CoordinateSystems keyword). The coordinates of the beam axis in this local plane are defined after the keyword AxisPosition (origin for evaluating properties for UserDefined cross sections). Then, the number of points that will be used to define the cross section external boundary is defined by the keyword NPoints. Finally, each point is defined in a list by the keyword Point, followed by an ascending identification number and the coordinates of the point.

## MultiCellSection

*Creates a multi-cell section details (for post-processing and visualization purposes).*

### Syntax:

MultiCellSection. SDID AxisPosition XAV YAV NPoints NP NWebs NW

Point PID XPV YPV

Web WID W1 W2

-   SDID: current section detail identification number

-   XAV and YAV: coordinates of the axis position in element local coordinate system (according to the properties employed to establish the element cross-section)

-   NP: number of points employed to define the cross-section

-   NW: number of webs employed to define the cross-section

-   PID: identification number of each point employed to define the cross-section external boundary

-   XPV and YPV: coordinates of each point employed to define the cross-section external boundary

-   WID: identification number of each web

-   W1 and W2: identification number of the points connecting web connections with the external cross-section boundary

### Example:

MultiCellSection 1 AxisPosition 0.0 0.0 NPoints 6 NWebs 1

Point 1 0.4 +0.2

Point 2 0.4 -0.2

Point 3 0 -0.2

Point 4 -0.4 -0.2

Point 5 -0.4 +0.2

Point 6 0 +0.2

Web 1 3 6

### Additional information:

This type of section detail is similar to the solid section details. The only difference is the presence of webs in visualization. Rendering also considers the structure as thin-walled and not solid cross-section.
