# Shell Sections

*Starts a command block for creation of shell sections.*

### Syntax:

ShellSections N

Name ID data

-   N: number of shell sections

-   Name: current shell section name

-   ID: current shell section identification number

-   data: current shell section data (depends on shell section resources and requirements)

### Example:

ShellSections 1

Homogeneous 1 Thickness 0.25

### Additional information:

Each shell section is defined by a specific keyword followed by the shell section identification number (must be an ascending sequence starting from number one) and additional data. Each shell section available and its input data is explained next.

Shell sections are to be used as parameters for the element Shell\_1.

## Homogeneous

*Creates a homogeneous shell section.*

### Syntax:

Homogeneous SID Thickness TV

-   SID: current shell section identification number

-   TV: shell section thickness value

### Example:

Homogeneous 1 Thickness 0.25

### Additional information:

It is assumed that the shell surface is located on the middle of the thickness height.

## Composite

*Creates a composite shell section.*

### Syntax:

Composite SID Laminas LN

//Material ID Thickness Orientation

table data

-   SID: current shell section identification number

-   LN: number of laminas that compose the composite shell section

-   table data: table containing each lamina information

column 1: material identification number

column 2: lamina thickness

column 3: lamina principal angle in degrees with respect to the local coordinate system

### Example:

Composite 1 Laminas 3

1 0.006 0

2 0.012 90

1 0.006 0

### Additional information:

It is assumed that the shell surface is located on the middle of the total thickness height. The principal angle orientation of each lamina with respect to the local coordinate system and the stacking order are defined as seen in Figure 18.

<table><thead><tr class="header"><th><img src="../../assets/manual/image20.png" style="width:2.21875in;height:1.49574in" alt="" /></th><th><img src="../../assets/manual/image21.png" style="width:1.80208in;height:1.22836in" alt="" /></th></tr></thead><tbody><tr class="odd"><td><ol type="a"><li><p>Principal angle orientation</p></li></ol></td><td><ol start="2" type="a"><li><p>Laminas stacking order</p></li></ol></td></tr><tr class="even"><td>Figure 18 – Composite shell section orientations</td><td></td></tr></tbody></table>
