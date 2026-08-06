# Tutorial 2: Modal analysis of a circular plate

!!! example "Input files"
    [`tutorial02.inp`](../assets/tutorials/inputs/tutorial02/tutorial02.inp)
    (also in the repository under `inputs/tutorial02/`)

Input file name: tutorial02.inp

*Developed by Paulo Roberto Refachinho & Gabriel Maki*

## Objective

This tutorial presents the modal analysis of a circular plate, clamped on its boundaries. The tutorial objective is to show how to perform a modal analysis and post-process results using Paraview<sup>TM</sup> software.

<img src="../../assets/tutorials/image9.png" style="width:3.61181in;height:2.4403in" alt="" /><img src="../../assets/tutorials/image10.png" style="width:2.0639in;height:2.03096in" alt="" />

Figure 2‑1 – Thin circular plate (a) Geometry (b) Mesh

## Description

The plate is meshed using 6197 nodes and 3036 elements, using the element Shell\_1. There is a clamp boundary condition on the edges of the plate. The section has a thickness of 2.5 mm. The plate is contained in the plane perpendicular to the direction z.

## Input File step-by-step

### Creation of nodes

A total of 6197 nodes are created to establish a finite element mesh. Note that here a precision of 14 digits is used to write the input file. This is recommended always when generating Giraffe input files. The scientific notation format is recommended, as shown in this example, unless when dealing with a very simple nodal sequence, such as described in tutorial 01.

<table><tbody><tr class="odd"><td><p>//Creation of nodes</p><p>Nodes 6197</p><p>Node 1 0.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 2 5.00000000000000E+00 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 3 1.00000000000000E+01 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 4 1.49999999999999E+01 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 5 1.99999999999999E+01 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 6 2.49999999999999E+01 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 7 2.99999999999999E+01 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 8 3.49999999999999E+01 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 9 4.00000000000000E+01 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 10 4.50000000000000E+01 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 11 5.00000000000000E+01 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 12 5.49999999999999E+01 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 13 6.00000000000000E+01 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 14 6.49999999999999E+01 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 15 6.99999999999999E+01 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 16 7.49999999999999E+01 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 17 7.99999999999999E+01 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 18 8.49999999999999E+01 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 19 9.00000000000000E+01 0.00000000000000E+00 0.00000000000000E+00</p><p>Node 20 9.50000000000000E+01 0.00000000000000E+00 0.00000000000000E+00</p></td></tr></tbody></table>

…

<table><tbody><tr class="odd"><td><p>Node 6165 1.53058982842033E+02 -1.14026258730851E+02 0.00000000000000E+00</p><p>Node 6166 1.47913129412780E+02 -1.21843677440843E+02 0.00000000000000E+00</p><p>Node 6167 1.42432380390826E+02 -1.29390594783211E+02 0.00000000000000E+00</p><p>Node 6168 1.74941348706255E+02 -4.84057409095254E+00 0.00000000000000E+00</p><p>Node 6169 1.74668214191052E+02 -1.43671212253454E+01 0.00000000000000E+00</p><p>Node 6170 1.74020872610914E+02 -2.38210467458590E+01 0.00000000000000E+00</p><p>Node 6171 1.72999502067091E+02 -3.31867762377589E+01 0.00000000000000E+00</p><p>Node 6172 1.71605951573866E+02 -4.24483843164989E+01 0.00000000000000E+00</p><p>Node 6173 1.69852358429996E+02 -5.15942395774600E+01 0.00000000000000E+00</p><p>Node 6174 1.67736588624359E+02 -6.06070008207730E+01 0.00000000000000E+00</p><p>Node 6175 1.65264625712240E+02 -6.94721590337491E+01 0.00000000000000E+00</p><p>Node 6176 1.62451464824384E+02 -7.81804125423577E+01 0.00000000000000E+00</p><p>Node 6177 1.59292155070117E+02 -8.67116171422745E+01 0.00000000000000E+00</p><p>Node 6178 1.55801013230708E+02 -9.50559369314019E+01 0.00000000000000E+00</p><p>Node 6179 1.51987536563878E+02 -1.03201993087607E+02 0.00000000000000E+00</p><p>Node 6180 1.47855378730407E+02 -1.11133271194406E+02 0.00000000000000E+00</p><p>Node 6181 1.43412535206243E+02 -1.18836169266648E+02 0.00000000000000E+00</p><p>Node 6182 1.38680440329816E+02 -1.26309596164411E+02 0.00000000000000E+00</p><p>Node 6183 1.64989172058302E+02 -4.71629909177040E+00 0.00000000000000E+00</p><p>Node 6184 1.64847824124615E+02 -1.39824570121293E+01 0.00000000000000E+00</p><p>Node 6185 1.64382576894756E+02 -2.31704236670644E+01 0.00000000000000E+00</p><p>Node 6186 1.63591871432435E+02 -3.22667440764171E+01 0.00000000000000E+00</p><p>Node 6187 1.62484093347860E+02 -4.12620306230941E+01 0.00000000000000E+00</p><p>Node 6188 1.61064694050891E+02 -5.01458000449104E+01 0.00000000000000E+00</p><p>Node 6189 1.59328661600592E+02 -5.89021671438646E+01 0.00000000000000E+00</p><p>Node 6190 1.57290898762299E+02 -6.75253544690905E+01 0.00000000000000E+00</p><p>Node 6191 1.54955010068830E+02 -7.60041636670082E+01 0.00000000000000E+00</p><p>Node 6192 1.52320974309620E+02 -8.43241501432042E+01 0.00000000000000E+00</p><p>Node 6193 1.49397692831169E+02 -9.24769180839118E+01 0.00000000000000E+00</p><p>Node 6194 1.46199126112956E+02 -1.00458125059544E+02 0.00000000000000E+00</p><p>Node 6195 1.42717259396612E+02 -1.08246358021446E+02 0.00000000000000E+00</p><p>Node 6196 1.39324988441771E+02 -1.15795202781407E+02 0.00000000000000E+00</p><p>Node 6197 1.35336888727668E+02 -1.23202163313770E+02 0.00000000000000E+00</p></td></tr></tbody></table>

### Creation of node sets

A node set, which has 248 nodes, is created to establish the constraints.

<table><tbody><tr class="odd"><td><p>//Creation of node sets (in this case, a single node set)</p><p>NodeSets 1</p><p>NodeSet 1 Nodes 248 List 41 42 43 44 45 46 47</p></td></tr></tbody></table>

…

<table><tbody><tr class="odd"><td><p>4739 4740 4741 4742 4743 4744 4745 4746 4747 4748 4749</p></td></tr></tbody></table>

### Creation of elements

3036 elements are defined. All the elements have the same material properties (identification number 1) and the same cross-section (identification number 1). Each Shell\_1 element has the shape of a triangle and is defined by six nodes, where three nodes are located at the vertices and the others are positioned at the midpoints of the sides of the triangle. The nodes within an element must lie in a plane at the given reference configuration mesh. After deforming, arbitrary displacements may occurs and the shell element will be no longer contained in a plane.

<table><tbody><tr class="odd"><td><p>//Creation of elements</p><p>Elements 3036</p><p>Shell_1 1 Mat 1 Sec 1 Nodes 143 157 186 156 185 1590</p><p>Shell_1 2 Mat 1 Sec 1 Nodes 150 186 188 192 187 1589</p><p>Shell_1 3 Mat 1 Sec 1 Nodes 152 188 190 193 189 1588</p><p>Shell_1 4 Mat 1 Sec 1 Nodes 154 190 71 194 1587 155</p><p>Shell_1 5 Mat 1 Sec 1 Nodes 157 159 196 158 195 1586</p><p>Shell_1 6 Mat 1 Sec 1 Nodes 186 196 198 202 197 1585</p><p>Shell_1 7 Mat 1 Sec 1 Nodes 188 198 190 203 1584 189</p><p>Shell_1 8 Mat 1 Sec 1 Nodes 190 200 69 204 1583 191</p><p>Shell_1 9 Mat 1 Sec 1 Nodes 159 161 206 160 205 1582</p><p>Shell_1 10 Mat 1 Sec 1 Nodes 196 206 198 212 1581 197</p><p>Shell_1 11 Mat 1 Sec 1 Nodes 198 208 200 213 1580 199</p><p>Shell_1 12 Mat 1 Sec 1 Nodes 200 210 67 214 1579 201</p><p>Shell_1 13 Mat 1 Sec 1 Nodes 161 163 216 162 215 1578</p><p>Shell_1 14 Mat 1 Sec 1 Nodes 206 216 208 222 1577 207</p></td></tr></tbody></table>

…

<table><tbody><tr class="odd"><td><p>Shell_1 3022 Mat 1 Sec 1 Nodes 4819 35 33 4820 34 6183</p><p>Shell_1 3023 Mat 1 Sec 1 Nodes 4817 4819 4783 4818 4834 6184</p><p>Shell_1 3024 Mat 1 Sec 1 Nodes 4815 4817 4781 4816 4833 6185</p><p>Shell_1 3025 Mat 1 Sec 1 Nodes 4813 4815 4779 4814 4832 6186</p><p>Shell_1 3026 Mat 1 Sec 1 Nodes 4811 4813 4777 4812 4831 6187</p><p>Shell_1 3027 Mat 1 Sec 1 Nodes 4809 4811 4775 4810 4830 6188</p><p>Shell_1 3028 Mat 1 Sec 1 Nodes 4807 4809 4773 4808 4829 6189</p><p>Shell_1 3029 Mat 1 Sec 1 Nodes 4805 4807 4771 4806 4828 6190</p><p>Shell_1 3030 Mat 1 Sec 1 Nodes 4803 4805 4769 4804 4827 6191</p><p>Shell_1 3031 Mat 1 Sec 1 Nodes 4801 4803 4767 4802 4826 6192</p><p>Shell_1 3032 Mat 1 Sec 1 Nodes 4799 4801 4765 4800 4825 6193</p><p>Shell_1 3033 Mat 1 Sec 1 Nodes 4797 4799 4763 4798 4824 6194</p><p>Shell_1 3034 Mat 1 Sec 1 Nodes 4795 4797 4761 4796 4823 6195</p><p>Shell_1 3035 Mat 1 Sec 1 Nodes 4795 4759 4757 4822 4758 6196</p><p>Shell_1 3036 Mat 1 Sec 1 Nodes 4793 4757 4750 4821 4756 6197</p></td></tr></tbody></table>

### Creation of materials

One only material is created, which identification number is 1.

<table><tbody><tr class="odd"><td><p>//Creation of materials</p><p>Materials 1</p><p>Hooke 1 E 200000 Nu 0.30 Rho 7.85e-9</p></td></tr></tbody></table>

### Creation of shell sections

One shell section is defined, which has a thickness of 2.5 cm and it is homogeneous. Its identification number is 1.

<table><tbody><tr class="odd"><td><p>//Creation of shell sections</p><p>ShellSections 1</p><p>Homogeneous 1 Thickness 2.5</p></td></tr></tbody></table>

### Creation of constraints

A nodal constraint is created to fix all degrees of freedom of the nodes contained in the node set defined previously. In this case, the nodal constraint is established for the solution step 1 (in current tutorial there is a single solution step), which will be shown in the next topic.

<table><tbody><tr class="odd"><td><p>//Creation of constraints</p><p>Constraints 1</p><p>NodalConstraint 1 NodeSet 1</p><p>UX BoolTable 1</p><p>UY BoolTable 1</p><p>UZ BoolTable 1</p><p>ROTX BoolTable 1</p><p>ROTY BoolTable 1</p><p>ROTZ BoolTable 1</p></td></tr></tbody></table>

### Creation of solution steps

This tutorial has a single solution step (modal). There is the possibility of choosing the number of modes to be evaluated, the tolerance for eigenvalues and the number of frames written for post-processing the mode shapes, using Paraview<sup>TM</sup>. Alternatively, the user may export matrices involved in modal analysis (mass and stiffness matrix). However, this will generate large-size text files. We chose here not to export such files.

<table><tbody><tr class="odd"><td><p>//Creation of solution steps</p><p>SolutionSteps 1</p><p>Modal 1</p><p>ExportMatrices 0</p><p>NumberModes 12</p><p>Tolerance 1e-6</p><p>ComputeEigenvectors 1</p><p>NumberFrames 12</p></td></tr></tbody></table>

### Creation of solver options

Here it is possible to choose the number of processors used (parallel processing). Besides, we must define the solver type for systems of linear equations. In this case, two threads are employed and the solver choice is “Direct”.

<table><tbody><tr class="odd"><td><p>//Creation of solver options</p><p>SolverOptions</p><p>Processors 2 LinSys Direct</p></td></tr></tbody></table>

In this tutorial, there is no need to create monitors and post-files.

## Results and post-processing

In the case of a single solution step analysis, Giraffe saves result files on “/post/solution\_1”. The file “solution\_1\_mesh.pvd” contain all the twelve requested modes.

<img src="../../assets/tutorials/image11.png" style="width:5.90556in;height:3.36181in" alt="" />

Figure 2‑2 – Opening the file in Paraview<sup>TM</sup>.

To view a specific vibration mode, first we must open in Paraview<sup>TM</sup> the file mentioned, as shown in Figure 2‑2, and click the preview icon and the Apply button on the Pipeline browser (Figure 2‑3).

<img src="../../assets/tutorials/image12.png" style="width:1.71642in;height:3.0742in" alt="" />

Figure 2‑3 – Pipeline Browser.

The following figure shows how Paraview<sup>TM</sup> displays the results.

<img src="../../assets/tutorials/image13.png" style="width:5.89583in;height:2.70139in" alt="" />

Figure 2‑4 – Results of modal analysis.

The second step is extracting the mode of interest, which can be done using the “Extract Block”. filter.

<img src="../../assets/tutorials/image14.png" style="width:4.25373in;height:2.26742in" alt="" />

Figure 2‑5 – Finding the “Extract Block”.

After clicking in the filter, a menu will appear to choose the desired block. Each block contains results for a mode. In this case, there will be twelve blocks (associated with the requested twelve first vibration modes). To animate results for a vibration mode, we must choose the desired block and click “Apply”.

<img src="../../assets/tutorials/image15.png" style="width:1.87313in;height:3.5482in" alt="" />

Figure 2‑6 – Choosing the sixth mode of vibration.

Figure 2‑7 shows an example of mode of vibration plot. It shows the expected shape when the structure experiences vibration patterns associated with a corresponding natural frequency.

<img src="../../assets/tutorials/image16.png" style="width:5.90278in;height:2.40278in" alt="" />

Figure 2‑7 – Animation of the chosen mode of vibration.

Finally, there is a text file “solution\_1\_eigenvalues” that contains the eigenvalues obtained from the modal analysis. It is located inside the folder “post/solution\_1. To obtain the natural frequencies, we must extract the square root of the eigenvalues and divide this by 2π.

<table><tbody><tr class="odd"><td><p>Eigenvalues - solution number 1</p><p>Mode Re Im</p><p>1 1009519.548036 0.000000</p><p>2 4316095.037999 0.000000</p><p>3 4316095.038019 0.000000</p><p>4 11630174.148244 0.000000</p><p>5 11561681.465562 0.000000</p><p>6 14897498.883645 0.000000</p><p>7 24801521.212492 0.000000</p><p>8 24801521.212497 0.000000</p><p>9 46155343.969949 0.000000</p><p>10 46177633.073695 0.000000</p><p>11 34857355.962447 0.000000</p><p>12 34857355.962431 0.000000</p></td></tr></tbody></table>
