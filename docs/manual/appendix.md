# Appendix

## Selection by element properties in Giraffe data using Paraview<sup>TM</sup>

### Problem statement:

Imagine that you have data results for a processed simulation in Giraffe, containing many kinds of elements, material ID’s (numbers), etc.

A hypothetic post-processing scenario is proposed: you would like to create a plot only containing cells, filtered by element type, material properties ID, or to a more complex extract based on data information. This is possible by using Paraview’s filters.

### Step-by-step procedure:

1.  Select in the model pipeline browser the data you would like to operate with. In this example, we choose a render mesh data.

    <img src="../../assets/manual/image37.png" style="width:4.68819in;height:2.62361in" alt="" />

2.  Create an “Extract Component” filter:

    <img src="../../assets/manual/image38.png" style="width:5.62361in;height:3.15556in" alt="" />

3.  Choose the input array to guide the filtering operation.

    Giraffe writes “ElementProperties” array associated with cells in a render mesh visualization, which contains numbers following the meaning:

    > 0 – element type number
    >
    > 1 – material ID
    >
    > 2 – section ID
    >
    > 3 – coordinate system ID
    >
    > In our example, we are interested in element type number, since we want to select only a given type of element. Choose a name for the output of the filter. In our example, we chose “ElementType”. Click “Apply” to make the filtering operation have effect.
    >
    > <img src="../../assets/manual/image39.png" style="width:2.01319in;height:4.22708in" alt="" />

    Now we need to extract only the “Beam\_1” elements from the just-filtered data. Element types numbering follow Table 7. This is done by the “Threshold” filter. Create it to operate in data just-filtered in previous step:

    <img src="../../assets/manual/image40.png" style="width:2.22708in;height:1.79236in" alt="" />

4.  In “Threshold” filter properties, choose the scalar to guide the new selection. In this example, “ElementType”, just created in previous step. Then, choose the Minimum/Maximum values to control the new selection range, based on the scalar selected. In our case, we are interested only in number 1 – associated to Beam\_1 elements.

    <img src="../../assets/manual/image41.png" style="width:2.10417in;height:2.61667in" alt="" />

5.  The new plot will filter only Beam\_1 elements:

    <img src="../../assets/manual/image42.png" style="width:4.27292in;height:2.15556in" alt="" />

Alternatively, one may follow another steps by running a Python routine in Paraview<sup>TM</sup>:

### Step-by-step procedure:

1.  Select in the model pipeline browser the data you would like to operate with. In this example, we choose a render mesh data.

2.  Click in Tools-&gt;Python Shell.

3.  Copy and Paste the following python routine in the Python Shell and press enter.

(script available in Giraffe Releases/Documentation/Giraffe&Paraview/selection\_script.py)

<table><tbody><tr class="odd"><td><p>##################################################################</p><p>#Python script for selection using Giraffe ElementProperties data#</p><p>##################################################################</p><p>src = GetActiveSource() #obtains the active source of data (selection in pipeline browser)</p><p>filt1 = ExtractComponent() #creates a filter in variable 'filt1'</p><p>filt1.InputArray = 'ElementProperties' #assigns 'ElementProperties' as the InputArray to 'filt1'</p><p>filt1.Component = 0 #assigns the index '0' (ElementType) as the Component to 'filt1'.</p><p>filt1.OutputArrayName = 'ElementTypeNumber' #assigns 'ElementType' as the OutputArrayName to 'filt1'</p><p>f = Threshold() #creates a filter</p><p>filt1.UpdatePipeline()</p><p>SetActiveSource(src) #sets the original source of data</p><p>filt2 = ExtractComponent() #creates a filter in variable 'filt2'</p><p>filt2.InputArray = 'ElementProperties' #assigns 'ElementProperties' as the InputArray to 'filt2'</p><p>filt2.Component = 1 #assigns the index '1' (MaterialNumber) as the Component to 'filt2'.</p><p>filt2.OutputArrayName = 'MaterialNumber' #assigns 'MaterialNumber' as the OutputArrayName to 'filt2'</p><p>f = Threshold() #creates a filter</p><p>filt2.UpdatePipeline()</p><p>SetActiveSource(src) #sets the original source of data</p><p>filt3 = ExtractComponent() #creates a filter in variable 'filt3'</p><p>filt3.InputArray = 'ElementProperties' #assigns 'ElementProperties' as the InputArray to 'filt3'</p><p>filt3.Component = 2 #assigns the index '2' (SectionNumber) as the Component to 'filt3'.</p><p>filt3.OutputArrayName = 'SectionNumber' #assigns 'SectionNumber' as the OutputArrayName to 'filt3'</p><p>f = Threshold() #creates a filter</p><p>filt3.UpdatePipeline()</p><p>SetActiveSource(src) #sets the original source of data</p><p>filt4 = ExtractComponent() #creates a filter in variable 'filt4'</p><p>filt4.InputArray = 'ElementProperties' #assigns 'ElementProperties' as the InputArray to 'filt4'</p><p>filt4.Component = 3 #assigns the index '3' (CSNumber) as the Component to 'filt4'.</p><p>filt4.OutputArrayName = 'CSNumber' #assigns 'CSNumber' as the OutputArrayName to 'filt4'</p><p>f = Threshold() #creates a filter</p><p>filt4.UpdatePipeline()</p><p>SetActiveSource(src) #sets the original source of data</p></td></tr></tbody></table>

## Post-processing modal analysis using Paraview<sup>TM</sup>

### Problem statement:

When performing a modal analysis as a solution step “i” the user may request the vibration modes evaluation. In this case, Giraffe saves result files on “/post/solution\_i”. These files are automatically loaded in Paraview<sup>TM</sup> by opening the related “solution\_i\_mesh.pvd” file, located in “/post” folder. Next, we show how to post-process vibration modes using Paraview<sup>TM</sup>.

### Step-by-step procedure:

1.  Open in Paraview<sup>TM</sup> the file “solution\_i\_mesh.pvd” (for a given “i”). Click Apply button on the Pipeline browser. With that, all modes will be opened simultaneously. Paraview<sup>TM</sup> will show results like this:

    <img src="../../assets/manual/image43.png" style="width:5.90556in;height:3.33264in" alt="" />

2.  We need to instruct Paraview<sup>TM</sup> to extract each mode of interest for plotting results. This is done by employing a filter named: “Extract Block”. It can be activated by selecting in the model tree the file “solution\_i\_mesh.pvd”. Then, go to:

    Filters-&gt;Alphabetical-&gt;Extract Block. Paraview<sup>TM</sup> will show a menu for the choice of the desired block, to post-process:

    <img src="../../assets/manual/image44.png" style="width:2.60531in;height:5.07156in" alt="" />

3.  Each block encompasses results for a given vibration mode. For example, by selecting “Part 4” block and clicking “Apply”, we are able to see and animate results for the vibration mode associated with the fourth eigenvalue found by Giraffe. This follows also for the other parts, associating the part number with the eigenvalue sequential number.

    <img src="../../assets/manual/image45.png" style="width:5.90556in;height:3.33264in" alt="" />
