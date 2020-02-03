# colEDGE
Analysis of location and composition of bacterial colony edges in microscopy data

## Introduction

When capturing the dynamics of bacterial colony expansion, it is often useful to use tiles of images to monitor the position and composition of the colony edge at high resolution over long periods of time. For example, the following image is taken from a 2 mm x 0.15 mm x 12 hr dataset, and has been stitched together from 20 adjacent fields of view:

<p align="center">
  <img src="https://raw.githubusercontent.com/Pseudomoaner/cellsOnEdge/master/Images/StitchedImg.png" alt="EdgeSchematic"/>
</p>

While this approach can provide a combination of excellent spatial resolution and long sampling times, the resulting datasets are typically too large to analyse efficiently by hand. colEDGE is a suite of tools designed to automatically find the position and composition of the colony edge within these datasets.

## Usage

### Part 0: Preprocessing
Begin in Fiji/ImageJ.

1. Brightfield images should be pre-processed to remove large-scale variations in intensity. They should also be split into separate frames, with each field of view stored in a separate folder. The script ProcessBulkBFSingleFrames.ijm does this job.

2. Fluorescence images should be split into separate frames, with each field of view stored in a separate folder. Despeckling is also advised. These jobs are performed by the script ProcessBulkFluo.ijm.

### Part 1: Finding the colony edge
Now switch to Matlab.

1. At this point you can use one of two scripts: PlotEdgeCoordsWithLinearToLinearModels.m or PlotEdgeCoordsWithExponentialToLinearModel.m. Both of these will extract the coordinates of the colony edge at all times, but they vary in the model that they fit for plotting:

   - PlotEdgeCoordsWithLinearToLinearModels.m fits a model consisting of two constant expansion rates, one early and one late.
   - PlotEdgeCoordsWithExponentialToLinearModel.m fits a model consisting of an exponentially increasing expansion rate, smoothly transitioning into a constant expansion rate later.
   
2. Edge extraction works through the following steps. Parameters involved at each stage are indicated in bold and in brackets:

   1. Binarise each frame in the tile using a combination of texture analysis and ridge detection (all fields of **segSets**).
   2. Stitch the resulting images to create a single fused binary image (all fields of **stitchSets**).
   3. Cut resulting binary image into strips running perpendicular to the direction of motion (**edgeSets.noBins**).
   4. Find the average coverage along each strip.
   5. The edge is found as a point at which the average coverage value increases above a chosen value (**edgeSets.colonyThresh**). If the current timepoint is before the user-defined time of confluence (**edgeSets.confluenceFrame**), the edge is chosen as the *first* instance this crossing happens (moving from the colony exterior to the colony interior). if the current timepoint is after confluence, the edge is chosen as the *final* instance of this crossing).
   6. The median edge position is then propogated to the next timepoint and used as the starting point to find the edge in the next set of strips. The algorithm can only look a fixed distance away from the previous timepoint's edge position (**edgeSets.sampleWindow**).
   
3. Assuming you have set the output to be verbose, you should see plots similar to the following appearing at this point:

<p align="center">
  <img src="https://raw.githubusercontent.com/Pseudomoaner/cellsOnEdge/master/Images/EdgeProfile2.PNG" alt="EdgeCapture"/>
</p>

   Each red circle indicates the detected position of the colony edge in each image strip, while the connecting lines indicate the estimated profile of the edge. It is best to monitor these plots continuously as the script runs to ensure that the edge detection is working properly. If not, cancel the script and adjust the analysis parameters.

4. The script now saves the edge coordinates as the variables edgeYs and edgeXs in the root directory, along with the segmented colony images. The file is called 'ExtractedProfiles.mat'. It also plots timecourses of the edge position and colony expansion rate, along with the fitted model:

<p align="center">
  <img src="https://raw.githubusercontent.com/Pseudomoaner/cellsOnEdge/master/Images/EdgePlots.PNG" alt="EdgePositionAndExpansionRate"/>
</p>

### Part 2: Finding the packing fraction of regions
Now that the edge position has been found, the density of the colony at different positions can be found. Run the script PlotEdgePackingFractions.m to do this.

This script will find the density of the colony at two positions: the front (the region just behind the leading edge), and the 'homeland' (the position just behind the edge at the beginning of imaging). In the image below, the front is represented by the green rectangles while the edge is represented by the purple rectangles:

<p align="center">
  <img src="https://raw.githubusercontent.com/Pseudomoaner/cellsOnEdge/master/Images/EdgeRegions.PNG" alt="EdgeRegions"/>
</p>

These rectanges indicate the region over which the packing fraction is calculated. Note that the 'front' window undulates to match the shape of the colony edge, preventing the inclusion of regions outside of the colony proper.

Once the script has finished running, plots indicating the density of the front and homeland will be generated:

<p align="center">
  <img src="https://raw.githubusercontent.com/Pseudomoaner/cellsOnEdge/master/Images/PackingFractionPlots.PNG" alt="EdgePacking"/>
</p>

On the left, the average packing fraction within each region at each timepoint is shown. On the right, the spatial composition of the 'front' region over time is shown.

The data underlying these plots is saved as the 'frontPackingFractions' and 'stationaryPackingFractions' variables in the file 'PackingFractions.mat'.

### Part 3: Finding the composition of monolayer regions
If your colony consists of two separate populations of cells marked with different fluorescent labels, you can also measure the relative number of each cell type in monolayer regions using the PlotEdgeMonolayerFluoProfiles.m script.

Similar to finding the packing fraction, once this script has finished running plots indicating the composition of the front and homeland will be generated:

<p align="center">
  <img src="https://raw.githubusercontent.com/Pseudomoaner/cellsOnEdge/master/Images/CompositionPlots.PNG" alt="EdgeComposition"/>
</p>

On the left, the average composition of the two regions within each region at each timepoint is shown. On the right, the spatial composition of the 'front' region over time is shown.

The data underlying these plots will also be saved as the 'YFPProfile', 'CFPProfile', 'YFPStationaryProfile' and 'CFPStationaryProfile' variables in the file 'FluoProfiles.mat'.

### Part 4: Finding the composition of multi-layered regions
Unfortunately, the algorithm of PlotEdgePackingFractions.m is unable to accurately assign populations if your system becomes multi-layered. This is because it fundamentally assumes that each pixel corresponds either to one cell from population 1, one cell from population 2, or no cell at all - it is unable to accurately quantify regions in which there are multiple layers of a single population, or a mixture of the two populations.

For datasets where this is the case, instead use the PlotEdgeMultilayerFluoProfiels.m script. This script estimates the fluorescence of individual cells in the two channels using the outermost rim of cells, which generally remains single-layered even in multi-layered datasets. These fluorescence values are then used to estimate the depth of the two cell types in each pixel of the image, allowing more accurate composition analysis.



## References

**Bacteria solve the problem of crowding by moving slowly**, Oliver J. Meacock, Amin Doostmohammadi, Kevin R. Foster, Julia M. Yeomans and William M. Durham.
