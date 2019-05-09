# cellsOnEdge
Analysis of location and composition of bacterial colony edges in microscopy data

## Introduction

When capturing the dynamics of bacterial colony expansion, it is often useful to use tiles of images to monitor the position and composition of the colony edge at high resolution over long periods of time. For example, the following image is taken from a 2 mm x 0.15 mm x 12 hr dataset, and has been stitched together from 20 adjacent fields of view:

While this approach can provide a combination of excellent spatial resolution and long sampling times, the resulting datasets are typically too large to analyse efficiently by hand. cellsOnEdge is a suite of tools designed to automatically find the position and composition of the colony edge within these datasets.

## Usage

### Part 0: Preprocessing
Begin in Fiji/ImageJ.

1. Brightfield images should be pre-processed to allow easy binarisation of the colony image (i.e. to allow cells and background to be distinguished with a simple threshold). They should also be split into separate frames, with each field of view stored in a separate folder. The scripts ProcessBulkBFSingleFrames63X.ijm and ProcessBulkBFSingleFrames20X.ijm perform these jobs. These two scripts contain settings for processing brightfield data from a Zeiss Observer microscope with 20X and 63X PlanApo objectives, but you will probably need to find parameters that fit your own dataset.

2. Fluorescence images should be split into separate frames, with each field of view stored in a separate folder. Despeckling is also advised. These jobs are performed by the script ProcessBulkFluo.ijm.

### Part 1: Finding the colony edge
Now switch to Matlab.



### Part 2: Finding the packing fraction

### Part 3: Finding the edge composition

## References

- Reference to my paper when it's finally out!
