function profile = getColonyPackingFractions(stitchPath,imStem,frameNos,pxSize,maxYfromEdge,edgeXs,edgeYs,verbose)
%GETCOLONYPACKINGFRACTIONS calculates the proportion of thresholded pixels 
%in a window just behind the leading edge of the spreading colony.
%
%   INPUTS:
%       -stitchPath: The path to the directory containing the pre-stitched,
%       pre-segmented colony images, generated during colony edge
%       localisation.
%       -imStem: The format of the filename each image is saved as.
%       -frameNos: The index of each field of view. Array should be ordered in such
%       a way that the element with the smallest array index corresonds to
%       the field of view index that is at the bottom of the stitched image
%       (closest to the coordinate (0,0) in whatever coordinate system you
%       wish to use).
%       -pxSize: The physical size of each pixel.
%       -maxYfromEdge: The distance back from the input edge position each
%       window should be drawn. In physical units.
%       -edgeYs: Y-coordinates of the colony edge at each timepoint. In
%       physical units.
%       -edgeXs: X-coordiantes of each strip at all timepoints. In physical
%       units.
%       -verbose: Whether to plot the current estimate of the edge position
%       as the code runs.
%
%   OUTPUTS:
%       -profile: Packing fraction at each position along the y-axis and at
%       each timepoint.
%
%   Author: Oliver J. Meacock, (c) 2019

pixelStripLength = abs(round(maxYfromEdge/pxSize));
profile = zeros(pixelStripLength,size(edgeYs,2));

for i = 1:size(frameNos,2)
    frameNo = frameNos(i);
    seg = double(imread([stitchPath,sprintf(imStem,frameNo)]))/255;
    
    %Do image opening to get rid of boundaries between densely packed
    %cells
    se = strel('disk',4);
    seg = imdilate(seg,se);
    seg = imerode(seg,se);
    
    %Extract strip, position 0 of the strip being the leading edge.
    stackMean = getLeadingEdgePixelProfile(seg,round(edgeXs/pxSize),round(edgeYs(:,i)/pxSize),pixelStripLength,verbose);
    
    if maxYfromEdge < 0
        stackMean = flip(stackMean);
    end
    
    profile(:,i) = stackMean;
    
    if verbose
        figure(1)
        hold off
        figure(2)
        plot(stackMean)
        pause(0.1)
    else
        progressbar(i/size(frameNos,2))
    end
end