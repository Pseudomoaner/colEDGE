function [edgeYs,edgeXs,times] = GetColonyEdgeCoords(BfImgPaths,imStem,stitchPath,frameNos,pxSize,dt,stitchSets,edgeSets,segSets)
%GETCOLONYEDGECOORDS calculates the location of the leading edge of a
%colony tile.
%
%Ensure that when looking at these images in Fiji, the colony starts from
%the bottom of the image and works its way upwards. BfImgPaths are assumed
%to be subfields of the entire image, from bottom up.
%
%   INPUTS:
%       -BfImgPaths: The path to each directory containing each image series.
%       Should be a Nx1 cell array of strings, where N is the total number
%       number of fields of view in the tile.
%       -ImStem: The format of the filename each image is saved as.
%       -stitchPath: The path to which the stitched colony segmentations
%       should be saved.
%       -frameNos: The index of each field of view. Array should be ordered in such
%       a way that the element with the smallest array index corresonds to
%       the field of view index that is at the bottom of the stitched image
%       (closest to the coordinate (0,0) in whatever coordinate system you
%       wish to use).
%       -pxSize: The physical size of each pixel.
%       -dt: The physical time gap between each pair of timepoints.
%       -stitchSets: Settings used to define the geometry of the stitching
%       of different fields of view.
%       -edgeSets: Settings used to determine the edge of the colony in
%       each strip.
%       -segSets: Settings used to perform the segmentation of the
%       spreading colony.
%
%   OUTPUTS:
%       -edgeYs: Y-coordinates of the colony edge at each timepoint. In
%       physical units.
%       -edgeXs: X-coordiantes of each strip at all timepoints. In physical
%       units.
%       -times: Times at which edge positions have been found. In physical
%       units.
%
%   Author: Oliver J. Meacock, (c) 2019

edgeYs = zeros(edgeSets.noBins,size(frameNos,2));

for i = 1:size(frameNos,2)
    frameNo = frameNos(i);
    
    %Start by stitching together the frame
    frameSeg = stitchFrame(BfImgPaths,imStem,[],frameNo,stitchSets,segSets);
    
    imwrite(frameSeg,[stitchPath,filesep,sprintf(imStem,frameNo)])
    
    binEdges = round(linspace(0,1,edgeSets.noBins+1)*(size(frameSeg,2)-1))+1; %Edges of the bins, defined perpendicular to edge of the colony.
    
    %Find the approximate location of the edge based on previous time point
    if frameNo == frameNos(1)
        minYCoord = 2*edgeSets.smoothWindow;
        maxYCoord = round(size(frameSeg,1)/3) - 2*edgeSets.smoothWindow;
        frameSeg = frameSeg(1:round(size(frameSeg,1)/3),:);
    else
        avgYCoord = median(edgeYs(:,i-1))/pxSize;
        minYCoord = max(avgYCoord - edgeSets.sampleWindow,edgeSets.smoothWindow*2);
        maxYCoord = min(avgYCoord + edgeSets.sampleWindow,size(frameSeg,1) - edgeSets.smoothWindow*2);
    end
    
    %Then calculate the profile of each strip of the image perpendicular to the leading edge (each bin)
    for j = 1:edgeSets.noBins
        bin = frameSeg(:,binEdges(j):binEdges(j+1));
        binProfile = mean(bin,2);
        binProfile = (binProfile-min(binProfile(:)))./(max(binProfile(:))-min(binProfile(:)));
        
        binProfile = smooth(binProfile,edgeSets.smoothWindow);
        binProfile = binProfile(minYCoord:maxYCoord);
        
        threshCrosses = find(diff(binProfile>edgeSets.colonyThresh) == -1);
        
        if ~isempty(threshCrosses)
            if i < edgeSets.confluenceFrame 
                edgeYs(j,i) = ((threshCrosses(end) + minYCoord)*pxSize)-edgeSets.edgeOffset;
            else
                edgeYs(j,i) = ((threshCrosses(1) + minYCoord)*pxSize)-edgeSets.edgeOffset;
            end
        elseif i > 1
            edgeYs(j,i) = edgeYs(j,i-1); %Use the previous timepoint's position if this one can't be found.
        else
            edgeYs(j,i) = edgeYs(j-1,i);
        end
    end

    edgeXs = ((binEdges(1:end-1) + binEdges(2:end))/2)*pxSize;
    
    %Plot results on image if want to verify output
    if edgeSets.verbose > 0
        displayWindow = 600;
        currMid = median(edgeYs(:,i))/pxSize; %Current Y location in pixels
        minYdisp = max(currMid - displayWindow,1);
        maxYdisp = min(currMid + displayWindow,size(frameSeg,1));
        
        f1 = figure(2);
        f1.Units = 'normalized';
        f1.Position = [0.2,0.2,0.5,0.7];
        
        imshow(frameSeg(minYdisp:maxYdisp,:),[])
        hold on
        plot(edgeXs/pxSize,(edgeYs(:,i)/pxSize) - currMid + displayWindow,'r-o','LineWidth',2,'MarkerSize',10)
        disp(frameNo*dt)
        pause(0.1);
        hold off
    else
        progressbar(i/size(frameNos,2))
    end
end

times = frameNos*dt;
