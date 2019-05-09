function [edgeYs,edgeXs,times] = GetColonyEdgeCoords(BfImgPaths,imStem,frameNos,pxSize,dt,segThresh,colonyThresh,sampleWindow,confluenceTime,imAngle,imOverlap,noBins,verbose)
%GETCOLONYEDGECOORDS calculates the location of the leading edge of a
%colony tile.
%
%Ensure that when looking at these images in Fiji, the colony starts from
%the bottom of the image and works its way upwards. BfImgPaths are assumed
%to be subfields of the entire image, from bottom up.
%
%The original image should be processed so that cells are white on a dark
%background. Make as binary as possible.
%
%   INPUTS:
%       -BfImgPaths: The path to each directory containing each image series.
%       Should be a Nx1 cell array of strings, where N is the total number
%       number of fields of view in the tile.
%       -ImStem: The format of the filename each image is saved as.
%       -frameNos: The index of each field of view. Array should be ordered in such
%       a way that the element with the smallest array index corresonds to
%       the field of view index that is at the bottom of the stitched image
%       (closest to the coordinate (0,0) in whatever coordinate system you
%       wish to use).
%       -pxSize: The physical size of each pixel.
%       -dt: The physical time gap between each pair of timepoints.
%       -segThresh: The global threshold value used to binarise the
%       stitched brightfield image.
%       -colonyThresh: The value above which the local packing fraction in
%       each strip should rise above to trigger detection of a potential colony edge.
%       -sampleWindow: The maximum distance above and below the previous 
%       position the algorithm should look relative to previous position to
%       find the next position of the colony edge.
%       -confluenceTime: Time at which confluence at the colony edge is
%       achieved.
%       -imAngle: The angle (in degrees) through which each sub-image
%       should be rotated to ensure the colony image is contiguous 
%       -imOverlap: The percentage overlap with which adjacent images in
%       the tile overlap.
%       -noBins: The number of strips each stitched image should be split
%       into.
%       -verbose: Whether to plot the current estimate of the edge position
%       as the code runs.
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

smoothWindow = 50;
edgeOffset = 0; %Number of um the actual colony edge is behind the maximum of the dark strip.

edgeYs = zeros(noBins,size(frameNos,2));

for i = 1:size(frameNos,2)
    frameNo = frameNos(i);
    
    %Start by stitching together the frame
    frame = stitchFrame(BfImgPaths,imStem,imAngle,imOverlap,frameNo,false);
    frame = medfilt2(frame,[2,2]);
    frameSeg = frame>segThresh;
    
    %Do median filtering to get rid of random spots
    frameSeg = medfilt2(frameSeg,[2,2]);
    
    %Do image opening to get rid of boundaries between densely packed
    %cells
    se = strel('disk',4);
    frameSeg = imdilate(frameSeg,se);
    frameSeg = imerode(frameSeg,se);
    
    binEdges = round(linspace(0,1,noBins+1)*(size(frame,2)-1))+1; %Edges of the bins, defined perpendicular to edge of the colony.
    
    %Find the approximate location of the edge based on previous time point
    if frameNo == frameNos(1)
        minYCoord = 2*smoothWindow;
        maxYCoord = round(size(frame,1)/3) - 2*smoothWindow;
        frameSeg = frameSeg(1:round(size(frameSeg,1)/3),:);
    else
        avgYCoord = median(edgeYs(:,i-1))/pxSize;
        minYCoord = max(avgYCoord - sampleWindow,smoothWindow*2);
        maxYCoord = min(avgYCoord + sampleWindow,size(frame,1) - smoothWindow*2);
    end
    
    %Then calculate the profile of each strip of the image perpendicular to the leading edge (each bin)
    for j = 1:noBins
        bin = frameSeg(:,binEdges(j):binEdges(j+1));
        binProfile = mean(bin,2);
        binProfile = (binProfile-min(binProfile(:)))./(max(binProfile(:))-min(binProfile(:)));
        
        binProfile = smooth(binProfile,smoothWindow);
        binProfile = binProfile(minYCoord:maxYCoord);
        
        threshCrosses = find(diff(binProfile>colonyThresh) == -1);
        
        if ~isempty(threshCrosses)
            if i < confluenceTime 
                edgeYs(j,i) = ((threshCrosses(end) + minYCoord)*pxSize)-edgeOffset;
            else
                edgeYs(j,i) = ((threshCrosses(1) + minYCoord)*pxSize)-edgeOffset;
            end
        elseif i > 1
            edgeYs(j,i) = edgeYs(j,i-1); %Use the previous timepoint's position if this one can't be found.
        else
            edgeYs(j,i) = edgeYs(j-1,i);
        end
    end

    edgeXs = ((binEdges(1:end-1) + binEdges(2:end))/2)*pxSize;
    
    %Plot results on image if want to verify output
    if verbose > 0
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
    end
end

times = frameNos*dt;