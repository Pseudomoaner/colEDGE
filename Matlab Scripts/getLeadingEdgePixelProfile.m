function profile = getLeadingEdgePixelProfile(imFrame,edgeX,edgeY,pixelStripLength,verbose)
%GETLEADINGEDGEPIXELPROFILE extracts the pixel profile in the given
%stitched image at the given coordinates by measuring the pixel intensities
%in a strip running behind each given coordinate.
%
%
%   INPUTS:
%       -imFrame: A stitched image.
%       -edgeYs: Y-coordinates of the colony edge. In pixels.
%       -edgeXs: X-coordiantes of each strip. In pixels.
%       -pixelStripLength: The distance each strip should run backwards
%       from the input coordinates. In pixels.
%       -verbose: Whether to plot output of code as it runs.
%
%   OUTPUTS:
%       -profile: The pixel profile of the image along the y-axis, averaged
%       along the x-axis.
%
%   Author: Oliver J. Meacock, (c) 2019

xSeparation = edgeX(2)-edgeX(1); %Size of bins in the x-direction

if verbose
    figure(1)
    imshow(imFrame(max(min(edgeY)-pixelStripLength,1):max(edgeY),:))
    hold on
end

stripStack = zeros(pixelStripLength,length(edgeX));
for j = 1:length(edgeX)
    maxY = min(edgeY(j),size(imFrame,1));
    minY = max(maxY - pixelStripLength + 1,1);
    
    maxX = (edgeX(j) + xSeparation/2)-1;
    minX = (edgeX(j) - xSeparation/2);
    
    imStrip = imFrame(minY:maxY,minX:maxX);
    
    meanStrip = mean(imStrip,2);
    
    %The strip might be shorter than the specified amount, so pad it with NaNs.
    if length(meanStrip) < pixelStripLength
        meanStrip = padarray(meanStrip,[pixelStripLength-length(meanStrip),0],NaN,'pre');
    end
    stripStack(:,j) = meanStrip;
    
    if verbose
        rectangle('Position',[minX,minY,maxX-minX,maxY-minY],'LineWidth',3,'EdgeColor','r')
    end
end
profile = nanmean(stripStack,2);