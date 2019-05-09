function frame = stitchFrame(ImRoots,ImStem,imAngle,imOverlap,frameIndex,normalize)
%STITCHFRAME stitches together the tiled image in the given directory.
%   
%   INPUTS:
%       -ImRoots: The path to each directory containing each image series.
%       Should be a Nx1 cell array of strings, where N is the total number
%       number of fields of view in the tile.
%       -ImStem: The format of the filename each image is saved as
%       -imAngle: The angle (in degrees) through which each sub-image
%       should be rotated to ensure the colony image is contiguous 
%       -imOverlap: The percentage overlap with which adjacent images in
%       the tile overlap.
%       -frameIndex: The current timepoint you wish to stitch the image
%       together for
%       -normalise: Whether to set the range of each sub-image in the tile
%       to vary between 0 and 1.
%
%   OUTPUTS:
%       -frame: The stitched image.
%
%   Author: Oliver J. Meacock, (c) 2019

frame = [];
for j = 1:length(ImRoots)
    framePath = [ImRoots{j},sprintf(ImStem,frameIndex)];
    frameTmp = double(imread(framePath));
    if normalize %If true, will set range of all sub-frames to be equal, varying between 0 and 1.
        maxFrameTmp = prctile(frameTmp(:),99.99);
        minFrameTmp = prctile(frameTmp(:),0.01);
        frameTmp = (frameTmp - minFrameTmp)./(maxFrameTmp-minFrameTmp);
        frameTmp(frameTmp > 1) = 1;
    end
    
    frameTmp = imrotate(frameTmp,imAngle);
    
    %Take into account the overlap between scenes
    halfOverlap = round((imOverlap/200)*size(frameTmp,1)); %ImageOverlap is in percent
    if j == length(ImRoots)
        frameTmp = frameTmp(halfOverlap:end,:);
    elseif j ==  1
        frameTmp = frameTmp(1:end-halfOverlap,:);
    else
        frameTmp = frameTmp(halfOverlap:end-halfOverlap,:);
    end
    
    frame = [frame;frameTmp];
end