function frame = stitchFrame(ImRoots,ImStem,flatField,frameIndex,stitchSets,segSets)
%STITCHFRAME stitches together the tiled image in the given directory.
%   
%   INPUTS:
%       -ImRoots: The path to each directory containing each image series.
%       Should be a Nx1 cell array of strings, where N is the total number
%       number of fields of view in the tile.
%       -ImStem: The format of the filename each image is saved as
%       -flatField: The flatfield (empty, time averaged) image for the
%       current channel. Can be empty.
%       -frameIndex: The current timepoint you wish to stitch the image
%       together for
%       -stitchSets: The settings to define the geometry of the image
%       stitching
%       -segSets: The settings used to define the segmentation of the
%       colony images.
%
%   OUTPUTS:
%       -frame: The stitched image.
%
%   Author: Oliver J. Meacock, (c) 2019

frame = [];
for j = 1:length(ImRoots)
    framePath = [ImRoots{j},sprintf(ImStem,frameIndex)];
    frameTmp = double(imread(framePath));
    
    if ~isempty(flatField)
        frameTmp = frameTmp./flatField;
    end
    
    if stitchSets.normalise %If true, will set range of all sub-frames to be equal, varying between 0 and 1.
        maxFrameTmp = prctile(frameTmp(:),99.99);
        minFrameTmp = prctile(frameTmp(:),0.01);
        frameTmp = (frameTmp - minFrameTmp)./(maxFrameTmp-minFrameTmp);
        frameTmp(frameTmp > 1) = 1;
    end
    
    if segSets.segment
        frameTmp = segmentColonyExpansionImgs(frameTmp,segSets.nHood,segSets.covThresh,segSets.ridgeScale,segSets.ridgeThresh,segSets.areaThresh);
    end
    
    frameTmp = imrotate(frameTmp,stitchSets.imAngle);
    
    %Take into account the overlap between scenes
    halfOverlap = round((stitchSets.imOverlap/200)*size(frameTmp,1)); %ImageOverlap is in percent
    if j == length(ImRoots)
        frameTmp = frameTmp(halfOverlap:end,:);
    elseif j ==  1
        frameTmp = frameTmp(1:end-halfOverlap,:);
    else
        frameTmp = frameTmp(halfOverlap:end-halfOverlap,:);
    end
    
    frame = [frame;frameTmp];
end