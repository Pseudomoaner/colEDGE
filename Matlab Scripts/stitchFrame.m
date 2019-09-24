function frame = stitchFrame(ImRoots,ImStem,frameIndex,stitchSets,segSets)
frame = [];
for j = 1:length(ImRoots)
    framePath = [ImRoots{j},sprintf(ImStem,frameIndex)];
    frameTmp = double(imread(framePath));
    
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