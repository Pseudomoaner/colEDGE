function outImg = segmentColonyExpansionImgs(inImg,seNhood,seThresh,ridgeScale,ridgeThresh,areaThresh)
%SEGMENTCOLONYEXPANSIONIMGS uses a combination of the local coefficient of
%variation (COV) and ridge detection to segment the given input image.
%
%   INPUTS:
%       -inImg: The raw (usually brightfield) image.
%       -seNhood: The neighbourhood within which the local COV should be
%       calculated.
%       -seThresh: The COV threshold below which pixels should not be
%       segmented as part of the colony.
%       -ridgeScale: The predicted scale (in pixels) of the ridges between
%       cells.
%       -ridgeThresh: The ridge threshold. Reduce to increase amount of
%       image classified as ridge.
%
%   OUTPUTS:
%       -outImg: The segmented image.
% 
%   Author: Oliver J. Meacock, (c) 2019

inImg = double(inImg); %Recast data type if not appropriate.

stdImg = stdfilt(inImg,ones(seNhood));
kernel = ones(seNhood) / seNhood^2; % Mean kernel
meanImg = conv2(inImg, kernel, 'same'); % Convolve keeping size of I
seImg = stdImg./(meanImg.^0.5); %Logic here is that dividing by meanImg gives COV. However, the COV itself scales as one over the square root of the image intensity (shot noise), so multiply by that number to get just the contribution from the cells.
seSeg = seImg > seThresh;

%Get rid of 'rims' around non-segmented patches
seStre = strel('disk',(seNhood-1)/2);
seSeg = imerode(seSeg,seStre);

ridgeSeg = bwRidgeCenterMod(inImg,ridgeScale,ridgeThresh);

outImg = and(seSeg,~ridgeSeg);

outImg = bwareaopen(outImg,areaThresh);