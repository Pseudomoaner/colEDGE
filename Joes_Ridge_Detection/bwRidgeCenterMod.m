function bwridge = bwRidgeCenterMod(I,  scales, valuethresh)

%% extract the centerline of scale-space valleys
[M_mag, M_deg] = im_scalablehess2( I, scales, 'dark');
% M_mag = nonmaxsup(M_mag, M_deg, 1.5);
if isempty(valuethresh)
    valuethresh = graythresh( M_mag );
end
bwridge = im2bw( M_mag, valuethresh );
% bwridge = bwmorph(bwridge, 'thin', Inf);
% bwridge = bwareaopen(bwridge, 100);
% bwridge = imdilate(bwridge, ones(2,2));
