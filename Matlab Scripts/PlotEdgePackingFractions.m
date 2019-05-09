%Finds and plots the spatial profile of system density
%the front and homeland of an expanding colony, imaged as a tile of
%individual fields of view. Ensure that you have run 
%PlotEdgeCoordsWithLinearToLinearModels.m or
%PlotEdgeCoordsWithExponentialToLinearModel.m beforehand to extract the
%edge position.

clear all
close all

set(0,'defaulttextinterpreter','latex')

Root = 'C:\Users\olijm\Desktop\cellOnEdgeTest\'; %Root directory where your images are stored

saveTracePath = [Root,'ExtractedProfiles.mat'];
load(saveTracePath,'edgeYs','edgeXs','Times','frameNos','imAngle','pxSize','dt','imOverlap','rootStem','imStem','imNos')

halfPaths = cell(length(imNos),1);

for i = 1:length(imNos)
    halfPaths{i} = [Root,sprintf(rootStem,imNos(i))];
end

savePackingPath = [Root,'PackingFractions.mat'];

reprocess = false;

%Second loop finds the packing fraction profile from the colony edge
if reprocess || ~exist(savePackingPath,'file')
    edgeYs = edgeYs(:,2:end); %Get rid of first time point (which is bound to be dodgy)
    Times = Times(2:end);
    frameNos = frameNos(2:end);
    noFrames = size(frameNos,2) - 1;
    
    maxYfromEdge = -50; %Image will be broken into strips from the colony edge to this position in the y-direction. In um.
    stationaryWindowEdge = median(edgeYs(:,1))-10; %In um - the stationary position of the front of the window for the stationary window (not the moving window)
   
    imageThresh = 0.35; %Once the edge has been found (in the previous step), need to keep this constant for results to be comparible between experiments
    bleachRate = 0;
    
    frontPackingFractions = getColonyPackingFractions(halfPaths,imStem,frameNos,imAngle,imageThresh,bleachRate,pxSize,maxYfromEdge,imOverlap,edgeXs,edgeYs,true);
    stationaryPackingFractions = getColonyPackingFractions(halfPaths,imStem,frameNos,imAngle,imageThresh,bleachRate,pxSize,maxYfromEdge,imOverlap,edgeXs,ones(size(edgeYs))*stationaryWindowEdge,false);
    
    save(savePackingPath)
else
    load(savePackingPath)
end

frontPackingFractions = packingFractions;

figure(3)
plot(Times,mean(frontPackingFractions,1),'r','LineWidth',2)
hold on
plot(Times,mean(stationaryPackingFractions,1),'c','LineWidth',2)
xlabel('\textsf{Time (mins)}','FontSize',15)
ylabel('\textsf{Areal packing fraction}','FontSize',15)
legend('Front','Homeland','Location','NorthWest')
ax = gca;
ax.LineWidth = 2;

figure(4)
imagesc(Times,[0,abs(maxYfromEdge)],medfilt2(frontPackingFractions,[3,3]))
caxis([0,1.0])
cbar = colorbar;
ylabel(cbar, 'Local packing fraction')
cbar.LineWidth = 2;
xlabel('\textsf{Time (mins)}','FontSize',15)
ylabel('\textsf{Distance from Edge ($\mu m$)}','FontSize',15)
ax = gca;
ax.LineWidth = 2;