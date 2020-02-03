%Finds the colony edge in a series of brightfield images corresponding to
%a tile covering the expansion of a bacterial colony. Then fits a model
%consisting of an exponentially increasing expansion rate, smoothly
%transitioning into a constant expansion rate later.

clear all
close all

set(0,'defaulttextinterpreter','latex')

Root = 'C:\Users\ph1ojm\Desktop\ColEDGE_Tester\20_10_17_PilH_YFP_WT_CFP_ColonyExpansion\'; %Root directory where your images are stored

rootStem = 'BF_Scene%i'; %Format of the folder name where your brightfield images are stored
imStem = '\\Frame_%i.tif'; %Format of the file name of each frame of your image sequence
imNos = 1:16; %Indices of the fields of view associated with this dataset. Note the order of input matters - the field of view to the right of this array will be the bottom most image in the stitching
reprocess = false; %Whether to reprocess the data from scratch, or load pre-processed data.
stitchStem = 'stitchedSegmentations';

halfPaths = cell(length(imNos),1);

for i = 1:length(imNos)
    halfPaths{i} = [Root,sprintf(rootStem,imNos(i))];
end

stitchPath = [Root,stitchStem];
if ~exist(stitchPath,'dir')
    mkdir(stitchPath)
end

saveTracePath = [Root,'ExtractedProfiles.mat'];

if ~exist(saveTracePath,'file') || reprocess
    frameNos = 1:47; %Frames that should be included in the analysis (can cut out frames in which, e.g., the colony edge has expended out of frame)
    
    pxSize = 0.1024; %In microns per pixel
    dt = 10; %In minutes per frame
    
    segSets.segment = true; %Whether to perform segmentation on the stitched images
    segSets.covThresh = 1.2; %Threshold local coefficient of variation to detect object - 2 for new Epi (63x w. 2x2 binning), 1.2 for old Epi (63x),
    segSets.nHood = 15; %Size of local neighbourhood over which COV should be calculated (in pixels) - 15 for 63x images
    segSets.ridgeScale = 18; %Scale of ridges to be detected - 18 for 63x images, 15 for 20x images
    segSets.ridgeThresh = 0.08; %Ridge detection threshold
    segSets.areaThresh = 15; %Minimum area (in pixels) an object needs to be to be included in segmentation. 50 for 63x images
    
    stitchSets.imAngle = 0; %Amount(in degrees) that the image should be rotated (counterclockwise) so the colony rises from the bottom of the frame.
    stitchSets.imOverlap = 7; %Percent overlap between frames
    stitchSets.normalise = false; %Whether each frame in the stitched tile should be scaled to cover the same range of pixel values.
    
    edgeSets.colonyThresh = 0.1; %Proportion of the binary image (within a strip) that must be white (on average) for it to be classified as part of the colony - 0.02 for new Epi (63x w. 2x2 binning), 0.2 for old Epi (63x),
    edgeSets.confluenceFrame = 58; %Frame at which the colony becomes confluent (packing fraction of the whole colony behind the edge becomes greater than colonyThresh). If unsure, make larger than the maximum frame index.
    edgeSets.sampleWindow = 1200; %Maximum distance forwards and backwards algorithm should look relative to previous position to find the next position of the colony edge.
    edgeSets.smoothWindow = 50; %The span of the smoothing function used to clean up the packing fraction profile for each strip.
    edgeSets.noBins = 20; %Number of strips that each image should be split up into.
    edgeSets.verbose = 1; %Amount of visualisation that should take place. 0 = none, 1 = show edge. Usually best to leave set to 1 to manually verify edge tracking.
    edgeSets.edgeOffset = 2; %Number of um the actual colony edge is behind the maximum of the dark strips - 4 for new Epi (63x w. 2x2 binning), 2 for old Epi (63x),
    
    [edgeYs,edgeXs,Times] = GetColonyEdgeCoords(halfPaths,imStem,stitchPath,frameNos,pxSize,dt,stitchSets,edgeSets,segSets);
    save(saveTracePath)
else
    load(saveTracePath)
end

edgeYs = edgeYs(:,2:end); %Get rid of first time point (which is bound to be dodgy)
Times = Times(2:end);
noFrames = size(frameNos,2) - 1;

medianLocs = median(edgeYs,1);
medianLocs = medianLocs - medianLocs(1);

ft = fittype(@(a,b,d,xTrans,x) exponentialToLinear(x,a,b,d,xTrans));
f = fit(Times',medianLocs',ft,'StartPoint',[1,0.03,2,edgeSets.confluenceFrame * dt],'Lower',[0,0,0,frameNos(1) * dt],'Upper',[Inf,1,Inf,frameNos(end) * dt]);

initX = Times(1:round(f.xTrans/dt));
lateX = Times(round(f.xTrans/dt)+1:noFrames);
initY = medianLocs(1:round(f.xTrans/dt));
lateY = medianLocs(round(f.xTrans/dt) + 1:noFrames);

modelExpY = f.a * exp(Times * f.b);
yTrans = f.a * exp(f.b*f.xTrans);
intercept = yTrans - f.d*f.xTrans;
modelLinY = intercept + (f.d * Times);

figure(1)
plot(Times,medianLocs,'LineWidth',2)
hold on
plot(Times,modelExpY,'r--','LineWidth',2)
plot(Times,modelLinY,'k--','LineWidth',2)

ylim([0,3500])
xlabel('\textsf{Time / mins}','FontSize',15)
ylabel({'\textsf{Distance from initial}' '\textsf{edge position / $\mu m$}'},'FontSize',15)
legend('Data','Initial Exponential','Late Linear','Location','NorthWest')

figure(2)
semilogy(Times,medianLocs,'LineWidth',2)
hold on
semilogy(Times,modelExpY,'r--','LineWidth',2)
semilogy(Times,modelLinY,'k--','LineWidth',2)

ylim([0,3500])
xlabel('\textsf{Time / mins}','FontSize',15)
ylabel({'\textsf{Distance from initial}' '\textsf{edge position / $\mu m$}'},'FontSize',15)
legend('Data','Initial Exponential','Late Linear','Location','NorthWest')

%Do plotting for colony expansion rate (the differential of the colony position)
modelDiffExpY = f.a * f.b * exp(Times(11:end-10) * f.b);
modelDiffLinY =  f.d * ones(size(Times(11:end-10)));
dxdt = diff(smooth(medianLocs,10))/dt;

figure(3)
plot(Times(11:end-10),dxdt(10:end-10),'LineWidth',2)
hold on
plot(Times(11:end-10),modelDiffExpY,'r--','LineWidth',2)
plot(Times(11:end-10),modelDiffLinY,'k--','LineWidth',2)

ylim([0,18.0])
ylabel('\textsf{Colony Expansion Rate ($\mu m$/min)}','FontSize',15)
xlabel('\textsf{Time (mins)}','FontSize',15)
legend('Data','Initial Exponential','Late Linear','Location','NorthWest')