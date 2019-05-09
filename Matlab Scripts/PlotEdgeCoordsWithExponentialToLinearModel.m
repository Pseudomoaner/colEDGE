%Finds the colony edge in a series of brightfield images corresponding to
%a tile covering the expansion of a bacterial colony. Then fits a model
%consisting of an exponentially increasing expansion rate, smoothly
%transitioning into a constant expansion rate later.

clear all
close all

set(0,'defaulttextinterpreter','latex')

Root = 'C:\Users\olijm\Desktop\cellOnEdgeTest\'; %Root directory where your images are stored

rootStem = 'BF_Scene%i'; %Format of the folder name where your brightfield images are stored
imStem = '\\Frame_%i.tif'; %Format of the file name of each frame of your image sequence
imNos = 1:20; %Indices of the fields of view associated with this dataset. Note the order of input matters - the field of view to the right of this array will be the bottom most image in the stitching
reprocess = false; %Whether to reprocess the data from scratch, or load pre-processed data.

halfPaths = cell(length(imNos),1);

for i = 1:length(imNos)
    halfPaths{i} = [Root,sprintf(rootStem,imNos(i))];
end

saveTracePath = [Root,'ExtractedProfiles.mat'];

if ~exist(saveTracePath,'file') || reprocess
    frameNos = 1:70; %Frames that should be included in the analysis (can cut out frames in which, e.g., the colony edge has expended out of frame)
    
    pxSize = 0.1024000; %In microns per pixel
    dt = 10; %In minutes per frame
    
    imageThresh = 220; %Global threshold to binarise image
    colonyThresh = 0.1; %Proportion of the binary image (within a strip) that must be white (on average) for it to be classified as part of the colony
    confluenceFrame = 34; %Frame at which the colony becomes confluent (packing fraction of the whole colony behind the edge becomes greater than colonyThresh). If unsure, make larger than the maximum frame index.
    sampleWindow = 800; %Maximum distance forwards and backwards algorithm should look relative to previous position to find the next position of the colony edge.
    imAngle = 0; %Amount(in degrees) that the image should be rotated (counterclockwise) so the colony rises from the bottom of the frame.
    imOverlap = 7; %Percent overlap between frames
    noBins = 10; %Number of strips that each image should be split up into.
    verbose = 1; %Amount of visualisation that should take place. 0 = none, 1 = show edge. Usually best to leave set to 1 to manually verify edge tracking.
    
    [edgeYs,edgeXs,Times] = GetColonyEdgeCoords(halfPaths,imStem,frameNos,pxSize,dt,imageThresh,colonyThresh,sampleWindow,confluenceFrame,imAngle,imOverlap,noBins,verbose);
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
f = fit(Times',medianLocs',ft,'StartPoint',[1,0.03,2,confluenceFrame * dt],'Lower',[0,0,0,frameNos(1) * dt],'Upper',[Inf,1,Inf,frameNos(end) * dt]);

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
modelDiffExpY = f.a * f.b * exp(Times(31:end-30) * f.b);
modelDiffLinY =  f.d * ones(size(Times(31:end-30)));
dxdt = diff(smooth(medianLocs,80))/dt;

figure(3)
plot(Times(31:end-30),dxdt(30:end-30),'LineWidth',2)
hold on
plot(Times(31:end-30),modelDiffExpY,'r--','LineWidth',2)
plot(Times(31:end-30),modelDiffLinY,'k--','LineWidth',2)

ylim([0,18.0])
ylabel('\textsf{Colony Expansion Rate ($\mu m$/min)}','FontSize',15)
xlabel('\textsf{Time (mins)}','FontSize',15)
legend('Data','Initial Exponential','Late Linear','Location','NorthWest')