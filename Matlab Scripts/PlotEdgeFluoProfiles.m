%Finds and plots the spatial profile of differentially labelled cells in
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

bfStem = rootStem;
YFPStem = 'Chan1_Scene%i';
CFPStem = 'Chan2_Scene%i';

halfBfPaths = cell(length(imNos),1);
halfYFPPaths = cell(length(imNos),1);
halfCFPPaths = cell(length(imNos),1);

for i = 1:length(imNos)
    halfBfPaths{i} = [Root,sprintf(bfStem,imNos(i))];
    halfYFPPaths{i} = [Root,sprintf(YFPStem,imNos(i))];
    halfCFPPaths{i} = [Root,sprintf(CFPStem,imNos(i))];
end

saveFluoPath = [Root,'FluoProfiles.mat'];

reprocess = false;
verbose = false;

if ~exist(saveFluoPath,'file') || reprocess
    edgeYs = edgeYs(:,2:end); %Get rid of first time point (which is bound to be dodgy)
    Times = Times(2:end);
    frameNos = frameNos(2:end);
    noFrames = size(frameNos,2);
    
    maxYfromEdge = -50; %in um - the distance back from the leading edge to include in the window
    stationaryWindowEdge = median(edgeYs(:,1))-10; %In um - the stationary position of the front of the window for the stationary window (not the moving window)
    
    BFimageThresh = 0.4; %Once the edge has been found (in the previous step), need to keep this constant for results to be comparible between experiments
    
    pixelStripLength = abs(round(maxYfromEdge/pxSize));
    CFPProfile = zeros(pixelStripLength,size(edgeYs,2));
    YFPProfile = zeros(pixelStripLength,size(edgeYs,2));
    
    CFPStationaryProfile = zeros(pixelStripLength,size(edgeYs,2));
    YFPStationaryProfile = zeros(pixelStripLength,size(edgeYs,2)); 
    
    %I was originally intending to use getColonyPixelFractions here, but we need to do binary operations on combinations of images so I've just copied the code here instead.
    for i = 1:size(edgeYs,2)
        BFframe = stitchFrame(halfBfPaths,imStem,imAngle,imOverlap,frameNos(i),true);
        BFframeSeg = BFframe>BFimageThresh;
        
        YFPframe = stitchFrame(halfYFPPaths,imStem,imAngle,imOverlap,frameNos(i),false);
        YFPframe = medfilt2(YFPframe,[2,2]);
        CFPframe = stitchFrame(halfCFPPaths,imStem,imAngle,imOverlap,frameNos(i),false);
        CFPframe = medfilt2(CFPframe,[2,2]);
        
        FluoRatio = YFPframe./CFPframe;
        
        %We wish to calculate a threshold based on the intensity of the
        %non-cell containing regions. 
        fluoThresh = nanmedian(FluoRatio(~BFframeSeg));
        YFPareas = FluoRatio>fluoThresh;
        CFPareas = FluoRatio<fluoThresh;
        
        YFPareas(~BFframeSeg) = 0;
        CFPareas(~BFframeSeg) = 0;
        
        %Extract strips, position 0 of the strip being the leading edge.
        YFPMean = getLeadingEdgePixelProfile(YFPareas,round(edgeXs/pxSize),round(edgeYs(:,i)/pxSize),pixelStripLength,false);
        CFPMean = getLeadingEdgePixelProfile(CFPareas,round(edgeXs/pxSize),round(edgeYs(:,i)/pxSize),pixelStripLength,false);
        
        %Extract strips, position 0 of the strip being the leading edge.
        YFPHomeMean = getLeadingEdgePixelProfile(YFPareas,round(edgeXs/pxSize),round((ones(size(edgeYs,1),1)*stationaryWindowEdge)/pxSize),pixelStripLength,false);
        CFPHomeMean = getLeadingEdgePixelProfile(CFPareas,round(edgeXs/pxSize),round((ones(size(edgeYs,1),1)*stationaryWindowEdge)/pxSize),pixelStripLength,false);
        
        if maxYfromEdge < 0
            YFPMean = flip(YFPMean);
            CFPMean = flip(CFPMean);
        end
        
        YFPProfile(:,i) = YFPMean;
        CFPProfile(:,i) = CFPMean;
        
        YFPStationaryProfile(:,i) = YFPHomeMean;
        CFPStationaryProfile(:,i) = CFPHomeMean;
        
        %Show image with sampling windows overlayed on top
        if verbose
            xSeparation = (edgeXs(2)-edgeXs(1))/pxSize; %Size of bins in the x-direction
            
            figure(1)
            cla
            hold on
            imshow(BFframe,[]);
            for j = 1:size(edgeYs,1)
                maxY = min((edgeYs(j,i)/pxSize),size(BFframe,1));
                minY = max(maxY - pixelStripLength + 1,1);
                
                maxX = ((edgeXs(j)/pxSize) + xSeparation/2)-1;
                minX = ((edgeXs(j)/pxSize) - xSeparation/2);
                rectangle('Position',[minX,minY,maxX-minX,maxY-minY],'LineWidth',1,'EdgeColor',[0,0.7,0])
                
                maxY = stationaryWindowEdge/pxSize;
                minY = (stationaryWindowEdge/pxSize) - pixelStripLength + 1;
                rectangle('Position',[minX,minY,maxX-minX,maxY-minY],'LineWidth',1,'EdgeColor',[0.5,0,1])
            end
            export_fig(sprintf('C:\\Users\\OliLocal\\Desktop\\TmpImg2\\Rand\\Frame_%04d',i),'-tif','-m3')
        end
    end
    save(saveFluoPath)
else
    load(saveFluoPath);
end

figure(4)
imagesc(Times,[0,abs(maxYfromEdge)],medfilt2(YFPProfile./(YFPProfile + CFPProfile),[3,3]))
cbar = colorbar;
ylabel(cbar, 'YFP-labelled to total cell ratio')
cbar.LineWidth = 2;
xlabel('\textsf{Time (mins)}','FontSize',15)
ylabel('\textsf{Distance from Edge ($\mu m$)}','FontSize',15)
ax = gca;
ax.LineWidth = 2;

figure(3)
hold on
plot(Times,mean(YFPProfile,1)./(mean(CFPProfile,1)+mean(YFPProfile,1)),'r','LineWidth',2)
plot(Times(1:26),mean(YFPStationaryProfile(:,1:26),1)./(mean(CFPStationaryProfile(:,1:26),1)+mean(YFPStationaryProfile(:,1:26),1)),'c','LineWidth',2)
xlabel('\textsf{Time (mins)}','FontSize',15)
ylabel('\textsf{YFP to total cell ratio}','FontSize',15)
legend('Front','Homeland','Location','NorthEast')
ax = gca;
ax.LineWidth = 2;
ax.Box = 'on';