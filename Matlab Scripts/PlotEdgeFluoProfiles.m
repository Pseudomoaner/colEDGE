clear all
% close all

Root = 'D:\Co-culture subsurface colonies\15_08_19_0pt8pct_Agr_CoCulture\Rep1\';

saveTracePath = [Root,'ExtractedProfiles.mat'];
load(saveTracePath,'edgeYs','imNos','edgeXs','stitchSets','segSets','imStem','Times','frameNos','stitchPath','pxSize','dt')
segSets.segment = false;

YFPStem = 'Chan1_Scene%i';
CFPStem = 'Chan2_Scene%i';

YFPsave = 'stitchedYFP';
CFPsave = 'stitchedCFP';

if ~exist([Root,YFPsave],'dir')
    mkdir([Root,YFPsave])
end
if ~exist([Root,CFPsave],'dir')
    mkdir([Root,CFPsave])
end

halfYFPPaths = cell(length(imNos),1);
halfCFPPaths = cell(length(imNos),1);

for i = 1:length(imNos)
    halfYFPPaths{i} = [Root,sprintf(YFPStem,imNos(i))];
    halfCFPPaths{i} = [Root,sprintf(CFPStem,imNos(i))];
end

saveFluoPath = [Root,'FluoProfiles.mat'];

reprocess = 1;
verbose = false;

if ~exist(saveFluoPath,'file') || reprocess > 0
    edgeYs = edgeYs(:,2:end); %Get rid of first time point (which is bound to be dodgy)
    Times = Times(2:end);
    frameNos = frameNos(2:end);
    noFrames = size(frameNos,2);
    
    maxYfromEdge = -100; %in um - the distance back from the leading edge to include in the window
    stationaryWindowEdge = median(edgeYs(:,1))-10; %In um - the stationary position of the front of the window for the stationary window (not the moving window)
    
    pixelStripLength = abs(round(maxYfromEdge/pxSize));
    CFPProfile = zeros(pixelStripLength,size(edgeYs,2));
    YFPProfile = zeros(pixelStripLength,size(edgeYs,2));
    
    CFPStationaryProfile = zeros(pixelStripLength,size(edgeYs,2));
    YFPStationaryProfile = zeros(pixelStripLength,size(edgeYs,2)); 
    
    %I was originally intending to use getColonyPixelFractions here, but we need to do binary operations on combinations of images so I've just copied the code here instead.
    for i = 1:size(edgeYs,2)
        if reprocess == 2 || ~exist([Root,YFPsave,filesep,sprintf(imStem,i)],'file')
            BFframeSeg = double(imread([stitchPath,sprintf(imStem,frameNos(i))])) > 0.5;
            
            YFPframe = stitchFrame(halfYFPPaths,imStem,frameNos(i),stitchSets,segSets);
            YFPframe = medfilt2(YFPframe,[2,2]);
            CFPframe = stitchFrame(halfCFPPaths,imStem,frameNos(i),stitchSets,segSets);
            CFPframe = medfilt2(CFPframe,[2,2]);
            
            FluoRatio = YFPframe./CFPframe;
            
            %Split the two populations according to a threshold calculated from
            %the ratio of the first frame
            if i == 1 %The ratio can be a bit dodgy for the first few frames if fluorescence is still accumulating in the cells.
                data = FluoRatio(BFframeSeg(:));
                model = fitgmdist(data,2); %Mixed Gaussian model of ratiometric pixel distribution in first frame
                idx = cluster(model,data);
                cluster1 = data(idx == 1);
                cluster2 = data(idx == 2);
                fluoThresh = min(max(cluster1),max(cluster2));
            end
            
            YFPareas = FluoRatio>fluoThresh;
            CFPareas = FluoRatio<fluoThresh;
            
            YFPareas(~BFframeSeg) = 0;
            CFPareas(~BFframeSeg) = 0;
            
            imwrite(YFPareas,[Root,YFPsave,filesep,sprintf(imStem,i)],'Compression','none')
            imwrite(CFPareas,[Root,CFPsave,filesep,sprintf(imStem,i)],'Compression','none')
        else
            YFPareas = imread([Root,YFPsave,filesep,sprintf(imStem,i)]);
            CFPareas = imread([Root,CFPsave,filesep,sprintf(imStem,i)]);
        end
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
%         
%         BFimageThresh = BFimageThresh * (1-BFbleachRate);
%         FluoImageThresh = FluoImageThresh * (1-FluoBleachRate);
%         
        %Show image with sampling windows overlayed on top
        if verbose
            xSeparation = (edgeXs(2)-edgeXs(1))/pxSize; %Size of bins in the x-direction
            
            figure(1)
            cla
            hold on
            imshow(BFframeSeg,[]);
            for j = 1:size(edgeYs,1)
                maxY = min((edgeYs(j,i)/pxSize),size(BFframeSeg,1));
                minY = max(maxY - pixelStripLength + 1,1);
                
                maxX = ((edgeXs(j)/pxSize) + xSeparation/2)-1;
                minX = ((edgeXs(j)/pxSize) - xSeparation/2);
                rectangle('Position',[minX,minY,maxX-minX,maxY-minY],'LineWidth',1,'EdgeColor',[0,0.7,0])
                
                maxY = stationaryWindowEdge/pxSize;
                minY = (stationaryWindowEdge/pxSize) - pixelStripLength + 1;
                rectangle('Position',[minX,minY,maxX-minX,maxY-minY],'LineWidth',1,'EdgeColor',[0.5,0,1])
            end
            export_fig(sprintf('C:\\Users\\olijm\\Desktop\\TmpImg2\\Rand\\Frame_%04d',i),'-tif','-m3')
        else
            progressbar(i/size(edgeYs,2))
        end
    end
    save(saveFluoPath)
else
    load(saveFluoPath);
end

figure(1)
cla
imagesc(Times,[0,abs(maxYfromEdge)],medfilt2(YFPProfile./(YFPProfile + CFPProfile),[3,3]))
colorbar
xlabel('\textsf{Time (mins)}','FontSize',15)
ylabel('\textsf{Distance from Edge ($\mu m$)}','FontSize',15)

figure(2)
cla
hold on
plot(Times,mean(YFPProfile,1)./(mean(CFPProfile,1)+mean(YFPProfile,1)),'r','LineWidth',2)
plot(Times(1:26),mean(YFPStationaryProfile(:,1:26),1)./(mean(CFPStationaryProfile(:,1:26),1)+mean(YFPStationaryProfile(:,1:26),1)),'c','LineWidth',2)
xlabel('\textsf{Time (mins)}','FontSize',15)
ylabel('\textsf{YFP to CFP cell proportion}','FontSize',15)
legend('Front','Homeland','Location','NorthEast')
ax = gca;
ax.LineWidth = 2;