clear all
close all

%File and path setup
Root = 'C:\Users\ph1ojm\Desktop\ColEDGE_Tester\20_10_17_PilH_YFP_WT_CFP_ColonyExpansion\';

saveTracePath = [Root,'ExtractedProfiles.mat'];
load(saveTracePath,'edgeYs','imNos','edgeXs','stitchSets','segSets','imStem','Times','frameNos','stitchPath','pxSize','dt')
segSets.segment = false;

YFPStem = 'Chan1_Scene%i';
CFPStem = 'Chan2_Scene%i';

YFPsegSave = 'stitchedYFPseg';
CFPsegSave = 'stitchedCFPseg';

YFPflatLoc = 'YFP_Flatfield.tif';
CFPflatLoc = 'CFP_Flatfield.tif';

if ~exist([Root,YFPsegSave],'dir')
    mkdir([Root,YFPsegSave])
end
if ~exist([Root,CFPsegSave],'dir')
    mkdir([Root,CFPsegSave])
end

halfYFPPaths = cell(length(imNos),1);
halfCFPPaths = cell(length(imNos),1);

for i = 1:length(imNos)
    halfYFPPaths{i} = [Root,sprintf(YFPStem,imNos(i))];
    halfCFPPaths{i} = [Root,sprintf(CFPStem,imNos(i))];
end

saveFluoPath = [Root,'FluoProfiles.mat'];

%Only need to deal with below if planning on exporting movie-quality stitchings
exportFluo = true; %Will export movie-quality stitchings if set to true, but requires additional parameters to be set below
correctionPath = [Root,'correctionCoordinates.mat']; %Sometimes the experiments can drift over time. If this is the case, can run tracking on static features (e.g. flecks of agar) and use the resulting coordinates to correct the drift.

YFPsave = 'stitchedYFP';
CFPsave = 'stitchedCFP';

if ~exist([Root,YFPsave],'dir')
    mkdir([Root,YFPsave])
end
if ~exist([Root,CFPsave],'dir')
    mkdir([Root,CFPsave])
end

reprocess = 0;
verbose = false;

%Actual data processing
if ~exist(saveFluoPath,'file') || reprocess > 0
    edgeYs = edgeYs(:,2:end); %Get rid of first time point (which is bound to be dodgy)
    Times = Times(2:end);
    frameNos = frameNos(2:end);
    noFrames = size(frameNos,2);
    
    maxYfromEdge = -50; %in um - the distance back from the leading edge to include in the window
    stationaryWindowEdge = median(edgeYs(:,1))-10; %In um - the stationary position of the front of the window for the stationary window (not the moving window)
    
    pixelStripLength = abs(round(maxYfromEdge/pxSize));
    CFPProfile = zeros(pixelStripLength,size(edgeYs,2));
    YFPProfile = zeros(pixelStripLength,size(edgeYs,2));
    
    CFPStationaryProfile = zeros(pixelStripLength,size(edgeYs,2));
    YFPStationaryProfile = zeros(pixelStripLength,size(edgeYs,2)); 
    
    YFPflat = double(imread([Root,YFPflatLoc]));
    CFPflat = double(imread([Root,CFPflatLoc]));
    
    if exist(correctionPath,'file')
        load(correctionPath)
    end
    
    %I was originally intending to use getColonyPixelFractions here, but we need to do binary operations on combinations of images so I've just copied the code here instead.
    for i = 1:size(edgeYs,2)
        if reprocess == 2 || ~exist([Root,YFPsave,filesep,sprintf(imStem,i)],'file')
            BFframeSeg = double(imread([stitchPath,sprintf(imStem,frameNos(i))])) > 0.5;
            
            YFPframe = stitchFrame(halfYFPPaths,imStem,YFPflat,frameNos(i),stitchSets,segSets);
            YFPframe = medfilt2(YFPframe,[2,2]);
            CFPframe = stitchFrame(halfCFPPaths,imStem,CFPflat,frameNos(i),stitchSets,segSets);
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
            
            %Write binary images
            imwrite(uint8(YFPareas*255),[Root,YFPsegSave,filesep,sprintf(imStem,i)])
            imwrite(uint8(CFPareas*255),[Root,CFPsegSave,filesep,sprintf(imStem,i)])
            
            if exportFluo
                %Prepare and write flatfield corrected fluoresence images (NOT
                %FOR QUANTIFICATION!)
                cropCoords = [1,1104,15468-7458-6144,15468-4144];
                
                YFPimg = YFPframe(cropCoords(3):cropCoords(4),cropCoords(1):cropCoords(2));
                CFPimg = CFPframe(cropCoords(3):cropCoords(4),cropCoords(1):cropCoords(2));
                
                YFPimg = (YFPimg-1)/max(YFPimg(:)-1);
                YFPimg(YFPimg<0) = 0;
                CFPimg = (CFPimg-1)/max(CFPimg(:)-1);
                CFPimg(CFPimg<0) = 0;
                
                if exist(correctionPath,'file')
                    %Subtract correction x-coordinate from image and pad with zeros
                    xCorr = corrCoords(i,1);
                    xCorrPx = round(xCorr/pxSize);
                    padding = zeros(abs(xCorrPx),size(YFPimg,2));
                    if xCorr > 0
                        YFPimg = [padding;YFPimg(1:end-xCorrPx,:)];
                        CFPimg = [padding;CFPimg(1:end-xCorrPx,:)];
                    elseif xCorr < 0
                        YFPimg = [YFPimg(-xCorrPx+1:end,:);padding];
                        CFPimg = [CFPimg(-xCorrPx+1:end,:);padding];
                    end
                end
                                               
                imwrite(imresize(YFPimg,0.5),[Root,YFPsave,filesep,sprintf(imStem,i)])
                imwrite(imresize(CFPimg,0.5),[Root,CFPsave,filesep,sprintf(imStem,i)])
            end
        else
            YFPareas = imread([Root,YFPsegSave,filesep,sprintf(imStem,i)]);
            CFPareas = imread([Root,CFPsegSave,filesep,sprintf(imStem,i)]);
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
    progressbar(1)
    save(saveFluoPath)
else
    load(saveFluoPath);
end

figure(1)
cla
imagesc(Times,[0,abs(maxYfromEdge)],medfilt2(YFPProfile./(YFPProfile + CFPProfile),[3,3]))
cBar = colorbar;
cBar.Label.String = 'YFP-labelled to total cell ratio';
cBar.LineWidth = 2;
xlabel('\textsf{Time (mins)}','FontSize',15)
ylabel('\textsf{Distance from Edge ($\mu m$)}','FontSize',15)
ax = gca;
ax.LineWidth = 2;

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
ax.Box = 'on';