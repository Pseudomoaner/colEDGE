clear all
close all

%File and path setup
Root = 'C:\Users\ph1ojm\Desktop\ColEDGE_Tester\20_10_17_PilH_YFP_WT_CFP_ColonyExpansion\';

saveTracePath = [Root,'ExtractedProfiles.mat'];
load(saveTracePath,'edgeYs','imNos','edgeXs','stitchSets','segSets','imStem','Times','frameNos','stitchPath','pxSize','dt')
stitchPath(1) = Root(1); %Just in case you changed directory mapping between calculating the expansions and this script...
segSets.segment = false;

reprocess = 3; %Set to 1 to reprocess based on existing YFPseg/CFPseg images. Set to 2 to reprocess those images with previously calculated temporal and spatial gradients. Set to 3 to also recalculate those gradients.
verbose = false;

%Only need to deal with below if planning on exporting movie-quality stitchings
exportFluo = true; %Will export movie-quality stitchings if set to true, but requires additional parameters to be set below
correctionPath = [Root,'correctionCoordinates.mat']; %Sometimes the experiments can drift over time. If this is the case, can run tracking on static features (e.g. flecks of agar) and use the resulting coordinates to correct the drift.

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

YFPsave = 'stitchedYFP';
CFPsave = 'stitchedCFP';

if ~exist([Root,YFPsave],'dir')
    mkdir([Root,YFPsave])
end
if ~exist([Root,CFPsave],'dir')
    mkdir([Root,CFPsave])
end

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
    
    %The CFP channel is a pain. Within a single dataset, there
    %can be two trends to deal with:
    %   -A spatial fluorescence gradient, where the end of the
    %   expansion path is somewhat brighter than the start
    %   -A temporal fluorescence gradient, where the final
    %   timepoint is somewhat brighter than the first.
    %Fortunately the YFP channel is less prone to these
    %problems, possibly because it has a higher signal to background.
    %But to be fair, we need to do the same corrections in that channel
    %as well.
    %
    %This section quantifies these two trends in both channels, so
    %they can be removed later. The operations used empirically give
    %the best results, i.e. the segmentations say there is one CFP cell
    %present where there should be one CFP cell, two where there
    %should be two etc.
    
    if reprocess == 3
        %First get the spatial gradient, based on the background of
        %the first frame in the dataset
        BFframeSeg = double(imread([stitchPath,sprintf(imStem,frameNos(1))])) > 0.5;
        SE = strel('disk',2);
        BFframeSeg = imdilate(BFframeSeg,SE); %Remove any fluorescence just surrounding cells
        
        CFPframe = stitchFrame(halfCFPPaths,imStem,CFPflat,frameNos(1),stitchSets,segSets);
        CFPframe = medfilt2(CFPframe,[2,2]);
        CFPframe(BFframeSeg) = NaN;
        CFPspatGrad = nanmean(CFPframe,2);
        
        YFPframe = stitchFrame(halfYFPPaths,imStem,YFPflat,frameNos(1),stitchSets,segSets);
        YFPframe = medfilt2(YFPframe,[2,2]);
        YFPframe(BFframeSeg) = NaN;
        YFPspatGrad = nanmean(YFPframe,2);
        
        %Then get the temporal gradient, based on the increase in
        %fluorescence of that portion of background that is NEVER
        %occupied by cells (i.e. is not segmented in the final
        %timepoint).
        BFframeSeg = double(imread([stitchPath,sprintf(imStem,frameNos(end))])) > 0.5;
        BFframeSeg = imdilate(BFframeSeg,SE); %Clear away any contribution from the region behind the leading edge
        CFPtempGrad = zeros(size(edgeYs,2),1);
        YFPtempGrad = zeros(size(edgeYs,2),1);
        for i = 1:size(edgeYs,2)
            CFPframe = stitchFrame(halfCFPPaths,imStem,CFPflat,frameNos(i),stitchSets,segSets);
            CFPframe(BFframeSeg) = NaN;
            CFPframe = CFPframe - repmat(CFPspatGrad,1,size(CFPframe,2)) + 1; %Correct for spatial gradient you just detected
            CFPtempGrad(i) = nanmean(CFPframe(:));
            
            YFPframe = stitchFrame(halfYFPPaths,imStem,YFPflat,frameNos(i),stitchSets,segSets);
            YFPframe(BFframeSeg) = NaN;
            YFPframe = YFPframe - repmat(YFPspatGrad,1,size(YFPframe,2)) + 1; %Correct for spatial gradient you just detected
            YFPtempGrad(i) = nanmean(YFPframe(:));
            
            progressbar(i/size(edgeYs,2))
        end
    end
    
    %Now run the main quantification steps
    for i = 1:size(edgeYs,2)
        if reprocess >= 2 || ~exist([Root,YFPsegSave,filesep,sprintf(imStem,i)],'file')
            BFframeSeg = double(imread([stitchPath,sprintf(imStem,frameNos(i))])) > 0.5;
            
            YFPframe = stitchFrame(halfYFPPaths,imStem,YFPflat,frameNos(i),stitchSets,segSets);
            YFPframe = medfilt2(YFPframe,[2,2]);
            CFPframe = stitchFrame(halfCFPPaths,imStem,CFPflat,frameNos(i),stitchSets,segSets);
            CFPframe = medfilt2(CFPframe,[2,2]);
            
            if reprocess == 2
                load(saveFluoPath,'CFPtempGrad','CFPspatGrad','YFPtempGrad','YFPspatGrad')
            end
            
            %Adjust the CFP frame to remove the trends detected in the
            %first stage
            CFPframe = CFPframe - repmat(CFPspatGrad,1,size(CFPframe,2)) + 1; %Spatial gradient correction
            CFPframe = CFPframe - CFPtempGrad(i) + 1; %Temporal gradient correction
            
            %Likewise for the YFP frame (shouldn't do much)
            YFPframe = YFPframe - repmat(YFPspatGrad,1,size(YFPframe,2)) + 1; %Spatial gradient correction
            YFPframe = YFPframe - YFPtempGrad(i) + 1; %Temporal gradient correction
            
            %Split the two populations according to an automatically calculated threshold
            
            %Threshold will be based on outermost 'skin' of cells,
            %which should be a single layer thick at all timepoints.
            %Isolate this skin from the segmentation:
            cellWidth = 5; %In pixels
            SE = strel('disk',7);
            thickSE = strel('disk',cellWidth);
            bulkColony = imclose(BFframeSeg,SE);
            outerSkin = bwmorph(bulkColony,'remove');
            thickSkin = imdilate(outerSkin,thickSE);
            realSkin = and(thickSkin,BFframeSeg);
            
            %Clear away the borders of the segmentated image
            realSkin(:,1:cellWidth*2) = false;
            realSkin(:,end-(cellWidth*2):end) = false;
            realSkin(1:cellWidth*2,:) = false;
            
            %Also remove anything that is far away from the leading
            %edge (debris etc.)
            currPos = round(median(edgeYs(:,i))/pxSize);
            realSkin(1:(currPos - 1000),:) = false;
            realSkin((currPos + 200):end,:) = false;
            
            FluoRatio = YFPframe./CFPframe;
            data = FluoRatio(realSkin(:));
            GMModels = cell(2,1);
            AIC = zeros(2,1);
            for k = 1:2
                GMModels{k} = fitgmdist(data,k);
                AIC(k)= GMModels{k}.AIC;
            end
            
            figure(1)
            histogram(data)
            pause(0.1)
            
            %Only run subsequent steps if the Gaussian mixture fitting
            %actually belives there are two populations present.
            %Otherwise, assume overtake by CFP is (near) complete and
            %only update its single-cell estimate
            model = GMModels{2};
            if AIC(2) < AIC(1)*1.05 && max(model.ComponentProportion) < 0.75 %Choice of parameters here is empirical, based on when the distribution 'looks' like the YFP cells have been lost.
                disp('opt1')
                
                idx = cluster(model,data);
                
                %Bleed through from YFP -> CFP is much greater than the
                %other way around. But we will treat both channels
                %equally in its correction.
                YFPdata = YFPframe(realSkin(:));
                YFPcluster1 = YFPdata(idx == 1);
                YFPcluster2 = YFPdata(idx == 2);
                YFPpeaks = [mean(YFPcluster1),mean(YFPcluster2)];
                
                YFP_bkd = mean(YFPframe(~bulkColony(:)));
                YFP_CFP = min(YFPpeaks) - YFP_bkd; %YFP fluorescence induced by presence of CFP-expressing cell
                YFP_YFP = max(YFPpeaks) - YFP_bkd; %YFP fluorescence induced by presence of YFP-expressing cell
                
                CFPdata = CFPframe(realSkin(:));
                CFPcluster1 = CFPdata(idx == 1);
                CFPcluster2 = CFPdata(idx == 2);
                CFPpeaks = [mean(CFPcluster1),mean(CFPcluster2)];
                
                CFP_bkd = mean(CFPframe(~bulkColony(:)));
                CFP_YFP = min(CFPpeaks) - CFP_bkd; %CFP fluorescence induced by presence of YFP-expressing cell
                CFP_CFP = max(CFPpeaks) - CFP_bkd; %CFP fluorescence induced by presence of CFP-expressing cell
            elseif mean(data) > 1.05 %Implies YFP is currently dominating front
                disp('opt2')
                
                YFP_bkd = mean(YFPframe(~bulkColony(:)));
                YFPdata = YFPframe(realSkin(:));
                YFP_YFP = mean(YFPdata) - YFP_bkd;
                
                CFP_bkd = mean(CFPframe(~bulkColony(:)));
            elseif mean(data) < 1 %Implies CFP is currently dominating front
                disp('opt3')
                
                CFP_bkd = mean(CFPframe(~bulkColony(:)));
                CFPdata = CFPframe(realSkin(:));
                CFP_CFP = mean(CFPdata) - CFP_bkd;
                
                YFP_bkd = mean(YFPframe(~bulkColony(:)));
            else %Ain't got no idea what's going on here. Best just to use the previous values and hope things improve
                YFP_bkd = mean(YFPframe(~bulkColony(:)));
                CFP_bkd = mean(CFPframe(~bulkColony(:)));
            end
            
            YFPareas = round((YFPframe - ((CFPframe-CFP_bkd)*(CFP_YFP/CFP_CFP)) - YFP_bkd)/YFP_YFP);
            CFPareas = round((CFPframe - ((YFPframe-YFP_bkd)*(YFP_CFP/YFP_YFP)) - CFP_bkd)/CFP_CFP);
            
            YFPareas(YFPareas<0) = 0;
            CFPareas(CFPareas<0) = 0;
            
            YFPareas(~BFframeSeg) = 0;
            CFPareas(~BFframeSeg) = 0;
            
            YFPareas = YFPareas .* bwareaopen(YFPareas,10);
            CFPareas = CFPareas .* bwareaopen(CFPareas,10);
            
            %Write binary images
            imwrite(uint8(YFPareas),[Root,YFPsegSave,filesep,sprintf(imStem,i)])
            imwrite(uint8(CFPareas),[Root,CFPsegSave,filesep,sprintf(imStem,i)])
            
            if exportFluo
                %Prepare and write flatfield corrected fluoresence images
                YFPimg = YFPframe;
                CFPimg = CFPframe;
                
                SE = strel('disk',15);
                
                YFPimg(~(imdilate(BFframeSeg,SE))) = 0;
                CFPimg(~(imdilate(BFframeSeg,SE))) = 0;
                
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
                
                imwrite(uint16(YFPimg*2^12),[Root,YFPsave,filesep,sprintf(imStem,i)])
                imwrite(uint16(CFPimg*2^12),[Root,CFPsave,filesep,sprintf(imStem,i)])
            end
        elseif reprocess == 1
            YFPareas = imread([Root,YFPsegSave,filesep,sprintf(imStem,i)]);
            CFPareas = imread([Root,CFPsegSave,filesep,sprintf(imStem,i)]);
        end
        %Extract strips, position 0 of the strip being the leading edge.
        YFPMean = getLeadingEdgePixelProfile(YFPareas,round(edgeXs/pxSize),round(edgeYs(:,i)/pxSize),pixelStripLength,false);
        CFPMean = getLeadingEdgePixelProfile(CFPareas,round(edgeXs/pxSize),round(edgeYs(:,i)/pxSize),pixelStripLength,false);
        
        if maxYfromEdge < 0
            YFPMean = flip(YFPMean);
            CFPMean = flip(CFPMean);
        end
        
        YFPProfile(:,i) = YFPMean;
        CFPProfile(:,i) = CFPMean;
        
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
xlabel('\textsf{Time (mins)}','FontSize',15)
ylabel('\textsf{YFP to CFP cell proportion}','FontSize',15)
legend('Front','Homeland','Location','NorthEast')
ax = gca;
ax.LineWidth = 2;
ax.Box = 'on';