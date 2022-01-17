%Copyright (C) 2013-2021  Martin Schmuck, Thomas Temme

%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU Affero General Public License as published
%by the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.

%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU Affero General Public License for more details.

%You should have received a copy of the GNU Affero General Public License
%along with this program.  If not, see <https://www.gnu.org/licenses/>.


function CompositeOligos_Callback(hObject, eventdata, handles)
% hObject    handle to CompositeOligos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
%foldername = uigetdir;
%Call EdgeFill Neurons and Skeleton Neurons and play with parameters in
%Option Handler
Oligo = 1;
ExecuteOptionExperimentOligo(handles,1,'Neurite');
guidata(handles.figure1, handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function ExecuteOptionExperimentOligo(handles,onlyCompFill,channelName)
h = waitbar(0.5,'Neuronal quantification in progress. This could take a few hours. Please be patient.', 'WindowStyle', 'modal');
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
optionHandler = handles.OptionHandler;

if(strcmp(channelName,'Oligo'))
    competingHandler = handles.NeuronCoordinates;
    if(isfield(handles,'OligoCoordinates'))
       neuronHandler = handles.OligoCoordinates;
    else
       neuronHandler = CSVCoordinates();
    end
else
    neuronHandler = handles.NeuronCoordinates;
    competingHandler = 0;
end
wellList = get(handles.lbWell, 'string');
foldername = imageHandler.Foldername;
[sizeY sizeX] = size(imageHandler.NucleusImage);
%Select STARTWELL

if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
if(numel(neuronHandler.AreaDictionary) ==0)
    neuronHandler.AreaDictionary = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
    neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFill) ==0)
    neuronHandler.NeuronPositionsEdgeFill = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFillSecond) ==0)
    neuronHandler.NeuronPositionsEdgeFillSecond = containers.Map();
end
if(numel(neuronHandler.NeuronStatMatrix) ==0)
    neuronHandler.NeuronStatMatrix = containers.Map();
end

%Create Matlab job
%Leave one Core out to ensure that GUI isn't freezing
maxCores =optionHandler.NumberOfCores;
if(maxCores > 1)
    c = parcluster();
    c.NumWorkers = maxCores;
    job = createJob(c,'JobData',optionHandler);
end

for(j=optionHandler.StartWell:numel(wellList))
    selectedWell=wellList(j);
    selectedWell=selectedWell{1};
    selectedWellLong=selectedWell;
    if(length(selectedWell) == 2)
              selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    end
    NucleusM = csvHandler.CellPosMatrix(selectedWell);    
    if(size(NucleusM,2) > sizeX || size(NucleusM,1) > sizeY)
        NucleusM = NucleusM(1:sizeY,1:sizeX);
        csvHandler.CellPosMatrix(selectedWell) = NucleusM;
        handles.CSVCoordinates = csvHandler;
        guidata(handles.figure1, handles);
    end
    if(maxCores > 1)
        %Check if watershed images are available. If not create them before
        %processing as ImageJ operation can't be done in parallel.
        %Check if there is already a watershed image. If not, create it.
        imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
        if(~exist(imagePathWatershed,'file'))
            WatershedNuclei(handles);
        end
        createTask(job, @ExecuteCompositeFillAndSkeletonOligo, 4, {selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler});
    else
       [neuronHandler, TP, FP, TN] = ExecuteCompositeFillAndSkeletonOligo(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler);
    end
end
if(maxCores > 1)
    submit(job);
    wait(job);
    %job = out(2);
    out=fetchOutputs(job);
    %ToDo: Merge NeuronHandler
    TP = sum(cell2mat(out(:,2)));
    FP = sum(cell2mat(out(:,3)));
    TN = sum(cell2mat(out(:,4)));
end
% TPPerWell = cell2mat(out(:,2));
% FPPerWell = cell2mat(out(:,3));
% TNPerWell = cell2mat(out(:,4));
% TPPerWell = TPPerWell ./ (TPPerWell+FPPerWell);
% FPPerWell = FPPerWell ./ (TPPerWell+FPPerWell);
i=0;

if(maxCores > 1)
    for(j=optionHandler.StartWell:numel(wellList))
     %Get current NeuronHandler result and write it in main Neuronhandler
        i=i+1;
        selectedWell=wellList(j);
        selectedWell=selectedWell{1};    
    
        currentNeuronHandler = out(i,1);
        currentNeuronHandler=currentNeuronHandler{1};
        %Merge Neuron positions as well as neurite length matrix.
        neuronHandler.NeuronPositionsSkeletonization(selectedWell) = currentNeuronHandler.NeuronPositionsSkeletonization(selectedWell);
        neuronHandler.AreaDictionary(selectedWell) = currentNeuronHandler.AreaDictionary(selectedWell);
        neuronHandler.NeuronPositionsEdgeFill(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFill(selectedWell);
        neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
        neuronHandler.NeuronStatMatrix(selectedWell) = currentNeuronHandler.NeuronStatMatrix(selectedWell);
        neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = currentNeuronHandler.NeuronPositionsSkelDeleted(selectedWell);
        %neuronHandler.NeuriteLengthMatrix(selectedWell) = currentNeuronHandler.NeuriteLengthMatrix;
    end
end
close(h);
if(strcmp(channelName,'Oligo'))
     handles.OligoCoordinates = neuronHandler;
else
    handles.NeuronCoordinates = neuronHandler;
end
guidata(handles.figure1, handles);


 function [neuronHandler, TP, FP, TN] = ExecuteCompositeFillAndSkeletonOligo(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler)    
        %Get first filter from selected Well
        firstY=1;
        firstX=1;
        lastY=sizeY;
        lastX=sizeX;      
         if(optionHandler.FilterMulPa==1)
             path = [foldername '/ConvertedCellomics/' selectedWell '.mat'];
             TP=0;
             FP=0;
             TN=0;
             if(exist(path,'file')) 
                 load(path);
                 %For every positive Filter
                 filterCount = numel(FMask.PositiveFilters);
                 for(i=1:filterCount);
                    f1 = FMask.PositiveFilters{i};
                    firstInd=find(f1,1, 'first');
                    lastInd=find(f1,1, 'last');
                    [firstY firstX] = ind2sub([sizeY sizeX], firstInd);
                    [lastY lastX] = ind2sub([sizeY sizeX], lastInd);
                    [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithmOligo(selectedWell,firstY,lastY, firstX, lastX,0,optionHandler,neuronHandler,NucleusM,foldername,channelName);                     
                end
             end
         else
             if(lastY - firstY > 10000 || lastX - firstX > 10000)        
                sizeXsmall = (lastX-firstX)/2;
                sizeYsmall = (lastY-firstY)/2;
                [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithmOligo(selectedWell,firstY,sizeYsmall,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithmOligo(selectedWell,firstY,sizeYsmall,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithmOligo(selectedWell,sizeYsmall+1,lastY,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithmOligo(selectedWell,sizeYsmall+1,lastY,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
              else
                [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithmOligo(selectedWell,firstY,lastY, firstX, lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);  
              end
              
         end    
         
        % system('taskkill /fi "WINDOWTITLE eq FBird"');
        %Calculate only for filter 
        
function [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithmOligo(selectedWell, yStartPos, yEndPos, xStartPos, xEndPos, save, optionHandler, neuronHandler, NucleusM, foldername, channelName)
[sizeYorig sizeXorig] = size(NucleusM);
NucleusM = NucleusM(yStartPos:yEndPos,xStartPos:xEndPos);

%[D IDX] = bwdist(full(NucleusM));
TN=0;
if(isfield(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
    NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
    NeuronManualM = NeuronManualM(yStartPos:yEndPos,xStartPos:xEndPos);
    TN = nnz(NeuronManualM);
else
    NeuronManualM = sparse(logical(zeros(sizeYorig,sizeXorig)));
    if(isa(neuronHandler.ManualNeuronPositionsSparse,'containers.Map') == 0)
        neuronHandler.ManualNeuronPositionsSparse = containers.Map();
    end
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
end
%Load neurite and nucleus image from file 
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end
imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
imagePathNeuriteBig = [foldername '/ConvertedCellomics/' selectedWellLong 'OligoBig' optionHandler.FileEnding];

%Check if NeuronStatMatrix exists
if(~isa(neuronHandler.NeuronStatMatrix,'containers.Map'))
    neuronHandler.NeuronStatMatrix = containers.Map();
end
if(~isa(neuronHandler.OmniNucleusPositions,'containers.Map'))
    neuronHandler.OmniNucleusPositions = containers.Map();
end
currentNeuronStat = containers.Map();
neuritePic = imread(imagePathNeuriteBig);
nucleusImage = imread(imagePathNucleusBig);
sizeY = yEndPos - yStartPos+1;
sizeX = xEndPos - xStartPos+1;
[sizeYOrig sizeXOrig]= size(neuritePic);

%Threshold Neurite Picture the known way
    
[sizeY sizeX] = size(neuritePic);
if(optionHandler.FilterMulPa~=1 && optionHandler.CutOutSphereCore == 1)
    neuritePic = CutOutCircles(neuritePic,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM, 0);
end
[neuritePic neuritePicStrong lMin] = ThresholdPic(neuritePic,optionHandler, neuronHandler.StretchlimResult, sizeY, sizeX,0,channelName);

SE = strel('disk', 3);
neuritePic = imdilate(neuritePic,SE);
neuritePic = neuritePic(yStartPos:yEndPos,xStartPos:xEndPos);
neuritePicStrong = neuritePicStrong(yStartPos:yEndPos,xStartPos:xEndPos);
nucleusPic=nucleusImage;
if(optionHandler.FilterMulPa~=1 && optionHandler.CutOutSphereCore == 1)
    nucleusPic = CutOutCircles(nucleusPic,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM,0, 0);
end
nucleusPic(nucleusPic<optionHandler.NucleusThreshold) = 0;
nucleusPic(nucleusPic>=optionHandler.NucleusThreshold) = 255;
nucleusPic=imcomplement(nucleusPic);
nucleusPic = uint8(nucleusPic);

%Do Watershed
%Check if there is already a watershed image. If not, create it.
imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
if(exist(imagePathWatershed,'file'))
    nucleusPic = imread(imagePathWatershed);
else
    javaaddpath('lib/ij.jar');
    javaaddpath('lib/mij.jar');
    MIJ.start(false);
    MIJ.createImage(nucleusPic);
    MIJ.run('Options...', 'iterations=1 count=1 edm=Overwrite do=Nothing');
    MIJ.run('Watershed');
    %nucleusPicBin = logical(MIJ.getCurrentImage);
    MIJ.run('Save', 'save=temp.tif');
    MIJ.closeAllWindows;
    MIJ.exit;
    nucleusPic = imread('temp.tif');
    delete temp.tif;
    nucleusPic = logical(1-nucleusPic);
    imwrite(nucleusPic,imagePathWatershed);
end
nucleusPic = nucleusPic(yStartPos:yEndPos,xStartPos:xEndPos);
[sizeYsmall sizeXsmall] = size(nucleusPic);
fusedPic = imfuse(nucleusPic, neuritePic);
%figure(1);
%imshow(fusedPic);
neuronPointIndices = fusedPic(:,:,1) > 250 & fusedPic(:,:,2) > 250 & fusedPic(:,:,3) > 250;
[sizeYneuronPoints sizeXneuronPoints] = size(neuronPointIndices);
%For every Nucleus:
%Check if Neurite with edges around.
%If yes: Count overlap of Nucleus and Neurites
%Count how many points of Nucleus are within the Neurite edges
EdgeFillNeuronPositionList = sparse(double(sizeY),double(sizeX));
EdgeFillNeuronPositionListSecond = sparse(double(sizeY),double(sizeX));
EdgeFillNeuronPositionList = EdgeFillNeuronPositionList(yStartPos:yEndPos,xStartPos:xEndPos);
EdgeFillNeuronPositionListSecond = EdgeFillNeuronPositionListSecond(yStartPos:yEndPos,xStartPos:xEndPos);
%stringCSVStatisticExport = sprintf('Position X;Position Y;Overlap (percent)\r\n');
[NonZeroY, NonZeroX] = find(NucleusM);
NonZeroY = uint16(NonZeroY);
NonZeroX = uint16(NonZeroX);

dbgYesCount = 0;
dbgNoCount = 0;
dbgNoStartFound = 0;
TP = 0;
FP = 0;
meanSize=zeros(0);
meanBrightness=zeros(0);
count=0;
nucleusPicBin = logical(nucleusPic);
% Get all connected components and label them
LabelComps=bwlabel(neuronPointIndices,8);
LabelCompsNeu = bwlabel(neuritePic);
%Check for each connected component, if also a place on NucleusPic is given
LabelCompsNuc = bwlabel(nucleusPicBin);
for i=1:numel(NonZeroY)    
    if(mod(i,10) == 0)
       disp([num2str(i) ' of ' num2str(numel(NonZeroY))]);
    end
    startX=0;
    startY=0;
    maxDistWhiteNuc = optionHandler.EdgeFillLookAround;
    maxDistWhiteNucIndex=-9;
    maxDistWhiteNucIndexSaved=-9;
    for j=-maxDistWhiteNuc:maxDistWhiteNuc
        for k=-maxDistWhiteNuc:maxDistWhiteNuc
            maxDistWhiteNucIndex = abs(abs(j)+abs(k));
            currentX = NonZeroX(i) + k;
            currentY = NonZeroY(i) + j;
            %If point or one neighbour is white in Overlaypic.            
            if(currentY < sizeY && currentX < sizeX && currentY>0 && currentX>0 && currentX<sizeXneuronPoints && currentY < sizeYneuronPoints && neuronPointIndices(currentY,currentX) > 0 && (startX ==0 || maxDistWhiteNucIndex < maxDistWhiteNucIndexSaved))
                startX= currentX;
                startY = currentY;                
            end
        end
    end
    if(startX > 0)
        %FloodFill white Area from current point. Check also nucleus points
        %stack = java.util.Stack();
        %Idea to be more efficient:
        %Cut out NucleusArea from startX - 50 until startX + 50
        %Do FloodFill Operation only on this subset
        if(startY-150>0 && startX-150>0 && startX+150 < sizeXsmall && startY + 150 < sizeYsmall)
            SubArea = neuronPointIndices(startY-150:startY+150,startX-150:startX+150);
            SubNucleusM = NucleusM(startY-150:startY+150,startX-150:startX+150);
            SubNucleusPic = nucleusPic(startY-150:startY+150,startX-150:startX+150);
            SubNucleusImage = nucleusImage(startY-150:startY+150,startX-150:startX+150);
            
            SubLabelComps = LabelComps(startY-150:startY+150,startX-150:startX+150);
            SubLabelCompsNuc = LabelCompsNuc(startY-150:startY+150,startX-150:startX+150);
            SubLabelCompsNeu = LabelCompsNeu(startY-150:startY+150,startX-150:startX+150);
            conCompIndex = SubLabelComps(150, 150);
            conCompIndexNuc = SubLabelCompsNuc(150, 150);
            conCompIndexNeu = SubLabelCompsNeu(150, 150);
            
            if(conCompIndex == 0)
               conCompIndex = SubLabelComps(151, 151); 
               if(conCompIndex == 0)
                   conCompIndex = SubLabelComps(149, 149); 
                   if(conCompIndex == 0)
                     conCompIndex = SubLabelComps(151, 149); 
                        if(conCompIndex == 0)
                            conCompIndex = SubLabelComps(149, 151); 
                        end
                   end
               end
            end
            if(conCompIndexNuc == 0)
               conCompIndexNuc = SubLabelCompsNuc(151, 151); 
               if(conCompIndexNuc == 0)
                   conCompIndexNuc = SubLabelCompsNuc(149, 149); 
                   if(conCompIndexNuc == 0)
                     conCompIndexNuc = SubLabelCompsNuc(151, 149); 
                        if(conCompIndexNuc == 0)
                            conCompIndexNuc = SubLabelCompsNuc(149, 151); 
                        end
                   end
               end
            end
            if(conCompIndexNeu == 0)
               conCompIndexNeu = SubLabelCompsNeu(151, 151); 
               if(conCompIndexNeu == 0)
                   conCompIndexNeu = SubLabelCompsNeu(149, 149); 
                   if(conCompIndexNeu == 0)
                     conCompIndexNeu = SubLabelCompsNeu(151, 149); 
                        if(conCompIndexNeu == 0)
                            conCompIndexNeu = SubLabelCompsNeu(149, 151); 
                        end
                   end
               end
            end
            
            
            %Extract area
            numberNeuritePixels = bwarea(SubArea == 1 & SubLabelCompsNuc == conCompIndexNuc);
            numberNucleusPixels = bwarea(SubLabelCompsNuc == conCompIndexNuc);            
            numberNeuritePixelsTotal = bwarea(conCompIndexNeu == SubLabelCompsNeu);
            
            count=count+1;
            meanSize = [meanSize numberNucleusPixels];
            meanBrightness = [meanBrightness CalculateNucleusBrightness(SubNucleusPic, SubNucleusImage)];
            [maxY, maxX] = size(SubNucleusPic);
            %Recalculate startX and startY            
         %   stack.push([50 50]);
            %FloodFill SubArea as well as Binary NucleusPic
            %Count number of pixels of both during FloodFill operation
      %      [numberNeuritePixels, numberNucleusPixels] = FloodFillForEdgeFill(SubNucleusPic,SubArea,maxX,maxY,stack);
            %WhiteArea2 = imfill(SubArea,[50 50]);         
            %[yIndices xIndices] = find(WhiteArea2);
            %yIndices = yIndices + (startY-50);
            %xIndices = xIndices + (startX-50);
        else
            %Extract connected component at given point
            conCompIndex = LabelComps(startY, startX);
            conCompIndexNuc = LabelCompsNuc(startY, startX);
            conCompIndexNeu = LabelCompsNeu(startY, startX);
            %Extract area
            numberNeuritePixels = bwarea(LabelCompsNuc == conCompIndexNuc & neuronPointIndices == 1);
            numberNucleusPixels = bwarea(LabelCompsNuc == conCompIndexNuc);
            numberNeuritePixelsTotal = bwarea(LabelCompsNeu == conCompIndexNeu);
%            stack.push([startX startY]);
            WhiteArea = zeros(sizeY,sizeX);
    %        [numberNeuritePixels, numberNucleusPixels] = FloodFillForEdgeFill(nucleusPic,neuronPointIndices,sizeX,sizeY,stack);
             %Schwerpunkt aller Punkte berechnen.
            %[yIndices xIndices] = find(WhiteArea);
        end     
        %Check if Area has not too less and not too much pixels

        %Check relation between numberNeuritePixels and numberNucleusPixels
        
        %Add neurite independent of it's neuritePixels per
        %NucleusPixels. Save that value.
        %Data structure: Matrix with positions AND relation
        
        %Mean calculation of brightness and NucleusArea
        
        area = optionHandler.EdgeFillNucleusAreaWithinNeurite;
        areaSecond = optionHandler.EdgeFillNucleusAreaWithinNeuriteSecond;
        areaAround = optionHandler.CompositeFillNeuriteAreaWithinNucleus;
        areaAroundMax = optionHandler.CompositeFillNeuriteAreaWithinNucleusMax;
        edgeFillMarked = 0;
        edgeFillMarkedSecond = 0;
        
        %Find next position of NucleusM from current NucleusOM
        %[curY curX] = ind2sub([sizeY,sizeX],IDX(NonZeroY(i),NonZeroX(i)));
        
        if(numberNeuritePixels/numberNucleusPixels >= area && numberNucleusPixels/numberNeuritePixelsTotal >= areaAround && numberNucleusPixels/numberNeuritePixelsTotal <= areaAroundMax && numberNeuritePixelsTotal > optionHandler.MinSizeNeuriteAreaCompFill) 
            EdgeFillNeuronPositionList(NonZeroY(i), NonZeroX(i)) = 1;
            EdgeFillNeuronPositionListSecond(NonZeroY(i), NonZeroX(i)) = 1;
            edgeFillMarked=1;
            edgeFillMarkedSecond=1;
        elseif(numberNeuritePixels/numberNucleusPixels >= areaSecond) 
            EdgeFillNeuronPositionListSecond(NonZeroY(i), NonZeroX(i)) = 1;   
            edgeFillMarkedSecond=1;
        end
        %stringCSVStatisticExport = [stringCSVStatisticExport num2str(NonZeroX(i)) ';' num2str(NonZeroY(i)) ';' num2str(numberNeuritePixels/numberNucleusPixels) ';' sprintf('\r\n')];
        %EdgeFillNeuronPositionList(NonZeroY(i), NonZeroX(i),2) = (numberNeuritePixels / numberNucleusPixels)*100;
        %Check if position is marked as Manual
        
        if(NeuronManualM(uint16(yStartPos)-1+uint16(NonZeroY(i)),uint16(xStartPos)-1+uint16(NonZeroX(i))))
            manualMarked=1;            
            if(edgeFillMarked)
                TP = TP + 1;
                TN = TN - 1;
            end
        else
            manualMarked=0;
            if(edgeFillMarked)
                FP = FP + 1;
            end
        end
        if(save)
            %Check if NeuronStatMatrix contains already key for selected Position
            key = mapPositionToKeyString(NonZeroY(i),NonZeroX(i));
            if(isKey(currentNeuronStat,key))
               list = currentNeuronStat(key);
            else
                list = zeros(8,1);            
            end
            list(4)=manualMarked;
            list(2) = edgeFillMarked;
            list(3) = edgeFillMarkedSecond;
            list(6) = numberNeuritePixels/numberNucleusPixels;
            currentNeuronStat(key) = list;
            eFillFirst = 0;
            eFillSecond=0;
        else
            eFillFirst = EdgeFillNeuronPositionList;
            eFillSecond=EdgeFillNeuronPositionListSecond;
        end
    else
        dbgNoStartFound = dbgNoStartFound +1;
    end    
end

meanSize = median(meanSize);
meanBrightness = median(meanBrightness);

if(numel(neuronHandler.NucleusBrightnessDict) ==0)% || ~isKey(neuronHandler.NucleusBrightnessDict, selectedWell)
    neuronHandler.NucleusBrightnessDict = containers.Map();
    neuronHandler.NucleusAreaDict = containers.Map();
    neuronHandler.NucleusBrightnessDict(selectedWell) = meanBrightness;
    neuronHandler.NucleusAreaDict(selectedWell) = meanSize;
else
    neuronHandler.NucleusAreaDict(selectedWell) = meanBrightness;
    neuronHandler.NucleusAreaDict(selectedWell) = meanSize;
end
if(numel(neuronHandler.NeuronPositionsEdgeFill) ==0)
    neuronHandler.NeuronPositionsEdgeFill = containers.Map();
end
if(~isKey(neuronHandler.NeuronPositionsEdgeFill, selectedWell)) 
    neuronHandler.NeuronPositionsEdgeFill(selectedWell) = sparse(logical(zeros(sizeYOrig,sizeXOrig)));
    fullMatrix = neuronHandler.NeuronPositionsEdgeFill(selectedWell);
    fullMatrix(yStartPos:yEndPos,xStartPos:xEndPos) = EdgeFillNeuronPositionList;
    neuronHandler.NeuronPositionsEdgeFill(selectedWell) = fullMatrix;
    %neuronHandler.NeuronPositionsEdgeFill(selectedWell) = EdgeFillNeuronPositionList;
else
    fullMatrix = neuronHandler.NeuronPositionsEdgeFill(selectedWell);
    fullMatrix(yStartPos:yEndPos,xStartPos:xEndPos) = EdgeFillNeuronPositionList;
    neuronHandler.NeuronPositionsEdgeFill(selectedWell) = fullMatrix;
end
if(numel(neuronHandler.NeuronPositionsEdgeFillSecond) ==0)
    neuronHandler.NeuronPositionsEdgeFillSecond = containers.Map();
end
if(~isKey(neuronHandler.NeuronPositionsEdgeFillSecond, selectedWell))
    neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = sparse(logical(zeros(sizeYOrig,sizeXOrig)));
    fullMatrix = neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
    fullMatrix(yStartPos:yEndPos,xStartPos:xEndPos) = EdgeFillNeuronPositionListSecond;
    neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = fullMatrix;
    %neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = EdgeFillNeuronPositionListSecond;
else
    fullMatrix = neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
    fullMatrix(yStartPos:yEndPos,xStartPos:xEndPos) = EdgeFillNeuronPositionListSecond;
    neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = fullMatrix;
end
neuronHandler.NeuronStatMatrix(selectedWell) = currentNeuronStat;
%handles.NeuronCoordinates = neuronHandler;
%handles.ImageHandler = imageHandler;
%guidata(handles.figure1, handles);    
%Check if Neuronhandler already has Dictionary of NeuronPositionsEdgeComposite        
        