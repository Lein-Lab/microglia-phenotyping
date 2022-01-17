function OmnispheroOnline(filnuc,filneu)
%Input parameters: 1. Filename nucleus as filnuc, 2. Filename neurons as filneu
%oh = OnlineHelper();
imageHandler = ImageHandler();
imageHandler.NucleusImage = imread(filnuc);
imageHandler.NeuriteImage = imread(filneu);
optionHandler = Option();
optionHandler.HistogramMaxNucleus = 4095;
optionHandler.HistogramMaxNeurite = 4095;
optionHandler.HistogramMinNucleus = 0;
optionHandler.HistogramMinNeurite = 0;
optionHandler.StartWell = 1;
optionHandler.MigDistLowerNucleusThreshold = 20;
optionHandler.MigDistLowerDensityThreshold = 250;
optionHandler.MigDistLowerFloodFillThreshold = 105;
optionHandler.DensityDistributionRingNumber = 10;
optionHandler.MigrationDistanceDensityImageXSize = 57;
optionHandler.MigrationDistanceDensityImageYSize = 57;
optionHandler.EdgeCompositNeuriteLowerThreshold = 20;
optionHandler.EdgeCompositMinArea = 70;
optionHandler.EdgeCompositeDistanceNucleusWhiteArea = 3;
optionHandler.EdgeFillNucleusAreaWithinNeurite = 0.55;
optionHandler.EdgeFillNucleusAreaWithinNeuriteSecond = 0.35;
optionHandler.FileEnding = '.png';
optionHandler.NucleusThreshold = 15;
optionHandler.SkeletonMinNeuriteLength = 40;
optionHandler.SkeletonNeuriteThresholdHardDistance = 0.4;
optionHandler.SkeletonNeuriteThresholdLow = 0.75;
optionHandler.MaxDistanceFromEndpoint = 33;
optionHandler.ToleranceAngleFromEndpoint = 40;
optionHandler.MaxAllowedFP = 15;
%medHop = zeros(numel(C));
%medHopArea = zeros(numel(C));
[sizeY sizeX]=size(imageHandler.NucleusImage);
[counts x] = imhist(imageHandler.NucleusImage);   
firstDer = diff(counts);
[maxval hop] = max(firstDer(10:60));  
hop = hop + 10 - 3;
medHop = hop;

[maxval hop] = max(firstDer(100:240));  
hop = hop + 100 - 15;
medHopArea = hop;
meanHop = medHop;
meanHopArea=medHopArea;
optionHandler.MigDistLowerNucleusThreshold = meanHop(1);
optionHandler.MigDistLowerFloodFillThreshold = meanHopArea(1);
nucImage = thresholdPic(imageHandler.NucleusImage);
%Watershed with Nucleus Image
nucImage = watershedNucleusImage(nucImage,optionHandler);
NucleusM = sparse(logical(zeros(sizeY,sizeX)));

%1. Create CSVHandler Object and save all Nucleus positions
%CC = bwlabel(nucImage);
CC = bwconncomp(nucImage);
for cou=1:CC.NumObjects
    %currentNeuriteArea = sparse(logical(zeros(sizeY, sizeX)));
    %currentNeuriteArea(CC.PixelIdxList{cou}) = 1;
    %Get balance point of currentNeuriteArea. 
    [Y,X] = ind2sub([sizeY sizeX],CC.PixelIdxList{cou});
    %Mean of Y and X is balance point -> Point of Nucleus
    meanY = uint16(mean(Y));
    meanX = uint16(mean(X));
    NucleusM(meanY,meanX) = 1;
end
figure(1);
imshow(nucImage);
[nucleusRows, nucleusCols] = find(NucleusM);
hold on;
plot((nucleusCols),(nucleusRows),'Linestyle','none','Marker','.','Markersize',2,'Color',[1 .5 0]);       
hold off;

%Calculate Migration Area
[filterDistance nonFilterDistance SphereArea markerPointCoordinates] = calculateMigrationDistance(32, optionHandler, NucleusM, sizeY, sizeX, imageHandler.NucleusImage, imageHandler.NeuriteImage)

%Now we can quantify Neurons With our two algorithms
%CompositeFill
%Skeletonization

function nucleusBinary = watershedNucleusImage(nucleusBinary,optionHandler)
    nucleusPicBinBigNuclei = xor(bwareaopen(nucleusBinary,250),  bwareaopen(nucleusBinary,9000));
    D = bwdist(~nucleusPicBinBigNuclei);
    D = -D;
    D(~nucleusBinary) = -Inf;
    L = watershed(D);
    nucleusBinary(find(~logical(L))) = 0;
    nucleusBinary = xor(bwareaopen(nucleusBinary,35),  bwareaopen(nucleusBinary,15000));          


function BW = thresholdPic(Pic)
stretchlimLow = stretchlim(Pic,[0.975 0.9999]);
brightNeurites = imadjust(Pic,stretchlimLow,[0 1]);
brightNeurites = medfilt2(brightNeurites);
unsharpFilter = fspecial('unsharp');
brightNeurites = imfilter(brightNeurites,unsharpFilter);
brightNeurites = imcomplement(brightNeurites);
level=isodata(brightNeurites);
BW = im2bw(brightNeurites, level);
BW = imcomplement(BW);



function [filterDistance nonFilterDistance SphereArea markerPointCoordinates result] = calculateDensityDistribution(SphereAreaSize, NucleusImage, NeuriteImage, NucleusM, NeuronM)
markerPointCoordinates=-1;
filterDistance = -1;
nonFilterDistance=-1;
neuronHandler = handles.NeuronCoordinates;
[sizeY sizeX] = size(imageHandler.NucleusImage);
SphereArea=-1;
optionHandler = handles.OptionHandler;
path = [imageHandler.Foldername '/ConvertedCellomics/MarkerPointCoordinates-' selectedWell '.mat'];
ringNumber = optionHandler.DensityDistributionRingNumber;
if(exist(path,'file'))        
    load(path);
    filterDistance = str2double(filterDistance);
    nonFilterDistance = str2double(nonFilterDistance);
    %hInner = impoly(handles.axes2,double(markerPointCoordinates('0')./10));
    %previousMask = createMask(hInner);
    %for i=1:ringNumber  
    %  hCurrent = impoly(handles.axes2,double(markerPointCoordinates(num2str(i*10))./10));
    %end
  else
     NucleusM = csvHandler.CellPosMatrix(selectedWell);
     NeuronM = neuronHandler.CellPosMatrix(selectedWell);
     foldername = imageHandler.Foldername;
     subfoldername = [foldername '/ConvertedCellomics'];
     [sizeY sizeX] = size(imageHandler.NeuriteImage);
     [filterDistance nonFilterDistance SphereArea markerPointCoordinates] = calculateMigrationDistance(SphereAreaSize, optionHandler, NucleusM, NeuronM, subfoldername, sizeY, sizeX); 
     %Save hInner to file
     
     %if(markerPointCoordinates~=0)
        save('-v7.3',strcat(subfoldername,'/MarkerPointCoordinates-',selectedWell),'markerPointCoordinates', 'filterDistance', 'nonFilterDistance');
     %end
  end
%[filterDistance nonFilterDistance SphereArea markerPointCoordinates] = calculateMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, handles);
%Migration Distance calculation gives distance from Sphere and Nucleus
%We need: LinePic
if(filterDistance == -1)
   result = zeros(ringNumber,9);
elseif (markerPointCoordinates==0)
   result = zeros(ringNumber,9);
elseif(numel(markerPointCoordinates('0')) == 0)
    result = zeros(ringNumber,9);
else

%NucleusM = csvHandler.CellPosMatrix;
%NeuronM = neuronHandler.CellPosMatrix;
NucleusM = csvHandler.CellPosMatrix(selectedWell);
NeuronM = neuronHandler.CellPosMatrix(selectedWell);
if(~isKey(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = sparse(sizeY, sizeX);
end
NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);

if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
if(~isKey(neuronHandler.NeuronPositionsSkeletonization,selectedWell))
    neuronHandler.NeuronPositionsSkeletonization(selectedWell) = sparse(sizeY, sizeX);
end
NeuronSkelM = neuronHandler.NeuronPositionsSkeletonization(selectedWell);

if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
    neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
end
if(~isKey(neuronHandler.NeuronPositionsSkelDeleted,selectedWell))
    neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = sparse(sizeY, sizeX);
end
NeuronsDeleted = neuronHandler.NeuronPositionsSkelDeleted(selectedWell);

    NeuronSkelM = NeuronSkelM - (NeuronsDeleted);
    ind = NeuronSkelM<0;
    NeuronSkelM(ind) = 0;


ringNumber = optionHandler.DensityDistributionRingNumber;
[sizeY, sizeX] = size(imageHandler.NucleusImage);
mat = sparse(sizeY, sizeX);
%Fix size of CSVCoordinate Matrices
[NucleusMSizeY NucleusMSizeX] = size(NucleusM);
[NeuronMSizeY NeuronMSizeX] = size(NeuronM);
[NeuronManualMSizeY NeuronManualMSizeX] = size(NeuronManualM);
[filterSizeY filterSizeX] = size(imageHandler.NucleusImage);
if(NucleusMSizeX > filterSizeX)
    NucleusM(:,filterSizeX+1:NucleusMSizeX) = [];
    csvHandler.CellPosMatrix(selectedWell) = NucleusM;
end
if(NeuronMSizeX > filterSizeX)
   NeuronM(:,filterSizeX+1:NeuronMSizeX) = [];
   neuronHandler.CellPosMatrix(selectedWell) = NeuronM;
end
if(NeuronManualMSizeX > filterSizeX)
                NeuronManualM(:,filterSizeX+1:NeuronManualMSizeX) = [];
                neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
end
if(NucleusMSizeY > filterSizeY)
                NucleusM(filterSizeY+1:NucleusMSizeY,:) = [];
                csvHandler.CellPosMatrix(selectedWell) = NucleusM;
end
if(NeuronMSizeY > filterSizeY)
                NeuronM(filterSizeY+1:NeuronMSizeY,:) = [];
                neuronHandler.CellPosMatrix(selectedWell) = NeuronM;
end
if(NeuronManualMSizeY > filterSizeY)
                NeuronManualM(filterSizeY+1:NeuronManualMSizeY,:) = [];
                neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
end
result = zeros(ringNumber,9);
hInner = impoly(handles.axes2,double(markerPointCoordinates('0')./10));
previousMask = createMask(hInner);
cbRings = get(handles.cbPolys ,'Value');
if(~cbRings)
    delete(hInner);
end
for i=1:ringNumber
%Create masks from LinePic and cut them from CSVs
%hTest = impoly(handles.axes2, [1 1; 500 500]);   
    hCurrent = impoly(handles.axes2,double(markerPointCoordinates(num2str(i*10))./10));    
    currentMask = logical(createMask(hCurrent));
    currentRing = logical(currentMask-previousMask); 
    currentRing = imresize(currentRing, [sizeY, sizeX]);
    previousMask = currentMask;
    %Save hCurrent to file
    if(~cbRings)
        delete(hCurrent);
    end
%Filter each ring on CSVs to calculate Density            
    Nucleus2 = logical(logical(NucleusM) .* currentRing);
    NeuronCell2 = logical(logical(NeuronM) .* currentRing);
    NeuronManual2 = logical(logical(NeuronManualM) .* currentRing);
    NeuronSkel2 = logical(logical(NeuronSkelM) .* currentRing);
    result(i,1) = nnz(currentRing);
    result(i,2) = nnz(Nucleus2);
    result(i,3) = nnz(NeuronCell2);
    result(i,4) = nnz(NeuronManual2);
    result(i,5) = nnz(NeuronSkel2);
    result(i,6) = result(i,3) / result(i,2);
    result(i,7) = result(i,4) / result(i,2);
    result(i,8) = result(i,5) / result(i,2);
    result(i,9) = result(i,2) / result(i,1);
end
end

function Imx = keepMaxObj(X)
%Function to keep only the maximum sized (biggest) object in an image
%SCd 11/30/2010
%
%Updates:
%   -02/03/2011: Added ability to handle an image directly
%
%Usage:
%   Imx = keepMaxObj(CC);
%   Imx = keepMaxObj(V);
%
%Input Arguments:
%   -CC: Connected components returned from bwconncomp
%   -V: Logical image with parts you want true
%   
%Output Arguments:
%   -Imx: Logical volume with only the biggest object left true.
%
%See Also: bwconncomp
%
    %Error checking:
    assert(islogical(X)||isstruct(X),'The first input argument is expected to be a struct or a logical');
    if isstruct(X)
        CC = X;
        parts = {'PixelIdxList','ImageSize'};
        assert(all(ismember(parts,fieldnames(CC))),'CC is expected to be the output from bwconncomp');
    else
        CC = bwconncomp(X);
    end  
    clear X;
    %Preallocate and find number of voxels/object
    Nvox = zeros(CC.NumObjects,1);
    for ii = 1:CC.NumObjects
        Nvox(ii) = numel(CC.PixelIdxList{ii});
    end
    %Find the biggest object's index, warn and save all if there are multiples
    [mx,midx] = max(Nvox);
    more_than1_max = sum(mx==Nvox);
    if more_than1_max > 1
        midx = find(mx == Nvox);
        warning('Multiple:Maxima', 'There were %i objects with the maximum size.\n  They are all left on!',more_than1_max);
    end    
    %Create the final image
    Imx = false(CC.ImageSize);
    Imx([CC.PixelIdxList{midx}]) = true;


function [filterDistance nonFilterDistance SphereArea markerPointCoordinates] = calculateMigrationDistance(SphereAreaSize, optionHandler, NucleusM, sizeY, sizeX, NucleusImage, NeuriteImage)
ringNumber = optionHandler.DensityDistributionRingNumber;
SphereArea=0;
NucleusArea=0;
filterDistance=-1;
nonFilterDistance=-1;
markerPointCoordinates = 0;
NucleusArea64=0;
[DensityM] = getDensityDistributionsFromCSV(NucleusM,32,32);

%Find circle in Density Edge
%Get point of Maximum Density
[C,rowmaxarray]=max(DensityM);
[maxvalue,xMaxSphere]=max(C);
yMaxSphere=rowmaxarray(xMaxSphere);
i=0;    
 thresholdNucImage = optionHandler.MigDistLowerNucleusThreshold;
 %New Calculation method for Nucleus Area
 [thresholdedRows thresholdedCols] = find(NucleusImage > thresholdNucImage);
 thresholdedRows = [thresholdedRows;sizeY];
 thresholdedCols = [thresholdedCols;sizeX];
 thresholdedRows = [thresholdedRows;0];
 thresholdedCols = [thresholdedCols;0];
 %1. Cut out border
 NucleusArea = NucleusImage(1025:sizeY-1024, 1025:sizeX-1024);
 NucleusArea = logical(im2bw(NucleusArea, thresholdNucImage/255));
 NucleusArea = logical(keepMaxObj(NucleusArea));
 vertical = logical(zeros(sizeY-(2*1024),1024));
 horizontal = logical(zeros((1024),sizeX));
 NucleusArea = logical([vertical NucleusArea vertical]);
 NucleusArea = logical([horizontal; NucleusArea; horizontal]);

CuttedIndicesNuc=0;
NucleusAreaTemp=0;

 [row,col] = find(NucleusArea);
 if(length(row) > 0 && length(col) > 0)

   sYMapped = sum(row)/length(row);
   sXMapped = sum(col)/length(col);
   %Map sY and sX to Big Picture
   [sX sY] = MapPoint(sXMapped,sYMapped,SphereAreaSize,SphereAreaSize,sizeX,sizeY,1);

   %Take out Nuclei from NucleusM, for Nuclei which are too far from sX
   %and sY
   %Dictionary Point -> Distance
   [yIndices,xIndices,values] = find(NucleusM);       
   %for every point in NucleusM
     %calculate pdist from sX, sY
   distanceList = zeros(numel(yIndices),1);
   for (i=1:numel(yIndices))
     distanceList(i) = pdist([yIndices(i) xIndices(i);sYMapped sXMapped]);
   end
   %Sort by Distance
   [sortedDistances, sortedDistancesIndices] = sort(distanceList);
   %for Top 90% in Sorted Dictionary
     %Copy point to NucleusMSorted
   NucleusMSorted = sparse(sizeY,sizeX);
   for (i=1:numel(yIndices) * 0.97)
     NucleusMSorted(yIndices(sortedDistancesIndices(i)), xIndices(sortedDistancesIndices(i))) = 1;
   end

   %Calculate Density with NucleusMSorted
   %Calculate SphereArea with that
   [DensityMSized] = getDensityDistributionsFromCSV(NucleusMSorted,SphereAreaSize, SphereAreaSize);
% NucleusArea64=imfill(NucleusArea64,'holes');
%Flood fill everything from max Density above Threshold
%Outer circle is area for Sphere
   stack = java.util.Stack();
%Validate point of max density: It's only valid, if Cellomics Neurons
%are also in region
  i=0;
  while( i<=5)
    DensityM(yMaxSphere,xMaxSphere) =0;
    [C,rowmaxarray]=max(DensityMSized);
    [maxvalue,xMaxSphere]=max(C);
    yMaxSphere=rowmaxarray(xMaxSphere);
    i=i+1;
    stack.push([xMaxSphere yMaxSphere]);
  end
  SphereArea = logical(zeros(SphereAreaSize,SphereAreaSize));
  SphereArea = FloodFillIterative(DensityMSized,SphereArea,1,1,SphereAreaSize,SphereAreaSize,stack,0);
  SphereArea=imfill(SphereArea,'holes');
  SE = [1 1 0;1 1 0; 0 0 0];
  % SE = strel('disk', 2);
  SphereArea = imdilate(SphereArea,SE);
  SphereArea = edge(SphereArea);
  %MigrationPic = SphereArea + NucleusArea64;
  %figure(1);
  %imshow(MigrationPic);
%hPixelInfo = impixelinfo;        
%imageHandler.MousePosition = hPixelInfo;
%F�r weitere Verarbeitung:
% Berechne den Schwerpunkt aller wei�en Pixel in NucleusArea

  minimizingFactorX = 32/SphereAreaSize;
  minimizingFactorY = 32/SphereAreaSize;
  mediumDistanceList = zeros(0,0);
  mediumDistanceListInInterval = zeros(0);

% Gehe jeweils im 10 Grad Winkel von dort aus los.
   successCount=0;
   %Declaration for density distribution calculation
   %Mapping: Distance in 5 steps
   markerPointCoordinates = containers.Map();
   for i=0:ringNumber
    markerPointCoordinates(num2str(i*10)) = zeros(0,0); 
   end
   counter=0;
   for i=pi/32:pi/32:2*pi

    deltaX = cos(i);
    deltaY = sin(i);
    xNucleus = 0;
    yNucleus = 0;
    xOuterSphere = 0;
    yOuterSphere = 0;    
    x=sX;
    y=sY;
    xBig=sXMapped;
    yBig=sYMapped;
    if(uint8(x)<1 | uint8(y)<1)
        continue;
    end
    SphereArea(uint8(y),uint8(x)) =0;
    while(xOuterSphere ==0 && uint8(x) <= SphereAreaSize && uint8(y) <= SphereAreaSize && uint8(x) > 0 && uint8(y) > 0 && uint8(xBig) > 0 && uint8(yBig) >0 && uint16(xBig) < sizeX && uint16(yBig) < sizeY)
      [x y] = MapPoint(xBig,yBig,SphereAreaSize,SphereAreaSize,sizeX,sizeY,1);
      if(uint8(x)<1 | uint8(y)<1)
          break;
      end
      %[xBig yBig] = MapPoint(x,y,sizeX,sizeY,SphereAreaSizeX,SphereAreaSizeY,0);
      if(NucleusArea(uint16(yBig),uint16(xBig)) == 1)
        xNucleus=xBig;
        yNucleus=yBig;
          %[xNucleus yNucleus]=[xBig yBig];%MapPoint(xBig,yBig,sizeX,sizeY,256,256,0);
      end

      if(SphereArea(uint8(y),uint8(x)) == 1)
        xOuterSphere=xBig;
        yOuterSphere=yBig;
          %[xOuterSphere yOuterSphere] = %MapPoint(x,y,sizeX,sizeY,SphereAreaSizeX,SphereAreaSizeY,1);
      end
      %x=x+deltaX;
      %y=y+deltaY;
      xBig=xBig+deltaX;
      yBig=yBig+deltaY;
    end
    if(xNucleus == 0)
        xNucleus=sXMapped;
        yNucleus=sYMapped;
        %[xNucleus yNucleus]= MapPoint(sX,sY,sizeX,sizeY,SphereAreaSizeX,SphereAreaSizeY,1);
    end

    if(xOuterSphere~=0 && yOuterSphere~=0 && xNucleus ~= 0 && yNucleus ~= 0)
      successCount = successCount+1;
      mediumDistanceList = [mediumDistanceList; pdist([xOuterSphere yOuterSphere;xNucleus yNucleus]) successCount];        
      %Create 5 points between xyOuterSphere and xyNucleus. Save these
      %points to markerPointCoordinates
      stepX = int16((int16(xOuterSphere) - int16(xNucleus)) / (ringNumber));
      stepY = int16((int16(yOuterSphere) - int16(yNucleus)) / (ringNumber));
      for k=0:ringNumber
        markerPointCoordinates(num2str(k*10)) = [markerPointCoordinates(num2str(k*10)); [int16(xNucleus)+k*stepX int16(yNucleus)+k*stepY]];
      end
    elseif(xNucleus ~=0 && yNucleus ~=0)
        %Set all points to Nucleus, if Outer Sphere is out of borders
        xOuterSphere = xNucleus;
        yOuterSphere = yNucleus;
        successCount = successCount+1;
        mediumDistanceList = [mediumDistanceList; pdist([xOuterSphere yOuterSphere;xNucleus yNucleus]) successCount];        
      %Create 5 points between xyOuterSphere and xyNucleus. Save these
      %points to markerPointCoordinates
      stepX = int16((int16(xOuterSphere) - int16(xNucleus)) / (ringNumber));
      stepY = int16((int16(yOuterSphere) - int16(yNucleus)) / (ringNumber));
      for k=0:ringNumber
        markerPointCoordinates(num2str(k*10)) = [markerPointCoordinates(num2str(k*10)); [int16(xNucleus)+k*stepX int16(yNucleus)+k*stepY]];
      end
    end
   end
 if(numel(mediumDistanceList)>0)
    distancestddev = std(mediumDistanceList(:,1));
    mediumDistance = mean(mediumDistanceList(:,1));
    maxDistance = max(mediumDistanceList(:,1));

 for i=1:numel(mediumDistanceList(:,1))
  currentDistance = mediumDistanceList(i,1);
  if(i-1>0)
    leftNeighbourDistance = mediumDistanceList(i-1,1);
  else
      leftNeighbourDistance = mediumDistanceList(numel(mediumDistanceList(:,1)),1);
  end
  if(i+1<numel(mediumDistanceList(:,1)))
    rightNeighbourDistance = mediumDistanceList(i+1,1);
  else
      rightNeighbourDistance = 1;
  end
  threshold = (maxDistance-(0.85*maxDistance));
  %threshold = (maxDistance-(2*distancestddev));
%Check if current distance is in Interval. Check also if one of both
%neighbourpoints is also in interval
  %if(currentDistance >= 150 && currentDistance >= (maxDistance-(2.0*distancestddev)))
  if(currentDistance >= threshold && (leftNeighbourDistance > threshold || rightNeighbourDistance > threshold))
    mediumDistanceListInInterval = [mediumDistanceListInInterval currentDistance];
  else
      %If not: Set points of MarkerPointCoordinates to inner point.
      currentCounter = mediumDistanceList(i,2);
      markerReferencePoint = markerPointCoordinates('0');
      markerReferencePoint = markerReferencePoint(currentCounter,:);
      for i=1:ringNumber
        markerPointList = markerPointCoordinates(num2str(i*10));
        markerPointList(currentCounter,:) = markerReferencePoint;
        markerPointCoordinates(num2str(i*10)) = markerPointList;
      end          
  end
  %Extra check: Spike prevention: At least two markerpoints on one ring
  %not set back.      
 end
 mediumDistance = num2str(mediumDistance,'%f');
 mediumDistanceSelected = num2str(mean(mediumDistanceListInInterval),'%f');
 disp('Medium distance total (in Pixel): ');
 disp(mediumDistance);
 disp('Medium distance selected (in Pixel): ');
 disp(mediumDistanceSelected);
 filterDistance = mediumDistanceSelected;
 nonFilterDistance = mediumDistance;

 else
    distancestddev = 0;
    mediumDistance = 0;
    mediumDistanceListInInterval=0;
    maxDistance = 0;
    filterDistance=-1;
    nonFilterDistance=-1;
 end


%mediumDistance = mediumDistance / successCount;
 end

      
function result = FloodFillIterative(DensityM,AreaResult,threshold,above,maxX,maxY,stack,image)
[sizeY, sizeX] = size(image);
while(stack.isEmpty() == 0)
    vector = stack.pop();
    x=vector(1);
    y=vector(2);
    if(above==0)
        [checkX, checkY]=MapPoint(x,y,sizeX,sizeY,64,64,1);
        if(DensityM(y,x) <= threshold && AreaResult(y,x) == 0 && x-1 > 0 && y-1 > 0 && y+1<= maxY && x+1 <= maxX && image(checkY,checkX) > 250)
            AreaResult(y,x) = 1;
            stack.push([x y+1]);
            stack.push([x-1 y]);
            stack.push([x y-1]);
            stack.push([x+1 y]);
            %8 Neighbourhood
            stack.push([x+1 y+1]);
            stack.push([x-1 y-1]);
            stack.push([x-1 y+1]);
            stack.push([x+1 y-1]);
        elseif(DensityM(y,x) <= threshold && AreaResult(y,x) == 0)
            AreaResult(y,x)=1;
        end
    else
        if(DensityM(y,x) >= threshold && AreaResult(y,x) == 0 && x-1 > 0 && y-1 > 0 && y+1<= maxY && x+1 <= maxX)
            AreaResult(y,x) = 1;
            stack.push([x y+1]);
            stack.push([x-1 y]);
            stack.push([x y-1]);
            stack.push([x+1 y]);
            %8 Neighbourhood
            stack.push([x+1 y+1]);
            stack.push([x-1 y-1]);
            stack.push([x-1 y+1]);
            stack.push([x+1 y-1]);
        elseif(DensityM(y,x) >= threshold && AreaResult(y,x) == 0)
            AreaResult(y,x)=1;
        end
    end
end
result=AreaResult;
 
 
function [x y] = MapPoint(xOriginal,yOriginal,sizeXOld,sizeYOld,sizeXNew,sizeYNew, calculateOffset)
xFactor = sizeXOld/sizeXNew;
yFactor=sizeYOld/sizeYNew;
if(calculateOffset==1)
    x=uint16((double(xOriginal)*xFactor)-xFactor/2);
    y=uint16((double(yOriginal)*yFactor)-yFactor/2);
else
    x=uint16((double(xOriginal)*xFactor));
    y=uint16((double(yOriginal)*yFactor));
end
 
function [DensityM] = getDensityDistributionsFromCSV(NucleusM,sizeXTarget,sizeYTarget)
    [sizeY, sizeX] = size(NucleusM);
    [row col] = find(NucleusM);
    row = [row;sizeY];
    col = [col;sizeX];
    row = [row;0];
    col = [col;0];
    DensityM = uint8(hist3([row col],[sizeYTarget sizeXTarget]));
    %DensityNeuron = uint8(hist3([neuronRow neuronCol],[sizeYTarget sizeXTarget]));
%Map density matrix to values between o and 1
    minBrightness=0;
    maxBrightness=max(DensityM(:));
    DensityM = double(DensityM)./255;   
    DensityM = imadjust(DensityM, [double(minBrightness)/255;double(maxBrightness)/255], [0;1]);
    DensityM = uint8(DensityM.*255); 

function level=isodata(I)
%   Source:http://www.mathworks.com/matlabcentral/fileexchange/3195-automatic-thresholding/content/isodata.m
%   ISODATA Compute global image threshold using iterative isodata method.
%   LEVEL = ISODATA(I) computes a global threshold (LEVEL) that can be
%   used to convert an intensity image to a binary image with IM2BW. LEVEL
%   is a normalized intensity value that lies in the range [0, 1].
%   This iterative technique for choosing a threshold was developed by Ridler and Calvard .
%   The histogram is initially segmented into two parts using a starting threshold value such as 0 = 2B-1, 
%   half the maximum dynamic range. 
%   The sample mean (mf,0) of the gray values associated with the foreground pixels and the sample mean (mb,0) 
%   of the gray values associated with the background pixels are computed. A new threshold value 1 is now computed 
%   as the average of these two sample means. The process is repeated, based upon the new threshold, 
%   until the threshold value does not change any more. 
  
%
%   Class Support
%   -------------
%   The input image I can be of class uint8, uint16, or double and it
%   must be nonsparse.  LEVEL is a double scalar.
%
%   Example
%   -------
%       I = imread('blood1.tif');
%       level = graythresh(I);
%       BW = im2bw(I,level);
%       imshow(BW)
%
%   See also IM2BW.
%
% Reference :T.W. Ridler, S. Calvard, Picture thresholding using an iterative selection method, 
%            IEEE Trans. System, Man and Cybernetics, SMC-8 (1978) 630-632.

% Convert all N-D arrays into a single column.  Convert to uint8 for
% fastest histogram computation.
I = im2uint8(I(:));

% STEP 1: Compute mean intensity of image from histogram, set T=mean(I)
[counts,N]=imhist(I);
i=1;
mu=cumsum(counts);
T(i)=(sum(N.*counts))/mu(end);
T(i)=round(T(i));

% STEP 2: compute Mean above T (MAT) and Mean below T (MBT) using T from
% step 1
mu2=cumsum(counts(1:T(i)));
MBT=sum(N(1:T(i)).*counts(1:T(i)))/mu2(end);

mu3=cumsum(counts(T(i):end));
MAT=sum(N(T(i):end).*counts(T(i):end))/mu3(end);
i=i+1;
% new T = (MAT+MBT)/2
T(i)=round((MAT+MBT)/2);
Threshold=T(i);
% STEP 3 to n: repeat step 2 if T(i)~=T(i-1)
while abs(T(i)-T(i-1))>=1
    mu2=cumsum(counts(1:T(i)));
    MBT=sum(N(1:T(i)).*counts(1:T(i)))/mu2(end);
    
    mu3=cumsum(counts(T(i):end));
    MAT=sum(N(T(i):end).*counts(T(i):end))/mu3(end);
    
    i=i+1;
    T(i)=round((MAT+MBT)/2); 
    Threshold=T(i);
end

 % Normalize the threshold to the range [i, 1].
level = (Threshold - 1) / (N(end) - 1);
