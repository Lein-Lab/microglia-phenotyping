%Copyright (C) 2017-2021  Martin Schmuck

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

function [summary] = MeasureGanglia (foldername, csvHandler, ExcelName, Macrophages)
foldername = [foldername '\ConvertedCellomics'];
allfiles = dir(foldername);
MacrophageSize =50;
n=0;
sSN=1;
tSO=1;
tSDP=1;
gg =1;
ll=1;

for(i=1:numel(allfiles))
    ind = strfind([foldername '\' allfiles(i).name],'NeuriteBig.png');
    % Check for nucleus images
    if(numel(ind) > 0 )%
    n = n+1;
    imgNeurite = imread([foldername '\' allfiles(i).name]);
    newfile1 = regexprep(allfiles(i).name, 'NeuriteBig','NucleusBigWatershed','ignorecase');
    imgNucleusBinary = imread([foldername '\' newfile1]);
    newfile2 = regexprep(allfiles(i).name, 'NeuriteBig','OligoBig','ignorecase');
    imgOligo = imread([foldername '\' newfile2]);
    newfile3 = regexprep(allfiles(i).name, 'NeuriteBig','NucleusBig','ignorecase');
    imgNucleus = imread([foldername '\' newfile3]);
    newfile4 = regexprep(allfiles(i).name, 'NeuriteBig','Mask','ignorecase');
    Mask = imread([foldername '\' newfile4]);
    Mask = logical(Mask);
   
   if(sum(sum(Mask))>50000) 
   %Get NucleusMatrix:
   imgNucleusBinary = logical(imgNucleusBinary); 
    Centroids = regionprops(imgNucleusBinary,'Centroid');
   sdimensions = size(Centroids);
   NucleusMatrix = zeros(size(imgNucleusBinary));
   for(k=1:sdimensions)
            x_centroid = Centroids(k).Centroid(1);
            y_centroid = Centroids(k).Centroid(2);
            x_centroid = round(x_centroid);
            y_centroid = round(y_centroid);
            NucleusMatrix(y_centroid,x_centroid)= 1;
   end
   NucleusMatrix = logical(NucleusMatrix);
   
   
   if(Macrophages==0)
   %Now we will preprocess Channel 1:
   [lehisto x] = imhist(imgNucleus(imgNucleus>0));
   level = triangle_th(lehisto,256);
   imgNNucleusBinary = im2bw(imgNucleus,level*1.2);
   imgNNucleusBinary = bwmorph(imgNNucleusBinary,'thicken',5);
   imgNNucleusMask = uint16(imgNNucleusBinary)*65535;
   imgNucleusMasked = imgNucleus - imgNNucleusMask;
   imgNucleusMasked = uint16(imgNucleusMasked);
   [lehisto x] = imhist(imgNucleusMasked(imgNucleusMasked>0));
   level = triangle_th(lehisto,256);
   imgNucleusMaskedBinary = im2bw(imgNucleusMasked,level*0.7);
   imgNucleusMaskedBinaryMask = uint16(imgNucleusMaskedBinary)*65535;
   imgNucleusMaskedBinaryMask = imcomplement(imgNucleusMaskedBinaryMask);
   imgNucleusMasked = imgNucleus - imgNucleusMaskedBinaryMask;
   imgNucleusMasked = uint16(imgNucleusMasked);
   %We will generate the binary of the Nuclei Image:
   imgNucleusMasked = double(imgNucleusMasked)./4095;
   imgNucleusMasked = uint8(imgNucleusMasked.*255);
   [lehisto x] = imhist(imgNucleusMasked(imgNucleusMasked>0));
   level = triangle_th(lehisto,256);
   imgNucleusMaskedBinary = im2bw(imgNucleusMasked,level*0.8);
   imgNucleusMaskedBinary = imfill(imgNucleusMaskedBinary,'holes');
   imgNucleusMaskedBinary = bwmorph(imgNucleusMaskedBinary,'open');
   imgNucleusMaskedBinary = bwareaopen(imgNucleusMaskedBinary,50);
   %Finally we exclude dim structures, by chekcing whether, the contain
   %bright pixels from imgNeuriteBinaryHard
   
   %Now we will preprocess Channel 2:
   [lehisto x] = imhist(imgNeurite(imgNeurite>0));
   level = triangle_th(lehisto,256);
   imgNeuriteBinary = im2bw(imgNeurite,level*1);
   %% add if statement, since some level samples are over 1, which cause conflit with im2bw. 
   if (level > 1/1.5)
       imgNeuriteBinaryHard = im2bw(imgNeurite,level);
   else
       imgNeuriteBinaryHard = im2bw(imgNeurite,level*1.5);   
   end
   imgNeuriteBinaryLow = im2bw(imgNeurite,level*0.75);
   
   imgNeuriteFilter = bwareaopen(imgNeuriteBinaryLow,1000);
   %Get of Fiber Structures:
   labeledImage = bwlabel(imgNeuriteBinary);
  measurements = regionprops(labeledImage,'Area','Perimeter');
   allAreas = [measurements.Area];
   allPerimeters = [measurements.Perimeter];
   circularities = allPerimeters.^2 ./ (4*pi*allAreas);
   keeperBlobs = circularities < 4;
   roundObjects = find(keeperBlobs);
   imgNeuriteBinaryFilterWoFiber = ismember(labeledImage, roundObjects) > 0;
   imgNeuriteBinaryRest = imgNeuriteBinary - imgNeuriteBinaryFilterWoFiber;
   %Now we try to seperate structures remaining and see if they pass the
   %circularity criterium afterwards
   % To rescue almost circles we will close these structures, but only the
   % larger particles to save computational time
   CircelRescue = bwareaopen(imgNeuriteBinaryRest,100);
   imgNeuriteBinaryRest = imgNeuriteBinaryRest - CircelRescue;
   imgNeuriteBinaryRest = uint8(imgNeuriteBinaryRest);
   imgNeuriteBinaryRest = logical(imgNeuriteBinaryRest);
   se = strel('disk',5);
   CircelRescue = imclose(CircelRescue,se);
   Holes = imfill(CircelRescue,'holes');
   Holes = logical(Holes);
   Holes = Holes - CircelRescue;
   Holes = logical(Holes);
   HolesBig = bwareaopen(Holes,1000);
   Holes = Holes - HolesBig;
   CircelRescue = CircelRescue + Holes;
   imgNeuriteBinaryRest = imgNeuriteBinaryRest + CircelRescue;
   FiberFind = bwmorph(imgNeuriteBinaryRest,'thin',16);
   FiberFindBody = bwmorph(FiberFind,'open');
   FiberFind = FiberFind - FiberFindBody;
   FiberFind = bwmorph(FiberFind,'spur',16);
   FiberFind = bwareaopen(FiberFind,30);
   se = strel('square',16);
   FiberFind = imdilate(FiberFind,se);
   imgNeuriteBinaryRest = imgNeuriteBinaryRest - FiberFind;
   imgNeuriteBinaryRest = uint8(imgNeuriteBinaryRest);
   imgNeuriteBinaryRest = logical(imgNeuriteBinaryRest);
   %Check again for circularities:
   labeledImage = bwlabel(imgNeuriteBinaryRest);
  measurements = regionprops(labeledImage,'Area','Perimeter');
   allAreas = [measurements.Area];
   allPerimeters = [measurements.Perimeter];
   circularities = allPerimeters.^2 ./ (4*pi*allAreas);
   keeperBlobs = circularities < 6;
   imgNeuriteBinaryRest = ismember(labeledImage, roundObjects) > 0;
   imgNeuriteBinary = imgNeuriteBinaryFilterWoFiber + imgNeuriteBinaryRest;
   imgNeuriteBinary = imfill(imgNeuriteBinary,'holes');
  imgNeuriteBinaryTest = bwareaopen(imgNeuriteBinary,2000);
   imgNeuriteBinarySmall =imgNeuriteBinary  - imgNeuriteBinaryTest;
   imgNeuriteBinarySmall = uint8(imgNeuriteBinarySmall);
   imgNeuriteBinarySmall = logical(imgNeuriteBinarySmall);
   imgNeuriteBinaryTest = bwareaopen(imgNeuriteBinary,100);
   %Check again for circularities:
   labeledImage = bwlabel( imgNeuriteBinarySmall);
   measurements = regionprops(labeledImage,'Area','Perimeter');
   allAreas = [measurements.Area];
   allPerimeters = [measurements.Perimeter];
   circularities = allPerimeters.^2 ./ (4*pi*allAreas);
   keeperBlobs = circularities < 6;
   roundObjects = find(keeperBlobs);
   imgNeuriteBinarySmall = ismember(labeledImage, roundObjects) > 0;
   imgNeuriteBinary = bwareaopen(imgNeuriteBinary,500);
   imgNeuriteBinary = imgNeuriteBinary + imgNeuriteBinarySmall;
   imgNeuriteBinary = bwareaopen(imgNeuriteBinary,500);
   se = strel('disk',3);
   Gaps = imdilate( imgNeuriteBinary,se);
   HolesGaps = imfill(Gaps,'holes');
   HolesGaps = HolesGaps - Gaps;
   se = strel('disk',5);
   HolesGaps =imdilate(HolesGaps,se);
   imgNeuriteBinary = imgNeuriteBinary + HolesGaps;
  imgNeuriteBinary = logical(imgNeuriteBinary);
   imgNeuriteBinary = imfill(imgNeuriteBinary,'holes');
   imgNeuriteBinary = bwmorph(imgNeuriteBinary,'open');
%    imgNeuriteSmall = bwareaopen(imgNeuriteBinary,500);
%    imgNeuriteSmall = logical(imgNeuriteBinary) - logical(imgNeuriteSmall);
%    imgNeuriteSmall = logical(imgNeuriteSmall);
%    Holes = imfill(imgNeuriteBinaryLow,'holes');
%    Holes = Holes - imgNeuriteBinaryLow;
%    Holes = bwareaopen(Holes,500);
%    HolesOrig = Holes;
%    HolesOrigL = bwlabel(HolesOrig);
%    HolesOrigN = max(max(HolesOrigL));
%    se = strel('disk',8,8);
%    Holes = imdilate(Holes,se);
%    Holes = imdilate(Holes,se);
%    Holes = imdilate(Holes,se);
%    HolesCut = imdilate(Holes,se);
%    HolesCut = imdilate(HolesCut,se);
%    HolesCut = imdilate(HolesCut,se);
%    imgNeuriteFilter = imgNeuriteBinary + HolesCut;
%    imgNeuriteFilter = imgNeuriteFilter-1;
%    imgNeuriteFilter = uint8(imgNeuriteFilter);
%    imgNeuriteFilter = logical(imgNeuriteFilter);
%    imgNeuriteBinaryLabel = bwlabel(imgNeuriteFilter);
%    imgNeuriteBinaryN = max(max(imgNeuriteBinaryLabel));
%    TestPunctaKeep = zeros(size(imgNeuriteBinary));
%    for(g=1:imgNeuriteBinaryN)
%        imgNeuriteBinaryOne = imgNeuriteBinaryLabel;
%        ii = imgNeuriteBinaryOne == g;
%        imgNeuriteBinaryOne(ii) = 60000;
%        imgNeuriteBinaryOne = imgNeuriteBinaryOne - 59995;
%        imgNeuriteBinaryOne = uint8(imgNeuriteBinaryOne);
%        imgNeuriteBinaryOne = logical(imgNeuriteBinaryOne);
%        TestOver = imgNeuriteBinaryOne+HolesCut;
%        TestOver = TestOver-1;
%        TestOver= uint8(TestOver);
%        TestOver=logical(TestOver);
%        if(sum(sum(TestOver))>0)
%            for(t=1:HolesOrigN)
%                HolesOrigOne = HolesOrigL;
%                ii = HolesOrigOne == t;
%                HolesOrigOne(ii) = 255;
%                HolesOrigOne = HolesOrigOne -250;
%                HolesOrigOne = uint8(HolesOrigOne);
%                HolesOrigOne = logical(HolesOrigOne);
%                SEE = strel('disk',2);
%                HolesOrigOneF = imdilate(HolesOrigOne,SEE);
%                Overlay = imgNeuriteBinaryOne + HolesOrigOneF;
%                Overlay = Overlay -1;
%                Overlay = uint8(Overlay);
%                Overlay = logical(Overlay);
%                if(sum(sum(Overlay))>0)
%                    HoleAdd = HolesOrigOne;
%                end   
%            end  
%            TestPunctaKeep = TestPunctaKeep + imgNeuriteBinaryOne+HolesOrigOne;
%        end
%    end    
%    
%    TestPuncta = imdilate(TestPunctaKeep,se);
%    TestPuncta = imfill(TestPuncta,'holes');
%    TestPuncta = bwmorph(TestPuncta,'thin',4);
%    TestPuncta = bwmorph(TestPuncta,'open');
%    TestPuncta = bwareaopen(TestPuncta,500);
%    Holes = imfill(imgNeuriteBinary,'holes');
%    Holes = Holes - imgNeuriteBinary;
%    imgNeuriteBinary = imgNeuriteBinary + Holes+TestPuncta;
%    imgNeuriteBinary= bwareaopen(imgNeuriteBinary,100);
%    imgNeuriteFilter = bwareaopen(imgNeuriteBinaryLow,4000);
%    se = strel('disk',8,8);
%    imgNeuriteSmall = bwareaopen( imgNeuriteBinaryHard,500);
%    imgNeuriteSmall = imgNeuriteBinaryHard - imgNeuriteSmall;
%    imgNeuriteFilter = imdilate(imgNeuriteFilter,se);
%    imgNeuriteFilter = imgNeuriteFilter+imgNeuriteSmall;
%    imgNeuriteFilter = imgNeuriteFilter-1;
%    imgNeuriteFilter = uint8(imgNeuriteFilter);
%    imgNeuriteFilter = logical(imgNeuriteFilter);
%    
%    TestPuncta = imdilate(imgNeuriteFilter ,se);
%    TestPuncta = imfill(TestPuncta,'holes');
%    TestPuncta = bwmorph(TestPuncta,'thin',4);
%    TestPuncta = bwmorph(TestPuncta,'open');
%    TestPuncta = bwareaopen(TestPuncta,2000);
%    imgNeuriteBinary = imgNeuriteBinary + TestPuncta;
%    imgNeuriteBinary = logical(imgNeuriteBinary);
%    
%    
%    se = strel('disk',5);
%    Gaps = imdilate( imgNeuriteBinary,se);
%    HolesGaps = imfill(Gaps,'holes');
%    HolesGaps = HolesGaps - Gaps;
%    se = strel('disk',5);
%    HolesGaps =imdilate(HolesGaps,se);
%    imgNeuriteBinary = imgNeuriteBinary + HolesGaps;
%    imgNeuriteBinary = logical(imgNeuriteBinary);
%    imgNeuriteBinary = imfill(imgNeuriteBinary,'holes');
%    imgNeuriteBinary = bwmorph(imgNeuriteBinary,'open');
%    %In order to complete small structures, we will replace them with
%    %imgBinaryNeuriteLow
%    imgNeuriteBinarySmall = bwareaopen(imgNeuriteBinary,1000);
%    imgNeuriteBinarySmall = imgNeuriteBinary - imgNeuriteBinarySmall;
%    imgNeuriteBinarySmall = uint8(imgNeuriteBinarySmall);
%    imgNeuriteBinarySmall = logical(imgNeuriteBinarySmall);
%    imgNeuriteBinarySmall = bwareaopen(imgNeuriteBinarySmall,500);
%    imgNeuriteBinaryLow = bwareaopen(imgNeuriteBinaryLow,1000);
%    imgNeuriteBinaryLowLabel = bwlabel(imgNeuriteBinaryLow);
%    imgNeuriteBinaryLowN = max(max(imgNeuriteBinaryLowLabel));
%    imgNeuriteBinarySmallKeep = zeros(size(imgNeuriteBinarySmall));
%    for(u=1:imgNeuriteBinaryLowN);
%        imgNeuriteBinaryLowOne = imgNeuriteBinaryLowLabel;
%        ii = imgNeuriteBinaryLowOne == u;
%        imgNeuriteBinaryLowOne(ii) = 60000;
%        imgNeuriteBinaryLowOne = imgNeuriteBinaryLowOne - 59995;
%        imgNeuriteBinaryLowOne = uint8(imgNeuriteBinaryLowOne);
%        imgNeuriteBinaryLowOne = logical(imgNeuriteBinaryLowOne);
%        TestFill = imgNeuriteBinaryLowOne + imgNeuriteBinarySmall;
%        TestFill = TestFill-1;
%        TestFill = uint8(TestFill);
%        TestFill = logical(TestFill);
%        if(sum(sum(TestFill))>0)
%           imgNeuriteBinarySmallKeep = imgNeuriteBinarySmallKeep +  imgNeuriteBinaryLowOne;
%        end
%    end
%    imgNeuriteBinarySmallKeep1 = bwareaopen(imgNeuriteBinarySmallKeep,10000);
%    imgNeuriteBinarySmallKeep = imgNeuriteBinarySmallKeep - imgNeuriteBinarySmallKeep1;
%    imgNeuriteBinarySmallKeep = uint8(imgNeuriteBinarySmallKeep);
%    imgNeuriteBinarySmallKeep = logical(imgNeuriteBinarySmallKeep);
%    imgNeuriteBinary = imgNeuriteBinary + imgNeuriteBinarySmallKeep;
    %Now we will perform a watershed:
   %Watershed for segmentation (from:
   %https://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/)
   %We will only do this on the non puncta objects to cricumvent
   %oversegmentation!
   bw = imgNeuriteBinary;
   bw2 = ~bwareaopen(~imgNeuriteBinary, 1000);
   D = -bwdist(~bw2);
   mask = imextendedmin(D,5);
   D2 = imimposemin(D,mask);
   Ld2 = watershed(D2);
   bw3 = bw;
   bw3(Ld2 == 0) = 0;
   imgNeuriteBinary = bw3;
   imgNeuriteBinary = bwareaopen(imgNeuriteBinary,1000);
   %Finally we exclude dim structures, by chekcing whether, the contain
   %bright pixels from imgNeuriteBinaryHard
   Keep =zeros(size(imgNeuriteBinary));
   imgNeuriteBinaryL = bwlabel(imgNeuriteBinary);
   imgNeuriteBinaryN = max(max(imgNeuriteBinaryL));
   for(o=1:imgNeuriteBinaryN);
       imgNeuriteBinaryOne = imgNeuriteBinaryL;
       ii = imgNeuriteBinaryOne == o;
       imgNeuriteBinaryOne(ii) = 255;
       imgNeuriteBinaryOne = imgNeuriteBinaryOne - 250;
       imgNeuriteBinaryOne = uint8(imgNeuriteBinaryOne);
       imgNeuriteBinaryOne = logical(imgNeuriteBinaryOne);
       TestBright = imgNeuriteBinaryOne + imgNeuriteBinaryHard;
       TestBright = TestBright - 1;
       TestBright = uint8(TestBright);
       TestBright = logical(TestBright);
       if(sum(sum(TestBright))>25)
           Keep = Keep + imgNeuriteBinaryOne;
       end
   end
   imgNeuriteBinary = Keep;
   
   
   
   
   
   
   
   
   %Now we will preprocess Channel 3:
   [lehisto x] = imhist(imgOligo(imgOligo>0));
   level = triangle_th(lehisto,256);
   imgOligoBinaryLow = im2bw(imgOligo,level*0.75);
   if(level>0.75)
   imgOligoBinary = im2bw(imgOligo,level*0.35);
   imgOligoBinaryLow = im2bw(imgOligo,level*0.25);
   else
   imgOligoBinary = im2bw(imgOligo,level*1);    
   end
   if(level>0.75)
       imgOligoBinaryHard = im2bw(imgOligo,level*0.5);
   else    
   imgOligoBinaryHard = im2bw(imgOligo,level*1.25);
   end
   imgOligoFilter = bwareaopen(imgOligoBinaryLow,1000);
   %Get of Fiber Structures:
   labeledImage = bwlabel(imgOligoBinary);
   measurements = regionprops(labeledImage,'Area','Perimeter');
   allAreas = [measurements.Area];
   allPerimeters = [measurements.Perimeter];
   circularities = allPerimeters.^2 ./ (4*pi*allAreas);
   keeperBlobs = circularities < 8;
   roundObjects = find(keeperBlobs);
   imgOligoFilterWoFiber = ismember(labeledImage, roundObjects) > 0;
   ImgOligoRest = imgOligoBinary - imgOligoFilterWoFiber;
   %Now we try to seperate structures remaining and see if they pass the
   %circularity criterium afterwards
   FiberFind = bwmorph(ImgOligoRest,'thin',16);
   FiberFindBody = bwmorph(FiberFind,'open');
   FiberFind = FiberFind - FiberFindBody;
   FiberFind = bwmorph(FiberFind,'spur',16);
   FiberFind = bwareaopen(FiberFind,30);
   se = strel('square',16);
   FiberFind = imdilate(FiberFind,se);
   ImgOligoRest = ImgOligoRest - FiberFind;
   ImgOligoRest = uint8(ImgOligoRest);
   ImgOligoRest = logical(ImgOligoRest);
   %Check again for circularities:
   labeledImage = bwlabel(ImgOligoRest);
   measurements = regionprops(labeledImage,'Area','Perimeter');
   allAreas = [measurements.Area];
   allPerimeters = [measurements.Perimeter];
   circularities = allPerimeters.^2 ./ (4*pi*allAreas);
   keeperBlobs = circularities < 8;
   roundObjects = find(keeperBlobs);
   ImgOligoRest = ismember(labeledImage, roundObjects) > 0;
   imgOligoBinary = imgOligoFilterWoFiber + ImgOligoRest;
   
%    
%    
%    
%    se = strel('disk',8,8);
%    imgOligoSmall = bwareaopen( imgOligoBinary,500);
%    imgOligoSmall = imgOligoBinary - imgOligoSmall;
%    imgOligoFilter = imdilate(imgOligoFilter,se);
%    imgOligoFilter = imgOligoFilter+imgOligoSmall;
%    imgOligoFilter = imgOligoFilter-1;
%    imgOligoFilter = uint8(imgOligoFilter);
%    imgOligoFilter = logical(imgOligoFilter);
%    
%    TestPuncta = imdilate(imgOligoFilter ,se);
%    TestPuncta = imfill(TestPuncta,'holes');
%    TestPuncta = bwmorph(TestPuncta,'thin',4);
%    TestPuncta = bwmorph(TestPuncta,'open');
%    TestPuncta = bwareaopen(TestPuncta,500);
%    
%    HolesClose = TestPuncta + imgOligoBinary ;
%    HolesClose = imfill(HolesClose,'holes');
%    HolesClose = HolesClose - imgOligoBinary ;
%    HolesClose = imfill(HolesClose,'holes');
%    
%    %Check again for circularities:
%    labeledImage = bwlabel(HolesClose);
%    measurements = regionprops(labeledImage,'MinorAxis','MajorAxis');
%    allAreas = [measurements.MinorAxisLength];
%    allPerimeters = [measurements.MajorAxisLength];
%    circularities = allPerimeters ./ allAreas;
%    keeperBlobs = circularities < 3;
%    roundObjects = find(keeperBlobs);
%    TestPunctaA = ismember(labeledImage, roundObjects) > 0;
%   
%    imgOligoBinary = bwareaopen(imgOligoBinary,1000);
%    imgOligoBinary = imgOligoBinary + TestPunctaA;
   imgOligoBinary = imfill(imgOligoBinary,'holes');
   imgOligoBinaryTest = bwareaopen(imgOligoBinary,2000);
   imgOligoBinarySmall = imgOligoBinary  - imgOligoBinaryTest;
   imgOligoBinarySmall = uint8(imgOligoBinarySmall);
   imgOligoBinarySmall = logical(imgOligoBinarySmall);
   imgOligoBinaryTest = bwareaopen(imgOligoBinary,100);
   %Check again for circularities:
   labeledImage = bwlabel( imgOligoBinarySmall);
   measurements = regionprops(labeledImage,'Area','Perimeter');
   allAreas = [measurements.Area];
   allPerimeters = [measurements.Perimeter];
   circularities = allPerimeters.^2 ./ (4*pi*allAreas);
   keeperBlobs = circularities < 2;
   roundObjects = find(keeperBlobs);
   imgOligoBinarySmall = ismember(labeledImage, roundObjects) > 0;
   imgOligoBinary = bwareaopen(imgOligoBinary,2000);
   imgOligoBinary = imgOligoBinary + imgOligoBinarySmall;
   imgOligoBinary = bwareaopen(imgOligoBinary,500);
      
%    imgOligoBinaryHard = imgOligoBinaryHard + TestPuncta;
%    imgOligoBinaryHard = bwmorph(imgOligoBinaryHard,'bridge');
%    Holes = imfill(imgOligoBinaryHard,'holes');
%    Holes = Holes - imgOligoBinaryHard;
%    imgOligoBinaryHard = imgOligoBinaryHard + Holes;
%    imgOligoBinaryHard = bwmorph(imgOligoBinaryHard,'bridge');
%    FiberFind = bwmorph(imgOligoBinaryHard,'thin',16);
%    FiberFindBody = bwmorph(FiberFind,'open');
%    FiberFind = FiberFind - FiberFindBody;
%    FiberFind = bwmorph(FiberFind,'spur',16);
%    FiberFind = bwareaopen(FiberFind,30);
%    se = strel('square',16);
%    FiberFind = imdilate(FiberFind,se);
%   
%    imgOligoBinary = imgOligoBinary + TestPuncta;
%    imgOligoBinary = bwareaopen(imgOligoBinary,500);
%    imgOligoBinary = bwmorph(imgOligoBinary,'bridge');
%    imgOligoBinary = logical(imgOligoBinary);
% %    TestPuncta = imgOligoBinary;
%    %Here we need to get rid of fibers again:
%    TestPunctaLabel = bwlabel(TestPuncta);
%    TestPunctaN = max(max(TestPunctaLabel));
%    TestPunctaKeep = zeros(size(imgOligoBinary));
%    for(t=1:TestPunctaN)
%        TestPunctaOne = TestPunctaLabel;
%        ii = TestPunctaOne == t;
%        TestPunctaOne(ii) = 60000;
%        TestPunctaOne = TestPunctaOne - 59995;
%        TestPunctaOne = uint8(TestPunctaOne);
%        TestPunctaOne = logical(TestPunctaOne);
%        FiberTestPuncta = bwmorph(TestPunctaOne,'thin',20);
%        FiberTestPuncta = bwmorph(FiberTestPuncta,'open');
%        if(sum(sum(FiberTestPuncta))>5)
%            TestPunctaKeep = TestPunctaKeep + TestPunctaOne;
%        else
%            tt=1;
%        end    
%    end    
%    TestPunctaKeep = logical(TestPunctaKeep);  
%    imgOligoBinary = TestPunctaKeep;
   
   %We will close some remaining gaps:
   se = strel('disk',3);
   Gaps = imdilate( imgOligoBinary,se);
   HolesGaps = imfill(Gaps,'holes');
   HolesGaps = HolesGaps - Gaps;
   se = strel('disk',5);
   HolesGaps =imdilate(HolesGaps,se);
   imgOligoBinary = imgOligoBinary + HolesGaps;
   imgOligoBinary = logical(imgOligoBinary);
   imgOligoBinary = imfill(imgOligoBinary,'holes');
   imgOligoBinary = bwmorph(imgOligoBinary,'open');
    %Now we will perform a watershed:
   %Watershed for segmentation (from:
   %https://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/)
   %We will only do this on the non puncta objects to cricumvent
   %oversegmentation!
   bw = imgOligoBinary;
   bw2 = ~bwareaopen(~imgOligoBinary, 500);
   D = -bwdist(~bw2);
   mask = imextendedmin(D,5);
   D2 = imimposemin(D,mask);
   Ld2 = watershed(D2);
   bw3 = bw;
   bw3(Ld2 == 0) = 0;
   imgOligoBinary = bw3;
   
   %Name for summary:
   if(n==1);
       SummaryName2{1,2}='Cell Area';
       SummaryName2{1,3}='Cell Intensity';
   end    
   
   %Here we try to seperate adjactend areas:
   opened= imgOligoBinary;
   openedLabel = bwlabel(opened);
   openedN = max(max(openedLabel));  
   NN_SO = 0;
   NN_SDP =0;
    Filtered = zeros(size(imgOligoBinary));

   for(k=1:openedN)
       openedOne = openedLabel;
       ii = openedOne == k;
       openedOne(ii) = 65000;
       openedOne = openedOne - 64000;
       openedOne = uint8(openedOne);
       openedOne = logical(openedOne);
       %Here we will eliminare fibers:
       FiberTest = bwmorph(openedOne,'thin',10);
       FiberTest = bwmorph(FiberTest,'open');
       FiberTest = sum(sum(FiberTest));
       if(FiberTest>0)
       %We will here also get the intensity value:
       %Generate inverse mask:
       openedInvers = imcomplement(openedOne);
       openedInvers = uint8(openedInvers)*255;
       Signal = imgOligo - openedInvers;
       Signal = uint8(Signal);
       Signal = sum(sum(Signal));
       %Now we will also check for overlapp with other channel
       CheckDoublePositive = openedOne + imgNeuriteBinary;
       CheckDoublePositive = CheckDoublePositive-1;
       CheckDoublePositive = uint8(CheckDoublePositive);
       CheckDoublePositive=logical(CheckDoublePositive);
       if(sum(sum(CheckDoublePositive))>100)
           CheckDoublePositive =1;
       else
           CheckDoublePositive =0;
       end    
      
           
           if(CheckDoublePositive==0)
               tSO=tSO+1;
               SummaryNameSO{tSO,1}=allfiles(i).name(1:3);
               SummarySO(tSO,2)= sum(sum(openedOne));
               SummarySO(tSO,3)= Signal;
%                ImageSO = ImageSO + openedOne;
               if(sum(sum(openedOne))>0)
                    NN_SO = NN_SO +1;
               end     
           else
               tSDP = tSDP+1;
               if(n<10)
                   SummaryNameSDP{tSDP,1}=allfiles(i).name(1:3);
               else
                   SummaryNameSDP{tSDP,1}=allfiles(i).name(1:3);
               end
               SummarySDP(tSDP,2)= sum(sum(openedOne));
               SummarySDP(tSDP,3)= Signal;
%                ImageSDP = ImageSDP + openedOne;
               if(sum(sum(openedOne))>0)
                    NN_SDP = NN_SDP +1;
               end
           end
      
       Filtered = Filtered + openedOne;
       end
  end    
  Filtered = logical(Filtered);  
  openedO= Filtered;
    
  %Now we have to do the same for the imgNeuriteBinary 
   opened= imgNeuriteBinary;
   openedLabel = bwlabel(opened);
   openedN = max(max(openedLabel));
   openedNeu = opened;
   NN_SN = 0;
   NN_SDP_Ctrl =0;
   
%    ImageSN = zeros(size(imgOligoBinary));
%    ImageSN_C= zeros(size(imgOligoBinary));
   for(k=1:openedN)
       openedOne = openedLabel;
       ii = openedOne == k;
       openedOne(ii) = 65000;
       openedOne = openedOne - 64000;
       openedOne = uint8(openedOne);
       openedOne = logical(openedOne);
       %We will here also get the intensity value:
       %Generate inverse mask:
       openedInvers = imcomplement(openedOne);
       openedInvers = uint8(openedInvers)*255;
       Signal = imgNeurite - openedInvers;
       Signal = uint8(Signal);
       Signal = sum(sum(Signal));
       
       %Now we will also check for overlapp with other channel
       CheckDoublePositive = openedOne + imgOligoBinary;
       CheckDoublePositive = CheckDoublePositive-1;
       CheckDoublePositive = uint8(CheckDoublePositive);
       CheckDoublePositive=logical(CheckDoublePositive);
       if(sum(sum(CheckDoublePositive))>100)
           CheckDoublePositive =1;
       else
           CheckDoublePositive =0;
       end    
       
           
           if(CheckDoublePositive==0)
               sSN=sSN+1;
               SummaryNameSN{sSN,1}=allfiles(i).name(1:3);
               SummarySN(sSN,2)= sum(sum(openedOne));
               SummarySN(sSN,3)= Signal;
%                ImageSN = ImageSN + openedOne;
               if(sum(sum(openedOne))>0)
                    NN_SN = NN_SN +1;
               end
           else
               
           end
   end     
   if(n<10)
        SummarySumName{n+1,1}=allfiles(i).name(1:3);
   else
        SummarySumName{n+1,1}=allfiles(i).name(1:3);
   end  
%    ImageSO = logical(ImageSO);
%    ImageSDP = logical(ImageSDP);
%    ImageSO_C = logical(ImageSO_C);
%    ImageSDP_C = logical(ImageSDP_C);
%    ImageSN = logical(ImageSN);
%    ImageSN_C = logical(ImageSN_C);
   
   SummarySumName2{1,1} = 'Well Name';
   SummarySumName2{1,2} = 'Number of only single substance p postive neurons';
   SummarySumName2{1,3} = 'Number of only single cgrp postive neurons';
   SummarySumName2{1,4} = 'Number of single double postive neurons';
   SummarySumName2{1,5} = 'Area ROI';
   
 
   
   
   SummarySum(n+1,2)= NN_SN;
   SummarySum(n+1,3)= NN_SO;
   SummarySum(n+1,4)= NN_SDP;
   SummarySum(n+1,5)= sum(sum(Mask));
   

   
   %Now we save the binary images
   newFile1 = regexprep(allfiles(i).name, 'NeuriteBig','Binary','ignorecase');
   imwrite(openedO,[foldername '/' newFile1]);
   newFile2 = regexprep(allfiles(i).name, 'NeuriteBig','Skeleton','ignorecase');
   imwrite(openedNeu,[foldername '/' newFile2]);
   imgNucleusMaskedBinary = uint8(imgNucleusMaskedBinary)*255;
   imgNucleusMaskedBinarySmall = imresize(imgNucleusMaskedBinary,0.1);
   newFile3 = regexprep(allfiles(i).name, 'NeuriteBig','AstroBig','ignorecase');
   imwrite(imgNucleusMaskedBinary,[foldername '/' newFile3]);
   newFile4 = regexprep(allfiles(i).name, 'NeuriteBig','AstroSmall','ignorecase');
   imwrite(imgNucleusMaskedBinarySmall,[foldername '/' newFile4]);
   
   
   
   
   
    summary=0;
   else
   %Here we will analyse macrophages:
    %Now we will preprocess Channel 1:
   [lehisto x] = imhist(imgNucleus(imgNucleus>0));
   level = triangle_th(lehisto,256);
   imgNNucleusBinary = im2bw(imgNucleus,level*0.75);
   imgNNucleusBinary = imfill(imgNNucleusBinary,'holes');
   imgNNucleusBinary = bwmorph(imgNNucleusBinary,'open');
    %Now we will perform a watershed:
   %Watershed for segmentation (from:
   %https://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/)
   %We will only do this on the non puncta objects to cricumvent
   %oversegmentation!
   bw = imgNNucleusBinary;
   D = -bwdist(~bw);
   mask = imextendedmin(D,0.7);
   D2 = imimposemin(D,mask);
   Ld2 = watershed(D2);
   bw3 = bw;
   bw3(Ld2 == 0) = 0;
   imgNNucleusBinary = bw3;
   %Now we will genearte a nuclei matrix:
   NucleusM = uint16(zeros(size(imgNNucleusBinary)));
   Nuclei = regionprops(imgNNucleusBinary,'centroid');
   
   for(y=1:length(Nuclei))
       x_centroid = Nuclei(y).Centroid(1);
       y_centroid = Nuclei(y).Centroid(2);
       NucleusM(uint16(y_centroid),uint16(x_centroid)) = 1;
   end   
   NucleusM = logical(NucleusM);
   
   
   %Now we will preprocess Channel 2:
   [lehisto x] = imhist(imgNeurite(imgNeurite>0));
   level = triangle_th(lehisto,256);
   if(level<0.3)
       imgNeuriteBinary = im2bw(imgNeurite,level*0.75);
   else    
       imgNeuriteBinary = im2bw(imgNeurite,level*0.5);
   end    
   imgNeuriteBinary = bwareaopen(imgNeuriteBinary,1000);
   se = strel('disk',3);
   imgNeuriteBinary= imopen(imgNeuriteBinary, se);
   %Watershed for segmentation (from:
   %https://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/)
   bw = imgNeuriteBinary;
   bw2 = ~bwareaopen(~imgNeuriteBinary, 1000);
   D = -bwdist(~bw2);
   mask = imextendedmin(D,2);
   D2 = imimposemin(D,mask);
   Ld2 = watershed(D2);
   bw3 = bw;
   bw3(Ld2 == 0) = 0;
   imgNeuriteBinary = bw3;
   imgNeuriteBinary = bwareaopen(imgNeuriteBinary,1000);
   imgNeuriteBinary = imfill(imgNeuriteBinary,'holes');
   
   %Now we get the NeuronNumber and area:
   imgNeuriteBinaryL = bwlabel(imgNeuriteBinary);
   imgNeuriteBinaryN = max(max(imgNeuriteBinaryL));
   %Here we will get the seperate values for cell area and intensity
   for(w=1:imgNeuriteBinaryN);
       gg = gg+1;
       imgNeuriteBinaryOne = imgNeuriteBinaryL;
       ii = imgNeuriteBinaryOne  == w;
       imgNeuriteBinaryOne(ii) = 60000;
       imgNeuriteBinaryOne = imgNeuriteBinaryOne  - 59995;
       imgNeuriteBinaryOne = uint8(imgNeuriteBinaryOne);
       imgNeuriteBinaryOne = logical(imgNeuriteBinaryOne);
       AreaNeuron = sum(sum(imgNeuriteBinaryOne));
       %Now we create a mask to get intensity of cell:
       imgNeuriteBinaryOneInvert = imcomplement(imgNeuriteBinaryOne);
       imgNeuriteBinaryOneInvert = uint8(imgNeuriteBinaryOneInvert)*255;
       Intensity = imgNeurite - imgNeuriteBinaryOneInvert;
       Intensity = uint8(Intensity);
       Intensity = sum(sum(Intensity));
       Summary2(gg,2) = AreaNeuron;
       Summary2(gg,3) = Intensity;
       Summary2_Name{gg,1} = allfiles(i).name(1:3);
       
   end    
   imgNeuriteBinaryArea = sum(sum(imgNeuriteBinary));
   
   
   
   
   
   %Generate Summaries:
   %Name for summary:
   if(n==1);
       SummaryName2{1,2}='Cell Area Neurons';
       SummaryName2{1,3}='Cell Intensities Neurons';
       SummaryName3{1,2}='Cell Area Macrophages';
       SummaryName3{1,3}='Cell Intensities Macrophages';
       SummaryNameGeneral{1,2} = 'Neuron Number';
       SummaryNameGeneral{1,3} = 'Number of Macrophages';
       SummaryNameGeneral{1,4} = 'ROI Area';
   end
   
   
   
   
   %Now we will preprocess Channel 3:
   [lehisto x] = imhist(imgOligo(imgOligo>0));
   level = triangle_th(lehisto,256);
   imgOligoBinary = im2bw(imgOligo,level);
   imgOligoBinaryHard = im2bw(imgOligo,level);
   % We will first eliminate to big structures;
   imgOligoBinaryBig = bwareaopen(imgOligoBinary,2000);
   imgOligoBinarySmall = imgOligoBinary- imgOligoBinaryBig;   
   se = strel('disk',3);
   imgOligoBinarySmall = imclose(imgOligoBinarySmall,se);
   %Now we need to get of very small structures:
   imgOligoBinary = bwareaopen(imgOligoBinarySmall,MacrophageSize);
   %Now we will close gaps:
   se = strel('disk',4);
   Holes = imdilate(imgOligoBinary,se);
   HolesFill = imfill(Holes,'holes');
   HolesFill = HolesFill - Holes;
   se = strel('disk',4);
   HolesFill = imdilate(HolesFill,se);
   
   imgOligoBinary = imgOligoBinary + HolesFill;
   imgOligoBinary = logical(imgOligoBinary);
   imgOligoBinary = imfill(imgOligoBinary,'holes');
   imgOligoBinarySmall = bwareaopen(imgOligoBinary,500);
   imgOligoBinarySmall = imgOligoBinary - imgOligoBinarySmall;
   imgOligoBinaryBig = imgOligoBinary- imgOligoBinarySmall;
   %Now we will check for overlapp and sort also estimate number of cells
   %per macrophage cluster
   imgOligoBinaryL = bwlabel(imgOligoBinaryBig);
   imgOligoBinaryN = max(max(imgOligoBinaryL));
   imgOligoBinaryKeep = zeros(size(imgOligoBinary));
   for(d=1:imgOligoBinaryN)
       imgOligoBinaryOne  = imgOligoBinaryL;
       ii = imgOligoBinaryOne == d;
       imgOligoBinaryOne(ii) = 60000;
       imgOligoBinaryOne = imgOligoBinaryOne -59995;
       imgOligoBinaryOne = uint8(imgOligoBinaryOne);
       imgOligoBinaryOne = logical(imgOligoBinaryOne);
       OverlayNuc = imgOligoBinaryOne + imgNNucleusBinary;
       OverlayNuc = OverlayNuc - 1;
       OverlayNuc = uint8(OverlayNuc);
       OverlayNuc = logical(OverlayNuc);
       if(sum(sum(OverlayNuc))>5)
           imgOligoBinaryKeep  = imgOligoBinaryOne+ imgOligoBinaryKeep;
       end    
   end  
   imgOligoBinary =    imgOligoBinarySmall+imgOligoBinaryKeep;
   %Get of Fiber Structures:
   labeledImage = bwlabel(imgOligoBinary);
   measurements = regionprops(labeledImage,'MajorAxis','MinorAxis');
   allAreas = [measurements.MinorAxisLength];
   allPerimeters = [measurements.MajorAxisLength];
   circularities = allPerimeters ./ allAreas;
   keeperBlobs = circularities < 4;
   roundObjects = find(keeperBlobs);
   imgOligoFilterWoFiber = ismember(labeledImage, roundObjects) > 0;
   imgOligoBinary = imgOligoFilterWoFiber;
   imgOligoBinaryOrig = imgOligoBinary;
   %Here we perform an opening to sperate adjacent cells in huge clusters:
   imgOligoBinaryBig = bwareaopen(imgOligoBinary,250);
   imgOligoBinary = imgOligoBinary - imgOligoBinaryBig;
   imgOligoBinaryBig = bwmorph(imgOligoBinaryBig,'open');
   imgOligoBinary = imgOligoBinary + imgOligoBinaryBig;
%    %Now we validate, that particles have bright pixels in them from
%    %imgOligoBinaryHard:
%    imgOligoKeep = zeros(size(imgOligoBinary));
%    imgOligoBinaryL = bwlabel(imgOligoBinary);
%    imgOligoBinaryN = max(max(imgOligoBinaryL));
%    for(e=1:imgOligoBinaryN)
%        imgOligoBinaryOne = imgOligoBinaryL;
%        ii = imgOligoBinaryOne == e;
%        imgOligoBinaryOne(ii) = 60000;
%        imgOligoBinaryOne = imgOligoBinaryOne - 59995;
%        imgOligoBinaryOne = uint8(imgOligoBinaryOne);
%        imgOligoBinaryOne = logical(imgOligoBinaryOne);
%        TestBright = imgOligoBinaryOne + imgOligoBinaryHard;
%        TestBright = TestBright -1;
%        TestBright = uint8(TestBright);
%        TestBright = logical(TestBright);
%        if(sum(sum(TestBright))>1)
%            imgOligoKeep = imgOligoKeep + imgOligoBinaryOne;
%        end
%    end
%    imgOligoBinary = imgOligoKeep;
     
   
   
   %Now we extract size, intensity and number of macrophages:
   imgOligoBinaryL = bwlabel(imgOligoBinary);
   imgOligoBinaryN = max(max(imgOligoBinaryL));
   NumberofCells = 0;
   %Here we will get the seperate values for cell area and intensity
   for(w=1:imgOligoBinaryN);
       ll = ll+1;
       imgOligoBinaryOne = imgOligoBinaryL;
       ii = imgOligoBinaryOne  == w;
       imgOligoBinaryOne(ii) = 60000;
       imgOligoBinaryOne = imgOligoBinaryOne  - 59995;
       imgOligoBinaryOne = uint8(imgOligoBinaryOne);
       imgOligoBinaryOne = logical(imgOligoBinaryOne);
       AreaMacrophage = sum(sum(imgOligoBinaryOne));
       %Now we create a mask to get intensity of cell:
       imgOligoBinaryOneInvert = imcomplement(imgOligoBinaryOne);
       imgOligoBinaryOneInvert = uint8(imgOligoBinaryOneInvert)*255;
       Intensity = imgOligo - imgOligoBinaryOneInvert;
       Intensity = uint8(Intensity);
       Intensity = sum(sum(Intensity));
       %Now we also get number of macrophages and will discriminate between
       %cluster and isolated ones
       OverlayNuc = imgOligoBinaryOne + imgNNucleusBinary;
       OverlayNuc = OverlayNuc - 1;
       OverlayNuc = uint8(OverlayNuc);
       OverlayNuc = logical(OverlayNuc); 
       OverlayNucL = bwlabel(OverlayNuc);
       OverlayNucN = max(max(OverlayNucL));
       
       if(OverlayNucN==0)
           OverlayNucN = 1;
       end
       if(OverlayNucN ==1)
       Summary3(ll,2) = AreaMacrophage;
       Summary3(ll,3) = Intensity;
       Summary3_Name{ll,1} = allfiles(i).name(1:3);
       else
           for(j=1:OverlayNucN)
           if(j>1)
               ll = ll+1;
           end    
           Summary3(ll,2) = AreaMacrophage/OverlayNucN;
           Summary3(ll,3) = Intensity/OverlayNucN;
           Summary3_Name{ll,1} = allfiles(i).name(1:3);
           
           end
       end  
       NumberofCells = NumberofCells + OverlayNucN;
   end    
   
   
   
   
   %Now we need to write the images:
   SummaryGeneral(n+1,2) =  imgNeuriteBinaryN;
   SummaryGeneral(n+1,3) = NumberofCells; 
   SummaryGeneral(n+1,4) = sum(sum(Mask));
   SummaryGeneral_Name{n+1,1} = allfiles(i).name(1:3);
   newFile1 = regexprep(allfiles(i).name, 'NeuriteBig','Binary','ignorecase');
   imwrite(imgNeuriteBinary,[foldername '/' newFile1]);
   newFile2 = regexprep(allfiles(i).name, 'NeuriteBig','Skeleton','ignorecase');
   imwrite(imgOligoBinary,[foldername '/' newFile2]);
  
   summary=0;
   end   
   end
    end
end    

%Now we write everything into an excell sheet:
if(Macrophages == 0)
filename = ExcelName;
sheet =1;
try
xlswrite(filename,SummarySN,sheet);
xlswrite(filename,SummaryNameSN,sheet)
xlswrite(filename,SummaryName2,sheet)
catch
end
try
sheet =2;
xlswrite(filename,SummarySO,sheet);
xlswrite(filename,SummaryNameSO,sheet)
xlswrite(filename,SummaryName2,sheet)
catch
end    
sheet =3;
try
xlswrite(filename,SummarySDP,sheet);
xlswrite(filename,SummaryNameSDP,sheet)
xlswrite(filename,SummaryName2,sheet)
catch
end    
sheet=4;
try
xlswrite(filename,SummarySum,sheet);
xlswrite(filename,SummarySumName,sheet)
xlswrite(filename,SummarySumName2,sheet)
catch
end    
else

sheet=1;    
try
xlswrite(filename,Summary2,sheet);
xlswrite(filename,Summary2_Name,sheet)
xlswrite(filename,SummaryName2,sheet)
catch
end
try
sheet =2;
xlswrite(filename,Summary3,sheet);
xlswrite(filename,Summary3_Name,sheet)
xlswrite(filename,SummaryName3,sheet)
catch
end    
sheet =3;
try
xlswrite(filename,SummaryGeneral,sheet);
xlswrite(filename,SummaryGeneral_Name,sheet)
xlswrite(filename,SummaryNameGeneral,sheet)
catch
end 

    
    
end
if(Macrophages ==1)
    filename = ExcelName;
    sheet = 1;
 xlswrite(filename,Summary2,sheet);
 xlswrite(filename,Summary2_Name,sheet)
 xlswrite(filename,SummaryName2,sheet)    
 sheet = 2;
 xlswrite(filename,Summary3,sheet);
 xlswrite(filename,Summary3_Name,sheet)
 xlswrite(filename,SummaryName3,sheet)  
 sheet = 3;
 xlswrite(filename,SummaryGeneral,sheet);
 xlswrite(filename,SummaryGeneral_Name,sheet)
 xlswrite(filename,SummaryNameGeneral,sheet) 
end 


end