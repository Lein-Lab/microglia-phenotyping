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

% Afer function are the variables obtained from Omnisphero in [], after =
% CorticalNeurons is variable which is returned to Omnisphero ()
function [wellList, foldername, ScalingFactor, csvHandler, neuronHandler,  Th] = CorticalNeuronsNucleusSkeletonownNeuronDetect (csvHandler, neuronHandler) 
% This is a quick script to convert 16-bit images of cortical neurons into
% binary image, to assess their center point and to generate a
% 'NucleusImage'.
n = 0;
u=0;
minOverlapp = 1;
A = 1000;
B = 4;
N = 50;
C = 250;
D=10;
E = 15000;
IPS = 0;
Cellomics = 0;
NeuronBody=200;
NeuriteThreshold = 0.7;
summary = cell(10,6);
NeuronNucleiArea = 500;
SE= strel('disk',8);
RatioFill = 0.6;
Maxlengthsquare = 100;
%Get Image folder
foldername = uigetdir;
mkdir(foldername,'ConvertedCellomics');
foldername1 = [foldername '/ConvertedCellomics'];
%Iterate over all files in this folder
allfiles = dir(foldername);
wellList = cell(0);
BackgroundBallFactor = 3;
SizeLeftoffBall = 5000;
TriangleBackground = 1;
TriangleBackgroundBall = 3;
answer = questdlg('Are these human cells?')
if(sum(answer(1:2)=='Ye')==2)
    IPS=1;
else
    IPS=0;
end


% Since we have different sized images we have to make them the same size:
% We will first find the biggest and add black to the others!
for(t=1:numel(allfiles))
   %Check if file ends with .tif or .tiff
    ind = strfind([foldername '/' allfiles(t).name],'.tif');
    ind2 = strfind([foldername '/' allfiles(t).name],'.jpg');
    ind3 = strfind([foldername '/' allfiles(t).name],'.jpeg');
    ind4 = strfind([foldername '/' allfiles(t).name],'.TIF');
    
    if(isempty(ind4)==1)
    if(numel(ind) > 0 | numel(ind2) > 0 | numel(ind3) > 0)
        
        u = u +1;
            
    %     Load RGB image
            imgRGB = imread([foldername '/' allfiles(t).name]);
            %Split Channels
            nucleusPic = imgRGB(:,:,3);
            nucleusImg = nucleusPic;
           
            if(u == 1)
                [r c] = size(nucleusImg);
            else
                if(u>1)
                  [r1 c1] = size(nucleusImg);
                  if(r<r1)
                      r = r1;
                  end
                  if(c<c1)
                      c = c1;
                  end
                end
            end  
        end
    end        
end

    for(i=1:numel(allfiles))
        %Check if file ends with .tif or .tiff
        ind = strfind([foldername '/' allfiles(i).name],'.tif');
        ind2 = strfind([foldername '/' allfiles(i).name],'.jpg');
        ind3 = strfind([foldername '/' allfiles(i).name],'.jpeg');
        ind4 = strfind([foldername '/' allfiles(i).name],'.TIF');
        ind5 = [];
        if(ind4>1)
             Extension = allfiles(i).name;
             Extension = Extension(end-5:end-4);
             if(strcmp(Extension,'d0')==1)
                 ind5 = 1;
             end
        else
            ind5 =[];
        end
        
        if(numel(ind) > 0 | numel(ind2) > 0 | numel(ind3) > 0 | numel(ind4) > 0)
        Found = 0;     
             
             
             
             if(isempty(ind5) == 1 && isempty(ind4) == 1)
                 n = n + 1;
                %Load RGB image
                imgRGB = imread([foldername '/' allfiles(i).name]);
                %Split Channels
                nucleusPic = imgRGB(:,:,3);
                nucleusImg = nucleusPic;
                imgSOX2 = imgRGB(:,:,2);
                imgNeurite = imgRGB(:,:,1); 
                Found = 1;
             else
                 %Check if extension is d0
                 if(ind5 ==1)
                     Cellomics =1;
                     n = n + 1;
                     Name = allfiles(i).name;
                     nucleusPic = imread([foldername '/' Name]);
                     newFile1 = regexprep(Name, 'd0','d1','ignorecase');
                     imgNeurite  = imread([foldername '/' newFile1]);
                     newFile2 = regexprep(Name, 'd0','d2','ignorecase');
                     imgSOX2  = imread([foldername '/' newFile2]);
                     Found = 1;
                 end
             end
          if(Found == 1)   
            if(IPS==1)
                if(Cellomics==1)
                    backgroundballNeuriteSmall =  imopen(uint8(imgNeurite),strel('ball',10,10));
                     imgNeuriteBallSmall = imgNeurite - backgroundballNeuriteSmall;
            %         imgNeuriteBallBig = imgNeurite - backgroundballNeuriteBig;
            %         imgNeuriteBallBigBig = imgNeurite - backgroundballNeuriteBigBig;
                    %imgNeuriteBall = imgNeuriteBallSmall; %+ imgNeuriteBallBig + imgNeuriteBallBigBig;
                    imgNeuriteBall = imgNeurite;
                else
                    backgroundballNeuriteSmall =  imopen(uint8(imgNeurite),strel('ball',10,10));
                     imgNeuriteBallSmall = imgNeurite - backgroundballNeuriteSmall;
            %         imgNeuriteBallBig = imgNeurite - backgroundballNeuriteBig;
            %         imgNeuriteBallBigBig = imgNeurite - backgroundballNeuriteBigBig;
                    imgNeuriteBall = imgNeuriteBallSmall; %+ imgNeuriteBallBig + imgNeuriteBallBigBig;
                end
            else
                backgroundballNeuriteSmall =  imopen(uint8(imgNeurite),strel('ball',10,10));
                backgroundballNeuriteBig =  imopen(uint8(imgNeurite),strel('ball',50,50));
                imgNeuriteBallSmall = imgNeurite - backgroundballNeuriteSmall;
                imgNeuriteBallBig = imgNeurite - backgroundballNeuriteBig;
                imgNeuriteBall = imgNeuriteBallSmall + imgNeuriteBallBig;
            end
    %         backgroundballNeuriteBig =  imopen(uint8(imgNeurite),strel('ball',20,20));
    %         backgroundballNeuriteBigBig =  imopen(uint8(imgNeurite),strel('ball',50,50));


    %         binaryimgNeurite=adaptivethreshold(imgNeurite,500,0.000001,0);
    %         binaryimgNeurite = bwareaopen(binaryimgNeurite,10);
    %         binaryimgNeurite = uint8(binaryimgNeurite);
    %         binaryimgNeurite = binaryimgNeurite*255;

    %        imgNeuriteMask = imgNeuriteBall - binaryimgNeurite;
    %        imgNeuriteMask = uint8(imgNeuriteMask);
    %        imgNeuriteBall = imgNeuriteBall - imgNeuriteMask;
    %         if(IPS==0)
    %           mu = mean(double(imgNeurite(:)));
    %           sigma = std(double(imgNeurite(:)));  
    %           snr= mu/sigma;
    %           
    %           
    %         %For rat:
    %         level = isodata(imgNeurite);
    %         
    %         if(snr<1)  
    %             if(level<0.1)
    %                 level =0.1;
    %             
    %             end
    %         else
    %             if(level<0.2)
    %                 level =0.25;
    %             end
    %         end
    %             
    %         HardMask = im2bw(imgNeurite,level);
    %         
    % %         imgRemain = imgNeurite - uint8(HardMask)*255;
    % %         imgRemain = uint8(imgRemain);
    % %         level = isodata(imgRemain(imgRemain>0));
    % %         if(level<0.1)
    % %             level = 0.1;
    % %         end
    % %         SoftMask = im2bw(imgRemain,level*2);
    %         imgNeuriteFill = HardMask; %+ SoftMask;
    %         imgNeuriteBall = imgNeurite;
    %         
    %         else

            try
        level = isodata(imgNeuriteBall);
        catch
            
            level = graythresh(imgNeuriteBall);
        end
        
        if(IPS ==1)
            imgNeuriteFill = im2bw(imgNeuriteBall,level);
        else
            level = graythresh(imgNeuriteBallSmall);
            imgNeuriteBallSmallFill = im2bw(imgNeuriteBallSmall,level);
            level = graythresh(imgNeuriteBallBig);
            imgNeuriteBallBigFill = im2bw(imgNeuriteBallBig,level);
            imgNeuriteFill = imgNeuriteBallSmallFill + imgNeuriteBallBigFill;
        end
            EdgeBinary = bwareaopen(imgNeuriteFill,10);
    %         end
            if(Cellomics ==0)
                   AreaCanny = 100;
                % Areas filtered from the edge processed image using prewitt method
                AreaPrewitt = 10;
                % Areas filtered from the edge processed image using sobel method
                AreaSobel = 10;
                % Areas filtered from the edge processed image using Roberts method
                AreaRoberts = 30;
                % Areas filtered from the edge processed image using Log method
                AreaLogHigh = 100;
                AreaLogLow = 30;
                % Edge Detectors:
                E1 = edge(imgNeuriteBall,'canny');

                % Filter image E1 for particles bigger AreaCanny
                E1old = bwareaopen(E1,AreaCanny);
                % Perform edge Prewitt
                E2 = edge(imgNeuriteBall,'prewitt');
                % Filter image E2 for particles bigger AreaPrewitt
                E2old = bwareaopen(E2,AreaPrewitt);
                % Perform edge Sobel
                E3 = edge(imgNeuriteBall,'sobel');
                % Filter image E3 for particles bigger AreaSobel
                E3old = bwareaopen(E3,AreaSobel);
                % Perform edge Roberts
                E4 = edge(imgNeuriteBall,'Roberts');
                % Filter image E4 for particles bigger AreaRoberts
                E4old = bwareaopen(E4,AreaRoberts);
                % Perform edge log
                E5 = edge(imgNeuriteBall,'log');
                % Filter image E5 for particles bigger AreaLogHigh
                E5old = bwareaopen(E5,AreaLogHigh);
                % Filter image E5 for particles bigger AreaLogLow
                E5 = bwareaopen(E5,AreaLogLow);

                EdgeBinary = E2+E3+E4+E5;


                EdgeBinary = logical(EdgeBinary);
                EdgeBinary = EdgeBinary + imgNeuriteFill;
                EdgeBinary = logical(EdgeBinary);
                EdgeBinary = bwmorph(EdgeBinary,'bridge');

                % Create image in which all holes are filled
                BW = imfill(EdgeBinary,'holes');
                % Substract original EdgeBinary from filled image to obtain only
                % areas of holes
                BW2 = BW - EdgeBinary;
                % Perform image opening
                if(IPS==1)
                BW2 = bwmorph(BW2,'open');
                end
                % Only filter Holes bigger than HoleArea
                HoleArea = 50;
                Holes = bwareaopen(BW2,HoleArea);
                if(IPS==0)
                    HolesElon= imopen(BW2,strel('disk',3));
                    HolesElon = Holes - HolesElon;
                    Holes = Holes - HolesElon;
                end
                % Substract Holes from BW to obtain only the small holes
                EdgeBinary = BW - Holes;
                if(IPS ==0)
                    Closed = imclose(EdgeBinary,strel('disk',5));
                    Closed = Closed - EdgeBinary;
                    ClosedElon = imopen(Closed,strel('disk',3));
                    Closed = Closed - ClosedElon;
                end
                % Perform brdige again to connect objects
                EdgeBinary = bwmorph(EdgeBinary,'bridge');    
                % Repeat a second time to fill new holes
                BW = imfill(EdgeBinary,'holes');
                BW2 = BW - EdgeBinary;
                BW2 = bwmorph(BW2,'open');
                Holes = bwareaopen(BW2,HoleArea);
                EdgeBinary = BW - Holes;
                % Perform majority to fill small one pixel missing holes
                EdgeBinary = bwmorph(EdgeBinary,'majority');
                % Sort again for big particles
                EdgeBinary = bwareaopen(EdgeBinary,50);

                % Gap filling algorithm:
                %1) Find all gaps!
                %1) dilatate image, so that gaps are closed
                SE = strel('disk',5);
                Gap = imdilate(EdgeBinary,SE);
                %2) Then erode again:
                Gap = imerode(Gap,SE);
                % Now get difference resulting in the actual 'gaps'
                Gap = Gap - EdgeBinary;
                % Convert image to 8-bit
                Gap = uint8(Gap);
                % Convert image to logical
                Gap = logical(Gap);
                % Dilatate Gaps
                SE = strel('disk',4);
                GapBig = imdilate(Gap,SE);
                % Identify big gaps, which are usually not real connections
                GapBig = bwareaopen(GapBig,500);
                % Substract big gaps to only obtain small true gaps
                Gap = Gap - GapBig;
                % Convert image to 8-bit
                Gap = uint8(Gap);
                % Convert image to logical
                Gap = logical(Gap);
                % Remove gaps attached to the image border
                Gap = imclearborder(Gap);
                % We should now further reduce reduntant gaps, by adding the
                % dialted image as a mask and only keep gaps within this mask
                % Imdialte image again, but this time even more to really close all
                % gaps
                SE2 = strel('disk',6);
                GapCorrected = imdilate(EdgeBinary,SE2);
                % Filter for big areas 
                GapCorrected = bwareaopen(GapCorrected,10000);
                % Add gaps to this mask. Gaps overlapping with the mask will have a
                % value of 2
                GapCorrected = GapCorrected + Gap;
                % Substraction of 1 results in gaps solely located on the mask
                GapCorrected = GapCorrected -1;
                % Convert image to 8-bit
                GapCorrected = uint8(GapCorrected);
                % Convert image to logical
                GapCorrected = logical(GapCorrected);
                % Rename GapCorrected to Gap 
                Gap = GapCorrected;
                %Now check expand each area for one pixel:
                GapFat = bwmorph(Gap,'thicken',1);
                % Now test for each area, wheather it is overlapping with two areas
                % from EdgeBinary which is only the case for a gap ar a pit:
                % Label all Gaps original
                GapLabel = bwlabel(Gap);
                % Label dilated Gaps
                GapLabelFat = bwlabel(GapFat);
                % Determine area with highest label
                sG = max(max(GapLabel));
                % Create image for GapAdd
                GapAdd = zeros(size(EdgeBinary));
                for(p=1:sG)
                    %Extract one of the gap areas
                    GapKeep = GapLabel;
                    % Index all pixles in GapKeep equal to p
                    ii = GapKeep == p;
                    % Since we could have more than 255 areas, we create a value
                    % depending on sG
                    ValueAdd = (100*sG);
                    % Set all indexed pixles to ValueAdd
                    GapKeep(ii) = ValueAdd;
                    ValueSubstract = ValueAdd -5; 
                    % Substract ValueSubstract form Gap Keep to only have positve
                    % pixel values for indexed pixels
                    GapKeep = GapKeep -ValueSubstract;
                    % Convert image to 8-bit
                    GapKeep = uint8(GapKeep);
                     % Convert image to logical
                    GapKeep = logical(GapKeep);
                    % Repeat for the dilated Gap imahe
                    GapOne = GapLabelFat;
                    ii = GapOne == p;
                    GapOne(ii) = ValueAdd;
                    GapOne = GapOne -ValueSubstract;
                    GapOne = uint8(GapOne);
                    GapOne = logical(GapOne);
                    %Now add labeled EdgeBinary and check if there are three
                    %or more different pixel values (means it connects two areas):
                    % Add GapOne*150 to have a high value which can be substracted
                    ConnectOne = (GapOne*150) + EdgeBinary;
                    % Substract 149 to only keep pixels from GapOne, with values of
                    % 1 (no overlapp) and 2 (overlapp with EdgeBinary)
                    ConnectOne = ConnectOne - 149;
                    % Convert image to 8-bit to get rid of negative values
                    ConnectOne = uint8(ConnectOne);
                    % Substract one to only keep pixell of overlapp
                    ConnectOne = ConnectOne - 1;
                    % Convert image to 8-bit to get rid of negative values
                    ConnectOne = uint8(ConnectOne);
                    % Perform bridge to connecte very near overlapp areas
                    ConnectOne = bwmorph(ConnectOne,'bridge');
                    % Look for number of areas within the Connected image
                    AreaOne = regionprops(ConnectOne,'Area');
                    AreaOne = size(AreaOne);
                    AreaOne = AreaOne(1);
                    % If the number of Areas is bigger than 1, the original gap is
                    % kept
                    if(AreaOne>1)
                        if(sum(sum(ConnectOne))<21)
                        GapAdd = GapAdd + GapKeep;
                        end
                    end
                end    
                % Converti image to 8-bit
                EdgeBinary = uint8(EdgeBinary);
                %Convert image to logical
                EdgeBinary = logical(EdgeBinary);
                % Add GapAdd to EdgeBinary to close Gaps
                EdgeBinary = EdgeBinary + GapAdd;
                % Converti image to 8-bit
                EdgeBinary = uint8(EdgeBinary);
                %Convert image to logical
                EdgeBinary = logical(EdgeBinary);

                % Fill holes again:
                EdgeBinary = bwmorph(EdgeBinary,'bridge');
            end



            % Count Synapsin!

            ballimgSOX2 = imopen(uint8(imgSOX2),strel('ball',10,10));
            imgSOX2Ball = imgSOX2 - ballimgSOX2;
            level = isodata(imgSOX2);
            if(Cellomics ==1)
                imgSOX2binary = im2bw(imgSOX2Ball,level);
            else
                imgSOX2binary = im2bw(imgSOX2Ball,level*0.25);
            end
            imgSOX2binaryBig = bwareaopen(imgSOX2binary,75);
            imgSOX2binary = imgSOX2binary - imgSOX2binaryBig;
            EdgeBinary = uint8(EdgeBinary);
            EdgeBinary = EdgeBinary*255;
            imgSOX2binary = uint8(imgSOX2binary);
            imgSOX2Mask = imgSOX2binary - EdgeBinary;
            imgSOX2Mask = uint8(imgSOX2Mask);
            imgSOX2new = imgSOX2binary - imgSOX2Mask;
            % Now we need also to further remove artifacts.
            % Delete very small particles:
            imgSOX2new = bwareaopen(imgSOX2new,3);
            %Now we check for non round particles. This can be done with the
            %bounding box and then the ratio between the minor and major axis!
            % First we label the particles
            imgSOX2newLabel = bwlabel(imgSOX2new);
            SynapseN = max(max(imgSOX2newLabel));
            imgSOX2new = zeros(size(EdgeBinary));
            for(y=1:SynapseN)
                SynapseOne = imgSOX2newLabel;
                ii = SynapseOne == y;
                Add =  SynapseN*10;
                SynapseOne(ii) = Add;
                %Sub value
                Sub = Add -5;
                SynapseOne = SynapseOne -Sub;
                SynapseOne = uint8(SynapseOne);
                SynapseOne = logical(SynapseOne);
                MajorSynapse = regionprops(SynapseOne,'MajorAxis');
                MinorSynapse = regionprops(SynapseOne,'MinorAxis');
                MajorSynapseAxis = MajorSynapse(1).MajorAxisLength(1);
                MinorSynapseAxis = MinorSynapse(1).MinorAxisLength(1);
                RatioAxis = MajorSynapseAxis/MinorSynapseAxis;
                Vector(y,1) = RatioAxis;
                if(RatioAxis<3)
                    imgSOX2new = imgSOX2new + SynapseOne;
                end    


            end    
            imgSOX2new = uint8(imgSOX2new*255);


            if(n<13)
                if(n<10)
                    a{1} = 'A0';
                else
                    a{1} = 'A';
                end
            a{2} = int2str(n);
            a{3} = 'space';
            end
            if(n>12&&n<25)
            p = n - 12;
            if(p<10)
                a{1} = 'B0';
            else
                a{1} = 'B';
            end
            a{2} = int2str(p);
            a{3} = 'space';
            end
            if(n>24&&n<37)

            p = n - 24;
            if(p<10)
                a{1} = 'C0';
            else
                a{1} = 'C';
            end
            a{2} = int2str(p);
            a{3} = 'space';
            end
            if(n>36&&n<49)

            p = n - 36;
            if(p<10)
                a{1} = 'D0';
            else
                a{1} = 'D';
            end
            a{2} = int2str(p);
            a{3} = 'space';
            end
            if(n>48&&n<61)

            p = n - 48;
            if(p<10)
                a{1} = 'E0';
            else
                a{1} = 'E';
            end
            a{2} = int2str(p);
            a{3} = 'space';
            end
            if(n>60&&n<73)

            p = n - 60;
            if(p<10)
                a{1} = 'F0';
            else
                a{1} = 'F';
            end
            a{2} = int2str(p);
            a{3} = 'space';
            end
            if(n>72&&n<85)

            p = n - 72;
            if(p<10)
                a{1} = 'G0';
            else
                a{1} = 'G';
            end
            a{2} = int2str(p);
            a{3} = 'space';
            end
            if(n>84&&n<96)
            p = n - 84;
            a{2} = int2str(p);
            a{3} = 'space';
            if(p<10)
                a{1} = 'H0';
            else
                a{1} = 'H';
            end
            end
            %Here we need to add the empty rim:
            if(ind5==0)
                EmptyImage = zeros(r,c);
                EmptyImage = uint8(EmptyImage);
                %Determine size of current Image:
                [r_s c_s] = size(nucleusImg);
                EmptyImage(1:r_s,1:c_s) = nucleusImg;
                nucleusImg = EmptyImage;
                EmptyImage(1:r_s,1:c_s) = imgNeurite;
                imgNeurite = EmptyImage;
                EmptyImage(1:r_s,1:c_s) = EdgeBinary;
                EdgeBinary = EmptyImage;
                EmptyImage(1:r_s,1:c_s) = imgSOX2new;
                imgSOX2new = EmptyImage;
                EmptyImage(1:r_s,1:c_s) = imgSOX2;
                imgSOX2 = EmptyImage;
            end
            newfile = [a{1} a{2} a{3}];
            newFile1 = regexprep(newfile, 'space','NucleusSmall.png','ignorecase');
            newFile2 = regexprep(newfile, 'space','NeuriteBig.png','ignorecase'); 
            newFile3 = regexprep(newfile, 'space','NeuriteSmall.png','ignorecase');
            newFile4 = regexprep(newfile, 'space','Binary.png','ignorecase');
            newFile5 = regexprep(newfile, 'space','NucleusBigWatershed.png','ignorecase');
            newFile7 = regexprep(newfile, 'space','AstroBig.png','ignorecase');
            newFile8 = regexprep(newfile, 'space','AstroSmall.png','ignorecase');
            newFile9 = regexprep(newfile, 'space','NucleusBig.png','ignorecase');
            newFile10 = regexprep(newfile, 'space','BinaryHard.png','ignorecase');
            newFile11 = regexprep(newfile, 'space','OligoBig.png','ignorecase');
            newFile12 = regexprep(newfile, 'space','OligoSmall.png','ignorecase');
            imwrite(nucleusPic,[foldername1 '/' newFile9]);
            imwrite(nucleusPic,[foldername1 '/' newFile1]);
            imwrite(imgNeurite,[foldername1 '/' newFile2]);
            imwrite(imgNeurite,[foldername1 '/' newFile3]);
            imwrite(EdgeBinary,[foldername1 '/' newFile4]);
            imwrite(imgSOX2new,[foldername1 '/' newFile7]);
            imwrite(imgSOX2new,[foldername1 '/' newFile8]);
            imwrite(imgSOX2,[foldername1 '/' newFile11]);
            imwrite(imgSOX2,[foldername1 '/' newFile12]);
            wellname = (newFile1(1:3));
            wellList{n} = wellname;
            NucleusM = sparse(zeros(size(imgNeurite)));
            csvHandler.CellPosMatrix(wellname) = NucleusM;
            neuronHandler.NeuronPositionsEdgeFill(wellname) = sparse(zeros(size(nucleusPic)));
            ScalingFactor = 1;
            clear newfile newfile1 newfile2 newfile3 newfile4 newfile5 newfile6 newfile7 newfile8 
        end
        end
end

% Comments: It would be good to include a suitable filter for bubble rims
