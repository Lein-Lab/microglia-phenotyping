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
function [wellList, foldername, ScalingFactor, csvHandler, neuronHandler,  Th] = SphereOutgrowth (csvHandler, neuronHandler) 
% This is a quick script to convert 16-bit images of cortical neurons into
% binary image, to assess their center point and to generate a
% 'NucleusImage'.
n = 0;
SphereCoreSize = 3000;
SE= strel('disk',50);
N = 10;
BallFactor = 1;
BoxFactor = 64;
%Get Image folder
foldername = uigetdir;
mkdir(foldername,'ConvertedCellomics');
foldername1 = [foldername '/ConvertedCellomics'];
%Iterate over all files in this folder
allfiles = dir(foldername);
wellList = cell(0);
for(i=1:numel(allfiles))
    %Check if file ends with .tif or .tiff
    ind = strfind([foldername '/' allfiles(i).name],'.jpg');
    if(numel(ind) > 0)%
         n = n + 1;
        %Load RGB image
        imgRGB = imread([foldername '/' allfiles(i).name]);
        %Convert to gray scale:
        imgNeurite = imgRGB(:,:,1); 
        
        %Fine sphere core:
        %Use low background to over sphere core and then try to erode
        %processes
        backgroundball = imopen(imgNeurite,SE);
        backgroundballimg = imgNeurite - backgroundball;
        
        
        
        %[lehisto x] = imhist(imgNeurite);
        %level = triangle_th(lehisto,256);
        level = isodata(imgNeurite);
        imgBinary = im2bw(imgNeurite,level);
        % Close small gaps originating from thresholding in the sphere
        % core!
        HolesSphereCore = imfill(imgBinary,'holes') - imgBinary;
        HolesSphereCore = bwareaopen(HolesSphereCore,50);
        imgBinary = imfill(imgBinary,'holes') - HolesSphereCore;
        % Create sort of a density image, by subdividing the image into a
        % grid and only keep grids with a certain number of white pixels:
        
        
        [r c] = size(imgNeurite);
        box = zeros(size(imgNeurite)/BoxFactor);
        boxNeighborRight = zeros(size(imgNeurite)/BoxFactor);
        boxNeighborLeft = zeros(size(imgNeurite)/BoxFactor);
        [q d] = size(box);
        deletebox = zeros(size(imgNeurite));
        nBoxesX = c/d;
        nBoxesY = r/q;
        f = 0;
        t = 0;
        for(g=1:nBoxesX)
          x1 = 1 +(f*d);
          x2 = d +(f*d);
          if(f<(floor(nBoxesX)-1))
          x3 = 1 +((f+1)*d);
          x4 = d +((f+1)*d);
          else
          x3 = 1 +((f-1)*d);
          x4 = d +((f-1)*d);    
          end
          if(f>0)
          x5 = 1 +((f-1)*d);
          x6 = d +((f-1)*d);
          else
          x5 = 1 +(f*d);
          x6 = d +(f*d);    
          end
                    
            for(o=1:nBoxesY)
             y1 = 1 +(t*q);
             y2 = q +(t*q);
              if(t<(floor(nBoxesY)-1))
              y3 = 1 +((t+1)*q);
              y4 = q +((t+1)*q);
              else
              y3 = 1 +((t-1)*q);
              y4 = q +((t-1)*q);    
              end
              if(t>0)
              y5 = 1 +((t-1)*q);
              y6 = q +((t-1)*q);
              else
              y5 = 1 +(t*q);
              y6 = q +(t*q);    
              end
             t = t + 1;
             box = imgBinary(y1:y2,x1:x2);
             boxNeighborRight = imgBinary(y3:y4,x3:x4);
             boxNeighborLeft = imgBinary(y5:y6,x5:x6);
                if(sum(sum(box))/(q*d) < 0.95)
                   if(sum(sum(boxNeighborRight))/(q*d) < 0.95 && sum(sum(boxNeighborLeft))/(q*d) < 0.95)
                       imgBinary(y1:y2,x1:x2) = deletebox(y1:y2,x1:x2);
                   end
                end 
                
            end 
            f = f + 1;
            t = 0;
        end
       
        imgBinary = bwareaopen(imgBinary,3000);
        imgBinaryThin = bwmorph(imgBinary,'thin',10);
        imgBinaryThin = bwmorph(imgBinaryThin,'open');
        sP = regionprops(imgBinaryThin,'area');
        LsP = length(sP);
        for(u=1:LsP)
            Area = sP(u).Area(1);
            V_Area(u) = Area;
        end
        maxParticle = max(V_Area);
        Particle = bwareaopen(imgBinaryThin,maxParticle-1);
        imgBinaryThin = imgBinaryThin - Particle;
        imgBinaryThin = uint8(imgBinaryThin);
        imgBinaryThin = logical(imgBinaryThin);
        imgBinaryThin = bwmorph(imgBinaryThin,'thicken',10);
        imgBinary = imgBinary - imgBinaryThin;
        imgBinary = bwareaopen(imgBinary,maxParticle-1);
        %Repeat to get rid of short branches
        imgBinaryThin = bwmorph(imgBinary,'thin',5);
        imgBinaryThin = bwmorph(imgBinaryThin,'open');
        sP = regionprops(imgBinaryThin,'area');
        LsP = length(sP);
        for(u=1:LsP)
            Area = sP(u).Area(1);
            V_Area(u) = Area;
        end
        maxParticle = max(V_Area);
        Particle = bwareaopen(imgBinaryThin,maxParticle-1);
        imgBinaryThin = imgBinaryThin - Particle;
        imgBinaryThin = uint8(imgBinaryThin);
        imgBinaryThin = logical(imgBinaryThin);
        imgBinaryThin = bwmorph(imgBinaryThin,'thicken',5);
        imgBinary = imgBinary - imgBinaryThin;
        imgBinary = bwareaopen(imgBinary,maxParticle-1);
        
        
        
        
        
        
        
        imgBinary= bwmorph(imgBinary,'thin');
        imgBinary = bwmorph(imgBinary,'majority');
        imgBinary = bwareaopen(imgBinary,SphereCoreSize);
        imgBinary = bwmorph(imgBinary,'thicken');
        SphereCore = imgBinary;
        %SphereCorePre = bwareaopen(imgBinary,SphereCoreSize);
        %SphereCorePre = bwmorph(SphereCorePre,'open');
        %SphereCorePre = bwareaopen(SphereCorePre,SphereCoreSize);
        %SphereCorePre = bwmorph(SphereCorePre,'thin',3);
        %SphereCorePre = bwmorph(SphereCorePre,'open');
        %SphereCorePre = bwareaopen(SphereCorePre,SphereCoreSize);
        %SphereCorePre = bwmorph(SphereCorePre,'thin',2);
        %SphereCorePre = bwmorph(SphereCorePre,'open');
        %SphereCore = bwareaopen(SphereCorePre,SphereCoreSize);
        %SphereCore = bwmorph(SphereCore,'thicken',20);
        
        %Neurite processis originating from the sphere core:
        
        
        % Cut out original image to only obtain areas within the mask:
            
        backgroundNeurite = imopen(uint8(imgNeurite),strel('disk',N));
        NeuriteBallBackground = imgNeurite - backgroundNeurite*BallFactor;
        [lehisto x] = imhist(NeuriteBallBackground);
        %Calculates new background for img1 using the triangle method
        level = triangle_th(lehisto,256);
        NeuriteBinary = im2bw(NeuriteBallBackground,level);
        NeuriteBinary = bwareaopen(NeuriteBinary,200);
        NeuriteOnly = NeuriteBinary - SphereCore;
        NeuriteOnly = uint8(NeuriteOnly);
        
        s = regionprops(SphereCore, 'Centroid');
        sdimensions = size(s);
        NucleusM = sparse(zeros(size(SphereCore)));
        if(sdimensions > 0)
            for(j=1:sdimensions)
            x_centroid = s(j).Centroid(1);
            y_centroid = s(j).Centroid(2);
            NucleusM(uint16(y_centroid),uint16(x_centroid)) = 1;
            end
        end
        
        
        
       SphereCore = uint8(SphereCore)*100; 
        
       
                
        
        if(n<6)
        a{1} = 'A0';
        a{2} = int2str(n);
        a{3} = 'space';
        end
        if(n>5&&n<11)
        a{1} = 'B0';
        p = n - 5;
        a{2} = int2str(p);
        a{3} = 'space';
        end
        if(n>10&&n<16)
        a{1} = 'C0';
        p = n - 10;
        a{2} = int2str(p);
        a{3} = 'space';
        end
        if(n>15&&n<21)
        a{1} = 'D0';
        p = n - 15;
        a{2} = int2str(p);
        a{3} = 'space';
        end
        if(n>20&&n<26)
        a{1} = 'E0';
        p = n - 20;
        a{2} = int2str(p);
        a{3} = 'space';
        end
        if(n>25&&n<31)
        a{1} = 'F0';
        p = n - 25;
        a{2} = int2str(p);
        a{3} = 'space';
        end
        if(n>30&&n<36)
        a{1} = 'G0';
        p = n - 30;
        a{2} = int2str(p);
        a{3} = 'space';
        end
        if(n>35&&n<41)
        p = n - 35;
        a{2} = int2str(p);
        a{3} = 'space';
        a{1} = 'H0';
        end
        
        newfile = [a{1} a{2} a{3}];
        newFile1 = regexprep(newfile, 'space','NucleusSmall.png','ignorecase');
        newFile2 = regexprep(newfile, 'space','NeuriteBig.png','ignorecase'); 
        newFile3 = regexprep(newfile, 'space','NeuriteSmall.png','ignorecase');
        newFile4 = regexprep(newfile, 'space','Binary.png','ignorecase');
        newFile5 = regexprep(newfile, 'space','NucleusBigWatershed.png','ignorecase');
        newFile7 = regexprep(newfile, 'space','OligoBig.png','ignorecase');
        newFile8 = regexprep(newfile, 'space','OligoSmall.png','ignorecase');
        newFile9 = regexprep(newfile, 'space','NucleusBig.png','ignorecase');
        newFile10 = regexprep(newfile, 'space','BinaryHard.png','ignorecase');
        imwrite(SphereCore,[foldername1 '/' newFile9]);
        imwrite(SphereCore,[foldername1 '/' newFile1]);
        imwrite(imgNeurite,[foldername1 '/' newFile2]);
        imwrite(imgNeurite,[foldername1 '/' newFile3]);
        %imwrite(binaryimgNeurite,[foldername1 '/' newFile4]);
        %Noch ImageJ fixen!
        imwrite(NeuriteOnly,[foldername1 '/' newFile4]);
        imwrite(SphereCore,[foldername1 '/' newFile5]);
        
        wellname = (newFile1(1:3));
        wellList{n} = wellname;
        csvHandler.CellPosMatrix(wellname) = NucleusM;
        neuronHandler.NeuronPositionsEdgeFill(wellname) = sparse(zeros(size(SphereCore)));
        ScalingFactor = 1;
        clear newfile newfile1 newfile2 newfile3 newfile4 newfile5 newfile6 newfile7 newfile8 
    end
    
end

% Comments: It would be good to include a suitable filter for bubble rims
