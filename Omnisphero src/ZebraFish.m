
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
function [wellList, foldername, ScalingFactor, csvHandler] = ZebraFish (csvHandler) 
% This is a quick script to convert 16-bit images of cortical neurons into
% binary image, to assess their center point and to generate a
% 'NucleusImage'.
n = 0;
A = 20;
N =10;
%Get Image folder
ScalingFactor = 1;
foldername = uigetdir;
foldername1 = [foldername '/ConvertedCellomics'];
mkdir(foldername1);
wellList = cell(0);
%Iterate over all files in this folder
allfiles = dir(foldername);
for(i=1:numel(allfiles))
    
    %Check if file ends with .tif or .tiff
    ind = strfind([foldername '/' allfiles(i).name],'w1');
    if(numel(ind) > 0)
        n = n+1
        %Load image and save it as .tiff
        img = imread([foldername '/' allfiles(i).name]);
        newFile1 = regexprep(allfiles(i).name, 'w1','w2','ignorecase');
        img2 = imread([foldername '/' newFile1]); 
        %Convert image to 8 bit
        img = double(img)./65535;
        img2 = double(img2)./65535;
        img2 = uint8(img2.*255);
        img = uint8(img.*255);
        backgroundball =  imopen(uint8(img),strel('ball',N,N));
        backgroundball2 =  imopen(uint8(img2),strel('ball',N,N));
        imgballbackground = img - backgroundball;
        imgballbackground2 = img2 - backgroundball2;
        imgc = img + (imgballbackground);
        imgc2 = img2 + (imgballbackground2);
        imgc = uint8(imgc);
        imgc2 = uint8(imgc2);
        backgroundball =  imopen(uint8(imgc),strel('ball',N,N));
        backgroundball2 =  imopen(uint8(imgc2),strel('ball',N,N));
        imgcballbackground = imgc - backgroundball;
        imgcballbackground2 = imgc2 - backgroundball2;
        
        %Calculates the histogram of img1
        [lehisto x] = imhist(imgcballbackground);
        %Calculates new background for img1 using the triangle method
        level = triangle_th(lehisto,256);
        %Creates a binary image of img1 using level
        binaryimg = im2bw(imgcballbackground,level*10);
        binaryimg = bwareaopen(binaryimg,A);
        s = regionprops(binaryimg, 'Centroid');
        sdimensions = size(s);
        NucleusM = sparse(zeros(size(img)));
        if(sdimensions > 0)
            for(j=1:sdimensions)
            x_centroid = s(j).Centroid(1);
            y_centroid = s(j).Centroid(2);
            NucleusM(uint16(y_centroid),uint16(x_centroid)) = 1;
            end
        end
        
        
        %Calculates the histogram of img1
        [lehisto x] = imhist(imgcballbackground2);
        %Calculates new background for img1 using the triangle method
        level2 = triangle_th(lehisto,256);
        %Creates a binary image of img1 using level
        binaryimg2 = im2bw(imgcballbackground2,level2*10);
        
        
        
        % Excludes all areas in 'imgnucleusbinary' smaller than 250000
        %Replace .tif(f) with .png in filename
        allfiles(i).name = strrep(allfiles(i).name,'.tiff','.png');
        allfiles(i).name = strrep(allfiles(i).name,'.tif','.png');
        if(n<10)
            newFileWell{1} = 'A0';
            newFileWell{2} = int2str(n);
            newFileWell{3} = 'space';
        else
            newFileWell{1} = 'A';
            newFileWell{2} = int2str(n);
            newFileWell{3} = 'space';
        end   
        newFileWell = [newFileWell{1} newFileWell{2} newFileWell{3}];
        newFile1 = regexprep(newFileWell, 'space','NucleusSmall.png','ignorecase'); 
        newFile2 = regexprep(newFileWell, 'space','NucleusBig.png','ignorecase');
        newFile3 = regexprep(newFileWell, 'space','NeuriteSmall.png','ignorecase');
        newFile4 = regexprep(newFileWell, 'space','NeuriteBig.png','ignorecase');
        newFile5 = regexprep(newFileWell, 'space','Binary.png','ignorecase');
        newFile6 = regexprep(newFileWell, 'space','Skeleton.png','ignorecase');
        imwrite(img,[foldername1 '/' newFile1]);
        imwrite(img,[foldername1 '/' newFile2]);
        imwrite(img2,[foldername1 '/' newFile3]);
        imwrite(img2,[foldername1 '/' newFile4]);
        imwrite(binaryimg,[foldername1 '/' newFile5]);
        imwrite(binaryimg2,[foldername1 '/' newFile6]);
        wellname = (newFile1(1:3));
        wellList{n} = wellname;
        csvHandler.CellPosMatrix(wellname) = NucleusM;
        ScalingFactor = 1;
    end
    
    
end
