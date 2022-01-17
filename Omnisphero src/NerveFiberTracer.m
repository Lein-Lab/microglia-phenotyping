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

% This script is designed to remove autofluorescence background from
% tissue, while maintaining nerve fibers. If possible also gaps should be
% closed

n = 0;
foldername = uigetdir;
mkdir(foldername,'ConvertedCellomics');
foldername1 = [foldername '/ConvertedCellomics'];
mkdir(foldername1,'GreyScale');
mkdir(foldername1,'Binary');
mkdir(foldername1,'Fibers');
foldername1 = [foldername '/ConvertedCellomics'];
foldername2 = [foldername1 '/GreyScale'];
foldername3 = [foldername1 '/Binary'];
foldername4 = [foldername1 '/Fibers'];
%Iterate over all files in this folder
allfiles = dir(foldername);

for(i=1:numel(allfiles))
    %Check if file ends with .tif or .tiff
    ind = strfind([foldername '/' allfiles(i).name],'.tif');
    if(numel(ind) > 0)%
         n = n + 1;
         img = imread([foldername '/' allfiles(i).name]);
         
         % In a first step we will perform a median filter in order to
         % smooth the initial image!
         imgPost = medfilt2(img, [3 3]);
         
         % Now we will perform the rolling ball to remove some background:
         SE = strel('ball',5,5);
         imgPostBall = imopen(imgPost,SE);
         imgPostBallCorr = imgPost - imgPostBall;
         imgPostBallCorr = double(imgPostBallCorr)./4095;
         imgPostBallCorr = uint8(imgPostBallCorr.*255);
         M = median(imgPostBallCorr(imgPostBallCorr>0));
         imgPostBallCorr = imgPostBallCorr - M;
         imgPostBallCorr = uint8(imgPostBallCorr);
         level = graythresh(imgPostBallCorr);
         if(level<0.1)
             level=0.1;
         end    
         imgPostBallBinary = im2bw(imgPostBallCorr,level);
         imgPostBallBinary = bwareaopen(imgPostBallBinary,100);
         imgPostBallBinary = bwmorph(imgPostBallBinary,'bridge');
         
         
         % Now we will apply a roundness filter to further narrow down the
         % fibers!
         imgPostBallBinaryLabel = bwlabel(imgPostBallBinary);
         LabelInd = max(max(imgPostBallBinaryLabel));
         Fibers = zeros(size(imgPostBallBinaryLabel));
         Fibers = logical(Fibers);
         for(f=1:LabelInd)
             BinaryAreaOne = imgPostBallBinaryLabel;
             ii = BinaryAreaOne == f;
             BinaryAreaOne(ii) = 1000000;
             BinaryAreaOne = BinaryAreaOne - (1000000-5);
             BinaryAreaOne = uint8(BinaryAreaOne);
             BinaryAreaOne = logical(BinaryAreaOne);
             Major = regionprops(BinaryAreaOne,'MajorAxisLength');
             Minor = regionprops(BinaryAreaOne,'MinorAxisLength');
             sLength = Major(1).MajorAxisLength(1);
             sHeight = Minor(1).MinorAxisLength(1);
             if(sLength>sHeight)
                 Roundness = sLength/sHeight;
             else
                 Roundness = sHeight/sLength;
             end    
             if(Roundness>2.5)
                 Fibers = Fibers + BinaryAreaOne;
             end    
            
         end   
           
         %We have to convert the Tifs to 8-bit!
         imgPostBallBinary = uint8(imgPostBallBinary*255);
         Fibers = uint8(Fibers*255);
         %Now we have to write the images in seperate folders!    
         imwrite(imgPostBallCorr,[foldername2 '/' allfiles(i).name]);    
         imwrite(imgPostBallBinary,[foldername3 '/' allfiles(i).name]); 
         imwrite(Fibers,[foldername4 '/' allfiles(i).name]); 
     
         
         
         
         
    end
end

%Now we perfrom a max-z-Projection, to get rid of all particles just
%occuring in in on plain! 
allfiles1 = dir(foldername3);
m = 0;
for(l=1:numel(allfiles1))
    %Check if file ends with .tif or .tiff
    ind = strfind([foldername3 '/' allfiles1(l).name],'.tif');
    if(numel(ind) > 0)%
         m = m + 1;
         img = imread([foldername3 '/' allfiles(l).name]);
         blub = 2;
    end
    

end
