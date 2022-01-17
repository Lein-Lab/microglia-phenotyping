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

%function [wellList, foldername, ScalingFactor, csvHandler] = GFAPonly (csvHandler)
% This function is designed to measure the intensity of a background
% corrected Nestin positive image

%Background = 13 ;

foldername = uigetdir;
foldername1 = uigetdir;
n=0;
wellList = cell(0);
%Iterate over all files in this folder
allfiles = dir(foldername);
summary = cell(25,2);
for(i=1:numel(allfiles))
            %Check if file ends with .tif or .tiff
    ind = strfind([foldername '/' allfiles(i).name],'AstroBig');
    if(numel(ind) > 0)
        n = n + 1
        %Load image and save it as .tiff
        img = imread([foldername '/' allfiles(i).name]);
        ind_2 = strfind([foldername '/' allfiles(i).name], '05AstroBig');
            if(numel(ind_2) > 0)
               Background = mean(mean(img));
            end   
        %Substract fixed background
        img = img - Background;
        Intensity = sum(sum(img));
        wellname = allfiles(i).name(1:3);
        summary{n,1} = wellname;
        summary(n,2) = {Intensity};
%        csvHandler.CellPosMatrix(wellname) = NucleusM;
%        ScalingFactor = 1;
        imwrite(img,[foldername1 '/' allfiles(i).name]);
    end
    
    
end