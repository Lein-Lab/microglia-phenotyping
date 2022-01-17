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

% This fucntion is counting the number of neurosphere cores on the
% Watershed image and should be used to set the proper size for sphere core
% area N:
N_1 = 75000;
N_2 = 180000;
foldername = uigetdir;
allfiles = dir(foldername);
n = 0;
summary = cell(8,2);
for(i=1:numel(allfiles))
    ind = strfind([foldername '/' allfiles(i).name],'NucleusBigWatershed.png');
        if(numel(ind) > 0)
            n=n+1;
            img = imread([foldername '/' allfiles(i).name]);
            imgsort = bwareaopen(img,N_1);
            imgsort = bwmorph(imgsort, 'bridge');
            imgsort = bwareaopen(imgsort,N_2);
            imgsort = imclearborder(imgsort);
            s = regionprops(imgsort,'centroid');
            s_L = length(s);
            wellname = allfiles(i).name(1:3);
            summary{n,1} = wellname;
            summary(n,2) = {s_L};
        end
end        

