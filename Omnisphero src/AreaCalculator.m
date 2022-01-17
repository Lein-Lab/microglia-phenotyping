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

%This function returns the area of white pixles in a binary image und is
%used to identfy changes made in Cortical Neurons:
foldername = uigetdir;
allfiles = dir(foldername);
n=0;
summary = 0;
for(i=1:numel(allfiles))
    %Check if file ends with Binary.png. This will ignore the original stack,
    %which is in tiff.
    ind = strfind([foldername '/' allfiles(i).name],'NucleusBigWatershed.png');
    % Check of image i is a png or not
    if(numel(ind) > 0)%
        n = n+1;
        imgNeurite = imread([foldername '/' allfiles(i).name]); 
        imgNeurite = logical(imgNeurite);
        Area = sum(sum(imgNeurite));
        summary(n,1) = Area;
        summary = summary
    end
    
end    