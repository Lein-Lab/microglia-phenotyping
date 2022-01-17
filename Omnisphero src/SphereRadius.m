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
function [summary] = SphereRadius (wellList, foldername, csvHandler) 
% This program is designed to quantify GFAP+ cells in brain slices and to
% determine their area in % of the image
n = 0;
foldername = [foldername '/ConvertedCellomics'];
SphereCoreSizeFrag = 2000;
SphereCore = 20000;
%Iterate over all files in this folder
allfiles = dir(foldername);


for(i=1:numel(allfiles))
    %Check if file ends with .tif or .tiff
    ind = strfind([foldername '/' allfiles(i).name],'NucleusBigWatershed.png');
    if(numel(ind) > 0)
        n = n + 1;
        %Load image and save it as .tiff
        img = imread([foldername '/' allfiles(i).name]);
        img = imclearborder(img);
        img = bwareaopen(img,SphereCoreSizeFrag);
        img = bwmorph(img,'bridge');
        img = bwareaopen(img,SphereCore);
        img = logical(img);
        
        
        AreaCores = regionprops(img,'Area');
        NAreas = size(AreaCores);
        NAreas = NAreas(1);
        VArea = 0;
        for(e=1:NAreas)
            AreaN = AreaCores(e).Area(1);
            VArea(e,1) = AreaN;
        end
        if(max(VArea>0))
        Cut = (max(VArea))-1;
        else
        Cut = 1;
        end
        img = bwareaopen(img,Cut);
        AreaCore = sum(sum(img));
        if(AreaCore >0)
        Major = regionprops(img,'MajorAxis');
        Minor = regionprops(img,'MinorAxis');
        
        MajorAxis = Major(1).MajorAxisLength(1);
        MinorAxis = Minor(1).MinorAxisLength(1);
                  
        
        SphereRadius = ((MajorAxis/2)+(MinorAxis/2))/2;
        else
        SphereRadius = 0;
        end
        
        wellname = allfiles(i).name(1:3);
        
        summary{1,1} = 'WellName';
        summary{1,2} = 'Sphere-Core Area';
        summary{1,3} = 'Sphere-Core Radius';
        
        p = n +1;
        
        summary{p,1} = wellname;
        summary{p,2} = AreaCore;
        summary{p,3} = SphereRadius;
        
        summary = summary
        
        filename = 'File Radius';
        xlswrite(filename,summary);
        
        newFile = regexprep(allfiles(i).name, 'NucleusBig','SphereCore','ignorecase');
        imwrite(img,[foldername '/' newFile]);
       
        
    end
    
    
end