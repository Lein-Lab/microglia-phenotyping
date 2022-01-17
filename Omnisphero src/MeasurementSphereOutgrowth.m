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

function [summary] = MeasurementSphereOutgrowth (foldername, csvHandler)
foldername = [foldername '/' 'ConvertedCellomics']
allfiles = dir(foldername);
wellList = cell(0);
n=0;

for(i=1:numel(allfiles))
    %Check if file ends with .tif or .tiff
    ind = strfind([foldername '/ConvertedCellomics' allfiles(i).name],'Binary');
    if(numel(ind) > 0)%
         n = n + 1;
         Binary = imread([foldername '/' allfiles(i).name]);
         newfile = regexprep(allfiles(i).name,'Binary','NucleusBig','ignorecase');
         Nuclei = imread([foldername '/' newfile]);
         
         Nuclei = logical(Nuclei);
         Binary = logical(Binary);
         
         EdgeNuclei = edge(Nuclei,'sobel');
         EdgeNuclei = bwareaopen(EdgeNuclei,10);
         
         gaus = fspecial('gaussian',4);
         Neurite = imfilter(Binary,gaus);
         NeuriteSkel = bwmorph(Neurite,'thin',Inf);
         
         
         %Extracted variables for CSV
         % 1) Nuclei features
         NucleusArea = sum(sum(Nuclei));
         NucleusEdge = sum(sum(EdgeNuclei));
         
         % 2) Neuron features
         NeuriteArea = sum(sum(Neurite));
         NeuriteLength = sum(sum(NeuriteSkel));
         
         a{1} = 'A0';
         a{2} = int2str(n);
         a{3} = 'space';
         newfile = [a{1} a{2} a{3}];
         newFile1 = regexprep(newfile, 'space','NucleusSmall.png','ignorecase');
         wellname = (newFile1(1:3));
         
        EdgeNuclei = uint8(EdgeNuclei); 
        Binary = logical(Binary);
        NeuriteArea = sum(sum(Binary));
        
        
        
        newFile1 = regexprep(allfiles(i).name,'Binary','Skeleton','ignorecase');
        imwrite(NeuriteSkel,[foldername '/' newFile1]);
        newFile2 = regexprep(allfiles(i).name,'Binary','AstroBig','ignorecase');
        newFile3 = regexprep(allfiles(i).name,'Binary','AstroSmall','ignorecase');
        imwrite(EdgeNuclei,[foldername '/' newFile2]);
        imwrite(EdgeNuclei,[foldername '/' newFile3]);
        
        
        %Create Matrix
        summary{1,1} =  'Wellname';
        summary{1,2} =  'CoreArea';
        summary{1,3} =  'CoreRim';
        summary{1,4} =  'NeuriteArea';
        summary{1,5} =  'NeuriteTotalLength';
        
        m=n+1;
        summary{m,1} =  wellname;
        summary{m,2} =  NucleusArea;
        summary{m,3} =  NucleusEdge;
        summary{m,4} =  NeuriteArea;
        summary{m,5} =  NeuriteLength;
        %summary = summary
        filename = 'Results.xlsx';
        xlswrite(filename,summary);
    end
end
         