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

%This function is designed to only load the different channels into the
%Omnisphero GUI to make the area drawings!
function [wellList, foldername, ScalingFactor, csvHandler,  Th] = GangliaImport (csvHandler)
n = 0;
foldername = uigetdir;
% Creates directory within the parent folder
allfiles = dir(foldername);
% Creates new folder within parent folder
mkdir(foldername,'ConvertedCellomics');
foldername1 = [foldername '/ConvertedCellomics'];
% Creeate directory in parent foldername
allfiles = dir(foldername);
t=0;
wellList = cell(0);
for(i=1:numel(allfiles))
    %Check if file ends with .png. This will ignore the original stack,
    %which is in tiff.
    ind = strfind([foldername '/' allfiles(i).name],'_w1');
   if(numel(ind) > 0)
     n = n + 1;   
     imgNucleus = imread([foldername '/' allfiles(i).name]);
     newFileNeuron = regexprep(allfiles(i).name, '_w1','_w2','ignorecase'); 
     imgNeuron = imread([foldername '/' newFileNeuron]);
     newFileChannel3 = regexprep(allfiles(i).name, '_w1','_w3','ignorecase');
     imgChannel3 = imread([foldername '/' newFileChannel3]);
     imgChannel3 = double(imgChannel3)./16383;
     imgChannel3 = uint8(imgChannel3.*255);
     imgNeuron = double(imgNeuron)./16383;
     imgNeuron = uint8(imgNeuron.*255);
%      imgNucleus = double(imgNucleus)./16383;
%      imgNucleus = uint8(imgNucleus.*255);
     
     
     
     %Generate imgNucleusWatershed and get Nuclei positions!
     [lehisto x] = imhist(imgNucleus);
    level = triangle_th(lehisto,256);
    NucleusBinary = im2bw(imgNucleus,level);
    %Remove huge Artefacts:
    Artefacts = bwareaopen(NucleusBinary,10000);
    NucleusBinary = NucleusBinary - Artefacts;
    %Fill small holes:
    Holes = imfill(NucleusBinary,'holes');
    HolesBig = bwareaopen(Holes,50);
    Holes = Holes - HolesBig;
    NucleusBinary = NucleusBinary + Holes;
    %From here:https://blogs.mathworks.com/steve/2006/06/02/cell-segmentation/
    PerimeterNuclei = bwperim(NucleusBinary);
    mask_em = imextendedmax(imgNucleus, 300);
    Nucleus_c = imcomplement(imgNucleus);
    I_mod = imimposemin(Nucleus_c, ~NucleusBinary | mask_em);
    L = watershed(I_mod);
    L = uint8(L);
    L = logical(L);
    L = ~L;
    NucleusBinary = NucleusBinary - L;
    Holes = imfill(NucleusBinary,'holes');
    HolesBig = bwareaopen(Holes,50);
    Holes = Holes - HolesBig;
    NucleusBinary = NucleusBinary + Holes;
    NucleusBinary = bwareaopen(NucleusBinary,50);
    sR = regionprops(NucleusBinary, 'Centroid');
    sRdimensions = size(sR);
    NucleusM = uint16(zeros(size(NucleusBinary)));
    for(k=1:sRdimensions)
            x_centroid = sR(k).Centroid(1);
            y_centroid = sR(k).Centroid(2);
            NucleusM(uint16(y_centroid),uint16(x_centroid)) = 1;
    end
        
     imgNucleusWatershed =  NucleusBinary;
     
     
     %This part is for saving the images
        %Wells are named in a linear order (max 99)
%         if(n>9)
%         a{1} = 'A';
%         else
%         a{1} = 'A0';
%         end
        WellName = allfiles(i).name(end-9:end-7);
        a{1} = WellName(3);
        a{2} ='_';
        a{3} = WellName(1);
        a{4} = 'space';
        newfile = [a{1} a{2} a{3} a{4}];
        imgNucleusSmall = imresize(imgNucleus,0.1);
        imgNeuronSmall = imresize(imgNeuron,0.1);
        imgChannel3Small = imresize(imgChannel3,0.1);
        
       % Names of images
        newFile1 = regexprep(newfile, 'space','NucleusSmall.png','ignorecase');
        newFile2 = regexprep(newfile, 'space','NeuriteBig.png','ignorecase'); 
        newFile3 = regexprep(newfile, 'space','NeuriteSmall.png','ignorecase');
        newFile5 = regexprep(newfile, 'space','NucleusBigWatershed.png','ignorecase');
        newFile7 = regexprep(newfile, 'space','OligoBig.png','ignorecase');
        newFile8 = regexprep(newfile, 'space','OligoSmall.png','ignorecase');
        newFile9 = regexprep(newfile, 'space','NucleusBig.png','ignorecase');
        % Saving of Images
        imwrite(imgNucleus,[foldername1 '/' newFile9]);
        imwrite(imgNucleusSmall,[foldername1 '/' newFile1]);
        imwrite(imgNeuron,[foldername1 '/' newFile2]);
        imwrite(imgNeuronSmall,[foldername1 '/' newFile3]);
        imwrite(imgNucleusWatershed,[foldername1 '/' newFile5]);
        imwrite(imgChannel3,[foldername1 '/' newFile7]);
        imwrite(imgChannel3Small,[foldername1 '/' newFile8]);
        % Well name gets name from filename
        wellname = (newFile1(1:3));
        % Welllist is filled with well names
        wellList{n} = wellname;
        % Cell centroids are send to csvHandler for Omnisphero GUI
        csvHandler.CellPosMatrix(wellname) = NucleusM;
        % ScalingFactor is set to 1
        ScalingFactor = 0.1;
        
        
        
        
   end
end
     