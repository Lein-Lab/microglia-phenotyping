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
function [wellList, foldername, ScalingFactor, csvHandler, Th] = CorticalNeuronsNucleus (csvHandler) 
% This is a quick script to convert 16-bit images of cortical neurons into
% binary image, to assess their center point and to generate a
% 'NucleusImage'.
n = 0;
A = 4000;
B = 4;
N = 50;
C = 250;

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

for(i=1:numel(allfiles))
    %Check if file ends with .tif or .tiff
    ind = strfind([foldername '/' allfiles(i).name],'.jpg');
    if(numel(ind) > 0)%
         n = n + 1;
        %Load RGB image
        imgRGB = imread([foldername '/' allfiles(i).name]);
        %Split Channels
        nucleusPic = imgRGB(:,:,3);
        nucleusImg = nucleusPic;
        imgSOX2 = imgRGB(:,:,1);
        imgNeurite = imgRGB(:,:,2); 
        %Indentification Nuclei
        backgroundNucleus = imopen(uint8(nucleusPic),strel('disk',N));
        NucleusCorrected = nucleusPic - backgroundNucleus;
        [lehisto x] = imhist(NucleusCorrected);
        %Calculates new background for img1 using the triangle method
        level = triangle_th(lehisto,256);
        NucleusBinary = im2bw(NucleusCorrected,level*1.5);
        %NucleusBinary = bwareaopen(NucleusBinary,C);
        nucleusPic = NucleusBinary;
        nucleusPic=uint8(nucleusPic);
        nucleusPic(nucleusPic==0)=255;
        nucleusPic(nucleusPic==1)=0;
        %Run Watershed from imageJ
        javaaddpath('lib/ij.jar');
        javaaddpath('lib/mij.jar');
        MIJ.start(false);
        MIJ.createImage(nucleusPic);
        MIJ.run('Options...', 'iterations=1 count=1 edm=Overwrite do=Nothing');
        MIJ.run('Watershed');
        %nucleusPicBin = logical(MIJ.getCurrentImage);
        MIJ.run('Save', 'save=temp.tif');
        MIJ.closeAllWindows;
        MIJ.exit;
        nucleusPicBinary = imread('temp.tif');
        delete temp.tif;
        nucleusPicBinary = logical(1-nucleusPicBinary);
        %Neuronal identification should be performed with Omnisphero.
        backgroundballNeurite =  imopen(uint8(imgNeurite),strel('disk',N));
        imgballbackgroundNeurite = imgNeurite - backgroundballNeurite;
        %Calculates the histogram of img1
        [lehisto x] = imhist(imgballbackgroundNeurite);
        %Calculates new background for img1 using the triangle method
        level = triangle_th(lehisto,256);
        %Creates a binary image of img1 using level
        binaryimgNeurite = im2bw(imgballbackgroundNeurite,level*0.8);
        % Excludes all areas in 'imgnucleusbinary' smaller than N
        binaryimgNeurite = bwareaopen(binaryimgNeurite,A);
        % The skeleton image is now created from the binary image
        %skelImg = bwmorph(binaryimgNeurite,'skel', Inf);
        %The raw skeleton image contains artifical branches and short processes. These are eliminated with the next step:
        s = regionprops(nucleusPicBinary, 'Centroid');
        sdimensions = size(s);
        NucleusM = sparse(zeros(size(nucleusPic)));
        if(sdimensions > 0)
            for(j=1:sdimensions)
            x_centroid = s(j).Centroid(1);
            y_centroid = s(j).Centroid(2);
            NucleusM(uint16(y_centroid),uint16(x_centroid)) = 1;
            end
        end
        a{1} = 'A0';
        a{2} = int2str(n);
        a{3} = 'space';
        newfile = [a{1} a{2} a{3}];
        newFile1 = regexprep(newfile, 'space','NucleusSmall.png','ignorecase');
        newFile2 = regexprep(newfile, 'space','NeuriteBig.png','ignorecase'); 
        newFile3 = regexprep(newfile, 'space','NeuriteSmall.png','ignorecase');
        newFile4 = regexprep(newfile, 'space','Binary.png','ignorecase');
        newFile5 = regexprep(newfile, 'space','NucleusBigWatershed.png','ignorecase');
        newFile6 = regexprep(newfile, 'space','Skeleton.png','ignorecase');
        newFile7 = regexprep(newfile, 'space','OligoBig.png','ignorecase');
        newFile8 = regexprep(newfile, 'space','OligoSmall.png','ignorecase');
        newFile9 = regexprep(newfile, 'space','NucleusBig.png','ignorecase');
        imwrite(nucleusImg,[foldername1 '/' newFile9]);
        imwrite(nucleusImg,[foldername1 '/' newFile1]);
        imwrite(imgNeurite,[foldername1 '/' newFile2]);
        imwrite(imgNeurite,[foldername1 '/' newFile3]);
        %imwrite(binaryimgNeurite,[foldername1 '/' newFile4]);
        %Noch ImageJ fixen!
        imwrite(nucleusPicBinary,[foldername1 '/' newFile5]);
        %imwrite(skelImg,[foldername1 '/' newFile6]);
        imwrite(imgSOX2,[foldername1 '/' newFile7]);
        imwrite(imgSOX2,[foldername1 '/' newFile8]);
        wellname = (newFile1(1:3));
        wellList{n} = wellname;
        csvHandler.CellPosMatrix(wellname) = NucleusM;
        ScalingFactor = 1;
        clear newfile newfile1 newfile2 newfile3 newfile4 newfile5 newfile6 newfile7 newfile8
    end
    
end

% Comments: It would be good to include a suitable filter for bubble rims
