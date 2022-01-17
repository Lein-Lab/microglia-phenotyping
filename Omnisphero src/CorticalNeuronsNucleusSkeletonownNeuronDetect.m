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
function [wellList, foldername, ScalingFactor, csvHandler, neuronHandler,  Th] = CorticalNeuronsNucleusSkeletonownNeuronDetect (csvHandler, neuronHandler) 
% This is a quick script to convert 16-bit images of cortical neurons into
% binary image, to assess their center point and to generate a
% 'NucleusImage'.
n = 0;
minOverlapp = 1;
A = 1000;
B = 4;
N = 50;
C = 250;
D=10;
E = 15000;
NeuronBody=200;
NeuriteThreshold = 0.7;
summary = cell(10,6);
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
BackgroundBallFactor = 2;
SizeLeftoffBall = 5000;
TriangleBackground = 2;
TriangleBackgroundBall = 2;
for(i=1:numel(allfiles))
    %Check if file ends with .tif or .tiff
    ind = strfind([foldername '/' allfiles(i).name],'.bmp');
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
        %nucleusPic = NucleusBinary;
        %nucleusPic=uint8(nucleusPic);
        %nucleusPic(nucleusPic==0)=255;
        %nucleusPic(nucleusPic==1)=0;
        %Run Watershed from imageJ
        %javaaddpath('lib/ij.jar');
        %javaaddpath('lib/mij.jar');
        %MIJ.start(false);
        %MIJ.createImage(nucleusPic);
        %MIJ.run('Options...', 'iterations=1 count=1 edm=Overwrite do=Nothing');
        %MIJ.run('Watershed');
        %nucleusPicBin = logical(MIJ.getCurrentImage);
        %MIJ.run('Save', 'save=temp.tif');
        %MIJ.closeAllWindows;
        %MIJ.exit;
        %nucleusPicBinary = imread('temp.tif');
        %delete temp.tif;
        %nucleusPicBinary = logical(1-nucleusPicBinary);
        %nucleusPicBinary = bwareaopen(nucleusPicBinary,150);
        %Neuronal identification should be performed with Omnisphero.
        D = bwdist(~NucleusBinary);
        D = -D;
        D(~NucleusBinary)= -Inf;
        nucleusPicBinary = watershed(D);
        backgroundballNeurite =  imopen(uint8(imgNeurite),strel('ball',10,10));
        % The strel or imageJ ball background correction is also killing
        % big bright areas such as the Neuron Cell body. Therefore it is
        % better to delete this fields from the background correction to
        % don't delete the cell soma:
        [lehisto x] = imhist(backgroundballNeurite);
        level = triangle_th(lehisto,256);
        % This generates a mask of all cut off soma parts
        MaskBackground = im2bw(backgroundballNeurite,level*TriangleBackgroundBall);
        Area = bwareaopen(MaskBackground,SizeLeftoffBall);
        MaskBackground = MaskBackground - Area;
        %MaskBackgroundHard = im2bw(backgroundballNeurite,level*4);
        %MaskBackground = logical(MaskBackground);
        %MaskBackgroundHard = logical(MaskBackgroundHard);
        %MaskLable = bwlabel(MaskBackground);
        %maxL = max(max(MaskLable));
        %[r c] = size(MaskLable);
        %for(j=1:maxL)
        %    Dummy = zeros(size(MaskLable));
        %    for(a=1:r)
        %        for(b=1:c)
        %        if(MaskLable(a,b)==j && MaskBackgroundHard(a,b) == 1)
        %          Dummy(a,b) = MaskLable(a,b);
        %        end
        %        end
        %    end
        %    Label = max(max(Dummy));
        %    MatrixLabel(j,1) = Label; 
        %end    
        
        %sM = size(MatrixLabel);
        %MaskComplete = zeros(size(MaskLable));
        %MaskComplete = logical(MaskComplete);
        %for(h=1:sM)
        %   OneMask = zeros(size(MaskComplete));
        %   for(a=1:r)
        %        for(b=1:c)
        %        if(MaskLable(a,b)==h && MatrixLabel(h,1) > 0)
        %          OneMask(a,b) = MaskLable(a,b);
        %        end
        %        end
        %   end
        %   OneMask = logical(OneMask);
        %   MaskComplete = MaskComplete + OneMask;
        %end   
        
        % Cut these areas out from the background correction
        BackgroundWoSoma = zeros(size(imgNeurite));
        BackgroundWoSoma = uint8(BackgroundWoSoma);
        [r c] = size(BackgroundWoSoma);
        for(g=1:r)
            for(h=1:c)
            if(MaskBackground(g,h)==0)
              BackgroundWoSoma(g,h) =  backgroundballNeurite(g,h);
            else
              BackgroundWoSoma(g,h) =0;
            end
            end
        end    
        
        imgballbackgroundNeurite = imgNeurite - BackgroundWoSoma*BackgroundBallFactor;
        [lehisto x] = imhist(imgballbackgroundNeurite);
        level = triangle_th(lehisto,256); 
        binaryimgNeurite = im2bw(imgballbackgroundNeurite,level*TriangleBackground);
        
        
        
        
        
        %javaaddpath('lib/ij.jar');
        %javaaddpath('lib/mij.jar');
        %MIJ.start(false);
        %MIJ.createImage(imgNeurite);
        %MIJ.run('Options...', 'iterations=1 count=1 edm=Overwrite do=Nothing');
        %MIJ.run('Subtract Background...', 'rolling=10');
        %nucleusPicBin = logical(MIJ.getCurrentImage);
        %MIJ.run('Save', 'save=temp.tif');
        %MIJ.closeAllWindows;
        %MIJ.exit;
        %imgNeuriteBall = imread('temp.tif');
        %delete temp.tif;
        
        % Version with adaptive threshold:
        
        %binaryimgNeurite=adaptivethreshold(imgballbackgroundNeurite,300,0.0001,0);
        binaryimgNeurite = bwareaopen(imgballbackgroundNeurite,A);
                        %Fixed Threshold approach:
                        %Calculates the histogram of img1
                        %[lehisto x] = imhist(imgNeuriteBall);
                        %Calculates new background for img1 using the triangle method
                        %level = triangle_th(lehisto,256);
                        %Creates a binary image of img1 using level
                        %binaryimgNeurite = im2bw(imgNeuriteBall,level*NeuriteThreshold);
                        % Excludes all areas in 'imgnucleusbinary' smaller than N
                        %binaryimgNeurite = bwareaopen(binaryimgNeurite,A);
                        % Perform gaussian blur:
        %H = fspecial('gaussian',8);
        %binaryimgNeurite = imfilter(binaryimgNeurite,H);
        
        % This generates a hard threshold for Axonal Tracing only:
        [lehisto x] = imhist(imgballbackgroundNeurite);
        level = triangle_th(lehisto,256);
        binaryimgNeuriteHard = im2bw(imgballbackgroundNeurite,level*3);
        binaryimgNeuriteHard = bwareaopen(binaryimgNeuriteHard, 1000);
        %binaryimgNeuriteHard = imfilter(binaryimgNeuriteHard,H);
        
                        % Create an image with only the diffuse background bubbles:
                        % 1) Create a binaryimage with hard streshold to delete most
                        % neurites
                        %imgThreshhard = im2bw(imgNeuriteBall,level*NeuriteThreshold*6);
                        % 2) Create normal threshold image
                        %binaryNormal = im2bw(imgNeurite,level*4);
                        % Substract hard threshold from normal threshold:
                        %binarySubstracted = binaryNormal - imgThreshhard;
                        %binarySubstracted = bwareaopen(binarySubstracted,A);
                        %binaryimgNeurite = binaryimgNeurite -binarySubstracted;
                        %binaryimgNeurite = bwareaopen(binaryimgNeurite,E);
                
                
        
        
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
        imwrite(nucleusImg,[foldername1 '/' newFile9]);
        imwrite(nucleusImg,[foldername1 '/' newFile1]);
        imwrite(imgNeurite,[foldername1 '/' newFile2]);
        imwrite(imgNeurite,[foldername1 '/' newFile3]);
        %imwrite(binaryimgNeurite,[foldername1 '/' newFile4]);
        %Noch ImageJ fixen!
        imwrite(binaryimgNeurite,[foldername1 '/' newFile4]);
        nucleusPicBinary = logical(nucleusPicBinary);
        imwrite(nucleusPicBinary,[foldername1 '/' newFile5]);
        
        imwrite(imgSOX2,[foldername1 '/' newFile7]);
        imwrite(imgSOX2,[foldername1 '/' newFile8]);
        imwrite(binaryimgNeuriteHard,[foldername1 '/' newFile10]);
        wellname = (newFile1(1:3));
        wellList{n} = wellname;
        csvHandler.CellPosMatrix(wellname) = NucleusM;
        neuronHandler.NeuronPositionsEdgeFill(wellname) = sparse(zeros(size(nucleusPic)));
        ScalingFactor = 1;
        clear newfile newfile1 newfile2 newfile3 newfile4 newfile5 newfile6 newfile7 newfile8 
    end
    
end

% Comments: It would be good to include a suitable filter for bubble rims
