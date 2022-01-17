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
minOverlapp = 0.3;
A = 15000;
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
        nucleusPicBinary = bwareaopen(nucleusPicBinary,150);
        %Neuronal identification should be performed with Omnisphero.
        
        %backgroundballNeurite =  imopen(uint8(imgNeurite),strel('disk',N));
        %imgballbackgroundNeurite = imgNeurite - backgroundballNeurite;
        %[lehisto x] = imhist(imgballbackgroundNeurite);
        %level = triangle_th(lehisto,256);
        %imgballbackgroundNeurite = im2bw(imgballbackgroundNeurite,level*NeuriteThreshold);
        %imgballbackgroundNeurite = bwareaopen(imgballbackgroundNeurite,A);
        %imgballbackgroundNeurite = uint8(imgballbackgroundNeurite);
        %Integration of the rolling ball method from imageJ into Matlab:
        
        
        
        
        javaaddpath('lib/ij.jar');
        javaaddpath('lib/mij.jar');
        MIJ.start(false);
        MIJ.createImage(imgNeurite);
        MIJ.run('Options...', 'iterations=1 count=1 edm=Overwrite do=Nothing');
        MIJ.run('Subtract Background...', 'rolling=10');
        %nucleusPicBin = logical(MIJ.getCurrentImage);
        MIJ.run('Save', 'save=temp.tif');
        MIJ.closeAllWindows;
        MIJ.exit;
        imgNeuriteBall = imread('temp.tif');
        delete temp.tif;
        
        % Version with adaptive threshold:
        
        binaryimgNeurite=adaptivethreshold(imgNeuriteBall,200,0.0001,0);
        binaryimgNeurite = bwareaopen(binaryimgNeurite,A);
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
        H = fspecial('gaussian',8);
        binaryimgNeurite = imfilter(binaryimgNeurite,H);
        
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
                
                
        % Creation of skeleton image with ImageJ:
        binaryimgNeurite=uint8(binaryimgNeurite);
        binaryimgNeurite(binaryimgNeurite==0)=255;
        binaryimgNeurite(binaryimgNeurite==1)=0;
        javaaddpath('lib/ij.jar');
        javaaddpath('lib/mij.jar');
        MIJ.start(false);
        MIJ.createImage(binaryimgNeurite);
        MIJ.run('Options...', 'iterations=1 count=1 edm=Overwrite do=Nothing');
        MIJ.run('Skeletonize');
        %nucleusPicBin = logical(MIJ.getCurrentImage);
        MIJ.run('Save', 'save=temp.tif');
        MIJ.closeAllWindows;
        MIJ.exit;
        SkeletonImage = imread('temp.tif');
        delete temp.tif;
        SkeletonImage = logical(1-SkeletonImage);
        binaryimgNeurite = logical(1-binaryimgNeurite);
        % Skeleton processing
        %1) Substracting the nucleus Area
        old = SkeletonImage;
        %BranchingPoints = bwmorph(SkeletonImage,'branchpoints');
        %skelWoBranches = SkeletonImage - BranchingPoints;
        %SaveBranches = BranchingPoints;
        %for(l=1:10)
        %if(sum(sum(BranchingPoints)) > 0)
        %    BranchingPoints = bwmorph(skelWoBranches,'branchpoints');
        %    SaveBranches = SaveBranches + BranchingPoints;
        %    skelWoBranches = skelWoBranches - BranchingPoints;
        %end
        %end
        %del = xor (bwareaopen(skelWoBranches,4),bwareaopen(skelWoBranches,10));
        %skelWoBranches = skelWoBranches -del;
        %skelWoBranches = skelWoBranches + SaveBranches;
        SkeletonImage = bwmorph(old,'spur',4);
        SkeletonImage = bwareaopen(SkeletonImage,50);
        %SkeletonImage = bwmorph(SkeletonImage,'spur',4);
        %SkeletonImage = SkeletonImage - nucleusPicBinary;
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
        
        %Neuronal Identification based on overlap criteria
        % Create overlap image:        
        binaryimgNeurite = uint8(binaryimgNeurite);
        % Look for close connected components for nucleusImage:
        %1) Label different areas:
        labelBinaryNucleus = bwlabel(nucleusPicBinary);
        %Count number of nuclei present
        k = max(max(labelBinaryNucleus));
        % Function to extract only neuron nuclei. 
        %Create empty image with size equal to original image
        [r,c] = size(labelBinaryNucleus);
        neuronNuclei = zeros(r,c);
        neuronNuclei = uint8(neuronNuclei);
        labelSingleNucleus = zeros(r,c);
        %count = 0;
        % Now extract each single lable (one separate Nuclei) and check
        % whether the overlapp with the neurite channel is sufficient to be
        % a neuron nuclei. If so keep it, else delete it.
        for(w=1:k)
        % This loop iterates over the image for a certain lable and extracts it. E.g. Nuclei with lable one is extracted in the first step.        
           for(u=1:r)
           for(v=1:c)
               if(labelBinaryNucleus(u,v)==w)
                 labelSingleNucleus(u,v)=1;
              else
                 labelSingleNucleus(u,v)=0; 
              end   
           end
           end
           % Add Neurite Image to the single nuclei image.
           labelSingleNucleus = uint8(labelSingleNucleus);
           CloseConnectedComponent = labelSingleNucleus +  binaryimgNeurite;
           % Count number of Pixels with value 2 (equals the close
           % connected component of nuclei and neurons staining!)
           numberPixelOverlay = length(CloseConnectedComponent(CloseConnectedComponent==2));
           % Count number of pixel of nuclei
           numberPixelNuclei = length(labelSingleNucleus(labelSingleNucleus==1));
           % Check if overlapp is bigger than minOverlapp:
           if(numberPixelOverlay/numberPixelNuclei>minOverlapp)
              neuronNuclei = neuronNuclei + labelSingleNucleus;
              %count = count +1;
           end
           tell = w;
        end
        neuronNuclei = logical(neuronNuclei);
        sn = regionprops(neuronNuclei,'centroid');
        sndimensions = size(sn);
        % Generating the Neuron matrix:
        NeuronM = sparse(zeros(size(neuronNuclei)));
        if(sndimensions > 0)
            for(b=1:sndimensions)
            xn_centroid = sn(b).Centroid(1);
            yn_centroid = sn(b).Centroid(2);
            NeuronM(uint16(yn_centroid),uint16(xn_centroid)) = 1;
            end
        end
        neuronNuclei = uint8(neuronNuclei);
        NumberofNeurons = sndimensions(1);
        %Substract Cell soma:
        SkeletonImage = uint8(SkeletonImage);
        SkeletonImage = SkeletonImage - neuronNuclei;
        
        %Extract morpholgoical characteristics
        %1) NeuriteMass:
        SkeletonImage = uint8(SkeletonImage);
        NeuriteMass = (sum(sum(binaryimgNeurite)))/NumberofNeurons;
        %2) Total Neurite length:
        TotalNeuritelength = ((sum(sum(SkeletonImage)))/NumberofNeurons);
        %3) Number of Neurites:
        NumberofNeurites =  max(max(bwlabel(SkeletonImage)));
        NumberofNeurites = NumberofNeurites/NumberofNeurons;
        %4) Average Neurte length:
        AverageNeuritelength = TotalNeuritelength/NumberofNeurites;
        %5) Number of branching points:
        NumberBranches = sum(sum(bwmorph(SkeletonImage,'branchpoints')))/NumberofNeurons;
        %Create Matrix
        summary{n,1} =  NeuriteMass;
        summary{n,2} =  TotalNeuritelength;
        summary{n,3} =  NumberofNeurites;
        summary{n,4} =  AverageNeuritelength;
        summary{n,5} =  NumberBranches;
        summary{n,6} =  NumberofNeurons;
        summary = summary
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
        imwrite(binaryimgNeurite,[foldername1 '/' newFile4]);
        imwrite(nucleusPicBinary,[foldername1 '/' newFile5]);
        imwrite(SkeletonImage,[foldername1 '/' newFile6]);
        imwrite(imgSOX2,[foldername1 '/' newFile7]);
        imwrite(imgSOX2,[foldername1 '/' newFile8]);
        wellname = (newFile1(1:3));
        wellList{n} = wellname;
        csvHandler.CellPosMatrix(wellname) = NucleusM;
        neuronHandler.NeuronPositionsEdgeFill(wellname) = NeuronM;
        ScalingFactor = 1;
        clear newfile newfile1 newfile2 newfile3 newfile4 newfile5 newfile6 newfile7 newfile8 
    end
    
end

% Comments: It would be good to include a suitable filter for bubble rims
