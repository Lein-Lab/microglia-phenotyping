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

function [summary] = Measurement (foldername, csvHandler, neuronHandler)
foldername = [foldername '/ConvertedCellomics'];
mkdir(foldername);
allfiles = dir(foldername);
wellList = cell(0);
n=0;
AxonMeasurement = 0;
for(i=1:numel(allfiles))
    %Check if file ends with .tif or .tiff
    ind = strfind([foldername '/ConvertedCellomics' allfiles(i).name],'Skeleton');
    if(numel(ind) > 0)%
         n = n + 1;
         SkeletonImage = imread([foldername '/' allfiles(i).name]);
         if(AxonMeasurement==1);
         newFileSkel = regexprep(allfiles(i).name,'Skeleton','AstroBig','ignorecase');
         SkeletonImageHard = imread([foldername '/' newFileSkel]);
         end
         newFileBinary = regexprep(allfiles(i).name,'Skeleton','Binary','ignorecase');  
         binaryimgNeurite = imread([foldername '/' newFileBinary]);
         if(AxonMeasurement==1);
         newFileBinaryHard = regexprep(allfiles(i).name,'Skeleton','BinaryHard','ignorecase');  
         binaryimgNeuriteHard = imread([foldername '/' newFileBinaryHard]);
         end
         %newFileNeuronNuc = regexprep(allfiles(i).name,'Skeleton','NeuronNuclei','ignorecase');
         %NeuronNuc = imread([foldername '/' newFileNeuronNuc]);
         %NeuronNuc = logical(NeuronNuc);
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
         wellname = (newFile1(1:3));
         NeuronM = neuronHandler.NeuronPositionsEdgeFill(wellname);
         NumberofNeurons=sum(sum(NeuronM));
        SkeletonImage = bwareaopen(SkeletonImage,200);
        SkeletonImage = uint8(SkeletonImage);
        %NeuronNuc = uint8(NeuronNuc);
        binaryimgNeurite = uint8(binaryimgNeurite);
        %binaryimgNeurite = binaryimgNeurite - NeuronNuc;
        %SkeletonWoNuc = SkeletonImage - NeuronNuc;
        %SkeletonWoNuc = logical(SkeletonWoNuc);
        binaryimgNeurite = logical(binaryimgNeurite);
        SkeletonImage = logical(SkeletonImage);
        %SkeletonWoNuc = logical(SkeletonWoNuc);
        if(AxonMeasurement==1);
        SkeletonImageHard = uint8(SkeletonImageHard);
        SkeletonImageHard = bwareaopen(SkeletonImageHard,200);
        NeuronNuc = uint8(NeuronNuc);
        binaryimgNeuriteHard = uint8(binaryimgNeuriteHard);
        %binaryimgNeuriteHard = binaryimgNeuriteHard - NeuronNuc;
        SkeletonWoNucHard = SkeletonImageHard - NeuronNuc;
        SkeletonWoNucHard = logical(SkeletonWoNucHard);
        binaryimgNeuriteHard = logical(binaryimgNeuriteHard);
        SkeletonImageHard = logical(SkeletonImageHard);
        SkeletonWoNucHard = logical(SkeletonWoNucHard);
        end
        %imwrite(SkeletonWoNuc,[foldername '/' allfiles(i).name]);
        %Morphological extracttion
        %1) Neurite Mass:
        NeuriteMass = sum(sum(binaryimgNeurite));
        %2) Total Neurite length:
        TotalNeuritelength = sum(sum(SkeletonImage));
        %3) Number of Neurites:
        NumberofNeurites = sum(sum(bwmorph(SkeletonImage,'endpoints')));
        %4) Average Neurte length:
        AverageNeuritelength = TotalNeuritelength;
        %5) Number of branching points:
        NumberBranches = sum(sum(bwmorph(SkeletonImage,'branchpoints')));
        %6) Number of endpoints:
        NumberofEndpoints = sum(sum(bwmorph(SkeletonImage,'endpoints')));
        %6) Neurite Mass Axons:
        if(AxonMeasurement==1);
        AxonMass = sum(sum(binaryimgNeuriteHard));
        %7) Total Neurite length:
        TotalAxonlength = sum(sum(SkeletonImageHard));
        %8) Number of Neurites:
        NumberofAxons = sum(sum(bwmorph(SkeletonImageHard,'endpoints')));
        %9) Average Neurte length:
        AverageAxonlength = TotalAxonlength
        %10) Number of branching points:
        NumberBranchesAxon = sum(sum(bwmorph(SkeletonWoNucHard,'branchpoints')));
        end
        
        
        
        
        
        
        %Create Matrix
        summary{1,1} =  '[NM]';
        summary{1,2} =  '[TNL]';
        %summary{1,3} =  '[NN]';
        %summary{1,4} =  '[ANL]';
        summary{1,3} =  '[NB]';
        summary{1,4} = '[NE]';
        if(AxonMeasurement==1);
        summary{1,6} =  '[AM]';
        summary{1,7} =  '[TAL]';
        summary{1,8} =  '[NA]';
        summary{1,9} =  '[AAL]';
        summary{1,10}=  '[NBA]';
        summary{1,11} = '[NN]';
        else
        %summary{1,6} = '[NN]'; 
        end
        m=n+1;
        summary{m,1} =  NeuriteMass;
        summary{m,2} =  TotalNeuritelength;
        %summary{m,3} =  NumberofNeurites;
        %summary{m,4} =  AverageNeuritelength;
        summary{m,3} =  NumberBranches;
        summary{m,4} = NumberofEndpoints;
        if(AxonMeasurement==1);
        summary{m,6} =  AxonMass;
        summary{m,7} =  TotalAxonlength;
        summary{m,8} =  NumberofAxons;
        summary{m,9} =  AverageAxonlength;
        summary{m,10}=  NumberBranchesAxon;
        summary{m,11} =  NumberofNeurons;
        else
        %summary{m,6} =  NumberofNeurons;
        end
        summary = summary
    end
end
         