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

function [wellList, foldername, ScalingFactor, csvHandler] = SCG_Morphology (csvHandler) 
n = 0;
foldername = uigetdir;
Name = foldername(end-10:end);
%foldername = [foldername '/TimePoint_1'];
mkdir(foldername,'ConvertedCellomics');
foldername1 = [foldername '/ConvertedCellomics'];
allfiles = dir(foldername);
wellList = cell(0);
for(i=1:numel(allfiles))
    ind = strfind([foldername '/' allfiles(i).name],'w1.TIF');
     if(numel(ind) > 0)
         n =n+1;    
         NucleusImage = imread([foldername '/' allfiles(i).name]);
         newFileNeuron = regexprep(allfiles(i).name, 'w1.TIF','w2.TIF','ignorecase');
         NeuriteImage = imread([foldername '/' newFileNeuron]);



        [lehisto x] = imhist(NeuriteImage);
        level = triangle_th(lehisto,256);
        if(level>0.01)
            level = level *4;
        end    
        CellBodies = im2bw(NeuriteImage,level);
        se = strel('disk',3);
        CellBodies = imopen(CellBodies,se);
        CellBodies = bwareaopen(CellBodies,5000);
        ProcessesThin = bwmorph(CellBodies, 'thin',5);
        ProcessesCellBodiesThin = bwmorph(ProcessesThin,'open');
        ProcessesThin = ProcessesThin -ProcessesCellBodiesThin; 
        Processes = bwmorph(CellBodies, 'thin',15);
        Processes = bwmorph(Processes,'spur',2);
        ProcessesCellBodies = bwmorph(Processes,'open');
        Processes = Processes- ProcessesCellBodies;
        Processes = uint8(Processes);
        Processes = logical(Processes);
        Processes = bwareaopen(Processes,20);
        %Here we will see if Processes and ProcessesThin overlapp the
        %combine them:
        ProcessesL = bwlabel(Processes);
        ProcessesN = max(max(ProcessesL));
        ProcessesKeep = zeros(size(CellBodies));
        for(r=1:ProcessesN)
            ProcessessOne = ProcessesL;
            ii = ProcessessOne == r;
            ProcessessOne(ii) = 1000;
            ProcessessOne = ProcessessOne - 995;
            ProcessessOne = uint8(ProcessessOne);
            ProcessessOne = logical(ProcessessOne);
            TestOver = ProcessessOne + ProcessesThin;
            TestOver = TestOver - 1;
            TestOver = uint8(TestOver);
            TestOver = logical(TestOver);
            if(sum(sum(TestOver))>5 && (sum(sum(TestOver))/sum(sum(ProcessessOne)))<0.95)
                ProcessesKeep = ProcessesKeep + ProcessessOne;
            end
        end    
        Processes = Processes + ProcessesKeep;
        Processes = bwmorph(Processes,'thin',inf);
        Processes = logical(Processes);
        Processes = bwareaopen(Processes,20);
        ProcessLength = sum(sum(Processes));
        CellBodiesL = bwlabel(CellBodies);
        CellN = max(max(CellBodiesL));
        %We can now also save the images, in particular only the Cellbodies
        %with processes:
        Binary = Processes + ProcessesCellBodies;
        Binary = logical(Binary);
        
        a{1} = allfiles(i).name(end-12:end-7);
        a{2} = 'space';
        newfile = [a{1} a{2}];
        newFile1 = regexprep(newfile, 'space','NucleusBig.png','ignorecase');
        imwrite(NucleusImage,[foldername1 '/' newFile1]);
        newFile2 = regexprep(newfile, 'space','NucleusSmall.png','ignorecase');
        imwrite(NucleusImage,[foldername1 '/' newFile2]);
        newFile3 = regexprep(newfile, 'space','NeuriteBig.png','ignorecase');
        imwrite(NeuriteImage,[foldername1 '/' newFile3]);
        newFile4 = regexprep(newfile, 'space','NeuriteSmall.png','ignorecase');
        imwrite(NeuriteImage,[foldername1 '/' newFile4]);
        newFile5 = regexprep(newfile, 'space','Binary.png','ignorecase');
        imwrite(Binary,[foldername1 '/' newFile5]);
        
        wellname = allfiles(i).name(end-12:end-7);
        wellList{n} = allfiles(i).name(end-12:end-7);
        csvHandler.CellPosMatrix(wellname) = sparse(zeros(size(Binary)));
        ScalingFactor = 1;
        
        ProcessLengthTabelHeader{1,1} = 'Image Name';
        ProcessLengthTabelHeader{1,2} = 'Length per cell';
        ProcessLengthTabelHeader{1,3} = 'Number of cells';
        ProcessLengthTabelName{n+1,1} = allfiles(i).name;
        ProcessLengthTabel(n+1,2) = ProcessLength/CellN;
        ProcessLengthTabel(n+1,3) = CellN;
        
     end
end     

