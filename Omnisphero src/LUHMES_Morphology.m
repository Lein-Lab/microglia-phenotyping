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

function [wellList, foldername, ScalingFactor, csvHandler,   Summary, SummaryName, Summary2, Summary2Name, Summary3Name, Summary3, Summary4 SummaryEnd SummaryEndTop Name Stats] = LUHMES_Morphology (csvHandler) 
% This function is designed to get the total neurite area, total neurite
% length, number of branches and number of neurites in Luhmes Cells!

n = 0;
ThresholdMultiplier = 1.5;
AreaCellbodieExclude = 300;
AreaLogLow = 50;
HoleArea = 50;
UpperLimitCellBody = 1129*1.5;
MaxAllowedArea = 250000;
NumberofSights = inputdlg('Number of sights');
NumberofSights = str2double(NumberofSights);
Question = questdlg('Do you want to analyze synapses?')
if(strcmp(Question,'No'))
    Stats = 0;
else
    Stats = 1;
end
foldername = uigetdir;
Name = foldername(end-10:end);
foldername = [foldername '/TimePoint_1'];
mkdir(foldername,'ConvertedCellomics');
foldername1 = [foldername '/ConvertedCellomics'];
allfiles = dir(foldername);
wellList = cell(0);
Multiplier = 0;
for(i=1:numel(allfiles))
    ind = strfind([foldername '/' allfiles(i).name],'w1.TIF');
     if(numel(ind) > 0)
     n =n+1;   
     % Reading in images:
     NucleusImage = imread([foldername '/' allfiles(i).name]);
     newFileNeuron = regexprep(allfiles(i).name, 'w1.TIF','w4.TIF','ignorecase');
     NeuriteImage = imread([foldername '/' newFileNeuron]);
     if(Stats ==1)
       newFileSyn = regexprep(allfiles(i).name, 'w1.TIF','w2.TIF','ignorecase');
       SynImage = imread([foldername '/' newFileSyn]);
       newFilePSD95 = regexprep(allfiles(i).name, 'w1.TIF','w3.TIF','ignorecase');
       PSD95Image = imread([foldername '/' newFilePSD95]);
     end    
     
     
     if(Stats==1)
     % First we will get the neurite area with normal thresholding
     [lehisto x] = imhist(NeuriteImage);
     level = triangle_th(lehisto,256);
     NeuriteArea = im2bw(NeuriteImage,level);
     %Now we will mask out this threshold area and subtract it from the
     %original image. Consequently we will threshold this for dim objects!
     NeuriteAreaMask = NeuriteArea*65535;
     NeuriteAreaMask = uint16(NeuriteAreaMask);
     NeuriteImageDim = NeuriteImage - NeuriteAreaMask;
     [lehisto x] = imhist(NeuriteImageDim);
     level = triangle_th(lehisto,256);
     if(level<0.01)
         [lehisto x] = imhist(NeuriteImage);
         level = triangle_th(lehisto,256);
         NeuriteImageDimBinary = im2bw(NeuriteImageDim,level);
     else        
     NeuriteImageDimBinary = im2bw(NeuriteImageDim,level*2);
     end
     NeuriteArea = NeuriteArea + NeuriteImageDimBinary;
     NeuriteMass = sum(sum(NeuriteArea));
     %Now we will threshold the Synapsin
     % We will use the rolling ball methd here first
     SE = strel('ball',5,5);
     backgroundball = imopen(SynImage,SE);
     SynImageBall = SynImage - backgroundball;
     [lehisto x] = imhist(SynImageBall);
     level = triangle_th(lehisto,256);
     SynImageBinary = im2bw(SynImageBall,level*2);
     SynImageBinary = bwareaopen(SynImageBinary,2);
     SynImageBinary = SynImageBinary + NeuriteArea;
     SynImageBinary = SynImageBinary -1;
     SynImageBinary = uint8(SynImageBinary);
     SynImageBinary = logical(SynImageBinary);
     SynArea = sum(sum(SynImageBinary));
     SynNumber = bwlabel(SynImageBinary);
     SynNumber = max(max(SynNumber));
     %Now we will threshold the PSD95
     % We will use the rolling ball methd here first
     SE = strel('ball',5,5);
     backgroundball = imopen(PSD95Image,SE);
     PSD95ImageBall = PSD95Image - backgroundball;
     [lehisto x] = imhist(PSD95ImageBall);
     level = triangle_th(lehisto,256);
     PSD95ImageBinary = im2bw(PSD95ImageBall,level*2);
     PSD95ImageBinary = bwareaopen(PSD95ImageBinary,2);
     PSD95ImageBinary = PSD95ImageBinary + NeuriteArea;
     PSD95ImageBinary = PSD95ImageBinary -1;
     PSD95ImageBinary = uint8(PSD95ImageBinary);
     PSD95ImageBinary = logical(PSD95ImageBinary);
     PSD95Area = sum(sum(PSD95ImageBinary));
     PSD95Number = bwlabel(PSD95ImageBinary);
     PSD95Number = max(max(PSD95Number));
     
     
     else    
     % Instead of using cell nuclei, we will try to use the cellbodies from
     % the Neurite channel for counting of neurons.
     % Idea we will take a very stringent threshold, to only  maintain the
     % cell bodies:
     try
     level = isodata(NeuriteImage);
     catch
         level = 0;
     end    
     if(level < 0.02)
         level = 0.02;
     end    
     CellBodies = im2bw(NeuriteImage,level*ThresholdMultiplier);
     CellBodies = bwareaopen(CellBodies,AreaCellbodieExclude);
     
     % Here we will delete artficat areas which are too big:
     Artefacts = bwareaopen(CellBodies,MaxAllowedArea);
     %Now we need to check whether it is a washing rim or a huge cell
     %cluster. Therefore, we will try to split up the area:
     
     ArtefactsLabel = bwlabel(Artefacts);
     NumberArtefacts = max(max(ArtefactsLabel));
     ArtefactsReal = zeros(size(NeuriteImage));
     
     for(u=1:NumberArtefacts)
         ArtefactOne = ArtefactsLabel;
         ii = ArtefactOne == u;
         ArtefactOne(ii) = 255;
         ArtefactOne = ArtefactOne - 250;
         ArtefactOne = uint8(ArtefactOne);
         ArtefactOne = logical(ArtefactOne);
         ArtefactOneTest = bwmorph(ArtefactOne,'thin',20);
         ArtefactOneTest = bwmorph(ArtefactOneTest,'open');
         ArtefactOneTest = bwareaopen(ArtefactOneTest,MaxAllowedArea);
         ArtefactOneTestHole = imfill(ArtefactOneTest,'holes');
         if(sum(sum(ArtefactOneTest))>0)
             if((sum(sum(ArtefactOneTestHole))/sum(sum(ArtefactOneTest)))>0.95);
             ArtefactsReal = ArtefactsReal + ArtefactOne;
             end
         end
     end    
         
     % The holes situation did not work on empty images which simply produce a very large area!    
     
     CellBodies = CellBodies - ArtefactsReal;
     CellBodies = uint8(CellBodies);
     CellBodies = logical(CellBodies);
 
     % Now we will save the area distribution of all objects for later
     % assessment of cell body count
     
     Area = regionprops(CellBodies,'Area');
     AreaVector = struct2table(Area);
     AreaVector = table2array(AreaVector);
     AreaNumber = length(AreaVector);
     
     % Here we will normalize the areas to the area of a single cell body
     % obtained from the median of all areas:
     
     AreaVectorNormalized = AreaVector;
     AreaNumber = length(AreaVectorNormalized);
     for(k = 1:AreaNumber)
         if(AreaVectorNormalized(k,1) > UpperLimitCellBody)
             AreaVectorNormalized(k,1) = AreaVectorNormalized(k,1)/UpperLimitCellBody;
             AreaVectorNormalized(k,1) = round(AreaVectorNormalized(k,1));
         else
             AreaVectorNormalized(k,1) = 1;
         end
     end   
     
     
     
     
     SummaryName{1,n} = allfiles(i).name;
     Summary(2:AreaNumber+1,n) = AreaVector;
     
     
     
     % Now we can compute the complete binary and subtract the cell bodies:
     % Compute complete image:
     
     % We can first subtract the Cellbodies:
     
     %CellBodies = bwmorph(CellBodies,'thicken',5);
     ArtefactsReal = bwmorph(ArtefactsReal,'thicken',30);
     ArtefactsReal = bwmorph(ArtefactsReal,'bridge');
     ArtefactsReal = uint16(ArtefactsReal);
     ArtefactsReal = ArtefactsReal * 65535;
     CellBodies = uint16(CellBodies);
     CellBodies = CellBodies *65535;
     NeuriteImageWoBody = NeuriteImage - CellBodies - ArtefactsReal;
     NeuriteImageWoBody = uint16(NeuriteImageWoBody);
     CellBodies = logical(CellBodies);
     try
     level = isodata(NeuriteImageWoBody);
     catch
     level = 1;    
     end   
     if(level<0.02)
         level = 0.1;
     end    
     NeuritesFill = im2bw(NeuriteImageWoBody,level);
     NeuritesFillOrig =  NeuritesFill;
     % Now we can threshold again the remaining part to catch dim
     % structures:
     %NeuritesFill = bwmorph(NeuritesFill,'thicken',5);
     if(sum(sum(NeuritesFill))>0)
     NeuritesFill = uint16(NeuritesFill);
     NeuritesFill = NeuritesFill *65535;
     DimNeurites = NeuriteImageWoBody -  NeuritesFill;
     DimNeurites = uint16(DimNeurites);
     try
     level = isodata(DimNeurites);
     if(level<0.01)
         level = 0.01;
     end    
         
     DimNeurites = im2bw(DimNeurites,level*1.5);
     catch
     end    
     %DimNeurites = bwareaopen(DimNeurites, 50);
     NeuritesFill = logical(NeuritesFill);
     DimNeurites = logical(DimNeurites);
     NeuritesComplete = NeuritesFillOrig + DimNeurites;
     
     s = regionprops(CellBodies, 'Centroid');
        sdimensions = size(s);
        NucleusM = sparse(zeros(size(CellBodies)));
        if(sdimensions > 0)
            for(j=1:sdimensions)
            x_centroid = s(j).Centroid(1);
            y_centroid = s(j).Centroid(2);
            NucleusM(uint16(y_centroid),uint16(x_centroid)) = 1;
            end
        end
     else
         NeuritesComplete = zeros(size(NeuritesFill));
         NucleusM = sparse(zeros(size(CellBodies)));
     end 
     
     
     % Here we will save the relevant endpoints:
     NumberofNeurons = sum(AreaVectorNormalized);
     NeuriteArea = sum(sum(NeuritesComplete));
     NeuriteAreaperNeuron = NeuriteArea/NumberofNeurons;
     
     if(n == 1);
     y = 2;
     t = 2;
     Summary2Name{1,y} =  allfiles(i).name;
     Summary3Name{t,1} =  allfiles(i).name(end-7);
     Summary2(t,y) =  NumberofNeurons;
     Summary3(t,y) =  NeuriteArea;
     Summary4(t,y) =  NeuriteAreaperNeuron;
     else
         if(sum(allfiles(i).name(end-11:end-10) == '01') ==2 && allfiles(i).name(end-7) == '1')
             Multiplier = Multiplier + 1;
             y = 1;
         else
             Multiplier = Multiplier;
         end    
         if(allfiles(i).name(end-7) == '1')
             t = 2+(NumberofSights*Multiplier);
             y = y+1;
             Summary3Name{t,1} =  allfiles(i).name(end-7);
             Summary2Name{1,y} =  allfiles(i).name;
             Summary2(t,y) = NumberofNeurons;
             Summary3(t,y) =  NeuriteArea;
             Summary4(t,y) =  NeuriteAreaperNeuron;
         else
             t = t + 1;
             Summary2(t,y) = NumberofNeurons;
             Summary3(t,y) =  NeuriteArea;
             Summary4(t,y) =  NeuriteAreaperNeuron;
             Summary3Name{t,1} =  allfiles(i).name(end-7);
         end
         
     end    
     
    Summary2 = Summary2
     
    
        % Here will get the well name which will be composed of three
        % parts:
        
        %1) Well letter
    
        if(allfiles(i).name(end-12) == 'A')
            a{1} = 'A';
        end
        if(allfiles(i).name(end-12) == 'B')
            a{1} = 'B';
        end
        if(allfiles(i).name(end-12) == 'C')
            a{1} = 'C';
        end
        if(allfiles(i).name(end-12) == 'D')
            a{1} = 'D';
        end
        if(allfiles(i).name(end-12) == 'E')
            a{1} = 'E';
        end
        if(allfiles(i).name(end-12) == 'F')
            a{1} = 'F';
        end
        if(allfiles(i).name(end-12) == 'G')
            a{1} = 'G';
        end
        if(allfiles(i).name(end-12) == 'H')
            a{1} = 'H';
        end
        
        
        %1) WellNumber
        
        if(allfiles(i).name(end-12) == 'A')
            a{1} = 'A';
        end
        if(allfiles(i).name(end-12) == 'B')
            a{1} = 'B';
        end
        if(allfiles(i).name(end-12) == 'C')
            a{1} = 'C';
        end
        if(allfiles(i).name(end-12) == 'D')
            a{1} = 'D';
        end
        if(allfiles(i).name(end-12) == 'E')
            a{1} = 'E';
        end
        if(allfiles(i).name(end-12) == 'F')
            a{1} = 'F';
        end
        if(allfiles(i).name(end-12) == 'G')
            a{1} = 'G';
        end
        if(allfiles(i).name(end-12) == 'H')
            a{1} = 'H';
        end
        
        %This list gives the row:
        
        if(allfiles(i).name(end-11:end-10) == '01')
            a{2} = '01';
        end
        if(allfiles(i).name(end-11:end-10) == '02')
            a{2} = '02';
        end
        if(allfiles(i).name(end-11:end-10) == '03')
            a{2} = '03';
        end
        if(allfiles(i).name(end-11:end-10) == '04')
            a{2} = '04';
        end
        if(allfiles(i).name(end-11:end-10) == '05')
            a{2} = '05';
        end
        if(allfiles(i).name(end-11:end-10) == '06')
            a{2} = '06';
        end
        if(allfiles(i).name(end-11:end-10) == '07')
            a{2} = '07';
        end
        if(allfiles(i).name(end-11:end-10) == '08')
            a{2} = '08';
        end
        if(allfiles(i).name(end-11:end-10) == '09')
            a{2} = '09';
        end
        if(allfiles(i).name(end-11:end-10) == '10')
            a{2} = '10';
        end
        if(allfiles(i).name(end-11:end-10) == '11')
            a{2} = '11';
        end
        if(allfiles(i).name(end-11:end-10) == '12')
            a{2} = '12';
        end
        
        %This list gives the sight:
           
       EndPhrase =    allfiles(i).name(end-12);
       if(EndPhrase == '0' | EndPhrase == '1' )
           a{3} = 's';
       else
           a{3} = 's0';
       end
       
       if(a{3} == 's')
           a{4} =allfiles(i).name(end-7:end-6);
       else
           a{4} = allfiles(i).name(end-7);
       end    
         
       a{5} = 'space';     
        
        newfile = [a{1} a{2} a{3} a{4} a{5}];
        
        newFileName{2,1} = 'A';
        j= 2+NumberofSights;
        newFileName{j,1} = 'B';
        j= j+NumberofSights;
        newFileName{j,1} = 'C';
        j= j+NumberofSights;
        newFileName{j,1} = 'D';
        j= j+NumberofSights;
        newFileName{j,1} = 'E';
        j= j+NumberofSights;
        newFileName{j,1} = 'F';
        j= j+NumberofSights;
        newFileName{j,1} = 'G';
        j= j+NumberofSights;
        newFileName{j,1} = 'H';
        
       Numberrow = [1 2 3 4 5 6 7 8 9 10 11 12];
       SummaryEndTop(1,2:13) = Numberrow;
       SummaryEnd = newFileName;
        
        % Names of images
     newFile1 = regexprep(newfile, 'space','NucleusSmall.png','ignorecase');
     newFile2 = regexprep(newfile, 'space','NucleusBig.png','ignorecase'); 
     newFile3 = regexprep(newfile, 'space','NeuriteSmall.png','ignorecase');
     newFile4 = regexprep(newfile, 'space','NeuriteBig.png','ignorecase');
     newFile5 = regexprep(newfile, 'space','Binary.png','ignorecase');
     imwrite(NucleusImage,[foldername1 '/' newFile1]);
     imwrite(NucleusImage,[foldername1 '/' newFile2]);
     imwrite(NeuriteImage,[foldername1 '/' newFile3]);
     imwrite(NeuriteImage,[foldername1 '/' newFile4]);
     imwrite(NeuritesComplete,[foldername1 '/' newFile5]);
     
     wellname = (newFile1(1:6));
     wellList{n} = wellname;
     csvHandler.CellPosMatrix(wellname) = NucleusM;
     ScalingFactor = 1;
     end
     if(Stats ==1)
        if(allfiles(i).name(end-12) == 'A' | allfiles(i).name(end-13) == 'A')
            a{1} = 'A';
        end
        if(allfiles(i).name(end-12) == 'B' | allfiles(i).name(end-13) == 'B')
            a{1} = 'B';
        end
        if(allfiles(i).name(end-12) == 'C' | allfiles(i).name(end-13) == 'C')
            a{1} = 'C';
        end
        if(allfiles(i).name(end-12) == 'D' | allfiles(i).name(end-13) == 'D')
            a{1} = 'D';
        end
        if(allfiles(i).name(end-12) == 'E' | allfiles(i).name(end-13) == 'E')
            a{1} = 'E';
        end
        if(allfiles(i).name(end-12) == 'F' | allfiles(i).name(end-13) == 'F')
            a{1} = 'F';
        end
        if(allfiles(i).name(end-12) == 'G'| allfiles(i).name(end-13) == 'G')
            a{1} = 'G';
        end
        if(allfiles(i).name(end-12) == 'H' | allfiles(i).name(end-13) == 'H')
            a{1} = 'H';
        end
        
         %This list gives the row:
        if(allfiles(i).name(end-10)=='_')
            Name2 = allfiles(i).name(end-12:end-11);
        else
            Name2 = allfiles(i).name(end-11:end-10);
        end    
        if(Name2 == '01')
            a{2} = '01';
        end
        if(Name2 == '02')
            a{2} = '02';
        end
        if(Name2 == '03')
            a{2} = '03';
        end
        if(Name2 == '04')
            a{2} = '04';
        end
        if(Name2 == '05')
            a{2} = '05';
        end
        if(Name2 == '06')
            a{2} = '06';
        end
        if(Name2 == '07')
            a{2} = '07';
        end
        if(Name2 == '08')
            a{2} = '08';
        end
        if(Name2 == '09')
            a{2} = '09';
        end
        if(Name2 == '10')
            a{2} = '10';
        end
        if(Name2 == '11')
            a{2} = '11';
        end
        if(Name2 == '12')
            a{2} = '12';
        end
        
        %Now we only need the sight number:
        if(allfiles(i).name(end-6)=='_')
            a{3} = allfiles(i).name(end-9:end-7);
        else
            a{3} = allfiles(i).name(end-8:end-7);
        end   
        a{4} = 'space';
        newfile = [a{1} a{2} a{3} a{4}];
         SynImageBinary = uint8(SynImageBinary);
         PSD95ImageBinary = uint8(PSD95ImageBinary);
         newFile1 = regexprep(newfile, 'space','NucleusSmall.png','ignorecase');
         newFile2 = regexprep(newfile, 'space','NucleusBig.png','ignorecase'); 
         newFile3 = regexprep(newfile, 'space','NeuriteSmall.png','ignorecase');
         newFile4 = regexprep(newfile, 'space','NeuriteBig.png','ignorecase');
         newFile5 = regexprep(newfile, 'space','Binary.png','ignorecase');
         newFile6 = regexprep(newfile, 'space','OligoSmall.png','ignorecase');
         newFile7 = regexprep(newfile, 'space','OligoBig.png','ignorecase');
         newFile8 = regexprep(newfile, 'space','AstroSmall.png','ignorecase');
         newFile9 = regexprep(newfile, 'space','AstroBig.png','ignorecase');
         imwrite(NucleusImage,[foldername1 '/' newFile1]);
         imwrite(NucleusImage,[foldername1 '/' newFile2]);
         imwrite(NeuriteImage,[foldername1 '/' newFile3]);
         imwrite(NeuriteImage,[foldername1 '/' newFile4]);
         imwrite(NeuriteArea,[foldername1 '/' newFile5]);
         imwrite(SynImageBinary,[foldername1 '/' newFile6]);
         imwrite(SynImageBinary,[foldername1 '/' newFile7]);
         imwrite(PSD95ImageBinary,[foldername1 '/' newFile8]);
         imwrite(PSD95ImageBinary,[foldername1 '/' newFile9]);
        
        
        NucleusM = sparse(zeros(size(SynImageBinary)));
        wellname = (newFile1(1:6));
        wellList{n} = wellname;
        csvHandler.CellPosMatrix(wellname) = NucleusM;
        ScalingFactor = 1;
        
        
        %Now we need to write the data:
        Summary{n+1,1} = regexprep(newfile, 'space','Name','ignorecase');
        SummaryName{1,1} = 'Image Name';
        SummaryName{1,2} = 'Neurite Mass';
        SummaryName{1,3} = 'Synaptophosin Area';
        SummaryName{1,4} = 'Synaptophosin Number';
        SummaryName{1,5} = 'PSD95 Area';
        SummaryName{1,6} = 'PSD95 Number';
        SummaryName{1,7} = 'Synaptophosin/NeuriteArea';
        SummaryName{1,8} = 'PSD95/NeuriteArea';
        
        
        Summary2(n+1,2) = NeuriteMass;
        Summary2(n+1,3) = SynArea ;
        Summary2(n+1,4) = SynNumber;
        Summary2(n+1,5) = PSD95Area;
        Summary2(n+1,6) = PSD95Number;
        Summary2(n+1,7) = SynNumber/NeuriteMass;
        Summary2(n+1,8) = PSD95Number/NeuriteMass;
        
        Summary2Name = 0;
        Summary3Name = 0;
        Summary3 = 0;
        Summary4 = 0;
        SummaryEnd = 0;
        SummaryEndTop =0;
        
     end    
     
    
    end
end   
   
     
     
     
      