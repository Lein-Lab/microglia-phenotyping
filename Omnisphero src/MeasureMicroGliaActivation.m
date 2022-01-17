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


% This algorithm is designed to cluster micro glia based on the morphology
% in resting and actiavtion states!
function [summary] = 	MeasureMicroGliaActivation (foldername, csvHandler, ExcelName, CD68,RAND)
foldernameOrig = foldername;
foldername = [foldername '/ConvertedCellomics'];
foldername1 = [foldername '/Boundary'];
foldername2 = [foldername '/Perimeter'];
mkdir(foldername1);
mkdir(foldername2);
allfiles = dir(foldername);
n = 0;
y=1;
g=1;
rng(1);
StepSize = 0.0025;
StepSize2 = 0.1;
SecondIntervall = 0.4;
Image = 0;
ImageCount = 0;
CD68M=0;
acc =1;
add =1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This parameters can be changed to change sample size of manual evaluation

RandomNumber = 300;
RealNumberofRandomCells = 100;
Activation = 0.06;
Inter = 0.3925;

% This parameters is responsible for distinguishing between cell body and
% processes. The higher the value, the more processes will be identified.

ThinningFactor = 4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Here we test, if transmission images are used!
%Check if transmission image:
%If only w1 exists, it has to be transmission!
ExtensionCheck = [foldernameOrig '/*_w2.tif'];
ExtensionCheck2 = [foldernameOrig '/*_w3.tif'];
CheckExtension = dir(ExtensionCheck);
CheckExtension2 = dir(ExtensionCheck2);
if(length(CheckExtension)==0 && length(CheckExtension2)==0 )
    Transmission =1;
    Extension = '_w1';
else
    Transmission = 0;
end

for(i=1:numel(allfiles))
    ind = strfind([foldername '/' allfiles(i).name],'NucleusBigWatershed.png');
    % Check for nucleus images
    if(numel(ind) > 0 )%
        n = n+1;
        ab=0;
        if(Transmission == 0)
            if(CD68==0)
                % Read in image files
                NucleusBinary = imread([foldername '/' allfiles(i).name]);
                newfile1 = regexprep(allfiles(i).name, 'NucleusBigWatershed','Binary','ignorecase');
                IBA1Binary = imread([foldername '/' newfile1]);
                newfile2 = regexprep(allfiles(i).name, 'NucleusBigWatershed','Mask','ignorecase');
                Mask = imread([foldername '/' newfile2]);
                Mask = logical(Mask);
                NucleusBinary = imclearborder(NucleusBinary);
                IBA1Binary = imclearborder(IBA1Binary);
                newfile3 = regexprep(allfiles(i).name, 'NucleusBigWatershed','AstroBig','ignorecase');
                try
                    CD68Img = imread([foldername '/' newfile3]);
                    if(max(max(CD68Img))>255)
                        CD68Img = uint8(double(CD68Img)/65535*255);
                    end
                    CD68M = 1;
                catch
                end
                newfile4 = regexprep(allfiles(i).name, 'NucleusBigWatershed','NeuriteBig','ignorecase');
                IBA1 = imread([foldername '/' newfile4]);
            else
                NucleusBinary = imread([foldername '/' allfiles(i).name]);
                newfile1 = regexprep(allfiles(i).name, 'NucleusBigWatershed','Binary','ignorecase');
                IBA1Binary = imread([foldername '/' newfile1]);
                newfile2 = regexprep(allfiles(i).name, 'NucleusBigWatershed','Mask','ignorecase');
                Mask = imread([foldername '/' newfile2]);
                newfile3 = regexprep(allfiles(i).name, 'NucleusBigWatershed','AstroBig','ignorecase');
                CD68Img = imread([foldername '/' newfile3]);
                if(max(max(CD68Img))>255)
                    CD68Img = uint8(double(CD68Img)/65535*255);
                end
                newfile4 = regexprep(allfiles(i).name, 'NucleusBigWatershed','NeuriteBig','ignorecase');
                IBA1 = imread([foldername '/' newfile4]);
                Mask = logical(Mask);
                NucleusBinary = imclearborder(NucleusBinary);
                IBA1Binary = imclearborder(IBA1Binary);
                CD68M = 1;
            end
        else
            NucleusBinary = imread([foldername '/' allfiles(i).name]);
            newfile1 = regexprep(allfiles(i).name, 'NucleusBigWatershed','Binary','ignorecase');
            IBA1Binary = imread([foldername '/' newfile1]);
            newfile2 = regexprep(allfiles(i).name, 'NucleusBigWatershed','Mask','ignorecase');
            Mask = imread([foldername '/' newfile2]);
            Mask = logical(Mask);
            NucleusBinary = imclearborder(NucleusBinary);
            IBA1Binary = imclearborder(IBA1Binary);
        end
        
        ImageName = regexprep(allfiles(i).name, 'NucleusBigWatershed.png','','ignorecase');
        ImageName2 = ImageName;
        %Get NucleusMatrix:
        NucleusBinary = logical(NucleusBinary);
        Centroids = regionprops(NucleusBinary,'Centroid');
        sdimensions = size(Centroids);
        NucleusMatrix = zeros(size(NucleusBinary));
        for(k=1:sdimensions)
            x_centroid = Centroids(k).Centroid(1);
            y_centroid = Centroids(k).Centroid(2);
            x_centroid = round(x_centroid);
            y_centroid = round(y_centroid);
            NucleusMatrix(y_centroid,x_centroid)= 1;
        end
        NucleusMatrix = logical(NucleusMatrix);
        
        
        Summary{1,1} = 'ImageName';
        Summary{1,2} = 'PercentFill';
        Summary{1,3} = 'Cellbody Area';
        Summary{1,4} = 'Process Length';
        Summary{1,5} = 'Ratio of Process Length and Cellbody Area';
        Summary{1,6} = 'Cell Number';
        Summary{1,7} = 'Area o f cell';
        if(CD68M==1)
            Summary{1,8} = 'IntensityCD68Raw';
            Summary{1,9} = 'Number of Endpoints';
            Summary{1,10} = 'Number of Branchingpoints';
            Summary{1,11} = 'Span Ratio';
            Summary{1,12} = 'Perimeter';
        else
            Summary{1,8} = 'Number of Endpoints';
            Summary{1,9} = 'Number of Branchingpoints';
            Summary{1,10} = 'Span Ratio';
            Summary{1,11} = 'Perimeter';
        end
        
        
        IBA1Binary = logical(IBA1Binary);
        SingleParticles = zeros(size(NucleusBinary));
        SingleParticlesCD68 = zeros(size(NucleusBinary));
        SingleParticlesCD68 = uint8(SingleParticlesCD68);
        if(CD68==0)
            IBA1Binary = bwmorph(IBA1Binary,'open');
            
            %We need to smoothen the image:
            
            %We will fill holes here originating from thresholding
            IntSummary = 0;
            IntSummarySingle=0;
            IntSummaryCluster=0;
            IBA1Binary = bwmorph(IBA1Binary,'close');
            Holes = imfill(IBA1Binary,'holes');
            Holes = Holes - IBA1Binary;
            HolesBig = bwareaopen(Holes,100);
            Holes = Holes - HolesBig;
            IBA1Binary = IBA1Binary + Holes;
            IBA1Binary = bwareaopen(IBA1Binary,250);
            IBA1BinaryLabel = bwlabel(IBA1Binary);
            NumberIBA1 = max(max(IBA1BinaryLabel));
            
            %Get exactly RealNumberofRandomCells!
            
            if(RAND==1)
                STOP = 0;
                Count  = 1;
                %Generate a random matrix in which each entry is only unique:
                RandomMatrix = [];
                for(j=1:30000)
                    if(STOP == 0)
                        RandomNumberOne = round((rand*NumberIBA1),0);
                        if(j==1)
                            if(RandomNumberOne>0)
                                RandomMatrix(j,1) = RandomNumberOne;
                            else
                                RandomMatrix(j,1) = 1;
                            end
                        else
                            CheckDouble = find(RandomMatrix == RandomNumberOne);
                            CheckDouble = length(CheckDouble);
                            if(CheckDouble==0 && RandomNumberOne>0)
                                Count = Count + 1;
                                RandomMatrix(Count,1) = RandomNumberOne;
                                
                                if(Count == NumberIBA1)
                                    STOP = 1;
                                end
                            end
                        end
                    end
                end
                
                
                IBA1BinaryRandom = zeros(size(IBA1Binary));
                CountCell = 0;
                STOP = 0;
                %Take only a cell with a centroid on it
                if(RandomNumber>NumberIBA1)
                    RandomNumber = NumberIBA1;
                end
                for(p=1:RandomNumber)
                    if(STOP == 0)
                        
                        IBA1BinaryRandomOne = IBA1BinaryLabel;
                        ii = IBA1BinaryRandomOne == RandomMatrix(p);
                        IBA1BinaryRandomOne(ii) = 50000;
                        IBA1BinaryRandomOne = IBA1BinaryRandomOne -49995;
                        IBA1BinaryRandomOne = uint8(IBA1BinaryRandomOne);
                        IBA1BinaryRandomOne = logical(IBA1BinaryRandomOne);
                        %Here we need to make sure, that the cells contain a
                        %centroid!
                        TestOver = IBA1BinaryRandomOne + NucleusMatrix;
                        TestOver = TestOver -1;
                        TestOver = uint8(TestOver);
                        TestOver = logical(TestOver);
                        if(sum(sum(TestOver))>0)
                            IBA1BinaryRandom = IBA1BinaryRandom + IBA1BinaryRandomOne;
                            CountCell = CountCell + 1;
                        end
                        if(CountCell == RealNumberofRandomCells)
                            %This ensures, that we will only get xxx cells!
                            STOP =1;
                        end
                    end
                    
                end
                IBA1Binary = IBA1BinaryRandom;
                IBA1BinaryLabel = bwlabel(IBA1Binary);
                NumberIBA1 = max(max(IBA1BinaryLabel));
                
            end
            
            
            
            %Here we will sperate big objects adjactend with one pixel:
            %If this is a manual evaluation, we will only keep objects with
            %centroids!
            IBA1BinaryCor = zeros(size(NucleusBinary));
            for(u=1:NumberIBA1)
                IBA1BinaryOne =  IBA1BinaryLabel;
                ii = IBA1BinaryOne == u;
                IBA1BinaryOne(ii) = 10*NumberIBA1;
                IBA1BinaryOne = IBA1BinaryOne - ((10*NumberIBA1)-5);
                IBA1BinaryOne = uint8(IBA1BinaryOne);
                IBA1BinaryOne = logical(IBA1BinaryOne);
                TestIBA1BinaryOne = bwmorph(IBA1BinaryOne,'open');
                TestIBA1BinaryOne = logical(TestIBA1BinaryOne);
                TestIBA1BinaryOneC = bwlabel(TestIBA1BinaryOne);
                TestIBA1BinaryOneC = max(max(TestIBA1BinaryOneC));
                if(TestIBA1BinaryOneC>1)
                    Count = 0;
                    AreaParticles = regionprops(TestIBA1BinaryOne,'Area');
                    NumberofParticles = length(AreaParticles);
                    for(t=1:NumberofParticles)
                        Area = AreaParticles(t).Area(1);
                        if(Area > 250)
                            Count = Count + 1;
                        end
                        
                    end
                    if(Count>1)
                        IBA1BinaryCor = IBA1BinaryCor +   TestIBA1BinaryOne;
                    else
                        IBA1BinaryCor = IBA1BinaryCor +   IBA1BinaryOne;
                    end
                else
                    IBA1BinaryCor = IBA1BinaryCor +   IBA1BinaryOne;
                end
            end
            IBA1Binary = IBA1BinaryCor;
            IBA1Binary = logical(IBA1Binary);
            
            IBA1BinaryLabel = bwlabel(IBA1Binary);
            NumberIBA1 = max(max(IBA1BinaryLabel));
            
            %Correct again, since splitting objects might give us objects than searched!
            if(RAND==1)
                %         STOP = 0;
                %         RandomMatrix = [];
                %         Count = 0;
                %         for(k=1:100000)
                %             if(STOP == 0)
                %             RandomNumberOne = round((rand*NumberIBA1),0);
                %             if(k==1)
                %                 RandomMatrix(k,1) = RandomNumberOne;
                %             else
                %                 CheckDouble = find(RandomMatrix == RandomNumberOne);
                %                 CheckDouble = length(CheckDouble);
                %                 if(CheckDouble==0)
                %                     Count = Count + 1;
                %                     RandomMatrix(Count,1) = RandomNumberOne;
                %
                %                     if(Count == NumberIBA1)
                %                         STOP = 1;
                %                     end
                %                 end
                %             end
                %             end
                %         end
                IBA1BinaryRandom = zeros(size(IBA1Binary));
                CountCell = 0;
                STOP = 0;
                for(p=1:NumberIBA1)
                    if(STOP == 0)
                        
                        IBA1BinaryRandomOne = IBA1BinaryLabel;
                        ii = IBA1BinaryRandomOne == p;
                        IBA1BinaryRandomOne(ii) = 50000;
                        IBA1BinaryRandomOne = IBA1BinaryRandomOne -49995;
                        IBA1BinaryRandomOne = uint8(IBA1BinaryRandomOne);
                        IBA1BinaryRandomOne = logical(IBA1BinaryRandomOne);
                        %Here we need to make sure, that the cells contain a
                        %centroid!
                        TestOver = IBA1BinaryRandomOne + NucleusMatrix;
                        TestOver = TestOver -1;
                        TestOver = uint8(TestOver);
                        TestOver = logical(TestOver);
                        if(sum(sum(TestOver))>0)
                            IBA1BinaryRandom = IBA1BinaryRandom + IBA1BinaryRandomOne;
                            CountCell = CountCell + 1;
                        end
                        if(CountCell == RealNumberofRandomCells)
                            %This ensures, that we will only get xxx cells!
                            STOP =1;
                        end
                    end
                    
                end
                IBA1Binary = IBA1BinaryRandom;
                IBA1BinaryLabel = bwlabel(IBA1Binary);
                NumberIBA1 = max(max(IBA1BinaryLabel));
                
            end
            
            a = 1;
            
            %    SeperateIBA11 = zeros(size(NucleusBinary));
            %    SeperateIBA12 = zeros(size(NucleusBinary));
            %    SeperateIBA13 = zeros(size(NucleusBinary));
            %    SeperateIBA14 = zeros(size(NucleusBinary));
            %    SeperateIBA15 = zeros(size(NucleusBinary));
            %    SeperateIBA16 = zeros(size(NucleusBinary));
            %    SeperateIBA17 = zeros(size(NucleusBinary));
            %    SeperateIBA18 = zeros(size(NucleusBinary));
            %    SeperateIBA19 = zeros(size(NucleusBinary));
            %    SeperateIBA110 = zeros(size(NucleusBinary));
            %    SeperateIBA111 = zeros(size(NucleusBinary));
            %    SeperateIBA112 = zeros(size(NucleusBinary));
            %    SeperateIBA113 = zeros(size(NucleusBinary));
            %    SeperateIBA114 = zeros(size(NucleusBinary));
            %    SeperateIBA115 = zeros(size(NucleusBinary));
            %    SeperateIBA116 = zeros(size(NucleusBinary));
            %    SeperateIBA117 = zeros(size(NucleusBinary));
            %    SeperateIBA118 = zeros(size(NucleusBinary));
            %    SeperateIBA119 = zeros(size(NucleusBinary));
            %    SeperateIBA120 = zeros(size(NucleusBinary));
            %Here we will also cout the number of objects within each cluster:
            %    Cluster_1 = 0;
            %    Cluster_2 = 0;
            %    Cluster_3 = 0;
            %    Cluster_4 = 0;
            %    Cluster_5 = 0;
            %    Cluster_6 = 0;
            %    Cluster_7 = 0;
            %    Cluster_8 = 0;
            %    Cluster_9 = 0;
            %    Cluster_10 = 0;
            %    Cluster_11 = 0;
            %    Cluster_12 = 0;
            %    Cluster_13 = 0;
            %    Cluster_14 = 0;
            %    Cluster_15 = 0;
            %    Cluster_16 = 0;
            %    Cluster_17 = 0;
            %    Cluster_18 = 0;
            %    Cluster_19 = 0;
            %    Cluster_20 = 0;
            %Here we will check on how many centroids of nuclei are located on a
            %given IBA1 particle
            Structureall = zeros(size(IBA1Binary));
            Activated = {};
            Activated_Cell = {};
            Intermediate = {};
            Intermediate_Cell = {};
            Resting = {};
            Resting_Cell = {};
            ArtefactCells = zeros(size(IBA1Binary));
            
            ac=0;
            as = 0;
            ae =0;
            af=0;
            IntermediateNumber =[];
            IntermediateType = [];
            ActivatedNumber =[];
            ActivatedType = [];
            RestingNumber =[];
            RestingType = [];
            
            
            
            for(u=1:NumberIBA1)
                % We will extract one particle of IBA1
                IBA1BinaryOne =  IBA1BinaryLabel;
                ii = IBA1BinaryOne == u;
                IBA1BinaryOne(ii) = 10*NumberIBA1;
                IBA1BinaryOne = IBA1BinaryOne - ((10*NumberIBA1)-5);
                IBA1BinaryOne = uint8(IBA1BinaryOne);
                IBA1BinaryOne = logical(IBA1BinaryOne);
                if(Transmission==0)
                    %Now we add the centroid matrix to the IBA1 particle and subtract one
                    %and count the number of remaining pixles, which will equal the
                    %number of nuclei within the particle!
                    TestIBA1One = IBA1BinaryOne + NucleusMatrix;
                    TestIBA1One  = TestIBA1One  -1 ;
                    TestIBA1One  =uint8(TestIBA1One );
                    TestIBA1One  = logical(TestIBA1One );
                    %Here we will inculde a rescue for those glia only overlapping
                    %partially with a second or third nucleus!
                    if(sum(sum(TestIBA1One ))>1)
                        testOver = IBA1BinaryOne + NucleusBinary;
                        testOver = testOver -1;
                        testOver = uint8(testOver);
                        testOver = logical(testOver);
                        AreaParticles = regionprops(testOver,'Area');
                        NParticles = length(AreaParticles);
                        Area = 0;
                        Count = 0;
                        for(b=1:NParticles)
                            Area(b,1) = AreaParticles(b).Area(1);
                            if(Area(b,1) > 199)
                                Count = Count + 1;
                            end
                        end
                        if(Count<2)
                            TestIBA1One = 1;
                        end
                        
                        
                    end
                    
                    if(sum(sum(TestIBA1One )) ==0)
                        AreaNucOver = NucleusBinary+IBA1BinaryOne;
                        AreaNucOver = AreaNucOver -1;
                        AreaNucOver = uint8(AreaNucOver);
                        AreaNucOver = logical(AreaNucOver);
                        %We area only checking for regions bigger than e.g. 200 and we
                        %will chekc if there is one or more area
                        AreaNucOver = bwareaopen(AreaNucOver,200);
                        NumberParticles = bwlabel(AreaNucOver);
                        NumberParticles = max(max(NumberParticles));
                        if(NumberParticles==1)
                            TestIBA1One = 1;
                        end
                        
                    end
                else
                    testOver = IBA1BinaryOne + NucleusBinary;
                    testOver = testOver -1;
                    testOver = uint8(testOver);
                    testOver = logical(testOver);
                    testOverL = bwlabel(testOver);
                    TestIBA1One = max(max(testOverL));
                end
                
                if(sum(sum(TestIBA1One ))>0)
                    if(sum(sum(TestIBA1One ))<2)
                        
                        Box = regionprops(IBA1BinaryOne,'BoundingBox');
                        Length = Box(1).BoundingBox(3);
                        Height = Box(1).BoundingBox(4);
                        AreaBox = Length*Height;
                        AreaBinary = sum(sum(IBA1BinaryOne));
                        PercentArea = AreaBinary/AreaBox;
                        IBA1Thin = bwmorph(IBA1BinaryOne,'thin',ThinningFactor);
                        %Now we will open the image to only obtain the cell body
                        %without processes:
                        IBA1Body = bwmorph(IBA1Thin,'open');
                        IBA1BodyBig = bwareaopen(IBA1Body,50);
                        IBA1BodySmall = IBA1Body - IBA1BodyBig;
                        IBA1BodyBig = bwmorph(IBA1BodyBig,'thicken',ThinningFactor);
                        IBA1BodyBig = bwmorph(IBA1BodyBig,'bridge');
                        IBA1Body =imfill(IBA1Body ,'holes');
                        IBA1BodyArea=sum(sum(IBA1BodyBig));
                        %Now we will also analyze the processis, by subtracting the
                        %Cell body from the pre thinned image
                        Skeleton = IBA1Thin - IBA1Body +IBA1BodySmall;
                        Skeleton = uint8(Skeleton);
                        Skeleton = logical(Skeleton);
                        Skeleton = bwmorph(Skeleton,'thin',inf);
                        Skeleton = bwareaopen(Skeleton,2);
                        Structure = Skeleton+IBA1BodyBig;
                        Structure=logical(Structure);
                        %Skeleton = bwmorph(Skeleton,'spur',5);
                        SkeletonLength = sum(sum(Skeleton));
                        Dummy = Skeleton + IBA1Body;
                        Endpoints = bwmorph(Dummy,'endpoints');
                        Endpoints = Endpoints - IBA1Body;
                        Endpoints = uint8(Endpoints);
                        Endpoints = logical(Endpoints);
                        Endpoints = sum(sum(Endpoints));
                        Branching = bwmorph(Skeleton,'branchpoints');
                        Branching = sum(sum(Branching));
                        ConvexHull = bwconvhull(IBA1BinaryOne,'objects');
                        majorAxis = regionprops(ConvexHull,'MajorAxisLength');
                        majorAxis = majorAxis.MajorAxisLength(1);
                        minorAxis = regionprops(ConvexHull,'MinorAxisLength');
                        minorAxis = minorAxis.MinorAxisLength(1);
                        SpanRatio = majorAxis/minorAxis;
                        %Compute edges:
                        Perimeter = regionprops(IBA1BinaryOne,'Perimeter');
                        Perimeter = Perimeter.Perimeter(1);
                        if(IBA1BodyArea>0)
                            y=y+1;
                            ab=ab+1;
                            ac=ac+1;
                            acc = acc+1;
                            Summary2{y,1} = ImageName;
                            Summary2_S{acc,1} = ImageName;
                            Summary3(y,2) = PercentArea;
                            Summary3(y,3) = IBA1BodyArea;
                            Summary3_S(acc,2) = PercentArea;
                            Summary3_S(acc,3) = IBA1BodyArea;
                            IntSummary(ab,1)= SkeletonLength/IBA1BodyArea;
                            IntSummarySingle(ac,1)= SkeletonLength/IBA1BodyArea;
                            Summary3(y,4) = SkeletonLength;
                            Summary3(y,5) = SkeletonLength/IBA1BodyArea;
                            Summary3(y,6) = 1;
                            Summary3(y,7) = sum(sum(IBA1BinaryOne));
                            Summary3_S(acc,4) = SkeletonLength;
                            Summary3_S(acc,5) = SkeletonLength/IBA1BodyArea;
                            Summary3_S(acc,6) = 1;
                            Summary3_S(acc,7) = sum(sum(IBA1BinaryOne));
                            if(CD68M==1)
                                MaskCD68 = imcomplement(IBA1BinaryOne)*255;
                                MaskCD68 = uint8(MaskCD68);
                                CD68Int = CD68Img - MaskCD68;
                                CD68Int = uint8(CD68Int);
                                CD68ImgS = CD68Int;
                                CD68Int = sum(sum(CD68Int));
                                Summary3(y,8) = sum(sum(CD68Int));
                                Summary3_S(acc,8) = sum(sum(CD68Int));
                            end
                            if(CD68M==0)
                                Summary3(y,8) =  Endpoints;
                                Summary3(y,9) = Branching;
                                Summary3(y,10) = SpanRatio;
                                Summary3(y,11) = Perimeter;
                                Summary3_S(acc,8) =  Endpoints;
                                Summary3_S(acc,9) = Branching;
                                Summary3_S(acc,10) = SpanRatio;
                                Summary3_S(acc,11) = Perimeter;
                            else
                                Summary3(y,9) =  Endpoints;
                                Summary3(y,10) = Branching;
                                Summary3(y,11) = SpanRatio;
                                Summary3(y,12) = Perimeter;
                                Summary3_S(acc,9) =  Endpoints;
                                Summary3_S(acc,10) = Branching;
                                Summary3_S(acc,11) = SpanRatio;
                                Summary3_S(acc,12) = Perimeter;
                            end
                            
                            
                            
                            
                            SingleParticles = SingleParticles+IBA1BinaryOne;
                            Structureall = Structureall+ Structure;
                            if (RAND == 1)
                                if((SkeletonLength/IBA1BodyArea)<=Activation)
                                    as = as + 1;
                                    BodyNew = bwmorph(Structure,'open');
                                    ProcessNew = Structure - BodyNew;
                                    BodyNew = edge(BodyNew,'Sobel');
                                    StructureNew = BodyNew + ProcessNew;
                                    BBox = regionprops(StructureNew,'BoundingBox');
                                    %Get Corrdinates:
                                    yStart = BBox(1).BoundingBox(1);
                                    xStart = BBox(1).BoundingBox(2);
                                    yEnd = yStart + BBox(1).BoundingBox(3);
                                    xEnd = xStart + BBox(1).BoundingBox(4);
                                    Activated_Coordinates{as,1} = yStart;
                                    Activated_Coordinates{as,2} = xStart;
                                    Activated_Coordinates{as,3} = yEnd;
                                    Activated_Coordinates{as,4} = xEnd;
                                    Activated_IBA1{as} = IBA1(xStart:xEnd,yStart:yEnd);
                                    StructureNew = StructureNew(xStart:xEnd,yStart:yEnd);
                                    CellIBA1 = IBA1(xStart:xEnd,yStart:yEnd);
                                    Activated_Cell{as} = CellIBA1;
                                    Activated{as} = StructureNew;
                                    ActivatedNumber(as,1) = SkeletonLength/IBA1BodyArea;
                                    ActivatedType{as,1} = 'Single';
                                    
                                end
                                if ((SkeletonLength/IBA1BodyArea)>Activation && (SkeletonLength/IBA1BodyArea)<=Inter)
                                    ae = ae +1;
                                    BodyNew = bwmorph(Structure,'open');
                                    ProcessNew = Structure - BodyNew;
                                    BodyNew = edge(BodyNew,'Sobel');
                                    StructureNew = BodyNew + ProcessNew;
                                    BBox = regionprops(StructureNew,'BoundingBox');
                                    %Get Corrdinates:
                                    yStart = BBox(1).BoundingBox(1);
                                    xStart = BBox(1).BoundingBox(2);
                                    yEnd = yStart + BBox(1).BoundingBox(3);
                                    xEnd = xStart + BBox(1).BoundingBox(4);
                                    Intermediate_IBA1{ae} = IBA1(xStart:xEnd,yStart:yEnd);
                                    Intermediate_Coordinates{ae,1} = yStart;
                                    Intermediate_Coordinates{ae,2} = xStart;
                                    Intermediate_Coordinates{ae,3} = yEnd;
                                    Intermediate_Coordinates{ae,4} = xEnd;
                                    
                                    StructureNew = StructureNew(xStart:xEnd,yStart:yEnd);
                                    CellIBA1 = IBA1(xStart:xEnd,yStart:yEnd);
                                    Intermediate_Cell{ae} = CellIBA1;
                                    Intermediate{ae} = StructureNew;
                                    IntermediateNumber(ae,1) = SkeletonLength/IBA1BodyArea;
                                    IntermediateType{ae,1} = 'Single';
                                    
                                end
                                if((SkeletonLength/IBA1BodyArea)>Inter)
                                    af = af +1;
                                    BodyNew = bwmorph(Structure,'open');
                                    ProcessNew = Structure - BodyNew;
                                    BodyNew = edge(BodyNew,'Sobel');
                                    StructureNew = BodyNew + ProcessNew;
                                    BBox = regionprops(StructureNew,'BoundingBox');
                                    %Get Corrdinates:
                                    yStart = BBox(1).BoundingBox(1);
                                    xStart = BBox(1).BoundingBox(2);
                                    yEnd = yStart + BBox(1).BoundingBox(3);
                                    xEnd = xStart + BBox(1).BoundingBox(4);
                                    CellIBA1 = IBA1(xStart:xEnd,yStart:yEnd);
                                    Resting_Cell{af} = CellIBA1;
                                    Resting_IBA1{af} = IBA1(xStart:xEnd,yStart:yEnd);
                                    StructureNew = StructureNew(xStart:xEnd,yStart:yEnd);
                                    Resting_Coordinates{af,1} = yStart;
                                    Resting_Coordinates{af,2} = xStart;
                                    Resting_Coordinates{af,3} = yEnd;
                                    Resting_Coordinates{af,4} = xEnd;
                                    Resting{af} = StructureNew;
                                    RestingNumber(af,1) = SkeletonLength/IBA1BodyArea;
                                    RestingType{af,1} = 'Single';
                                end
                            end
                            
                            if(CD68M ==1)
                                SingleParticlesCD68 = SingleParticlesCD68 + CD68ImgS;
                            end
                            
                        else
                            %Remove cells without cell body:
                            ArtefactCells = ArtefactCells + IBA1BinaryOne;
                        end
                        %              if((SkeletonLength/IBA1BodyArea)==0)
                        %                 SeperateIBA11 = SeperateIBA11 +  IBA1BinaryOne;
                        %                 Cluster_1 = Cluster_1 + 1;
                        %              end
                        %              if((SkeletonLength/IBA1BodyArea)<0.0025 && (SkeletonLength/IBA1BodyArea)>0)
                        %                 SeperateIBA12 = SeperateIBA12 +  IBA1BinaryOne;
                        %                 Cluster_2 = Cluster_2 + 1;
                        %              end
                        %              if((SkeletonLength/IBA1BodyArea)<0.005 && (SkeletonLength/IBA1BodyArea)>=0.0025)
                        %                 SeperateIBA13 = SeperateIBA13 +  IBA1BinaryOne;
                        %                 Cluster_3 = Cluster_3 + 1;
                        %              end
                        %               if((SkeletonLength/IBA1BodyArea)<0.0075 && (SkeletonLength/IBA1BodyArea)>=0.005)
                        %                 SeperateIBA14 = SeperateIBA14 +  IBA1BinaryOne;
                        %                 Cluster_4 = Cluster_4 + 1;
                        %              end
                        %               if((SkeletonLength/IBA1BodyArea)<0.01 && (SkeletonLength/IBA1BodyArea)>=0.0075)
                        %                 SeperateIBA15 = SeperateIBA15 +  IBA1BinaryOne;
                        %                 Cluster_5 = Cluster_5 + 1;
                        %               end
                        %              if((SkeletonLength/IBA1BodyArea)<0.0125 && (SkeletonLength/IBA1BodyArea)>=0.01)
                        %                 SeperateIBA16 = SeperateIBA16 +  IBA1BinaryOne;
                        %                 Cluster_6 = Cluster_6 + 1;
                        %              end
                        %              if((SkeletonLength/IBA1BodyArea)<0.015 && (SkeletonLength/IBA1BodyArea)>=0.0125)
                        %                 SeperateIBA17 = SeperateIBA17 +  IBA1BinaryOne;
                        %                 Cluster_7 = Cluster_7 + 1;
                        %              end
                        %              if((SkeletonLength/IBA1BodyArea)<0.0175 && (SkeletonLength/IBA1BodyArea)>=0.015)
                        %                 SeperateIBA18 = SeperateIBA18 +  IBA1BinaryOne;
                        %                 Cluster_8 = Cluster_8 + 1;
                        %              end
                        %              if((SkeletonLength/IBA1BodyArea)<0.45 && (SkeletonLength/IBA1BodyArea)>=0.4)
                        %                 SeperateIBA19 = SeperateIBA19 +  IBA1BinaryOne;
                        %                 Cluster_9 = Cluster_9 + 1;
                        %              end
                        %              if((SkeletonLength/IBA1BodyArea)<0.5 && (SkeletonLength/IBA1BodyArea)>=0.45)
                        %                 SeperateIBA110 = SeperateIBA110 +  IBA1BinaryOne;
                        %                 Cluster_10 = Cluster_10 + 1;
                        %              end
                        %              if((SkeletonLength/IBA1BodyArea)<0.55 && (SkeletonLength/IBA1BodyArea)>=0.5)
                        %                 SeperateIBA111 = SeperateIBA111 +  IBA1BinaryOne;
                        %                 Cluster_11 = Cluster_11 + 1;
                        %              end
                        %              if((SkeletonLength/IBA1BodyArea)<0.6 && (SkeletonLength/IBA1BodyArea)>=0.55)
                        %                 SeperateIBA112 = SeperateIBA112 +  IBA1BinaryOne;
                        %                 Cluster_12 = Cluster_12 + 1;
                        %              end
                        %              if((SkeletonLength/IBA1BodyArea)<0.65 && (SkeletonLength/IBA1BodyArea)>=0.6)
                        %                 SeperateIBA113 = SeperateIBA113 +  IBA1BinaryOne;
                        %                 Cluster_13 = Cluster_13 + 1;
                        %              end
                        %              if((SkeletonLength/IBA1BodyArea)<0.7 && (SkeletonLength/IBA1BodyArea)>=0.65)
                        %                 SeperateIBA114 = SeperateIBA114 +  IBA1BinaryOne;
                        %                 Cluster_14 = Cluster_14 + 1;
                        %              end
                        %              if((SkeletonLength/IBA1BodyArea)<0.75 && (SkeletonLength/IBA1BodyArea)>=0.7)
                        %                 SeperateIBA115 = SeperateIBA115 +  IBA1BinaryOne;
                        %                 Cluster_15 = Cluster_15 + 1;
                        %              end
                        %              if((SkeletonLength/IBA1BodyArea)<0.8 && (SkeletonLength/IBA1BodyArea)>=0.75)
                        %                 SeperateIBA116 = SeperateIBA116 +  IBA1BinaryOne;
                        %                 Cluster_16 = Cluster_16 + 1;
                        %              end
                        %              if((SkeletonLength/IBA1BodyArea)<0.85 && (SkeletonLength/IBA1BodyArea)>=0.8)
                        %                 SeperateIBA117 = SeperateIBA117 +  IBA1BinaryOne;
                        %                 Cluster_17 = Cluster_17 + 1;
                        %              end
                        %              if((SkeletonLength/IBA1BodyArea)<0.9 && (SkeletonLength/IBA1BodyArea)>=0.85)
                        %                 SeperateIBA118 = SeperateIBA118 +  IBA1BinaryOne;
                        %                 Cluster_18 = Cluster_18 + 1;
                        %              end
                        %              if((SkeletonLength/IBA1BodyArea)<0.95 && (SkeletonLength/IBA1BodyArea)>=0.9)
                        %                 SeperateIBA119 = SeperateIBA119 +  IBA1BinaryOne;
                        %                 Cluster_19 = Cluster_19 + 1;
                        %              end
                        %              if((SkeletonLength/IBA1BodyArea)>= 1)
                        %                 SeperateIBA120 = SeperateIBA120 +  IBA1BinaryOne;
                        %                 Cluster_20 = Cluster_20 + 1;
                        %              end
                    end
                end
                
            end
            
            %Now we will analyze the clusters:
            
        else
            SingleParticles = zeros(size(NucleusBinary));
        end
        ClusteredCells = IBA1Binary - SingleParticles - ArtefactCells;
        ClusteredCells = uint8(ClusteredCells);
        ClusteredCells = logical(ClusteredCells);
        if(Transmission ==0)
            ClusteredCells = bwareaopen(ClusteredCells,250);
        end
        
        ClusteredCellsLabeld = bwlabel(ClusteredCells);
        NumberClusters = max(max(ClusteredCellsLabeld));
        ad=0;
        ArtefactCells  = zeros(size(IBA1Binary));
        for(q=1:NumberClusters);
            IBA1BinaryClusterOne =  ClusteredCellsLabeld;
            ii = IBA1BinaryClusterOne == q;
            IBA1BinaryClusterOne(ii) = 10*NumberClusters;
            IBA1BinaryClusterOne = IBA1BinaryClusterOne - ((10*NumberClusters)-5);
            IBA1BinaryClusterOne = uint8(IBA1BinaryClusterOne);
            IBA1BinaryClusterOne = logical(IBA1BinaryClusterOne);
            
            Box = regionprops(IBA1BinaryClusterOne,'BoundingBox');
            Length = Box(1).BoundingBox(3);
            Height = Box(1).BoundingBox(4);
            AreaBox = Length*Height;
            AreaBinary = sum(sum(IBA1BinaryClusterOne));
            PercentArea = AreaBinary/AreaBox;
            IBA1ClusterThin = bwmorph(IBA1BinaryClusterOne,'thin',ThinningFactor);
            %Now we will open the image to only obtain the cell body
            %without processes:
            IBA1ClusterBody = bwmorph(IBA1ClusterThin,'open');
            IBA1ClusterBodyBig = bwareaopen(IBA1ClusterBody,50);
            IBA1ClusterBodySmall = IBA1ClusterBody - IBA1ClusterBodyBig;
            IBA1ClusterBodyBig = bwmorph(IBA1ClusterBodyBig,'thicken',ThinningFactor);
            IBA1ClusterBodyBig = bwmorph(IBA1ClusterBodyBig,'bridge');
            IBA1ClusterBodyBig = imfill(IBA1ClusterBodyBig ,'holes');
            IBA1ClusterBodyArea = sum(sum(IBA1ClusterBodyBig));
            %Now we will also analyze the processis, by subtracting the
            %Cell body from the pre thinned image
            Skeleton = IBA1ClusterThin - IBA1ClusterBody+IBA1ClusterBodySmall;
            Skeleton = uint8(Skeleton);
            Skeleton = logical(Skeleton);
            Skeleton = bwmorph(Skeleton,'thin',inf);
            Skeleton = bwareaopen(Skeleton,2);
            Structure = Skeleton+IBA1ClusterBodyBig;
            Structure=logical(Structure);
            %Skeleton = bwmorph(Skeleton,'spur',5);
            SkeletonLength = sum(sum(Skeleton));
            Dummy = Skeleton + IBA1ClusterBody;
            Endpoints = bwmorph(Dummy,'endpoints');
            Endpoints = Endpoints - IBA1ClusterBody;
            Endpoints = uint8(Endpoints);
            Endpoints = logical(Endpoints);
            Endpoints = sum(sum(Endpoints));
            Branching = bwmorph(Skeleton,'branchpoints');
            Branching = sum(sum(Branching));
            ConvexHull = bwconvhull(IBA1BinaryClusterOne,'objects');
            majorAxis = regionprops(ConvexHull,'MajorAxisLength');
            majorAxis = majorAxis.MajorAxisLength(1);
            minorAxis = regionprops(ConvexHull,'MinorAxisLength');
            minorAxis = minorAxis.MinorAxisLength(1);
            SpanRatio = majorAxis/minorAxis;
            %Compute edges:
            Perimeter = regionprops(IBA1BinaryClusterOne,'Perimeter');
            Perimeter = Perimeter.Perimeter(1);
            
            
            
            
            
            
            %Depending on the number of nuclei we will write the same cell
            %multiple times:
            NNucleus = IBA1BinaryClusterOne + NucleusMatrix;
            NNucleus = NNucleus -1;
            NNucleus = uint8(NNucleus);
            NNucleus = logical(NNucleus);
            NNucleus = sum(sum(NNucleus));
            if(CD68==0)
                if(NNucleus <2)
                    NNucleus = IBA1BinaryClusterOne + NucleusBinary;
                    NNucleus = NNucleus -1;
                    NNucleus = uint8(NNucleus);
                    NNucleus = logical(NNucleus);
                    NNucleus = sum(sum(NNucleus))/400;
                    NNucleus = round(NNucleus);
                end
            else
                NNucleus = 1;
            end
            if(Transmission==1)
                NNucleus = IBA1BinaryClusterOne + NucleusBinary;
                NNucleus = NNucleus -1;
                NNucleus = uint8(NNucleus);
                NNucleus = logical(NNucleus);
                NNucleusL = bwlabel(NNucleus);
                NNucleus = max(max(NNucleusL));
            end
            if(IBA1ClusterBodyArea>0)
                if(NNucleus==1)
                    y=y+1;
                    ab=ab+1;
                    ac=ac+1;
                    acc = acc+1;
                    Summary2{y,1} = ImageName;
                    Summary2_S{acc,1} = ImageName;
                    Summary3(y,2) = PercentArea;
                    Summary3(y,3) = IBA1ClusterBodyArea;
                    Summary3_S(acc,2) = PercentArea;
                    Summary3_S(acc,3) = IBA1ClusterBodyArea;
                    IntSummary(ab,1)= SkeletonLength/IBA1ClusterBodyArea;
                    IntSummarySingle(ac,1)= SkeletonLength/IBA1ClusterBodyArea;
                    Summary3(y,4) = SkeletonLength;
                    Summary3(y,5) = SkeletonLength/IBA1ClusterBodyArea;
                    Summary3(y,6) = 1;
                    Summary3(y,7) = sum(sum(IBA1BinaryOne));
                    Summary3_S(acc,4) = SkeletonLength;
                    Summary3_S(acc,5) = SkeletonLength/IBA1ClusterBodyArea;
                    Summary3_S(acc,6) = 1;
                    Summary3_S(acc,7) = sum(sum(IBA1BinaryOne));
                    if (RAND == 1)
                        if((SkeletonLength/IBA1ClusterBodyArea)<=Activation)
                            as = as + 1;
                            BodyNew = bwmorph(Structure,'open');
                            ProcessNew = Structure - BodyNew;
                            BodyNew = edge(BodyNew,'Sobel');
                            StructureNew = BodyNew + ProcessNew;
                            %Get Corrdinates:
                            yStart = BBox(1).BoundingBox(1);
                            xStart = BBox(1).BoundingBox(2);
                            yEnd = yStart + BBox(1).BoundingBox(3);
                            xEnd = xStart + BBox(1).BoundingBox(4);
                            Activated_Coordinates{as,1} = yStart;
                            Activated_Coordinates{as,2} = xStart;
                            Activated_Coordinates{as,3} = yEnd;
                            Activated_Coordinates{as,4} = xEnd;
                            Activated_IBA1{as} = IBA1(xStart:xEnd,yStart:yEnd);
                            CellIBA1 = IBA1(xStart:xEnd,yStart:yEnd);
                            Activated_Cell{as} = CellIBA1;
                            Activated{as} = StructureNew;
                            ActivatedNumber(as,1) = SkeletonLength/IBA1ClusterBodyArea;
                            ActivatedType{as,1} = 'Cluster';
                            
                        end
                        if ((SkeletonLength/IBA1ClusterBodyArea)>Activation && (SkeletonLength/IBA1BodyArea)<=Inter)
                            ae = ae + 1;
                            BodyNew = bwmorph(Structure,'open');
                            ProcessNew = Structure - BodyNew;
                            BodyNew = edge(BodyNew,'Sobel');
                            StructureNew = BodyNew + ProcessNew;
                            BBox = regionprops(StructureNew,'BoundingBox');
                            %Get Corrdinates:
                            yStart = BBox(1).BoundingBox(1);
                            xStart = BBox(1).BoundingBox(2);
                            yEnd = yStart + BBox(1).BoundingBox(3);
                            xEnd = xStart + BBox(1).BoundingBox(4);
                            CellIBA1 = IBA1(xStart:xEnd,yStart:yEnd);
                            Intermediate_Coordinates{ae,1} = yStart;
                            Intermediate_Coordinates{ae,2} = xStart;
                            Intermediate_Coordinates{ae,3} = yEnd;
                            Intermediate_Coordinates{ae,4} = xEnd;
                            Intermediate_Cell{ae} = CellIBA1;
                            Intermediate_IBA1{ae} = IBA1(xStart:xEnd,yStart:yEnd);
                            StructureNew = StructureNew(xStart:xEnd,yStart:yEnd);
                            Intermediate{ae} = StructureNew;
                            IntermediateNumber(ae,1) = SkeletonLength/IBA1ClusterBodyArea;
                            IntermediateType{ae,1} = 'Cluster';
                        end
                        if((SkeletonLength/IBA1ClusterBodyArea)>Inter)
                            af = af + 1;
                            BodyNew = bwmorph(Structure,'open');
                            ProcessNew = Structure - BodyNew;
                            BodyNew = edge(BodyNew,'Sobel');
                            StructureNew = BodyNew + ProcessNew;
                            BBox = regionprops(StructureNew,'BoundingBox');
                            %Get Corrdinates:
                            yStart = BBox(1).BoundingBox(1);
                            xStart = BBox(1).BoundingBox(2);
                            yEnd = yStart + BBox(1).BoundingBox(3);
                            xEnd = xStart + BBox(1).BoundingBox(4);
                            Resting_Coordinates{af,1} = yStart;
                            Resting_Coordinates{af,2} = xStart;
                            Resting_Coordinates{af,3} = yEnd;
                            Resting_Coordinates{af,4} = xEnd;
                            Resting_IBA1{af} = IBA1(xStart:xEnd,yStart:yEnd);
                            CellIBA1 = IBA1(xStart:xEnd,yStart:yEnd);
                            Resting_Cell{af} = CellIBA1;
                            StructureNew = StructureNew(xStart:xEnd,yStart:yEnd);
                            Resting{af} = StructureNew;
                            RestingNumber(af,1) = SkeletonLength/IBA1ClusterBodyArea;
                            RestingType{af,1} = 'Cluster';
                        end
                    end
                    if(CD68M==1)
                        MaskCD68 = imcomplement(IBA1BinaryOne)*255;
                        MaskCD68 = uint8(MaskCD68);
                        CD68Int = CD68Img - MaskCD68;
                        CD68Int = uint8(CD68Int);
                        CD68ImgS = CD68Int;
                        CD68Int = sum(sum(CD68Int));
                        Summary3(y,8) = sum(sum(CD68Int));
                        Summary3_S(acc,8) = sum(sum(CD68Int));
                    end
                    if(CD68M==0)
                        Summary3(y,8) =  Endpoints;
                        Summary3(y,9) = Branching;
                        Summary3(y,10) = SpanRatio;
                        Summary3(y,11) = Perimeter;
                        Summary3_S(acc,8) =  Endpoints;
                        Summary3_S(acc,9) = Branching;
                        Summary3_S(acc,10) = SpanRatio;
                        Summary3_S(acc,11) = Perimeter;
                    else
                        Summary3(y,9) =  Endpoints;
                        Summary3(y,10) = Branching;
                        Summary3(y,11) = SpanRatio;
                        Summary3(y,12) = Perimeter;
                        Summary3_S(acc,9) =  Endpoints;
                        Summary3_S(acc,10) = Branching;
                        Summary3_S(acc,11) = SpanRatio;
                        Summary3_S(acc,12) = Perimeter;
                    end
                end
                
                
                if(NNucleus>1)
                    for(z=1:NNucleus)
                        ab=ab+1;
                        y=y+1;
                        Summary2{y,1} = ImageName;
                        Summary3(y,2) = PercentArea;
                        Summary3(y,3) = IBA1ClusterBodyArea;
                        Summary3(y,4) = SkeletonLength;
                        Summary3(y,5) = SkeletonLength/IBA1ClusterBodyArea;
                        if(NNucleus>1)
                            ad=ad+1;
                            add=add+1;
                            Summary2_C{add,1} = ImageName;
                            Summary3_C(add,2) = PercentArea;
                            Summary3_C(add,3) = IBA1ClusterBodyArea;
                            Summary3_C(add,4) = SkeletonLength;
                            Summary3_C(add,5) = SkeletonLength/IBA1ClusterBodyArea;
                            Summary3_C(add,6) = NNucleus;
                            Summary3_C(add,7) = sum(sum(IBA1BinaryClusterOne))/NNucleus;
                            IntSummaryCluster(ad,1)= SkeletonLength/IBA1ClusterBodyArea;
                        end
                        IntSummary(ab,1)=SkeletonLength/IBA1ClusterBodyArea;
                        Summary3(y,6) = NNucleus;
                        Summary3(y,7) = sum(sum(IBA1BinaryClusterOne))/NNucleus;
                        
                        if(CD68M==1)
                            MaskCD68 = imcomplement(IBA1BinaryClusterOne)*255;
                            MaskCD68 = uint8(MaskCD68);
                            CD68Int = CD68Img - MaskCD68;
                            CD68Int = uint8(CD68Int);
                            CD68ImgS = CD68Int;
                            CD68Int = sum(sum(CD68Int));
                            Summary3(y,8) = sum(sum(CD68Int))/NNucleus;
                            if(NNucleus>1)
                                Summary3_C(add,8) = sum(sum(CD68Int))/NNucleus;
                            end
                        end
                        if(CD68M==0)
                            Summary3(y,8) =  Endpoints/NNucleus;
                            Summary3(y,9) = Branching/NNucleus;
                            Summary3(y,10) = SpanRatio;
                            Summary3(y,11) = Perimeter/NNucleus;
                            if(NNucleus>1)
                                Summary3_C(add,8) =  Endpoints/NNucleus;
                                Summary3_C(add,9) = Branching/NNucleus;
                                Summary3_C(add,10) = SpanRatio;
                                Summary3_C(add,11) = Perimeter/NNucleus;
                            end
                        else
                            
                            Summary3(y,9) =  Endpoints/NNucleus;
                            Summary3(y,10) = Branching/NNucleus;
                            Summary3(y,11) = SpanRatio;
                            Summary3(y,12) = Perimeter/NNucleus;
                            if(NNucleus>1)
                                Summary3_C(add,9) =  Endpoints/NNucleus;
                                Summary3_C(add,10) = Branching/NNucleus;
                                Summary3_C(add,11) = SpanRatio;
                                Summary3_C(add,12) = Perimeter/NNucleus;
                            end
                        end
                    end
                    %              if(IBA1ClusterBodyArea==0)
                    %                  Summary3(y,5) = 0;
                    %                  IBA1ClusterBodyArea = 1;
                    %              end
                    
                end
            else
                ArtefactCells = ArtefactCells + IBA1BinaryClusterOne;
                
            end
            if(NNucleus>0 && IBA1ClusterBodyArea >0 )
                SingleParticles = SingleParticles+IBA1BinaryClusterOne;
                Structureall = Structureall+ Structure;
                if (RAND == 1)
                    if((SkeletonLength/IBA1ClusterBodyArea)<=Activation)
                        as = as + 1;
                        BodyNew = bwmorph(Structure,'open');
                        ProcessNew = Structure - BodyNew;
                        BodyNew = edge(BodyNew,'Sobel');
                        StructureNew = BodyNew + ProcessNew;
                        BBox = regionprops(StructureNew,'BoundingBox');
                        %Get Corrdinates:
                        yStart = BBox(1).BoundingBox(1);
                        xStart = BBox(1).BoundingBox(2);
                        yEnd = yStart + BBox(1).BoundingBox(3);
                        xEnd = xStart + BBox(1).BoundingBox(4);
                        Activated_Coordinates{as,1} = yStart;
                        Activated_Coordinates{as,2} = xStart;
                        Activated_Coordinates{as,3} = yEnd;
                        Activated_Coordinates{as,4} = xEnd;
                        CellIBA1 = IBA1(xStart:xEnd,yStart:yEnd);
                        Activated_IBA1{as} = IBA1(xStart:xEnd,yStart:yEnd);
                        Activated_Cell{as} = CellIBA1;
                        StructureNew = StructureNew(xStart:xEnd,yStart:yEnd);
                        Activated{as} = StructureNew;
                        if(IBA1ClusterBodyArea>0)
                            ActivatedNumber(as,1) = SkeletonLength/IBA1ClusterBodyArea;
                        else
                            ActivatedNumber(as,1) = 1;
                        end
                        ActivatedType{as,1} = 'Single';
                        
                    end
                    if ((SkeletonLength/IBA1ClusterBodyArea)>Activation && (SkeletonLength/IBA1ClusterBodyArea)<=Inter)
                        ae = ae+1;
                        BodyNew = bwmorph(Structure,'open');
                        ProcessNew = Structure - BodyNew;
                        BodyNew = edge(BodyNew,'Sobel');
                        StructureNew = BodyNew + ProcessNew;
                        BBox = regionprops(StructureNew,'BoundingBox');
                        %Get Corrdinates:
                        yStart = BBox(1).BoundingBox(1);
                        xStart = BBox(1).BoundingBox(2);
                        yEnd = yStart + BBox(1).BoundingBox(3);
                        xEnd = xStart + BBox(1).BoundingBox(4);
                        Intermediate_Coordinates{ae,1} = yStart;
                        Intermediate_Coordinates{ae,2} = xStart;
                        Intermediate_Coordinates{ae,3} = yEnd;
                        Intermediate_Coordinates{ae,4} = xEnd;
                        Intermediate_IBA1{ae} = IBA1(xStart:xEnd,yStart:yEnd);
                        CellIBA1 = IBA1(xStart:xEnd,yStart:yEnd);
                        Intermediate_Cell{ae} = CellIBA1;
                        StructureNew = StructureNew(xStart:xEnd,yStart:yEnd);
                        Intermediate{ae} = StructureNew;
                        if(IBA1ClusterBodyArea > 0)
                            IntermediateNumber(ae,1) = SkeletonLength/IBA1ClusterBodyArea;
                        else
                            IntermediateNumber(ae,1) = 1;
                        end
                        IntermediateType{ae,1} = 'Cluster';
                    end
                    if((SkeletonLength/IBA1ClusterBodyArea)>Inter)
                        af = af + 1;
                        BodyNew = bwmorph(Structure,'open');
                        ProcessNew = Structure - BodyNew;
                        BodyNew = edge(BodyNew,'Sobel');
                        StructureNew = BodyNew + ProcessNew;
                        BBox = regionprops(StructureNew,'BoundingBox');
                        %Get Corrdinates:
                        yStart = BBox(1).BoundingBox(1);
                        xStart = BBox(1).BoundingBox(2);
                        yEnd = yStart + BBox(1).BoundingBox(3);
                        xEnd = xStart + BBox(1).BoundingBox(4);
                        Resting_Coordinates{af,1} = yStart;
                        Resting_Coordinates{af,2} = xStart;
                        Resting_Coordinates{af,3} = yEnd;
                        Resting_Coordinates{af,4} = xEnd;
                        Resting_IBA1{af} = IBA1(xStart:xEnd,yStart:yEnd);
                        CellIBA1 = IBA1(xStart:xEnd,yStart:yEnd);
                        Resting_Cell{af} = CellIBA1;
                        StructureNew = StructureNew(xStart:xEnd,yStart:yEnd);
                        Resting{af} = StructureNew;
                        if(IBA1ClusterBodyArea >0)
                            RestingNumber(af,1) = SkeletonLength/IBA1ClusterBodyArea;
                        else
                            RestingNumber(af,1) = 1;
                        end
                        RestingType{af,1} = 'Cluster';
                    end
                end
                if(CD68M ==1)
                    SingleParticlesCD68 = SingleParticlesCD68 + CD68ImgS;
                end
            end
            
            
        end
        
        
        
        %Now we will analyze Table One for features:
        
        %   if(Image ==0)
        %   IntSummary = Summary3(2:end,5);
        %   ImageCount = ImageCount + length(IntSummary);
        %   Image =1;
        %   else
        %   IntSummary = Summary3(ImageCount+2:end,5);
        %   ImageCount = ImageCount + length(IntSummary);
        %   end
        
        
        TotalCellNumber = length( IntSummary);
        Cluster_1 = sum(IntSummary==0);
        Cluster_1S = sum(IntSummarySingle==0);
        Cluster_1C = sum(IntSummaryCluster==0);
        if(length(IntSummary)==1 && IntSummary(1)==0);
            Cluster(2) =0;
        else
            Cluster(2) = Cluster_1;
        end
        if(length(IntSummarySingle==1) && IntSummarySingle(1)==0)
            Cluster_S(2) = 0;
        else
            Cluster_S(2) = Cluster_1S;
        end
        if(length(IntSummaryCluster==1) && IntSummaryCluster(1)==0)
            Cluster_C(2) = 0;
        else
            Cluster_C(2) = Cluster_1C;
        end
        %Now we make Clusters in a loop:
        if(n==1)
            Cluster_Name = [];
            Cluster_Name{1} = ' Image Name';
            Cluster_Name{2} = ' Ratio = 0';
        end
        for(h=1:160)
            Large = StepSize*h;
            if(h==1)
                Low = 0;
            else
                Low = StepSize*(h-1);
            end
            if(h==1)
                Cluster(h+2) = length(IntSummary(IntSummary<Large & IntSummary>Low ));
                Cluster_S(h+2) = length(IntSummarySingle(IntSummarySingle<Large & IntSummarySingle>Low ));
                Cluster_C(h+2) = length(IntSummaryCluster(IntSummaryCluster<Large & IntSummaryCluster>Low ));
                if(n==1)
                    Cluster_Name{h+2} = [num2str(Low) '<Ratio<' num2str(Large)];
                end
            else
                Cluster(h+2) = length(IntSummary(IntSummary<Large & IntSummary>=Low ));
                Cluster_S(h+2) = length(IntSummarySingle(IntSummarySingle<Large & IntSummarySingle>=Low ));
                Cluster_C(h+2) = length(IntSummaryCluster(IntSummaryCluster<Large & IntSummaryCluster>=Low ));
                if(n==1)
                    Cluster_Name{h+2} = [num2str(Low) '<=Ratio<' num2str(Large)];
                end
            end
            
        end
        %till 1:
        for(r=1:6)
            h = h+1;
            Large = SecondIntervall + (StepSize2*r);
            Low =   SecondIntervall +   (StepSize2*(r-1));
            Cluster(h+2) = length(IntSummary(IntSummary<Large & IntSummary>=Low ));
            Cluster_S(h+2) = length(IntSummarySingle(IntSummarySingle<Large & IntSummarySingle>=Low ));
            Cluster_C(h+2) = length(IntSummaryCluster(IntSummaryCluster<Large & IntSummaryCluster>=Low ));
            Cluster_Name{h+2} = [num2str(Low) '<=Ratio<' num2str(Large)];
        end
        h=h+1;
        Cluster_Name{h+2} = ' Ratio >= 1';
        Cluster(h+2) = length(IntSummary(IntSummary>=1));
        Cluster_S(h+2) = length(IntSummarySingle(IntSummarySingle>=1));
        Cluster_C(h+2) = length(IntSummaryCluster(IntSummaryCluster>=1));
        %Now we can calculate the total number of cells:
        TotalCellNumber = sum(Cluster);
        TotalCellNumberSingle = sum(Cluster_S);
        TotalCellNumberCluster = sum(Cluster_C);
        %Now we can also calculate % directly from cluster:
        Cluster_Percent = Cluster/TotalCellNumber *100;
        Overview1 = Cluster_Name;
        Overview2{n+1,1} = ImageName2;
        Overview3(n+1,:) = Cluster;
        Overview4(n+1,:) = Cluster_S;
        Overview5(n+1,:) = Cluster_C;
        OverviewPercent1 = Cluster_Name;
        OverviewPercent3(n+1,:) = Cluster_Percent;
        OverviewPercent2{n+1,1}=ImageName2;
        
        
        
        AreaSignal = sum(sum(SingleParticles));
        %2) The area of IBA1:
        AreaIBA1 = sum(sum(IBA1Binary));
        
        d = n+1;
        OverviewTotal{1,1} = 'ImageName';
        OverviewTotal{1,2} = 'Total number of cells';
        OverviewTotal{1,3} = 'Total Area';
        OverviewTotal{1,4} = 'ROI Area';
        OverviewTotal{1,5} = 'IBA1 Area';
        OverviewTotal{1,7} = 'Total number of single cells';
        OverviewTotal{1,8} = 'Total number of clustered cells';
        OverviewTotal1{d,1} = ImageName;
        OverviewTotal2(d,2) = TotalCellNumber;
        OverviewTotal2(d,3) = AreaSignal;
        OverviewTotal2(d,4) = sum(sum(Mask));
        OverviewTotal2(d,5) = AreaIBA1;
        OverviewTotal2(d,7) = TotalCellNumberSingle;
        OverviewTotal2(d,8) = TotalCellNumberCluster;
        
        %    ExcelName = foldername(1:end-19);
        %   [token,remain] = strtok(ExcelName, '\');
        %    while(length(remain) > 0)
        %         [token,remain] = strtok(remain, '\');
        %    end
        if (RAND == 1)
            %             for i = 1:length(Activated)
            %                 Filled = imfill(Activated,'holes');
            %                 Filled = bwmorph(Filled,'open');
            %                 Filled = bwmorph(Filled,'thin',1);
            %             ActivatedLabel = bwlabel(Activated);
            %             ActivatedCorr = zeros(size(Activated));
            %             for(o=1:max(max(ActivatedLabel)))
            %                 ActivatedSingle = ActivatedLabel;
            %                 ii = ActivatedSingle == o;
            %                 ActivatedSingle(ii) =255;
            %                 ActivatedSingle = ActivatedSingle -250;
            %                 ActivatedSingle = uint8(ActivatedSingle);
            %                 ActivatedSingle = logical(ActivatedSingle);
            %                 Overlay = ActivatedSingle + Filled;
            %                 Overlay = Overlay -1;
            %                 Overlay = uint8(Overlay);
            %                 if(sum(sum(ActivatedSingle))>sum(sum(Overlay)))
            %                     ActivatedCorr = ActivatedCorr + ActivatedSingle;
            %                 end
            %             end
            
            %             Activated = logical(ActivatedCorr);
            
            
            %             Filled = imfill(Intermediate,'holes');
            %             Filled = bwmorph(Filled,'open');
            %             Filled = bwmorph(Filled,'thin',1);
            %             IntermediateLabel = bwlabel(Intermediate);
            %             IntermediateCorr = zeros(size(Intermediate));
            %             for(o=1:max(max(IntermediateLabel)))
            %                 IntermediateSingle = IntermediateLabel;
            %                 ii = IntermediateSingle == o;
            %                 IntermediateSingle(ii) =255;
            %                 IntermediateSingle = IntermediateSingle -250;
            %                 IntermediateSingle = uint8(IntermediateSingle);
            %                 IntermediateSingle = logical(IntermediateSingle);
            %                 Overlay = IntermediateSingle + Filled;
            %                 Overlay = Overlay -1;
            %                 Overlay = uint8(Overlay);
            %                 if(sum(sum(IntermediateSingle))>sum(sum(Overlay)))
            %                     IntermediateCorr = IntermediateCorr + IntermediateSingle;
            %                 end
            %             end
            %
            %             Intermediate = logical(IntermediateCorr);
            
            %             Filled = imfill(Resting,'holes');
            %             Filled = bwmorph(Filled,'open');
            %             Filled = bwmorph(Filled,'thin',1);
            %             RestingLabel = bwlabel(Resting);
            %             RestingCorr = zeros(size(Resting));
            %             for(o=1:max(max(RestingLabel)))
            %                 RestingSingle = RestingLabel;
            %                 ii = RestingSingle == o;
            %                 RestingSingle(ii) =255;
            %                 RestingSingle = RestingSingle -250;
            %                 RestingSingle = uint8(RestingSingle);
            %                 RestingSingle = logical(RestingSingle);
            %                 Overlay = RestingSingle + Filled;
            %                 Overlay = Overlay -1;
            %                 Overlay = uint8(Overlay);
            %                 if(sum(sum(RestingSingle))>sum(sum(Overlay)))
            %                     RestingCorr = RestingCorr + RestingSingle;
            %                 end
            %             end
            %
            %             Resting = logical(RestingCorr);
            
            
            
            Count = 0;
            %            BBoxActivated = zeros(size(SingleParticles));
            %            BBoxIntermediate = zeros(size(SingleParticles));
            %            BBoxResting = zeros(size(SingleParticles));
            %Activated = logical(Activated);
            %Intermediate = logical(Intermediate);
            %Resting = logical(Resting);
            %ActivatedLabel = bwlabel(Activated);
            %IntermediateLabel = bwlabel(Intermediate);
            %RestingLabel = bwlabel(Resting);
            %            NumberofCellsActivated = max(max(ActivatedLabel));
            %            NumberofCellsIntermediate = max(max(IntermediateLabel));
            %            NumberofCellsResting = max(max(RestingLabel));
            NumberofCellsActivated = length(Activated);
            NumberofCellsIntermediate = length(Intermediate);
            NumberofCellsResting = length(Resting);
            SurroundingActivated = zeros(size(IBA1Binary));
            SurroundingIntermediate = zeros(size(IBA1Binary));
            SurroundingResting = zeros(size(IBA1Binary));
            CutActivated = uint8(zeros(size(IBA1Binary)));
            CutActivated = cat(3,CutActivated ,CutActivated ,CutActivated );
            CutIntermediate = uint8(zeros(size(CutActivated)));
            CutResting = uint8(zeros(size(CutActivated)));
            for (u=1:NumberofCellsActivated)
                ActivatedOne = Activated{u};
                ActivatedOne_IBA1 = Activated_Cell{u};
                %Now we will draw a bounding box arround the croped image;
                Box = ones(size(ActivatedOne));
                Box = bwmorph(Box,'thin',2);
                Box = imcomplement(Box);
                %Reconstruct image montage:
                SurroundingActivated(Activated_Coordinates{u,2}: Activated_Coordinates{u,4},Activated_Coordinates{u,1}: Activated_Coordinates{u,3}) = Box;
                
                
                R = ActivatedOne_IBA1 + uint8(ActivatedOne)*255;
                G = ActivatedOne_IBA1 ;
                B = ActivatedOne_IBA1 ;
                
                Cell = cat(3, R, G, B);
                
                CutActivated(Activated_Coordinates{u,2}: Activated_Coordinates{u,4},Activated_Coordinates{u,1}: Activated_Coordinates{u,3},:) = Cell;
                
                Count = Count + 1;
                
                
                
                
                %Generate overlapps:
                R = ActivatedOne_IBA1+ uint8(Box)*255;
                G = ActivatedOne_IBA1 ;
                B = ActivatedOne_IBA1 ;
                
                RGBActivatedB = cat(3, R, G, B);
                
                R = ActivatedOne_IBA1 + uint8(ActivatedOne)*255;
                G = ActivatedOne_IBA1 ;
                B = ActivatedOne_IBA1 ;
                
                RGBActivatedM = cat(3, R, G, B);
                
                
                
                
                %Write images out
                try
                    prt{1} = num2str(ActivatedNumber(u,1));
                catch
                    STOP=1;
                end
                %prt{2} = ActivatedType{u,1};
                prt{2} = num2str(Count)
                String = [ '_' prt{2} '_'  prt{1}];
                newFile1 = regexprep(allfiles(i).name, 'NucleusBigWatershed',String,'ignorecase');
                imwrite(RGBActivatedB,[foldername1 '/' newFile1]);
                newFile2 = regexprep(allfiles(i).name, 'NucleusBigWatershed',String,'ignorecase');
                imwrite(RGBActivatedM,[foldername2 '/' newFile2]);
            end
            for (u=1:NumberofCellsIntermediate)
                IntermediateOne = Intermediate{u};
                IntermediateOne_IBA1 = Intermediate_Cell{u};
                %Now we will draw a bounding box arround the croped image;
                Box = ones(size(IntermediateOne));
                Box = bwmorph(Box,'thin',2);
                Box = imcomplement(Box);
                
                
                SurroundingIntermediate(Intermediate_Coordinates{u,2}: Intermediate_Coordinates{u,4},Intermediate_Coordinates{u,1}: Intermediate_Coordinates{u,3}) = Box;
                
                R = IntermediateOne_IBA1 ;
                G = IntermediateOne_IBA1 + uint8(IntermediateOne)*255;
                B = IntermediateOne_IBA1 ;
                
                Cell = cat(3, R, G, B);
                
                CutIntermediate(Intermediate_Coordinates{u,2}: Intermediate_Coordinates{u,4},Intermediate_Coordinates{u,1}: Intermediate_Coordinates{u,3},:) = Cell;
                
                
                Count = Count + 1;
                
                
                
                
                %Generate overlapps:
                R = IntermediateOne_IBA1;
                G = IntermediateOne_IBA1 + uint8(Box)*255;
                B = IntermediateOne_IBA1 ;
                
                RGBIntermediateB = cat(3, R, G, B);
                
                R = IntermediateOne_IBA1 ;
                G = IntermediateOne_IBA1 + uint8(IntermediateOne)*255;
                B = IntermediateOne_IBA1 ;
                
                RGBIntermediateM = cat(3, R, G, B);
                
                
                
                
                %Write images out
                try
                    prt{1} = num2str(IntermediateNumber(u,1));
                catch
                    STOP=1;
                end
                %prt{2} = IntermediateType{u,1};
                prt{2} = num2str(Count)
                String = [ '_' prt{2} '_'  prt{1}];
                newFile1 = regexprep(allfiles(i).name, 'NucleusBigWatershed',String,'ignorecase');
                imwrite(RGBIntermediateB,[foldername1 '/' newFile1]);
                newFile2 = regexprep(allfiles(i).name, 'NucleusBigWatershed',String,'ignorecase');
                imwrite(RGBIntermediateM,[foldername2 '/' newFile2]);
            end
            for (u=1:NumberofCellsResting)
                RestingOne = Resting{u};
                RestingOne_IBA1 = Resting_Cell{u};
                %Now we will draw a bounding box arround the croped image;
                Box = ones(size(RestingOne));
                Box = bwmorph(Box,'thin',2);
                Box = imcomplement(Box);
                
                SurroundingResting(Resting_Coordinates{u,2}: Resting_Coordinates{u,4},Resting_Coordinates{u,1}: Resting_Coordinates{u,3}) = Box;
                Count = Count + 1;
                
                R = RestingOne_IBA1 ;
                G = RestingOne_IBA1 ;
                B = RestingOne_IBA1 + uint8(RestingOne)*255;
                
                Cell = cat(3, R, G, B);
                
                CutResting(Resting_Coordinates{u,2}: Resting_Coordinates{u,4},Resting_Coordinates{u,1}: Resting_Coordinates{u,3},:) = Cell;
                
                
                %Generate overlapps:
                R = RestingOne_IBA1;
                G = RestingOne_IBA1 ;
                B = RestingOne_IBA1 + uint8(Box)*255;
                
                RGBRestingB = cat(3, R, G, B);
                
                R = RestingOne_IBA1 ;
                G = RestingOne_IBA1 ;
                B = RestingOne_IBA1 + uint8(RestingOne)*255;
                
                RGBRestingM = cat(3, R, G, B);
                
                
                
                
                %Write images out
                try
                    prt{1} = num2str(RestingNumber(u,1));
                catch
                    STOP=1;
                end
                %prt{2} = RestingType{u,1};
                prt{2} = num2str(Count)
                String = [ '_' prt{2} '_'  prt{1}];
                newFile1 = regexprep(allfiles(i).name, 'NucleusBigWatershed',String,'ignorecase');
                imwrite(RGBRestingB,[foldername1 '/' newFile1]);
                newFile2 = regexprep(allfiles(i).name, 'NucleusBigWatershed',String,'ignorecase');
                imwrite(RGBRestingM,[foldername2 '/' newFile2]);
            end
            %Get Only edges
            %        MaskInverse = imcomplement(logical(BBoxActivated+BBoxIntermediate+BBoxResting));
            %        MaskInverse = uint8(MaskInverse)*255;
            %        CutIba1 = IBA1 - MaskInverse;
            %        SurroundingActivated = edge(BBoxActivated,'Canny');
            %        SurroundingActivated = bwmorph(SurroundingActivated,'thicken',2);
            %        SurroundingIntermediate = edge(BBoxIntermediate,'Canny');
            %        SurroundingIntermediate = bwmorph(SurroundingIntermediate,'thicken',2);
            %        SurroundingResting = edge(BBoxResting,'Canny');
            %        SurroundingResting = bwmorph(SurroundingResting,'thicken',2);
            
            %Now we generate RGB overlay images:
            R = IBA1 + (uint8(SurroundingActivated)*255);
            G = IBA1 + (uint8(SurroundingIntermediate)*255);
            B = IBA1 + (uint8(SurroundingResting)*255);
            
            RGBBox = cat(3, R, G, B);
            
            RGBMask = CutResting + CutActivated + CutIntermediate;
            
            %        R = CutIba1 + (uint8(Activated)*255);
            %        G = CutIba1 + (uint8(Intermediate)*255);
            %        B = CutIba1 + (uint8(Resting)*255);
            %
            %        RGBMask = cat(3, R, G, B);
            
            newFile1 = regexprep(allfiles(i).name, 'NucleusBigWatershed','BoxOverlay','ignorecase');
            imwrite(RGBBox,[foldername '/' newFile1]);
            newFile2 = regexprep(allfiles(i).name, 'NucleusBigWatershed','BoundaryOverlay','ignorecase');
            imwrite(RGBMask,[foldername '/' newFile2]);
            
            
        end
        
        
        
        SingleParticles_Thin = bwmorph(SingleParticles,'thin',ThinningFactor);
        SingleParticles_Body = bwmorph(SingleParticles_Thin,'open');
        SingleParticles_BodyBig = bwareaopen(SingleParticles_Body,50);
        SingleParticles_BodySmall = SingleParticles_Body- SingleParticles_BodyBig;
        SingleParticles_BodyBig = bwmorph(SingleParticles_BodyBig,'thicken',ThinningFactor);
        SingleParticles_BodyBig = bwmorph(SingleParticles_BodyBig,'bridge');
        SingleParticles_Body = imfill(SingleParticles_Body,'holes');
        SingleParticles_Process = SingleParticles_Thin - SingleParticles_Body+SingleParticles_BodySmall;
        SingleParticles_Process = uint8(SingleParticles_Process);
        SingleParticles_Process = logical(SingleParticles_Process);
        SingleParticles_Process = bwmorph(SingleParticles_Process,'thin',inf);
        SingleParticles_Process = bwareaopen(SingleParticles_Process,2);
        SingleParticles_ProcessBody = SingleParticles_Process+SingleParticles_BodyBig;
        SingleParticles_ProcessBody = logical(SingleParticles_ProcessBody);
        SingleParticles_ProcessBody = uint8(SingleParticles_ProcessBody);
        if(CD68==0)
            if(size(SingleParticles_ProcessBody,1)<3000&&size(SingleParticles_ProcessBody,2)<3000)
                SingleParticles_ProcessBodySmall = imresize(SingleParticles_ProcessBody,0.5);
            else
                SingleParticles_ProcessBodySmall = imresize(SingleParticles_ProcessBody,0.1);
            end
            SingleParticles_ProcessBodySmall = uint8(SingleParticles_ProcessBodySmall);
        else
            SingleParticles_ProcessBodySmall =SingleParticles_ProcessBody;
            SingleParticles_ProcessBodySmall =uint8(SingleParticles_ProcessBodySmall);
        end
        %Here we add some additional endpoints for CD68:
        if(CD68M==1)
            %1) The intensity of CD68:
            MaskCD68 = imcomplement(IBA1Binary)*255;
            MaskCD68 = uint8(MaskCD68);
            CD68Int = CD68Img - MaskCD68;
            CD68Int = uint8(CD68Int);
            CD68Int = sum(sum(CD68Int));
            
            OverviewTotal2(d,6) = CD68Int;
            OverviewTotal{1,6} = 'CD68 Intensity';
            
        end
        
        
        imwrite(NucleusBinary,[foldername '/' allfiles(i).name]);
        newFile1 = regexprep(allfiles(i).name, 'NucleusBigWatershed','Binary','ignorecase');
        imwrite(SingleParticles,[foldername '/' newFile1]);
        newFile2 = regexprep(allfiles(i).name, 'NucleusBigWatershed','OligoBig','ignorecase');
        imwrite(SingleParticles_ProcessBody,[foldername '/' newFile2]);
        if(CD68==0)
            newFile3 = regexprep(allfiles(i).name, 'NucleusBigWatershed','OligoSmall','ignorecase');
            imwrite(SingleParticles_ProcessBodySmall,[foldername '/' newFile3]);
        else
            newFile3 = regexprep(allfiles(i).name, 'NucleusBigWatershed','OligoSmall','ignorecase');
            imwrite(SingleParticles_ProcessBody,[foldername '/' newFile3]);
            Summary{1,1} = 'ImageName';
            Summary{1,2} = 'PercentFill';
            Summary{1,3} = 'Cellbody Area';
            Summary{1,4} = 'Process Length';
            Summary{1,5} = 'Ratio of Process Length and Cellbody Area';
        end
        summary = 0;
    end
    
end


%Here we write results in Excel files:
try
    sheet = 1;
    filename = ExcelName;
    xlswrite(filename,Summary3,sheet);
    xlswrite(filename,Summary2,sheet);
    xlswrite(filename,Summary,sheet);
catch
end
try
    sheet = 2;
    xlswrite(filename,Overview3,sheet);
    xlswrite(filename,Overview2,sheet);
    xlswrite(filename,Overview1,sheet);
catch
end
try
    sheet = 3;
    xlswrite(filename,OverviewPercent3,sheet);
    xlswrite(filename,OverviewPercent2,sheet);
    xlswrite(filename,OverviewPercent1,sheet);
catch
end
try
    sheet = 4;
    xlswrite(filename,OverviewTotal2,sheet);
    xlswrite(filename,OverviewTotal1,sheet);
    xlswrite(filename,OverviewTotal,sheet);
catch
end
try
    sheet = 5;
    filename = ExcelName;
    xlswrite(filename,Summary3_S,sheet);
    xlswrite(filename,Summary2_S,sheet);
    xlswrite(filename,Summary,sheet);
catch
end
try
    sheet = 6;
    filename = ExcelName;
    xlswrite(filename,Summary3_C,sheet);
    xlswrite(filename,Summary2_C,sheet);
    xlswrite(filename,Summary,sheet);
catch
end
try
    sheet = 7;
    xlswrite(filename,Overview4,sheet);
    xlswrite(filename,Overview2,sheet);
    xlswrite(filename,Overview1,sheet);
catch
end
try
    sheet = 8;
    xlswrite(filename,Overview5,sheet);
    xlswrite(filename,Overview2,sheet);
    xlswrite(filename,Overview1,sheet);
catch
end
end
