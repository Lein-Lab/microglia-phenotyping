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

%This function is designed to generate a matlab save file for Omnisphero,
%solely containing images and NucleusMatrix!
function [wellList, foldername, ScalingFactor, csvHandler,   Th] = MicrogliaMorphology (csvHandler)
foldername = uigetdir;
foldername1 = [foldername '/ConvertedCellomics'];
mkdir(foldername1);
allfiles = dir(foldername);
n = 0;
wellList = cell(0);
st = 0;
%Max size confocal images:
x_Size = 6929;
y_Size = 4622;
%Check if transmission image:
%If only w1 exists, it has to be transmission!
ExtensionCheck = [foldername '/*_w3.tif'];
CheckExtension = dir(ExtensionCheck);
if(length(CheckExtension)==0)
    Transmission =1;
    Extension = '_w1';
else
    Transmission = 0;
end
if(Transmission ==1)
else
    
    %Check if confocal!
    Question = questdlg('Do you have confocal images?')
    if(strcmp(Question,'Yes'))
        Confocal = 1;
    else
        Confocal = 0;
    end
    %Check if we have CD68!
    Question = questdlg('Do you want to analyze CD68?')
    if(strcmp(Question,'No'))
        CD68 = 0;
        Extension = '_w1';
        CD68M = 1;
        %We will also need a loop to generate the maximum intensity stacks of
        %the different images and to save those as the
    else
        CD68 = 1;
        Question2 = questdlg('Are you working with image stacks?')
        if(strcmp(Question2,'No'))
            Extension = '_w1';
            CD68M = 1;
        else
            Extension = '_c1';
            CD68M = 0;
            %We will also need a loop to generate the maximum intensity stacks of
            %the different images and to save those as pngs
            for(i=1:numel(allfiles))
                ind = strfind([foldername '/' allfiles(i).name],'.tif');
                ind2 = strfind([foldername '/' allfiles(i).name],'_10');
                if(numel(ind)>0)
                    st = st + 1;
                    RGB = imread([foldername '/' allfiles(i).name]);
                    NucleusPicSingle = RGB(:,:,1);
                    IBA1Single = RGB(:,:,2);
                    CD68Single = RGB(:,:,3);
                    %Now we will create 3 stacks to generate the z-projection:
                    
                    if(numel(ind2)==1)
                        NucleusStack(:,:,st) = NucleusPicSingle;
                        IBA1Stack(:,:,st) = IBA1Single;
                        CD68Stack(:,:,st) = CD68Single;
                        %Now we need to generate the max projection image:
                        NucleusZ = max(NucleusStack,[],3);
                        IBA1Z = max(IBA1Stack,[],3);
                        CD68Z = max(CD68Stack,[],3);
                        newFile = regexprep(allfiles(i).name, '_10','_c1','ignorecase');
                        imwrite(NucleusZ,[foldername '/' newFile]);
                        newFile1 = regexprep(allfiles(i).name, '_10','_c2','ignorecase');
                        imwrite(IBA1Z,[foldername '/' newFile1]);
                        newFile2 = regexprep(allfiles(i).name, '_10','_c3','ignorecase');
                        imwrite(CD68Z,[foldername '/' newFile2]);
                        %Now we need to reset the parameters:
                        st = 0;
                        NucleusStack  = zeros(size(NucleusZ));
                        NucleusStack  = uint16(NucleusStack);
                        IBA1Stack = zeros(size(NucleusZ));
                        IBA1Stack = uint16(IBA1Stack);
                        CD68Stack = zeros(size(NucleusZ));
                        CD68Stack = uint16(CD68Stack);
                    else
                        NucleusStack(:,:,st) = NucleusPicSingle;
                        IBA1Stack(:,:,st) = IBA1Single;
                        CD68Stack(:,:,st) = CD68Single;
                    end
                    
                end
            end
        end
    end
end
allfiles = dir(foldername);

for(i=1:numel(allfiles))
    ind = strfind([foldername '/' allfiles(i).name],Extension);
    % Check for nucleus images
    if(numel(ind) > 0 )%
        n = n+1;
        % Read in image files
        if(Transmission==1)
            IBA1 = imread([foldername '/' allfiles(i).name]);
        else
            if(CD68==0)
                Nucleus = imread([foldername '/' allfiles(i).name]);
                newFileIBA1 = regexprep(allfiles(i).name, '_w1','_w3','ignorecase');
                IBA1 = imread([foldername '/' newFileIBA1]);
                ImageName = allfiles(i).name(1:end-7);
                if(strcmp(class(IBA1),'uint8')==0)
                    Nucleus = double(Nucleus)./4095;
                    Nucleus = uint8(Nucleus.*255);
                    IBA1 = double(IBA1)./4095;
                    IBA1 = uint8(IBA1.*255);
                end
                if(size(Nucleus,1)<3000 && size(Nucleus,2)<3000)
                    NucleusSmall = imresize(Nucleus,0.5);
                    IBA1Small = imresize(IBA1,0.5);
                else
                    NucleusSmall = imresize(Nucleus,0.1);
                    IBA1Small = imresize(IBA1,0.1);
                end
                
                try
                    newFileCD68 = regexprep(allfiles(i).name, '_w1','_w2','ignorecase');
                    CD68Img = imread([foldername '/' newFileCD68]);
                    if(strcmp(class(IBA1),'uint8')==0)
                        CD68Img = double(CD68Img)./4095;
                        CD68Img = uint8(CD68Img.*255);
                    end
                    if(size(Nucleus,1)<3000 && size(Nucleus,2)<3000)
                        CD68ImgSmall = imresize(CD68Img,0.5);
                    else
                        CD68ImgSmall = imresize(CD68Img,0.1);
                    end
                catch
                end
            else
                if(CD68M == 0)
                    Nucleus = imread([foldername '/' allfiles(i).name]);
                    newFileIBA1 = regexprep(allfiles(i).name, '_c1','_c2','ignorecase');
                    newFileCD68 = regexprep(allfiles(i).name, '_c1','_c3','ignorecase');
                    IBA1 = imread([foldername '/' newFileIBA1]);
                    CD68Img = imread([foldername '/' newFileCD68]);
                    ImageName = allfiles(i).name(1:end-7);
                    if(strcmp(class(IBA1),'uint8')==0)
                        Nucleus = double(Nucleus)./8191;
                        Nucleus = uint8(Nucleus.*255);
                    end
                    NucleusSmall = Nucleus;
                    if(strcmp(class(IBA1),'uint8')==0)
                        IBA1 = double(IBA1)./8191;
                        IBA1 = uint8(IBA1.*255);
                    end
                    IBA1Small = IBA1;
                    if(strcmp(class(IBA1),'uint8')==0)
                        CD68Img = double(CD68Img)./8191;
                        CD68Img = uint8(CD68Img.*255);
                    end
                    CD68ImgSmall = CD68Img;
                else
                    Nucleus = imread([foldername '/' allfiles(i).name]);
                    newFileIBA1 = regexprep(allfiles(i).name, '_w1','_w3','ignorecase');
                    IBA1 = imread([foldername '/' newFileIBA1]);
                    ImageName = allfiles(i).name(1:end-7);
                    if(strcmp(class(IBA1),'uint8')==0)
                        Nucleus = double(Nucleus)./4095;
                        Nucleus = uint8(Nucleus.*255);
                    end
                    
                    if(strcmp(class(IBA1),'uint8')==0)
                        IBA1 = double(IBA1)./4095;
                        IBA1 = uint8(IBA1.*255);
                    end
                    if(size(Nucleus,1)<3000 && size(Nucleus,2)<3000)
                        IBA1Small = imresize(IBA1,0.5);
                        NucleusSmall = imresize(Nucleus,0.5);
                    else
                        
                        IBA1Small = imresize(IBA1,0.1);
                        NucleusSmall = imresize(Nucleus,0.1);
                    end
                    newFileCD68 = regexprep(allfiles(i).name, '_w1','_w2','ignorecase');
                    CD68Img = imread([foldername '/' newFileCD68]);
                    if(strcmp(class(IBA1),'uint8')==0)
                        CD68Img = double(CD68Img)./4095;
                        CD68Img = uint8(CD68Img.*255);
                    end
                    if(size(Nucleus,1)<3000 && size(Nucleus,2)<3000)
                        CD68ImgSmall = imresize(CD68Img,0.5);
                    else
                        CD68ImgSmall = imresize(CD68Img,0.1);
                    end
                end
            end
        end
        % Now we preprocess the images to obtain binary images
        if(Transmission == 1)
            IBA1 = imcomplement(IBA1);
            IBA1 = im2uint8(IBA1);
            %Trying to reduce background with rolling ball method
            se = strel('ball',10,10);
            backgroundBall = imopen(IBA1,se);
            IBA1Cor = IBA1 - backgroundBall;
            %We will threshold twice:
            [lehisto x] = imhist(IBA1Cor);
            level = triangle_th(lehisto,256);
            HardThresh = im2bw(IBA1Cor,level);
            SoftThresh = im2bw(IBA1Cor,level*0.75);
            %Repair with edge detection:
            Edge1 = edge(IBA1Cor,'log');
            Edge2 = edge(IBA1Cor,'Sobel');
            Edge3 = edge(IBA1Cor,'Roberts');
            Edge4 = edge(IBA1Cor,'zerocross');
            Edge5 = edge(IBA1Cor,'prewitt');
            Edge = Edge1 + Edge4 + Edge3 +Edge2+Edge5;;
            Edge = logical(Edge);
            SoftThresh = SoftThresh + Edge;
            SoftThresh = logical(SoftThresh);
            %Now fill holes, but only small one and compare!
            Holes = imfill(SoftThresh,'holes');
            Holes = Holes - SoftThresh;
            HolesBig = bwareaopen(Holes,50);
            Holes = Holes- HolesBig;
            Holes = logical(Holes);
            SoftThresh = SoftThresh + Holes;
            SoftThresh = logical(SoftThresh);
            SoftThresh = bwmorph(SoftThresh,'bridge');
            Holes = imfill(SoftThresh,'holes');
            Holes = Holes - SoftThresh;
            HolesBig = bwareaopen(Holes,50);
            Holes = Holes- HolesBig;
            SoftThresh = SoftThresh + Holes;
            SoftThresh = logical(SoftThresh);
            SoftThresh = bwmorph(SoftThresh,'spur',50);
            %Now we keep only structures, which have pixels from HardThresh
            %First delete small objects to improve speed:
            SoftThresh = bwmorph(SoftThresh,'bridge');
            %Close small holes:
            Holes = imfill(SoftThresh,'holes');
            HolesBig = bwareaopen(Holes,10);
            Holes = Holes- HolesBig;
            SoftThresh = SoftThresh + Holes;
            SoftThresh = bwareaopen(SoftThresh,250);
            %         SoftThreshKeep = zeros(size(IBA1Cor));
            %         SoftThreshL = bwlabel(SoftThresh);
            %         SoftThreshN = max(max(SoftThreshL));
            %         for(y=1:SoftThreshN)
            %             SoftThreshOne = SoftThreshL;
            %             ii = SoftThreshOne ==y;
            %             SoftThreshOne(ii) = 1000000;
            %             SoftThreshOne = SoftThreshOne - 999995;
            %             SoftThreshOne = uint8(SoftThreshOne);
            %             SoftThreshOne = logical(SoftThreshOne);
            %             TestOver = logical(SoftThreshOne) + logical(HardThresh);
            %             TestOver = TestOver - 1;
            %             TestOver = uint8(TestOver);
            %             TestOver = logical(TestOver);
            %             if(sum(sum(TestOver))>0)
            %                 SoftThreshKeep = SoftThreshKeep + SoftThreshOne;
            %             end
            %         end
            SoftThreshKeep = logical(SoftThresh);
            
            %Now we try to get the cell somata!
            %We will do this with thinning!
            CellSoma = bwmorph(SoftThreshKeep,'thin',4);
            CellSoma = bwmorph(CellSoma,'open');
            CellSoma = bwareaopen(CellSoma,10);
            CellSoma = bwmorph(CellSoma,'thicken',4);
            CellSoma = bwmorph(CellSoma,'bridge');
            CellSoma = imfill(CellSoma,'holes');
            
            %We will use the CellSoma image to get a cut nucleus Image like in
            %Cortical Neurons!
            CellSomaMask = imcomplement(CellSoma);
            CellSomaMask = uint8(CellSomaMask)*255;
            NucleusImage = IBA1-CellSomaMask;
            NucleusImage = uint8(NucleusImage);
            Nucleus = NucleusImage;
            
            %Get Centroids:
            Centroids = regionprops(CellSoma,'Centroid');
            sdimensions = size(Centroids);
            NucleusMatrix = zeros(size(CellSoma));
            for(k=1:sdimensions)
                x_centroid = Centroids(k).Centroid(1);
                y_centroid = Centroids(k).Centroid(2);
                x_centroid = round(x_centroid);
                y_centroid = round(y_centroid);
                NucleusMatrix(y_centroid,x_centroid)= 1;
            end
            NucleusM = logical(NucleusMatrix);
            
            %Write images:
            CellSoma = uint8(CellSoma)*255;
            token=allfiles(i).name(1:end-4);
            a{2} = token;
            a{3} = 'space';
            newfile = [a{2} a{3}];
            newFile1 = regexprep(newfile, 'space','NucleusBig.png','ignorecase');
            imwrite(NucleusImage,[foldername1 '/' newFile1]);
            if(size(Nucleus,1)<3000 && size(Nucleus,2)<3000)
                NucleusImageSmall = imresize(CellSoma,0.5);
            else
                NucleusImageSmall = imresize(CellSoma,0.1);
            end
            newFile2 = regexprep(newfile, 'space','NucleusSmall.png','ignorecase');
            imwrite(NucleusImageSmall,[foldername1 '/' newFile2]);
            IBA1 = imcomplement(IBA1);
            newFile3 = regexprep(newfile, 'space','NeuriteBig.png','ignorecase');
            imwrite(IBA1,[foldername1 '/' newFile3]);
            if(size(Nucleus,1)<3000 && size(Nucleus,2)<3000)
                IBA1Small = imresize(IBA1,0.5);
            else
                IBA1Small = imresize(IBA1,0.1);
            end
            newFile4 = regexprep(newfile, 'space','NeuriteSmall.png','ignorecase');
            imwrite(IBA1Small,[foldername1 '/' newFile4]);
            newFile5 = regexprep(newfile, 'space','Binary.png','ignorecase');
            imwrite(SoftThreshKeep,[foldername1 '/' newFile5]);
            newFile6 = regexprep(newfile, 'space','NucleusBigWatershed.png','ignorecase');
            imwrite(CellSoma,[foldername1 '/' newFile6]);
        else
            %1) Nuclei:
            % We will try a fixed triangle threshold:
            [lehisto x] = imhist(Nucleus);
            if(CD68==0)
                level = triangle_th(lehisto,256);
                NucleusBinary=im2bw(Nucleus,level);
            else
                if(CD68M == 0)
                    level = isodata(Nucleus);
                else
                    level = triangle_th(lehisto,256);
                end
                NucleusBinary=im2bw(Nucleus,level);
            end
            %Filter of small particles
            NucleusBinary = bwareaopen(NucleusBinary,50);
            %Close small holes
            Holes = imfill(NucleusBinary,'holes');
            Holes = Holes - NucleusBinary;
            HolesBig = bwareaopen(Holes,5);
            Holes = Holes - HolesBig;
            NucleusBinary = NucleusBinary + Holes;
            
            %1) Generate centroids:
            NucleusBinary = uint8(NucleusBinary);
            NucleusBinary = logical(NucleusBinary);
            %NucleusBinary = imclearborder(NucleusBinary);
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
            NucleusM = logical(NucleusMatrix);
            %2) IBA1
            % We will try a rolling ball first!
            if(CD68==0)
                if(Confocal ==1)
                    level = graythresh(IBA1);
                    if(level<0.02)
                        level = 0.025;
                    end
                    IBA1Binary = im2bw(IBA1,level);
                    IBA1Binary = bwareaopen(IBA1Binary,5);
                    windowSize = 2;
                    kernel = ones(windowSize) / windowSize ^ 2;
                    blurryImage = conv2(single(IBA1Binary), kernel, 'same');
                    binaryImage = blurryImage > 0.5;
                    IBA1Binary = imclose(binaryImage,strel('disk',2));
                    level = graythresh(Nucleus);
                    if(level<0.03)
                        level = 0.03;
                    end
                    NucleusBinary = im2bw(Nucleus,level);
                    NucleusBinary = bwareaopen(NucleusBinary,5);
                    NucleusBinary= imclose(NucleusBinary,strel('disk',2));
                    Holes = imfill(NucleusBinary,'holes');
                    Holes = Holes - NucleusBinary;
                    HolesBig = bwareaopen(Holes,50);
                    Holes = Holes - HolesBig;
                    NucleusBinary = NucleusBinary + Holes;
                    NucleusBinary = bwareaopen(NucleusBinary,25);
                else
                    SE = strel('ball',20,20);
                    BackgroundBallBig = imopen(IBA1,SE);
                    IBA1BallBig =  IBA1 - BackgroundBallBig;
                    SE2 = strel('ball',5,5);
                    BackgroundBall = imopen(IBA1,SE2);
                    IBA1Ball =  IBA1 - BackgroundBall;
                    [lehisto x] = imhist(BackgroundBall);
                    level = triangle_th(lehisto,256);
		    if(level*1.5<1)
                    	BackgroundBinary = im2bw(BackgroundBall,level*1.5);
	            else
                        BackgroundBinary = im2bw(BackgroundBall,level);
                    end
                    [lehisto x] = imhist(IBA1BallBig);
                    level = triangle_th(lehisto,256);
                    IBA1Binary=im2bw(IBA1BallBig,level*2);
                    [lehisto x] = imhist(IBA1Ball);
                    level = triangle_th(lehisto,256);
                    IBA1BinarySmall=im2bw(IBA1Ball,level);
                    %Adding up images to close gaps in cell body and to maintain cell processis:
                    IBA1Binary = IBA1Binary + IBA1BinarySmall+BackgroundBinary;
                    IBA1Binary = uint8(IBA1Binary);
                    IBA1Binary = logical(IBA1Binary);
                    Holes = imfill(IBA1Binary,'holes');
                    Holes = Holes - IBA1Binary;
                    HolesBig = bwareaopen(Holes,5);
                    Holes = Holes - HolesBig;
                    IBA1Binary = IBA1Binary + Holes;
                    IBA1Binary = imclearborder(IBA1Binary);
                end
            else
                if(CD68M == 0)
                    level = isodata(IBA1) ;
                    IBA1Binary = im2bw(IBA1,level/2);
                    IBA1Binary = uint8(IBA1Binary);
                    IBA1Binary = logical(IBA1Binary);
                    Holes = imfill(IBA1Binary,'holes');
                    Holes = Holes - IBA1Binary;
                    HolesBig = bwareaopen(Holes,200);
                    Holes = Holes - HolesBig;
                    IBA1Binary = IBA1Binary + Holes;
                    IBA1Binary = bwareaopen(IBA1Binary,250);
                else
                    if(Confocal ==1)
                        level = graythresh(IBA1);
                        IBA1Binary = im2bw(IBA1,level);
                        IBA1Binary = bwareaopen(IBA1Binary,5);
                        windowSize = 2;
                        kernel = ones(windowSize) / windowSize ^ 2;
                        blurryImage = conv2(single(IBA1Binary), kernel, 'same');
                        binaryImage = blurryImage > 0.5;
                        IBA1Binary = imclose(binaryImage,strel('disk',2));
                    else
                        SE = strel('ball',20,20);
                        BackgroundBallBig = imopen(IBA1,SE);
                        IBA1BallBig =  IBA1 - BackgroundBallBig;
                        SE2 = strel('ball',5,5);
                        BackgroundBall = imopen(IBA1,SE2);
                        IBA1Ball =  IBA1 - BackgroundBall;
                        [lehisto x] = imhist(BackgroundBall);
                        level = triangle_th(lehisto,256);
                        level_Back = level;
			if(level*1.5<1)
                        	BackgroundBinary = im2bw(BackgroundBall,level*1.5);
		        else
				BackgroundBinary = im2bw(BackgroundBall,level);
			end
                        [lehisto x] = imhist(IBA1BallBig);
                        level = triangle_th(lehisto,256);
                        IBA1Binary=im2bw(IBA1BallBig,level*2);
                        [lehisto x] = imhist(IBA1Ball);
                        level = triangle_th(lehisto,256);
                        IBA1BinarySmall=im2bw(IBA1Ball,level);
                        %Adding up images to close gaps in cell body and to maintain cell processis:
                        if(level_Back<0.1)
                            IBA1Binary = IBA1Binary + IBA1BinarySmall;
                        else
                            IBA1Binary = IBA1Binary + IBA1BinarySmall+BackgroundBinary;
                        end
                        IBA1Binary = uint8(IBA1Binary);
                        IBA1Binary = logical(IBA1Binary);
                        Holes = imfill(IBA1Binary,'holes');
                        Holes = Holes - IBA1Binary;
                        HolesBig = bwareaopen(Holes,5);
                        Holes = Holes - HolesBig;
                        IBA1Binary = IBA1Binary + Holes;
                        IBA1Binary = imclearborder(IBA1Binary);
                    end
                end
            end
            
            
            if(Confocal == 1)
                %Resize images for confocal images:
                Dummy = uint8(zeros(y_Size,x_Size));
                [r c] = size(Nucleus);
                Dummy(1:r,1:c) = Nucleus;
                Nucleus = Dummy;
                NucleusSmall = imresize(Nucleus,0.1);
                Dummy = uint8(zeros(y_Size,x_Size));
                Dummy(1:r,1:c) = IBA1;
                IBA1 = Dummy;
                IBA1Small = imresize(IBA1,0.1);
                DummyBinary = logical(zeros(y_Size,x_Size));
                DummyBinary(1:r,1:c) = NucleusBinary;
                NucleusBinary = DummyBinary;
                DummyBinary(1:r,1:c) = IBA1Binary;
                IBA1Binary = DummyBinary;
                
                %Get Nucleus matrix again:
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
                NucleusM = logical(NucleusMatrix);
            end
            if(CD68M==1)
                [token,remain] = strtok(ImageName, '_');
                a{1} = token;
                while(length(remain) > 0)
                    [token,remain] = strtok(remain, '_');
                end
                a{2} = '_';
                a{3} = token;
                a{4} = 'space';
                newfile = [a{1} a{2} a{3} a{4}];
            else
                token=allfiles(i).name(9:end-4);
                a{2} = token;
                a{3} = 'space';
                newfile = [a{1} a{2} a{3}];
            end
            
            if(CD68M==0)
                IBA1Binary = bwmorph(IBA1Binary,'majority');
            end
            
            newFile1 = regexprep(newfile, 'space','NucleusBig.png','ignorecase');
            imwrite(Nucleus,[foldername1 '/' newFile1]);
            newFile2 = regexprep(newfile, 'space','NucleusSmall.png','ignorecase');
            imwrite(NucleusSmall,[foldername1 '/' newFile2]);
            newFile3 = regexprep(newfile, 'space','NeuriteBig.png','ignorecase');
            imwrite(IBA1,[foldername1 '/' newFile3]);
            newFile4 = regexprep(newfile, 'space','NeuriteSmall.png','ignorecase');
            imwrite(IBA1Small,[foldername1 '/' newFile4]);
            newFile5 = regexprep(newfile, 'space','NucleusBigWatershed.png','ignorecase');
            imwrite(NucleusBinary,[foldername1 '/' newFile5]);
            newFile6 = regexprep(newfile, 'space','Binary.png','ignorecase');
            imwrite(IBA1Binary,[foldername1 '/' newFile6]);
            if(CD68==1)
                newFile7 = regexprep(newfile, 'space','AstroBig.png','ignorecase');
                imwrite(CD68Img,[foldername1 '/' newFile7]);
                newFile8 = regexprep(newfile, 'space','AstroSmall.png','ignorecase');
                imwrite(CD68ImgSmall,[foldername1 '/' newFile8]);
            end
        end
        
        if(Transmission==1)
            ScalingFactor = 0.1;
            wellname = allfiles(i).name(1:end-4);
            wellList{n} = wellname;
            csvHandler.CellPosMatrix(wellname) = NucleusM ;
        else
            if(CD68M==1)
                ScalingFactor = 0.1;
                wellname = [a{1} a{2} a{3}];
                wellList{n} = wellname;
                csvHandler.CellPosMatrix(wellname) = NucleusM;
            else
                ScalingFactor = 1;
                wellname = [a{1} a{2}];
                wellList{n} = wellname;
                csvHandler.CellPosMatrix(wellname) = NucleusM;
            end
        end
        
        
    end
end