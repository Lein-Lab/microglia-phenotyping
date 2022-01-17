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

function [wellList, foldername, ScalingFactor, csvHandler, neuronHandler,  Th] = RainbowCells (csvHandler, neuronHandler)
% This script is designed to monateg images from the 'Rainbow Cell"
% protocol and to count RyR1 and to analyze the intensity of mTor and
% FKBP12
%Running Variables:
n= 0;
o=0;
NumberofImages = 6;
SizeImage = 2160;
ImageMontageNuclei = 0;
ImageMontageNeuron = 0;
ImageMontageTexasRed = 0;
ImageMontageCy5 = 0;
ImageMontageBox = 0;
level = 0;
% Transfection performed y or n:
Transfection = 1;


foldername = uigetdir;
mkdir(foldername,'ConvertedCellomics');
foldername1 = [foldername '/ConvertedCellomics'];
allfiles = dir(foldername);
wellList = cell(0);
MaxIntTexasRyR1 = 0;
MaxIntTexasmTor = 0;
MaxIntCy5FKBP12 = 0;
MaxIntCy5mTor =0;

for(i=1:numel(allfiles))
    ind = strfind([foldername '/' allfiles(i).name],'w1.TIF');
     if(numel(ind) > 0)
     n =n+1;   
     % Reading in images:
     NucleusImage = imread([foldername '/' allfiles(i).name]);
     if(Transfection == 0)
         newFileNeuron = regexprep(allfiles(i).name, 'w1.TIF','w2.TIF','ignorecase');
         NeuriteImage = imread([foldername '/' newFileNeuron]);
         newFileTexasRed = regexprep(allfiles(i).name, 'w1.TIF','w3.TIF','ignorecase');
         TexasRedImage = imread([foldername '/' newFileTexasRed]);
         newFileCy5 = regexprep(allfiles(i).name, 'w1.TIF','w4.TIF','ignorecase');
         Cy5Image = imread([foldername '/' newFileCy5]);
     else
         newFileTexasRed = regexprep(allfiles(i).name, 'w1.TIF','w2.TIF','ignorecase');
         TexasRedImage = imread([foldername '/' newFileTexasRed]);
         newFileCy5 = regexprep(allfiles(i).name, 'w1.TIF','w3.TIF','ignorecase');
         Cy5Image = imread([foldername '/' newFileCy5]);
     end    
     %Before we montage images, we will test whether they are 'empty'
     %Therefore we simply calculate the mean, max and min value:
     MEANNuc = mean(mean(NucleusImage));
     MINNuc = min(min(NucleusImage));
     MAXNuc = max(max(NucleusImage));
     MEDIANNuc = median(median(NucleusImage));
     % Get threshold value for nuclei:
     IntNucSmall = max(max(NucleusImage))*1.2;
     IntNucSmall = double(IntNucSmall);
     NucleusImageT = double(NucleusImage)./IntNucSmall;
     NucleusImageT = uint8(NucleusImageT.*255);
     if(level == 0)
        level = isodata(NucleusImageT);
     end
     
     if(MEANNuc <400 && MEDIANNuc <400)
        NucleusImage = zeros(size(NucleusImage));
        NeuriteImage = zeros(size(NucleusImage));
        TexasRedImage =zeros(size(NucleusImage));
        Cy5Image =zeros(size(NucleusImage));
        Box1 = ones(size(NucleusImage));
     else
        Box1 =  zeros(size(NucleusImage));
     end   
     
     if(Transfection == 0)
     MEANNeu = mean(mean(NeuriteImage));
     MINNeu = min(min(NeuriteImage));
     MAXNeu = max(max(NeuriteImage));
     MEDIANNeu = median(median(NeuriteImage));
     
     if(MEANNeu <400 && MEDIANNeu <400)
        NucleusImage = zeros(size(NeuriteImage));
        NeuriteImage = zeros(size(NeuriteImage));
        TexasRedImage =zeros(size(NeuriteImage));
        Cy5Image =zeros(size(NeuriteImage));
        Box2 = ones(size(NeuriteImage));
     else
        Box2 =  zeros(size(NeuriteImage));
     end   
     end
     
     MEANTex = mean(mean(TexasRedImage));
     MINTex = min(min(TexasRedImage));
     MAXTex = max(max(TexasRedImage));
     MEDIANTex = median(median(TexasRedImage));
     
     if(MEANTex <400 && MEDIANTex <400)
        NucleusImage = zeros(size(NucleusImage));
        if(Transfection == 0)
        NeuriteImage = zeros(size(NucleusImage));
        end
        TexasRedImage =zeros(size(NucleusImage));
        Cy5Image =zeros(size(NucleusImage));
        Box3 = ones(size(NucleusImage));
     else
        Box3 =  zeros(size(NucleusImage));
     end 
     
     
     MEANCy = mean(mean(Cy5Image));
     MINCy = min(min(Cy5Image));
     MAXCy = max(max(Cy5Image));
     MEDIANCy = median(median(Cy5Image));
     
     if(allfiles(i).name(8) == '7' | allfiles(i).name(8) == '8')
         if(MEANCy <200 && MEDIANCy <200)
             NucleusImage = zeros(size(NucleusImage));
             if(Transfection == 0)
                NeuriteImage = zeros(size(NeuriteImage));
             end
        TexasRedImage =zeros(size(NucleusImage));
        Cy5Image =zeros(size(NucleusImage));
        Box4 = ones(size(NucleusImage));
         else
         Box4 =  zeros(size(NucleusImage));
         end
     else    
     if(MEANCy <400 && MEDIANCy <400)
        NucleusImage = zeros(size(NucleusImage));
        NeuriteImage = zeros(size(NucleusImage));
        TexasRedImage =zeros(size(NucleusImage));
        Cy5Image =zeros(size(NucleusImage));
        Box4 = ones(size(NucleusImage));
     else
        Box4 =  zeros(size(NucleusImage));
     end
     end
     %Now we can add all boxes and make them logical again:
     if(Transfection == 0)
     Box = Box1 + Box2 + Box3 + Box4;
     else
     Box = Box1 + Box3 + Box4;
     end
     Box = logical(Box);
     
     
     %Now we have to write the images depending of their relative position!
     %First row will be 1-6 or max NumberofImages!
     %Therefore we need to find the sight number! We have to be careful
     %again since we have single and double digit numbers
     if(allfiles(i).name(end-8) == 's')
        SightNumber = allfiles(i).name(end-7);
     else
        SightNumber = allfiles(i).name(end-8:end-7);
     end    
     SightNumber = str2num(SightNumber);
     % Now we have to identify the row number the Image will be in!
     STOP = 0;
     for(v=1:NumberofImages)
     if(STOP == 0)    
     if(SightNumber<(v*NumberofImages)+1)
        Row = v; 
        STOP =1;
     end 
     end
     end
     
     % Now we need to determine the coloumn number in the row!
     %Uneven rows are in the correct order and even are reversed!
     %Uneven row will be 1!
        if(Row == 1)
            Colomn = SightNumber;
        else
            Colomn = SightNumber - ((Row -1) * NumberofImages)
        end
     
     % Now we will montage the images:
     %Create an empty image:
     if(Colomn == 4 && Row == 2)
     Imagesize = NumberofImages*SizeImage;
     ImageMontageNuclei = zeros(Imagesize,Imagesize);
     ImageMontageNuclei = uint16(ImageMontageNuclei);
     if(Transfection == 0)
     ImageMontageNeuron = zeros(Imagesize,Imagesize);
     ImageMontageNeuron = uint16(ImageMontageNeuron);
     end     
     ImageMontageTexasRed = zeros(Imagesize,Imagesize);
     ImageMontageTexasRed = uint16(ImageMontageTexasRed);
     ImageMontageCy5 = zeros(Imagesize,Imagesize);
     ImageMontageCy5 = uint16(ImageMontageCy5);
     ImageMontageBox = zeros(Imagesize,Imagesize);
     end
     %Now we can calculate the x,y coordinates:
     x_Start = 1 + SizeImage*(Colomn -1);
     x_End = SizeImage*Colomn;
     y_Start= 1+ SizeImage*(Row -1);
     y_End = SizeImage*Row;
     
     %Now we write the image into the Montage:
     %NucleusImage = NucleusImage*255;
     ImageMontageNuclei(y_Start:y_End,x_Start:x_End) = NucleusImage;
     if(Transfection == 0)
     ImageMontageNeuron(y_Start:y_End,x_Start:x_End) = NeuriteImage;
     end
     ImageMontageTexasRed(y_Start:y_End,x_Start:x_End) = TexasRedImage;
     ImageMontageCy5(y_Start:y_End,x_Start:x_End) = Cy5Image;
     ImageMontageBox(y_Start:y_End,x_Start:x_End) = Box;
     
     %Before we write the images, we need to convert them to 8-bit!
     %Since a normal conversion would lead to an imense stretch, we will
     %calculate the maximum intensity, add 20% of each single image (Thereby we will ignore images which are exceptional bright). We will ignore picture with max intensity of 65535:
     %We generate an empty vector containing the intensity values:
     IntNuc = max(max(NucleusImage));
     IntNuc = double(IntNuc);
     if(Transfection == 0)
     IntNeurite = max(max(NeuriteImage));
     IntNeurite = double(IntNeurite);
     end
     IntTexas = max(max(TexasRedImage));
     IntTexas = double(IntTexas);
     IntCy5 = max(max(Cy5Image));
     IntCy5 = double(IntCy5);
     if(IntNuc<65535)
     IntNucAll(n,1) =  IntNuc;
     end
     if(Transfection == 0)
     if(IntNeurite<65535)
     IntNeuriteAll(n,1) =  IntNeurite;
     end
     end
     if(IntTexas<65535)
     IntTexasAll(n,1) =  IntTexas;
     end
     if(IntCy5<65535)
     IntCy5All(n,1) =  IntCy5;
     end
     
     %Now we write the images. However we need to write after the last image is placed and then replace images with empties again!:
     if(SightNumber == 9)
     o = o +1;
     if(o < 10)
         a{1} = 'A0';
         a{2} = int2str(o);
         a{3} = 'space';
     else
         a{1} = 'A';
         a{2} = int2str(o);
         a{3} = 'space';
     end
     newfile = [a{1} a{2} a{3}];
     Name = [a{1} a{2}];
     %Here we calculate the final threshold. However, since we used
     %different markers, we need to calculate different thresholds!
     %1) Dapi and FICI is alwas the same:
     MaxIntNuc = max(IntNucAll)*1.2;
     if(MaxIntNuc>12000)
         MaxIntNuc = 12000;
     end 
     if(Transfection == 0)
     MaxIntNeurite = max(IntNeuriteAll)*1.2;
     if(MaxIntNeurite>65535)
         MaxIntNeurite = 65535;
     end  
     end
     ImageMontageNuclei = double(ImageMontageNuclei)./MaxIntNuc;
     ImageMontageNuclei = uint8(ImageMontageNuclei.*255);
     if(Transfection == 0)
     ImageMontageNeuron = double(ImageMontageNeuron)./MaxIntNeurite;
     ImageMontageNeuron = uint8(ImageMontageNeuron.*255);
     end
     %2) Texas Red
     %For Well A01, A02, A04 and A05, the Texas Red threshold is the same for RyR1!
     %Since we need to compare the images, we need to use the same stretch
     %on all of them. Therefore we will calculate it in the first well for a respective marker:
     %2 a) RyR1
     if(MaxIntTexasRyR1 == 0)
     if(Name(3) == '1')
     MaxIntTexasRyR1 = max(IntTexasAll)*1.2;
     if(MaxIntTexasRyR1>65535)
         MaxIntTexasRyR1 = 65535;
     end    
     end
     end

     if(Name(3) == '1' | Name(3) == '2'| Name(3) == '4' | Name(3) == '5')
     ImageMontageTexasRed = double(ImageMontageTexasRed)./MaxIntTexasRyR1;
     ImageMontageTexasRed = uint8(ImageMontageTexasRed.*255);
     %Now we can also clean up the image, so that only RYR1 receptors
     %remain!
%      SE = strel('ball',10,10);
%      backgroundball = imopen(ImageMontageTexasRed,SE);
%      ImageMontageTexasRed = ImageMontageTexasRed -  backgroundball;
     % Now we have to readjust the histogramm!
     MaxIntTexasRyR1Ball = max(max(ImageMontageTexasRed));
     MaxIntTexasRyR1Ball = double(MaxIntTexasRyR1Ball);
     ImageMontageTexasRed = double(ImageMontageTexasRed)./MaxIntTexasRyR1Ball;
     ImageMontageTexasRed = uint8(ImageMontageTexasRed.*255);
     end
     
     %2 b) mTor
     if(MaxIntTexasmTor ==0)
     if(Name(3) == '3')
         MaxIntTexasmTor = max(IntTexasAll)*1.2;
         if(MaxIntTexasmTor>65535)
             MaxIntTexasmTor = 65535;
         end    
     end
     end
     
     if(Name(3) == '3' | Name(3) == '6')
        ImageMontageTexasRed = double(ImageMontageTexasRed)./MaxIntTexasmTor;
        ImageMontageTexasRed = uint8(ImageMontageTexasRed.*255);
     end
     
     %3) For Cy5:
     %3 a) FKBP12
     if(MaxIntCy5FKBP12 == 0)
     if(Name(3) == '1')
     MaxIntCy5FKBP12 = max(IntCy5All)*1.2;
     if(MaxIntCy5FKBP12>65535)
         MaxIntCy5FKBP12=65535;
     end    
     end
     end
     
     if(Name(3) == '1' | Name(3) == '3' | Name(3) == '4' | Name(3) == '6')
     ImageMontageCy5 = double(ImageMontageCy5)./MaxIntCy5FKBP12;
     ImageMontageCy5 = uint8(ImageMontageCy5.*255);
     end
     
     %3) For Cy5:
     %3 a) mTor
     if(MaxIntCy5mTor == 0)
     if(Name(3) == '2')
     MaxIntCy5mTor = max(IntCy5All)*1.2;
     if(MaxIntCy5mTor>65535)
         MaxIntCy5mTor = 65535;
     end    
     end
     end
     
     if(Name(3) == '2' | Name(3) == '5' )
     ImageMontageCy5 = double(ImageMontageCy5)./MaxIntCy5mTor;
     ImageMontageCy5 = uint8(ImageMontageCy5.*255);
     end
     
     
     %4) Secondary AK Controls!
     if(Name(3) == '7')
         ImageMontageTexasRed = double(ImageMontageTexasRed)./MaxIntTexasRyR1;
         ImageMontageTexasRed = uint8(ImageMontageTexasRed.*255);
         %Now we can also clean up the image, so that only RYR1 receptors
         %remain!
         SE = strel('ball',10,10);
         backgroundball = imopen(ImageMontageTexasRed,SE);
         ImageMontageTexasRed = ImageMontageTexasRed -  backgroundball;
         % Now we have to readjust the histogramm!
         ImageMontageTexasRed = double(ImageMontageTexasRed)./MaxIntTexasRyR1Ball;
         ImageMontageTexasRed = uint8(ImageMontageTexasRed.*255);
         ImageMontageCy5 = double(ImageMontageCy5)./MaxIntCy5mTor;
         ImageMontageCy5 = uint8(ImageMontageCy5.*255);
     end
     
     if(Name(3) == '8')
        ImageMontageTexasRed = double(ImageMontageTexasRed)./MaxIntTexasmTor;
        ImageMontageTexasRed = uint8(ImageMontageTexasRed.*255);
        ImageMontageCy5 = double(ImageMontageCy5)./MaxIntCy5FKBP12;
        ImageMontageCy5 = uint8(ImageMontageCy5.*255);
     end   
     
     
     
     % Now we can also perform nucleus identification:
     binaryImgNuc = im2bw(ImageMontageNuclei,level);
     binaryImgNuc = bwareaopen(binaryImgNuc,5000);
     level = 0;
     binaryFill = imfill(binaryImgNuc,'holes');
     Holes = binaryFill - binaryImgNuc;
     Holes = bwareaopen(Holes,50);
     binaryImgNuc = binaryImgNuc + Holes;
     binaryImgNuc = imclose(binaryImgNuc,strel('disk',5));
     binaryImgNuc = logical(binaryImgNuc);
     s = regionprops(binaryImgNuc, 'Centroid');
        sdimensions = size(s);
        NucleusM = sparse(zeros(size(ImageMontageNuclei)));
        if(sdimensions > 0)
            for(j=1:sdimensions)
            x_centroid = s(j).Centroid(1);
            y_centroid = s(j).Centroid(2);
            NucleusM(uint16(y_centroid),uint16(x_centroid)) = 1;
            end
        end
     
     
     
     
     newFile1 = regexprep(newfile, 'space','NucleusSmall.png','ignorecase');
     newFile2 = regexprep(newfile, 'space','NucleusBig.png','ignorecase');
     if(Transfection == 0)
     newFile3 = regexprep(newfile, 'space','NeuriteSmall.png','ignorecase');
     newFile4 = regexprep(newfile, 'space','NeuriteBig.png','ignorecase');
     end
     newFile5 = regexprep(newfile, 'space','OligoSmall.png','ignorecase');
     newFile6 = regexprep(newfile, 'space','OligoBig.png','ignorecase');
     newFile7 = regexprep(newfile, 'space','AstroSmall.png','ignorecase');
     newFile8 = regexprep(newfile, 'space','AstroBig.png','ignorecase');
     newFile9 = regexprep(newfile, 'space','NucleusBigWatershed.png','ignorecase');
     newFile10 = regexprep(newfile, 'space','Box.png','ignorecase');
     ScalingFactor = 0.1;
     ImageMontageNucleiSmall = imresize(ImageMontageNuclei,ScalingFactor);
     if(Transfection == 0)
     ImageMontageNeuronSmall = imresize(ImageMontageNeuron,ScalingFactor);
     end
     ImageMontageTexasRedSmall = imresize(ImageMontageTexasRed,ScalingFactor);
     ImageMontageCy5Small = imresize(ImageMontageCy5,ScalingFactor);
     imwrite(ImageMontageNucleiSmall,[foldername1 '/' newFile1]);
     imwrite(ImageMontageNuclei,[foldername1 '/' newFile2]);
     if(Transfection == 0)
     imwrite(ImageMontageNeuronSmall,[foldername1 '/' newFile3]);
     imwrite(ImageMontageNeuron,[foldername1 '/' newFile4]);
     end
     imwrite(ImageMontageTexasRedSmall,[foldername1 '/' newFile5]);
     imwrite(ImageMontageTexasRed,[foldername1 '/' newFile6]);
     imwrite(ImageMontageCy5Small,[foldername1 '/' newFile7]);
     imwrite(ImageMontageCy5,[foldername1 '/' newFile8]);
     imwrite(binaryImgNuc,[foldername1 '/' newFile9]);
     imwrite(ImageMontageBox,[foldername1 '/' newFile10]);
     wellname = (newFile1(1:3));
     wellList{o} = wellname;
     csvHandler.CellPosMatrix(wellname) = NucleusM;
     neuronHandler.NeuronPositionsEdgeFill(wellname) = sparse(zeros(size(ImageMontageNuclei)));
     ImageMontageNuclei = 0;
     if(Transfection == 0)
     ImageMontageNeuron = 0;
     end
     ImageMontageTexasRed =0;
     ImageMontageCy5 = 0;
     ImageMontageBox = 0;
     end
     
     end
end

blub = 2;