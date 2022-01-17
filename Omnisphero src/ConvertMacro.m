%Copyright (C) 2013-2021  Martin Schmuck, Thomas Temme

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

%Convert Tiff from imageJ to png and create small images:

n= 0;
foldername = uigetdir;
mkdir(foldername,'ConvertedCellomics');
foldername1 = [foldername '/ConvertedCellomics'];
allfiles = dir(foldername);
Transfection = 1;
for(i=1:numel(allfiles))
    ind = strfind([foldername '/' allfiles(i).name],'NucleusBig.tif')
     if(numel(ind) > 0)
     NucleusImage = imread([foldername '/' allfiles(i).name]);
     if(Transfection ==0)
     newFileNeuron = regexprep(allfiles(i).name, 'NucleusBig','NeuriteBig','ignorecase');
     NeuriteImage = imread([foldername '/' newFileNeuron]);
     end
     newFileTexasRed = regexprep(allfiles(i).name, 'NucleusBig','OligoBig','ignorecase');
     TexasRedImage = imread([foldername '/' newFileTexasRed]);
     newFileCy5 = regexprep(allfiles(i).name, 'NucleusBig','AstroBig','ignorecase');
     Cy5Image = imread([foldername '/' newFileCy5]);
     newFileWater = regexprep(allfiles(i).name, 'NucleusBig','NucleusBigWatershed','ignorecase');
     WatershedImage = imread([foldername '/' newFileWater]);
     
     Name = allfiles(i).name(1:3);
     b{1} = Name;
     b{2} = 'space';
     newfile = [b{1} b{2}];
     
     newFile1 = regexprep(newfile, 'space','NucleusSmall.png','ignorecase');
     newFile2 = regexprep(newfile, 'space','NucleusBig.png','ignorecase');
     if(Transfection ==0)
     newFile3 = regexprep(newfile, 'space','NeuriteSmall.png','ignorecase');
     newFile4 = regexprep(newfile, 'space','NeuriteBig.png','ignorecase');
     end
     newFile5 = regexprep(newfile, 'space','OligoSmall.png','ignorecase');
     newFile6 = regexprep(newfile, 'space','OligoBig.png','ignorecase');
     newFile7 = regexprep(newfile, 'space','AstroSmall.png','ignorecase');
     newFile8 = regexprep(newfile, 'space','AstroBig.png','ignorecase');
     newFile9 = regexprep(newfile, 'space','NucleusBigWatershed.png','ignorecase');
     ScalingFactor = 0.1;
     ImageMontageNucleiSmall = imresize(NucleusImage,ScalingFactor);
     if(Transfection ==0)
     ImageMontageNeuronSmall = imresize(NeuriteImage,ScalingFactor);
     end
     ImageMontageTexasRedSmall = imresize(TexasRedImage,ScalingFactor);
     ImageMontageCy5Small = imresize(Cy5Image,ScalingFactor);
     imwrite(ImageMontageNucleiSmall,[foldername1 '/' newFile1]);
     imwrite(NucleusImage,[foldername1 '/' newFile2]);
     if(Transfection ==0)
     imwrite(ImageMontageNeuronSmall,[foldername1 '/' newFile3]);
     imwrite(NeuriteImage,[foldername1 '/' newFile4]);
     end
     imwrite(ImageMontageTexasRedSmall,[foldername1 '/' newFile5]);
     imwrite(TexasRedImage,[foldername1 '/' newFile6]);
     imwrite(ImageMontageCy5Small,[foldername1 '/' newFile7]);
     imwrite(Cy5Image,[foldername1 '/' newFile8]);
     imwrite(WatershedImage,[foldername1 '/' newFile9]);
    
     
     end
end     