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

% This script is designed to analyze the number of RYR1 receptors within
% the preprocessed images of the rainbow cells, as well as the total
% area/intensity of the mTor and FKBP12 markers!
n= 0;
level = 0;
OneCell = 10000;
foldername = uigetdir;
allfiles = dir(foldername);
Transfection = 1;

for(i=1:numel(allfiles))
    ind = strfind([foldername '/' allfiles(i).name],'NucleusBig.');
     if(numel(ind) > 0)
     n =n+1;   
     if(foldername(end-23) == '\')
         FileName = foldername(end-22:end-19);
     else
         FileName = foldername(end-23:end-19);
     end
     % Reading in images:
     NucleusImage = imread([foldername '/' allfiles(i).name]);
     if(Transfection == 0)
     newFileNeuron = regexprep(allfiles(i).name, 'NucleusBig','NeuriteBig','ignorecase');
     NeuriteImage = imread([foldername '/' newFileNeuron]);
     end
     newFileTexasRed = regexprep(allfiles(i).name, 'NucleusBig','OligoBig','ignorecase');
     TexasRedImage = imread([foldername '/' newFileTexasRed]);
     newFileCy5 = regexprep(allfiles(i).name, 'NucleusBig','AstroBig','ignorecase');
     Cy5Image = imread([foldername '/' newFileCy5]);
     newFileWater = regexprep(allfiles(i).name, 'NucleusBig','NucleusBigWatershed','ignorecase');
     WatershedImage = imread([foldername '/' newFileWater]);
     
     WatershedImage = logical(WatershedImage);
     WaterIn = uint8(WatershedImage);
     WaterIn(WaterIn==0)=255;
     WaterIn(WaterIn==1)=0;
     
     % Now we precossess TexasRedImage and CyImages. We have to be again
     % carefull to check for different markers in different channels:
     %1) RyR1
     Name = allfiles(i).name(1:3);
     a{1} = Name;
     a{2} = 'space';
     newfile = [a{1} a{2}];
     if(Name == 'A01' | Name == 'A02' | Name == 'A04' | Name == 'A05' | Name == 'A07')
         %Do rolling ball:
         SE = strel('ball',10,10);
         BackgroundBall = imopen(TexasRedImage,SE);
         RyR1Spots =  TexasRedImage - BackgroundBall;
         if(level == 0);
         [lehisto x] = imhist(RyR1Spots);
         level = triangle_th(lehisto,256);
         end
         RyR1Binary=im2bw(RyR1Spots,level*3);
         RyR1Binary = bwareaopen(RyR1Binary,2);
         NumberofRyR = regionprops(RyR1Binary,'centroid');
         NumberofRyR = length(NumberofRyR);
         newfile1 = regexprep(newfile, 'space','Binary.png','ignorecase');
         imwrite(RyR1Binary,[foldername '/' newfile1]);
         TexasRedImageInt = TexasRedImage - 23;
         IntensityRyR1 = sum(sum(TexasRedImageInt));
         RyR1Nuc = TexasRedImageInt- WaterIn;
         RyR1Nuc = uint8(RyR1Nuc);
         IntensityRyR1Nuc = sum(sum(RyR1Nuc));
         %Now we can also get intranuclear markes
         WatershedImage = logical(WatershedImage);
         NumberofRyRinNuc = RyR1Binary + WatershedImage;
         NumberofRyRinNuc = NumberofRyRinNuc -1;
         NumberofRyRinNuc = uint8(NumberofRyRinNuc);
         NumberofRyRinNuc = logical(NumberofRyRinNuc);
         NumberofRyRinNuc = regionprops(NumberofRyRinNuc,'centroid');
         NumberofRyRinNuc = length(NumberofRyRinNuc);
     else
         NumberofRyR = 0;
         NumberofRyRinNuc = 0;
         IntensityRyR1 = 0;
         IntensityRyR1Nuc = 0;
         RyR1Binary = zeros(size(NucleusImage));
         newfile1 = regexprep(newfile, 'space','Binary.png','ignorecase');
         imwrite(RyR1Binary,[foldername '/' newfile1]);
     end
     
    
     
     
     %2) FKBP12 in Cy5
     if(Name(1:3) == 'A01' | Name(3) == 'A03' | Name(3) == 'A04' | Name(3) == 'A06' | Name(3) == 'A08')
         Cy5ImageInt = Cy5Image -28;
         IntensityFKBP12Cy5 = sum(sum(Cy5ImageInt)); 
         %Intranuclear:
         FKBP12Cy5Nuc = Cy5ImageInt- WaterIn;
         FKBP12Cy5Nuc = uint8(FKBP12Cy5Nuc);
         IntensityFKBP12Cy5Nuc = sum(sum(FKBP12Cy5Nuc));
     else
         IntensityFKBP12Cy5 =0;
         IntensityFKBP12Cy5Nuc = 0;
     end
     
     %3) mTOR in  Texas:
     if(Name(1:3) == 'A03' | Name(3) == 'A06' | Name(3) == 'A08')
        TexasRedImageInt =  TexasRedImage - 19;
        IntensitymTorTexas = sum(sum(TexasRedImageInt));  
        mTorTexasNuc = TexasRedImageInt - WaterIn;
        mTorTexasNuc = uint8(mTorTexasNuc);
        IntensitymTorTexasNuc = sum(sum(mTorTexasNuc));
     else
        IntensitymTorTexas = 0;
        IntensitymTorTexasNuc =0;
     end
     
     %4) mTOR in  Cy5:
     if(Name(1:3) == 'A02' | Name(3) == 'A05' | Name(3) == 'A07')
         Cy5ImageInt = Cy5Image - 15;
         IntensitymTorCy5 = sum(sum(Cy5ImageInt)); 
         Cy5ImageNuc = Cy5ImageInt - WaterIn;
         Cy5ImageNuc = uint8(Cy5ImageNuc);
         IntensitymTorCy5Nuc = sum(sum(Cy5ImageNuc));
         
     else
         IntensitymTorCy5 =0;
         IntensitymTorCy5Nuc = 0;
     end
     
     
    
     %5) Nuclei:
%      WatershedImage = uint8(WatershedImage);
%      WatershedImage(WatershedImage==0)=255;
%      WatershedImage(WatershedImage==1)=0;
     WatershedImage = logical(WatershedImage);
     AreaNuc = sum(sum(WatershedImage));
     %Now we will also calculate the number of Nuclei:
     AreaNumberVector = regionprops(WatershedImage,'Area');
     AreaNumberVector = struct2table(AreaNumberVector);
     AreaNumberVector = table2array(AreaNumberVector);
     % Now we also 
     NumberofCells = AreaNuc/OneCell;
     
     
     %Now we create a matrix in which we will save all endpoints to later
     %save it in excel:
     NameTabel{1,2} = 'Cell Number';
     NameTabel{1,3} = 'Number of RyR1-Receptors';
     NameTabel{1,4} = 'Intensity RyR1';
     NameTabel{1,5} = 'Intensity mTor';
     NameTabel{1,6} = 'Intensity FKBP12';
     NameTabel{1,7} = 'Number of RyR1-Receptors per Nucleus';
     NameTabel{1,8} = 'Intensity RyR1 per Nucleus';
     NameTabel{1,9} = 'Intensity mTor per Nucleus';
     NameTabel{1,10} = 'Intensity FKBP12 per Nucleus';
     NameTabel{1,11} = 'Intensity RyR1 within Nuclei';
     NameTabel{1,12} = 'Intensity FKBP12 within Nuclei';
     NameTabel{1,13} = 'Intensity mTor within Nuclei';
     NameTable2{3,1} = 'Well Number';
     NameTable2{n+3,1} = Name;
     
     % Now we will in Data:
     %1) Number of Nuclei:
     Table(n+3,2) = NumberofCells;
     %2) Number of RyR1 receptors:
     Table(n+3,3) = NumberofRyR;
     Table(n+3,4) = IntensityRyR1;
     %3) Intensity mTor:
     if(IntensitymTorTexas>0)
         IntensitymTor = IntensitymTorTexas;
     else
         if(IntensitymTorCy5>0)
            IntensitymTor = IntensitymTorCy5;
         else
            IntensitymTor = 0;
         end
     end     
     Table(n+3,5) = IntensitymTor;
     %4) Intensity FKBP12:
     Table(n+3,6) = IntensityFKBP12Cy5;
     %5) Number of RyR1 receptors per Nucleus:
     NumberofRyR1perNuc = NumberofRyR/NumberofCells;
     Table(n+3,7) = NumberofRyR1perNuc;
     IntensityRyR1perNuc = IntensityRyR1/NumberofCells;
     Table(n+3,8) = IntensityRyR1perNuc;
     
     %6) ntensity mTor per Nucleus:
     IntensitymTorNuc = IntensitymTor/NumberofCells;
     Table(n+3,9) = IntensitymTorNuc;
     %6) ntensity FKBP12 per Nucleus:
     IntensityFKBP12Nuc = IntensityFKBP12Cy5/NumberofCells;
     Table(n+3,10) = IntensityFKBP12Nuc;
     Table(n+3,11) = IntensityRyR1Nuc;
     Table(n+3,12) = IntensityFKBP12Cy5Nuc;
     if(IntensitymTorTexasNuc>0)
         IntensitymTorIN = IntensitymTorTexasNuc;
     else
         if(IntensitymTorCy5Nuc>0)
            IntensitymTorIN = IntensitymTorCy5Nuc;
         else
            IntensitymTorIN = 0;
         end
     end
     Table(n+3,13) = IntensitymTorIN;
     
     end
end     
b{1} = FileName;
b{2} = '.xlsx';
File = [b{1} b{2}];
xlswrite(File,Table);
xlswrite(File,NameTable2);
xlswrite(File,NameTabel);