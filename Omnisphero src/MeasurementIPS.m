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

function [summary] = MeasurementIPS (foldername, csvHandler, neuronHandler)
foldername2 = foldername;
foldername = [foldername '/ConvertedCellomics'];
mkdir(foldername);
allfiles = dir(foldername);
allfiles2 = dir(foldername2);

%Get all file names from allfiles2:
NameList = struct2table(allfiles2);
NameList = NameList(:,1);
%Keep only tiff, jpg, jpeg:
NameList = table2cell(NameList);
Index = strfind(NameList,'tif');
Index2 = strfind(NameList,'jpg');
Index3 = strfind(NameList,'jpeg');
NewList = {};
Count = 0;
for(y=1:length(Index))
    Index1 = Index{y};
    if(Index1>0)
        Count = Count + 1;
        NewList{Count,1} = allfiles2(y).name;
    end
end



wellList = cell(0);
n=0;
AxonMeasurement = 0;
for(i=1:numel(allfiles))
    %Check if file ends with .tif or .tiff
    ind = strfind([foldername '/ConvertedCellomics' allfiles(i).name],'NeuriteBig');
    if(numel(ind) > 0)%
         n = n + 1;
         NeuronImage = imread([foldername '/' allfiles(i).name]);
         newFileSynapsen = regexprep(allfiles(i).name,'NeuriteBig','AstroBig','ignorecase');
         SynapseImage = imread([foldername '/' newFileSynapsen]);
         
                                                                            %         if(n<13)
                                                                            %             a{1} = 'A0';
                                                                            %             a{2} = int2str(n);
                                                                            %             a{3} = 'space';
                                                                            %         end
                                                                            %         if(n>12&&n<25)
                                                                            %         a{1} = 'B0';
                                                                            %         a{2} = int2str(p);
                                                                            %         a{3} = 'space';
                                                                            %         end
                                                                            %         if(n>10&&n<16)
                                                                            %         a{1} = 'C0';
                                                                            %         p = n - 10;
                                                                            %         a{2} = int2str(p);
                                                                            %         a{3} = 'space';
                                                                            %         end
                                                                            %         if(n>15&&n<21)
                                                                            %         a{1} = 'D0';
                                                                            %         p = n - 15;
                                                                            %         a{2} = int2str(p);
                                                                            %         a{3} = 'space';
                                                                            %         end
                                                                            %         if(n>20&&n<26)
                                                                            %         a{1} = 'E0';
                                                                            %         p = n - 20;
                                                                            %         a{2} = int2str(p);
                                                                            %         a{3} = 'space';
                                                                            %         end
                                                                            %         if(n>25&&n<31)
                                                                            %         a{1} = 'F0';
                                                                            %         p = n - 25;
                                                                            %         a{2} = int2str(p);
                                                                            %         a{3} = 'space';
                                                                            %         end
                                                                            %         if(n>30&&n<36)
                                                                            %         a{1} = 'G0';
                                                                            %         p = n - 30;
                                                                            %         a{2} = int2str(p);
                                                                            %         a{3} = 'space';
                                                                            %         end
                                                                            %         if(n>35&&n<41)
                                                                            %         p = n - 35;
                                                                            %         a{2} = int2str(p);
                                                                            %         a{3} = 'space';
                                                                            %         a{1} = 'H0';
                                                                            %         end
                                                                            
        if(n<13)
            if(n<10)
                a{1} = 'A0';
            else
                a{1} = 'A';
            end
        a{2} = int2str(n);
        a{3} = 'space';
        end
        if(12>5&&n<25)
        p = n - 12;
        if(p<10)
            a{1} = 'B0';
        else
            a{1} = 'B';
        end
        a{2} = int2str(p);
        a{3} = 'space';
        end
        if(n>24&&n<37)
        
        p = n - 24;
        if(p<10)
            a{1} = 'C0';
        else
            a{1} = 'C';
        end
        a{2} = int2str(p);
        a{3} = 'space';
        end
        if(n>36&&n<49)
        
        p = n - 36;
        if(p<10)
            a{1} = 'D0';
        else
            a{1} = 'D';
        end
        a{2} = int2str(p);
        a{3} = 'space';
        end
        if(n>48&&n<61)
        
        p = n - 48;
        if(p<10)
            a{1} = 'E0';
        else
            a{1} = 'E';
        end
        a{2} = int2str(p);
        a{3} = 'space';
        end
        if(n>60&&n<73)
        
        p = n - 60;
        if(p<10)
            a{1} = 'F0';
        else
            a{1} = 'F';
        end
        a{2} = int2str(p);
        a{3} = 'space';
        end
        if(n>72&&n<85)
        
        p = n - 72;
        if(p<10)
            a{1} = 'G0';
        else
            a{1} = 'G';
        end
        a{2} = int2str(p);
        a{3} = 'space';
        end
        if(n>84&&n<96)
        p = n - 84;
        a{2} = int2str(p);
        a{3} = 'space';
        if(p<10)
            a{1} = 'H0';
        else
            a{1} = 'H';
        end
        end

        NeuriteMass = sum(sum(NeuronImage));
        NumberofSynapses = bwlabel(SynapseImage);
        NumberofSynapses = max(max(NumberofSynapses));
        
        
        
        %Create Matrix
        summary{1,1} =  'Image name';
        summary{1,2} =  '[Neurite Mass]';
        summary{1,3} =  '[Number of Synapses]';
        summary{1,4} =  '[Synapse per Area]';
        
        m=n+1;
        summary{m,1} =  NewList{n};
        summary{m,2} =  NeuriteMass;
        summary{m,3} =  NumberofSynapses;
        summary{m,4} =  NumberofSynapses/NeuriteMass;
        
        summary = summary
    end
end
         