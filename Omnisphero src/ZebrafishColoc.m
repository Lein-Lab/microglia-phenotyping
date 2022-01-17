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

% This script is designed to read in the CSV files from the coordinates of
% Imaris and to calculate the shortest distance bewteen red labeled cells
% and gree. If the distance falls below a threshold, the two markers are
% defined as co-localized:

foldername = uigetdir;
foldername1 = [foldername '/Red'];
foldername2 = [foldername '/Green'];
allfiles = dir(foldername1);
allfiles1 = dir(foldername2);
n=0;
minDist = 5;
c{1} = foldername(end-2:end);
c{2} = '.xlsx';
% c{3} = '_Red';
% c{4} = '_Green';
File = [c{1} c{2}];
% File1 = [c{1} c{3} c{2}];
% File2 = [c{1} c{4} c{2}];
for(i=1:numel(allfiles))
     ind = strfind([foldername1 '/' allfiles(i).name],'.xls');
     if(numel(ind) > 0)
     n= n+1;
     sheet = 17;
     xlRange = 'A3:C1000';
     a{1} = foldername1;
     a{2} = '/';
     a{3} = allfiles(i).name;
     b{1} = foldername2;
     b{2} = '/';
     b{3} = allfiles1(i).name;
     Red = [a{1} a{2} a{3}];
     Green = [b{1} b{2} b{3}];
     RedList = xlsread(Red,sheet,xlRange);
     GreenList = xlsread(Green,sheet,xlRange);
         
     %Now we will load chck for each coordinate, whether there is a
     %coordinate with a user defined distance (This is dependant of the
     %size of the spheroid, which was set to 10um)
     LengthRed = length(RedList);
     LengthGreen = length(GreenList);
     RedCount = 0;
     y = 0;
     row = 0;
     DistanceV = 0;
     V_CoordinatesRG_Red = 0;
     for(o=1:LengthRed)
         xRed = RedList(o,1);
         yRed = RedList(o,2);
         zRed = RedList(o,3);
         for(t=1:LengthGreen);
             xGreen = GreenList(t,1);
             yGreen = GreenList(t,2);
             zGreen = GreenList(t,3);
             %Now we calculate the euclidean distance to the red spot:
             distGR = ((xGreen-xRed)^2+(yGreen-yRed)^2+(zGreen-zRed)^2)^0.5;
             DistanceV(t,1) = distGR;
         end
         minDistV = min(DistanceV);
         inds = find(DistanceV==minDistV);
         row = ind2sub(size(DistanceV),inds);
         xGreen_min = GreenList(row,1);
         yGreen_min = GreenList(row,2);
         zGreen_min = GreenList(row,3);
         if(minDistV < minDist)
             y =y+1;
             r = y +1;  
             sheet = r;
                 V_CoordinatesRG_Red(y,1) = xGreen_min;
                 V_CoordinatesRG_Red(y,2) = yGreen_min;
                 V_CoordinatesRG_Red(y,3) = zGreen_min;
%                V_CoordinatesRG_Red(y,4) = minDistV;
                 %xlswrite(File1,V_CoordinatesRG_Red,sheet);  
                 RedCount  = RedCount + 1;
             
         end    
     end    
     
     % Do the same for green:
     GreenCount = 0;
     DistanceV1 = 0;
     row = 0;
     DistanceV1 = 0;
     V_CoordinatesRG_Red1 = 0;
     y=0;
     for(f=1:LengthGreen)
         xGreen = GreenList(f,1);
         yGreen = GreenList(f,2);
         zGreen = GreenList(f,3);
         for(p=1:LengthRed);
             xRed = RedList(p,1);
             yRed = RedList(p,2);
             zRed = RedList(p,3);
             %Now we calculate the euclidean distance to the red spot:
             distGR_1 = ((xGreen-xRed)^2+(yGreen-yRed)^2+(zGreen-zRed)^2)^0.5;
             DistanceV1(p,1) = distGR_1;
         end
         minDistV_1 = min(DistanceV1);
         inds = find(DistanceV1==minDistV_1);
         row = ind2sub(size(DistanceV1),inds);
         xRed_min = RedList(row,1);
         yRed_min = RedList(row,2);
         zRed_min = RedList(row,3);
         if(minDistV_1 < minDist)
             y =y+1;
             r = y +1;  
             sheet = r;
             V_CoordinatesRG_Red1(y,1) = xRed_min;
             V_CoordinatesRG_Red1(y,2) = yRed_min;
             V_CoordinatesRG_Red1(y,3) = zRed_min;
%            V_CoordinatesRG_Red1(y,4) = minDistV;
             %xlswrite(File2,V_CoordinatesRG_Red1,sheet);  
             GreenCount = GreenCount + 1;
         end    
     end
    % Now we write the values in a final table:
    Name = allfiles(i).name(1:end-4);
    m = n+1;
    TableName{m,1} = Name;
    TableName2{1,2} = 'Number of red cells';
    TableName2{1,3} = 'Number of green cells';
    TableName2{1,4} = 'Number of co-localized cells (red with green)';
    TableName2{1,5} = 'Number of co-localized cells (green with red)';
    Table(m,2) = LengthRed;
    Table(m,3) = LengthGreen;
    Table(m,4) =RedCount;
    Table(m,5) =GreenCount;
    RedList = 0;
    GreenList = 0;
    
     
     
     
     end
end     


xlswrite(File,Table);
xlswrite(File,TableName);
xlswrite(File,TableName2);
clear all