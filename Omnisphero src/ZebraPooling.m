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
% This matlab script is designed for pooling data from fish beahivior
% assays in the Lein Lab:

%Running variable:
n = 0;
Multiplier = 0;

%Directory;
foldername = uigetdir;
allfiles = dir(foldername);

% Loop over all the files in the folder uigetdir;

for(i=1:numel(allfiles))
     ind = strfind([foldername '/' allfiles(i).name],'4 dpf');
     if(numel(ind) > 0)
     n =n+1; 
     ExpName = 'PooledStatistics.xlsx'
     % Now we need to define the Names of the sheets:
     % 1) Movement of individual fish:
     SummaryName1{1,1} = 'Minute';
     for(t=1:30)
         SummaryName1{t+2,1} = t;
     end    
     
     SummaryName2{1,2} = 'DMSO';
     SummaryName2{1,2} = 'A1';
     
     
     
     SummaryName2{1,2} = 'EM';
     SummaryName2{1,3} = 'PCB95';
     SummaryName2{1,4} = '0.3';
     SummaryName2{1,5} = '1';
     SummaryName2{1,6} = '3';
     SummaryName2{1,7} = '10';
     SummaryName2{1,8} = '30';
     
     % 2) 
     
     
     
     
     filename = allfiles(i).name;
     sheet = 1;
     xlRange = 'H2890:H2919';
     DMSO_A1 = xlsread(filename,sheet,xlRange);
     
     
     
     
     
     
     
     
     blub = 1;
     
     
     end
end
