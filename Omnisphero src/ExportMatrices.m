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

%This script exports coordinates of all matrices:
function [check] = ExportMatrices ( foldername, csvHandler, neuronHandler)
allfilesOrig = dir(foldername);
n=0;
o=0;
for(j=1:numel(allfilesOrig))
    % Check of parent folder contains png imageso
    indO = strfind([foldername '/' allfilesOrig(j).name],'.png');
    if(numel(indO) > 0)
      n= n+1;
      %Now we have to generate the well name:
      if(n<10)
        a{1} = 'A0';
      else
        a{1} = 'A';
      end 
      a{2} = num2str(n);
      ExpName = allfilesOrig(j).name(5:end-4);
      wellName = [a{1} a{2}];
      NucleusMatrix = csvHandler.CellPosMatrix(wellName);
      Dummy = zeros(size(NucleusMatrix));
      NucleusMatrix = NucleusMatrix + Dummy;
      NucleusMatrix = logical(NucleusMatrix);
       % Here we need to include a try function since it could be a manual
      % or aan automatic experiment!
      try
        AutomaticMatrix = csvHandler.NeuronPositionsEdgeFill(wellName);
      catch
        AutomaticMatrix = 0;  
      end
      try
        ManualMatrix = csvHandler.ManualPositions1(wellName);
      catch
        ManualMatrix = 0;  
      end
      %Name or header
      Header{1,1} = 'Well Name';
      Header{1,2} = 'x All';
      Header{1,3} = 'y All';
      Header{1,4} = 'x Auto';
      Header{1,5} = 'y Auto';
      Header{1,6} = 'x Manual';
      Header{1,7} = 'y Manual';
      
      %Now we have to extract the coordinates an write them in a new
      %comparison matrix. It makes sense to use the NucleusMatrix as
      %reference since it will contain all possibe coordinates:
      NNeurons = sum(sum(NucleusMatrix));
      NeuronsList = regionprops(NucleusMatrix ,'Centroid');
      %Here we will remove the nuclei matrix from the previous well
      MatrixNewX = 0;
      MatrixNewY = 0;
      for(t=1:NNeurons)
          o=o+1;
          xValue = NeuronsList(t).Centroid(1);
          yValue = NeuronsList(t).Centroid(2);
          NeuronList(o+1,1) = xValue;
          NeuronList(o+1,2) = yValue;
          MatrixNewX(t,1) = xValue;
          MatrixNewY(t,1) = yValue;
          NeuronName{o+1,1} = wellName;
      end
      %Now we need to check whether coordinates of manual or automatic
      %algorithms are the same:
      if(size(AutomaticMatrix>0))
            AutomaticMatrix = logical(AutomaticMatrix);
            NNeuronsAuto = sum(sum(AutomaticMatrix));
            try
                NeuronsListAuto = regionprops(AutomaticMatrix,'Centroid');
            catch
                Dummy = zeros(size(NucleusMatrix));
                AutomaticMatrix = Dummy + AutomaticMatrix;
                AutomaticMatrix = logical(AutomaticMatrix);
                NeuronsListAuto = regionprops(AutomaticMatrix,'Centroid');
            end    
            if(n==1)
            [r c] = size(NeuronList);
            end
            for(u=1:NNeuronsAuto);
                xValueAuto = NeuronsListAuto(u).Centroid(1);
                yValueAuto = NeuronsListAuto(u).Centroid(2);
                FindX = find(MatrixNewX==xValueAuto);
                FindY = find(MatrixNewY==yValueAuto);
                if(length(FindX)>1)
                    lengthX = length(FindX);
                    for(e=1:lengthX)
                        if(FindX(e)==FindY)
                            FindXNew = FindX(e);
                        end
                    end 
                    FindX = FindXNew;
                end
                
                if(length(FindY)>1)
                    lengthY = length(FindY);
                    for(e=1:lengthY)
                        if(FindY(e)==FindX)
                            FindYNew = FindY(e);
                        end
                    end 
                    FindY = FindYNew;
                end
                
                if(FindX >0 && FindY>0)
                    if(n==1)
                        NeuronList(FindX+1,3)=xValueAuto;
                        NeuronList(FindX+1,4)=yValueAuto;
                    else
                        NeuronList(FindX+r,3)=xValueAuto;
                        NeuronList(FindX+r,4)=yValueAuto;
                    end
                end
                   
                
            end    
      end
      
      
      if(size(ManualMatrix>0))
            ManualMatrix = ManualMatrix+Dummy;
            ManualMatrix = logical(ManualMatrix);
            NNeuronsManual = sum(sum(ManualMatrix));
            NeuronsListManual = regionprops(ManualMatrix,'Centroid');
            for(u=1:NNeuronsManual);
                xValueManual = NeuronsListManual(u).Centroid(1);
                yValueManual = NeuronsListManual(u).Centroid(2);
                FindXM = find(MatrixNewX==xValueManual);
                FindYM = find(MatrixNewY==yValueManual);
                if(length(FindXM)>1)
                    lengthXM = length(FindXM);
                    for(e=1:lengthXM)
                        if(FindXM(e)==FindYM)
                            FindXMNew = FindXM(e);
                        end
                    end 
                FindXM = FindXMNew;    
                end
                
                if(length(FindYM)>1)
                    lengthYM = length(FindYM);
                    for(e=1:lengthYM)
                        if(FindYM(e)==FindXM)
                            FindYMNew = FindYM(e);
                        end
                    end   
                FindYM = FindYMNew;    
                end  
                
                if(FindXM >0 && FindYM>0)
                    if(n==1)
                       
                        NeuronList(FindXM+1,5)=xValueManual;
                        NeuronList(FindXM+1,6)=yValueManual;    
                    else   
                        
                        NeuronList(FindXM+r,5)=xValueManual;
                        NeuronList(FindXM+r,6)=yValueManual;
                    end
                end
                
            end    
      end
      [r c] = size(NeuronList);
      SummaryHeader{1,1} = 'WellName';
      SummaryHeader {1,2} = 'Number Neurons';
      SummaryHeader  {1,3} = 'Number automatic Neurons';
      SummaryHeader  {1,4} = 'Number manual Neurons';
      SummaryName{n+1,1} = wellName;
      Summary(n+1,2) = NNeurons;
      Summary(n+1,3) = NNeuronsAuto;
      Summary(n+1,4) = NNeuronsManual;
    end
end    
% Now we have to write matrices:
% Add empty coloumn!
[r c] = size(NeuronList);
Dummy = zeros(r,c+1);
Dummy(1:r,2:c+1) = NeuronList;
NeuronList = Dummy;
b{1} = ExpName;
b{2} = '_Coordinates.xlsx';
filename = [b{1} b{2}]
xlswrite(filename,NeuronList);
xlswrite(filename,NeuronName);
xlswrite(filename,Header);
sheet = 2;
xlswrite(filename,Summary,sheet);
xlswrite(filename,SummaryName,sheet);
xlswrite(filename,SummaryHeader,sheet);
check =1;