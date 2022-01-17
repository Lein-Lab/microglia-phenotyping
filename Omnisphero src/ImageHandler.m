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

classdef ImageHandler < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        NucleusActive=1;
        NeuriteActive=1;
        OligoActive=0;
        AstroActive=0;
        NucleusPicArray;
        NeuritePicArray;
        C3PicArray;
        C4PicArray;
        NeuriteImage;
        NucleusImage;
        OligoImage;
        AstroImage;
        ResizedNeuriteImage;
        ResizedNucleusImage;
        ResizedOligoImage;
        ResizedAstroImage;
        BinaryImage;
        ResizedBinaryImage;
        SkeletonImage;
        ResizedSkeletonImage;
        ShowNucleusImage;
        ShowNeuriteImage;
        MousePosition;
        SphereMasks;
        SphereArray;
        ExtractSphereMode;
        %CurrentSphere;
        Foldername;
        ZoomState;
        PositiveFilters;
        NegativeFilters;
        MediumMigDistances;
        CheckSquareMode;
        MeasureMigDistMode;
        MigDistPoints;
        MigDistWellMapping;
        CurrentSquareMask;
        PreviousWell;
        ExamplePicturesXIndices;
        ExamplePicturesYIndices;
        ExamplePictureMap;
        MultiParamStepContainer;
        
    end
     %Singleton Pattern
    methods (Access = public)    
        %Constructor
        function obj = ImageHandler()
            obj.NucleusPicArray = cell(2500);
            obj.NeuritePicArray = cell(2500);
            obj.C3PicArray = cell(2500);
            obj.SphereArray = cell(1);
            obj.CheckSquareMode=0;
            obj.MeasureMigDistMode=0;
            obj.MigDistWellMapping = containers.Map();
            
            %obj.NucleusPicMap = containers.Map();
            %obj.NeuritePicMap = containers.Map();
        end %Constructor
    end
    
    methods (Static)
        function inst = getInstance()
            persistent myObj;
            if(isempty(myObj))
                myObj = ImageHandler();
                inst = myObj;
            else
                inst = myObj;
            end
        end
    end
    
    methods
        function GetSize(this)
        props = properties(this);
        totSize = 0;
        for ii=1:length(props)
            currentProperty = getfield(this, char(props(ii)));
            s = whos('currentProperty');
            totSize = totSize + s.bytes;
        end
        fprintf(1, '%d bytes\n', totSize);
        end
    end
end



