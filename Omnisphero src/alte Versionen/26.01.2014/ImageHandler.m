classdef ImageHandler < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        NucleusPicArray;
        NeuritePicArray;
        NeuriteImage;
        NucleusImage;
        ResizedNeuriteImage;
        ResizedNucleusImage;
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
    end
     %Singleton Pattern
    methods (Access = public)    
        %Constructor
        function obj = ImageHandler()
            obj.NucleusPicArray = cell(2500);
            obj.NeuritePicArray = cell(2500);
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



