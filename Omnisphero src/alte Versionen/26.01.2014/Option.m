classdef Option < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        HistogramMaxNucleus;
        HistogramMaxNeurite;
        HistogramMinNucleus;
        HistogramMinNeurite;
        MigDistLowerNucleusThreshold;
        MigDistLowerDensityThreshold;
        MigDistLowerFloodFillThreshold;
        DensityDistributionRingNumber;
        MigrationDistanceDensityImageXSize;
        MigrationDistanceDensityImageYSize;
        NucleusThreshold;
        
        EdgeCompositNeuriteLowerThreshold;
        EdgeCompositMinArea;
        EdgeCompositeDistanceNucleusWhiteArea;
        
        EdgeFillNucleusAreaWithinNeurite;
        EdgeFillNucleusAreaWithinNeuriteSecond;
                
        SkeletonMinNeuriteLength;
        SkeletonThresholdMethod;
        SkeletonNeuriteThresholdHardDistance;
        SkeletonNeuriteThresholdLow;
    end
    
    methods (Access = public)    
        %Constructor
        function obj = Option()            
            %obj.NucleusPicMap = containers.Map();
            %obj.NeuritePicMap = containers.Map();
        end %Constructor
    end
    
end