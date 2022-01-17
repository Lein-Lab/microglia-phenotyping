classdef Option < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        FileEnding;
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
        SkeletonThresholdMethod=3;
        SkeletonNeuriteThresholdHardDistance;
        SkeletonNeuriteThresholdLow;
        
        %Maximum distance between endpoint of Skeleton and Neuronnucleus,
        %if Nucleus is not on Neurite.
        MaxDistanceFromEndpoint = 33;
        ToleranceAngleFromEndpoint = 40;
        MinNumberStrongThresholdPixels = 18;
        MinDistanceBetweenNeurons = 23;
        MinSizeNeuriteArea = 2000;
        
        %New Parameters for Multi Parameter Analysis
        MaxAllowedFP = 15;
        RectangleSizeA = 1000;
        RectangleSizeB = 500;
        %SkeletonEdgesInSimplifiedSkeletonHigh;        
        FuzzFilterActivated = 1;
        
        StartWell = 1;
        ExcludedWells='';
    end
    
    methods (Access = public)    
        %Constructor
        function F = Option(rhs)            
          if(nargin > 0)
            % copy constructor
            fns = properties(rhs);
            for i=1:length(fns)
               F.(fns{i}) = rhs.(fns{i});
            end
          end
        end %Constructor
        
        
        function value = GetValue(self, optionNumber)
            switch optionNumber
                case 1        
                    value=self.EdgeFillNucleusAreaWithinNeurite;
                case 2
                    value=self.EdgeFillNucleusAreaWithinNeuriteSecond;                
                case 3
                    value = self.SkeletonMinNeuriteLength;                
                case 4
                    value = self.MaxDistanceFromEndpoint;                
                case 5
                    value = self.ToleranceAngleFromEndpoint;     
                case 6
                    value = self.MinNumberStrongThresholdPixels;     
                case 7
                    value = self.MinDistanceBetweenNeurons;     
                case 8
                    value = self.MinSizeNeuriteArea;   
                case 9
                    value = self.SkeletonNeuriteThresholdLow;
                case 10
                    value = self.FuzzFilterActivated;
                case 11
                    value = self.NucleusThreshold;
            end
        end %function
        
        function SetValue(self, optionNumber,value)
            switch optionNumber
                case 1        
                    self.EdgeFillNucleusAreaWithinNeurite=value;
                case 2
                    self.EdgeFillNucleusAreaWithinNeuriteSecond=value;                
                case 3
                    self.SkeletonMinNeuriteLength = value;                
                case 4
                    self.MaxDistanceFromEndpoint = value;                
                case 5
                    self.ToleranceAngleFromEndpoint = value;
                case 6
                    self.MinNumberStrongThresholdPixels = value;
                case 7
                    self.MinDistanceBetweenNeurons = value;
                case 8
                    self.MinSizeNeuriteArea = value;
                case 9
                    self.SkeletonNeuriteThresholdLow = value;
                case 10
                    self.FuzzFilterActivated = value;
                case 11
                    self.NucleusThreshold = value;
            end
        end %function       

    end    
end