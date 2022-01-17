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

classdef Option < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        FileEnding = '.png';
        HistogramMaxNucleus=255;
        HistogramMaxNeurite=100;
        HistogramMinNucleus=0;
        HistogramMinNeurite=0;
        HistogramMinOligo=0;
        HistogramMinAstro=0;
        HistogramMaxOligo=100;
        HistogramMaxAstro=100;
        MigDistLowerNucleusThreshold;
        MigDistLowerDensityThreshold;
        MigDistLowerFloodFillThreshold;
        DensityDistributionRingNumber;
        MigrationDistanceDensityImageXSize;
        MigrationDistanceDensityImageYSize;
        NucleusThreshold=1;
        ScalingFactor=1;
        NeuronType=0;
        EdgeCompositNeuriteLowerThreshold;
        EdgeCompositMinArea;
        EdgeCompositeDistanceNucleusWhiteArea;
        
        EdgeFillNucleusAreaWithinNeurite=0.55;
        EdgeFillNucleusAreaWithinNeuriteSecond=0.35;
        CompositeFillNeuriteAreaWithinNucleus=0;
        CompositeFillNeuriteAreaWithinNucleusMax=1;
                
        SkeletonMinNeuriteLength=40;
        SkeletonThresholdMethod=3;
        SkeletonNeuriteThresholdHardDistance=0.01;
        SkeletonNeuriteThresholdLow=0.97;
        SkeletonOligoThresholdLow=0.85;
        
        %Maximum distance between endpoint of Skeleton and Neuronnucleus,
        %if Nucleus is not on Neurite.
        MaxDistanceFromEndpoint = 33;
        ToleranceAngleFromEndpoint = 40;
        MinNumberStrongThresholdPixels = 18;
        MinDistanceBetweenNeurons = 23;
        MinSizeNeuriteArea = 2000;
        MinSizeNeuriteAreaCompFill = 100;
        
        %New Parameters for Multi Parameter Analysis
        MaxAllowedFPFirst = 10;
        MaxAllowedFPSecond = 15;
        FilterMulPa=0;
        RectangleSizeA = 1200;
        RectangleSizeB = 800;
        %SkeletonEdgesInSimplifiedSkeletonHigh;        
        FuzzFilterActivated = 1;
        FuzzFilterSecondActivated = 0;
        EdgeFillLookAround=1;
        MaxOverlaySize=999999;
        NumberOfCores=1;
        
        StartWell = 1;
        ExcludedWells='';
        ExcludedWellsMulPa='';
        CutOutSphereCore=0;
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
                case 12
                    value = self.FuzzFilterSecondActivated;
                case 13
                    value = self.NucleusThreshold;
                case 14
                    value = self.EdgeFillLookAround;
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
                case 12
                    self.FuzzFilterSecondActivated = value;
                case 13
                    self.NucleusThreshold = value;
                case 14
                    self.EdgeFillLookAround = value;
            end
        end %function     
        
        
       
        function SetStringValue(self, optionName,optionValueString)
            if(~isempty(str2num(optionValueString)))
                optionValueString=str2num(optionValueString);
            end
            if(strcmp(optionName,'FileEnding'))
                self.FileEnding = optionValueString;
            elseif(strcmp(optionName,'HistogramMaxNucleus'))
                self.HistogramMaxNucleus = optionValueString;
            elseif(strcmp(optionName,'HistogramMinNucleus'))
                self.HistogramMinNucleus = optionValueString;
            elseif(strcmp(optionName,'HistogramMinNeurite'))
                self.HistogramMinNeurite = optionValueString;
            elseif(strcmp(optionName,'HistogramMinOligo'))
                self.HistogramMinOligo = optionValueString;
            elseif(strcmp(optionName,'HistogramMinAstro'))
                self.HistogramMinAstro = optionValueString;
            elseif(strcmp(optionName,'HistogramMaxOligo'))
                self.HistogramMaxOligo = optionValueString;
            elseif(strcmp(optionName,'HistogramMaxAstro'))
                self.HistogramMaxAstro = optionValueString;
            elseif(strcmp(optionName,'MigDistLowerNucleusThreshold'))
                self.MigDistLowerNucleusThreshold = optionValueString;
            elseif(strcmp(optionName,'MigDistLowerDensityThreshold'))
                self.MigDistLowerDensityThreshold = optionValueString;
            elseif(strcmp(optionName,'MigDistLowerFloodFillThreshold'))
                self.MigDistLowerFloodFillThreshold = optionValueString;
            elseif(strcmp(optionName,'DensityDistributionRingNumber'))
                self.DensityDistributionRingNumber = optionValueString;
            elseif(strcmp(optionName,'MigrationDistanceDensityImageXSize'))
                self.MigrationDistanceDensityImageXSize = optionValueString;
            elseif(strcmp(optionName,'MigrationDistanceDensityImageYSize'))
                self.MigrationDistanceDensityImageYSize = optionValueString;
            elseif(strcmp(optionName,'NucleusThreshold'))
                self.NucleusThreshold = optionValueString;
            elseif(strcmp(optionName,'EdgeCompositNeuriteLowerThreshold'))
                self.EdgeCompositNeuriteLowerThreshold = optionValueString;
            elseif(strcmp(optionName,'EdgeCompositMinArea'))
                self.EdgeCompositMinArea = optionValueString;
            elseif(strcmp(optionName,'EdgeCompositeDistanceNucleusWhiteArea'))
                self.EdgeCompositeDistanceNucleusWhiteArea = optionValueString;
            elseif(strcmp(optionName,'EdgeFillNucleusAreaWithinNeurite'))
                self.EdgeFillNucleusAreaWithinNeurite = optionValueString;                
            elseif(strcmp(optionName,'EdgeFillNucleusAreaWithinNeuriteSecond'))
                self.EdgeFillNucleusAreaWithinNeuriteSecond = optionValueString;
            elseif(strcmp(optionName,'CompositeFillNeuriteAreaWithinNucleus'))
                self.CompositeFillNeuriteAreaWithinNucleus = optionValueString;
            elseif(strcmp(optionName,'CompositeFillNeuriteAreaWithinNucleusMax'))
                self.CompositeFillNeuriteAreaWithinNucleusMax = optionValueString;
            elseif(strcmp(optionName,'SkeletonNeuriteThresholdHardDistance'))
                self.SkeletonNeuriteThresholdHardDistance = optionValueString;
            elseif(strcmp(optionName,'SkeletonThresholdMethod'))
                self.SkeletonThresholdMethod = optionValueString;
            elseif(strcmp(optionName,'SkeletonNeuriteThresholdLow'))
                self.SkeletonNeuriteThresholdLow = optionValueString;
            elseif(strcmp(optionName,'SkeletonOligoThresholdLow'))
                self.SkeletonOligoThresholdLow = optionValueString;   
             
            elseif(strcmp(optionName,'MaxDistanceFromEndpoint'))
                self.MaxDistanceFromEndpoint = optionValueString;
            elseif(strcmp(optionName,'ToleranceAngleFromEndpoint'))
                self.ToleranceAngleFromEndpoint = optionValueString;
            elseif(strcmp(optionName,'MinNumberStrongThresholdPixels'))
                self.MinNumberStrongThresholdPixels = optionValueString;
            elseif(strcmp(optionName,'MinDistanceBetweenNeurons'))
                self.MinDistanceBetweenNeurons = optionValueString;
            elseif(strcmp(optionName,'MinSizeNeuriteArea'))
                self.MinSizeNeuriteArea = optionValueString;
            elseif(strcmp(optionName,'MinSizeNeuriteAreaCompFill'))
                self.MinSizeNeuriteAreaCompFill = optionValueString;
            elseif(strcmp(optionName,'MaxAllowedFPFirst'))
                self.MaxAllowedFPFirst = optionValueString;
            elseif(strcmp(optionName,'MaxAllowedFPSecond'))
                self.MaxAllowedFPSecond = optionValueString;
            elseif(strcmp(optionName,'FilterMulPa'))
                self.FilterMulPa = optionValueString;
            elseif(strcmp(optionName,'RectangleSizeA'))
                self.RectangleSizeA = optionValueString;
            elseif(strcmp(optionName,'RectangleSizeB'))
                self.RectangleSizeB = optionValueString;
            elseif(strcmp(optionName,'FuzzFilterActivated'))
                self.FuzzFilterActivated = optionValueString;
            elseif(strcmp(optionName,'FuzzFilterSecondActivated'))
                self.FuzzFilterSecondActivated = optionValueString;
            elseif(strcmp(optionName,'EdgeFillLookAround'))
                self.EdgeFillLookAround = optionValueString;
            elseif(strcmp(optionName,'MaxOverlaySize'))
                self.MaxOverlaySize = optionValueString;
            elseif(strcmp(optionName,'NumberOfCores'))
                self.NumberOfCores = optionValueString;
            elseif(strcmp(optionName,'StartWell'))
                self.StartWell = optionValueString;
            elseif(strcmp(optionName,'ExcludedWells'))
                self.ExcludedWells = optionValueString;
            elseif(strcmp(optionName,'ExcludedWellsMulPa'))
                self.ExcludedWellsMulPa = optionValueString;
            elseif(strcmp(optionName,'CutOutSphereCore'))
                self.CutOutSphereCore = optionValueString;            
            end
        end%function
                


    end    
end