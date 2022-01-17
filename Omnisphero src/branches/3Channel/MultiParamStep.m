classdef MultiParamStep < handle
    properties
        OptionHandler;
        TPWells;
        FPWells;
        TNWells;
        ChangedDim;
        T;
    end
    
    
    methods (Access = public)    
        %Constructor
        function obj = MultiParamStep()                                   
        end %Constructor
        
        
        function [TP, FP] = GetOverallValues(self)
            TP = 0;
            FP = 0;
            for(i=1:numel(self.TPWells.keys))
                TP = TP + TPWells(TPWells.keys(i));
                FP = FP + FPWells(FPWells.keys(i));
            end
        end
    end
end