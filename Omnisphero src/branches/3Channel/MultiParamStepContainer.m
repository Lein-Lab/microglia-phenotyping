classdef MultiParamStepContainer < handle
    
    properties
        MultiParamStepsDict;
    end
    
    
    methods (Access = public)    
        %Constructor
        function obj = MultiParamStepContainer()            
            obj.MultiParamStepsDict = containers.Map(0,0,'uniformValues',false);
        end %Constructor
    end
end