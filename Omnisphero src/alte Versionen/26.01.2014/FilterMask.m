classdef FilterMask < handle
    properties
        PositiveFilters;
        NegativeFilters;
    end
    
    methods (Access = public)    
        %Constructor
        function obj = FilterMask()
            obj.PositiveFilters = containers.Map();
            obj.NegativeFilters = containers.Map();
        end %Constructor
    end
end