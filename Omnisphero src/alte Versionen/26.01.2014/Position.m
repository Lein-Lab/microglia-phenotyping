classdef Position
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x;
        y;
    end
    
    methods
        function [obj,idx] = sort(obj,varargin{2})
            [~,idx]=sort([obj.x],varargin{:});
            obj = obj(idx);
        end
    end
    
end

