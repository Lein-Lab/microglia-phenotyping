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

classdef Nucleus
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ID;
        CompositeFillOverlap;
        AlgoNeuron=0;
        TP=0;
        TN=0;
        FP=0;
        DistanceNextNeuron=0;
        ConnectionList;
    end
    
    methods        
    end
    
end

