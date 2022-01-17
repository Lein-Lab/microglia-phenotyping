%Copyright (C) 2017-2021  Martin Schmuck, Thomas Temme

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