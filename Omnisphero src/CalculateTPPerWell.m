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


function [TP FP] = CalculateTPPerWell(out)
    TPPerWell = cell2mat(out(:,2));
    FPPerWell = cell2mat(out(:,3));
    TNPerWell = cell2mat(out(:,4));
    %TPPerWell = TPPerWell ./ (TPPerWell+FPPerWell);
    %FPPerWell = FPPerWell ./ (TPPerWell+FPPerWell);
    for(j=1:numel(TPPerWell))
        if(TPPerWell(j)+FPPerWell(j)==0)
            FPPerWell(j)=FPPerWell(j);
            TPPerWell(j)=1;
        else
            FPPerWell(j)=(FPPerWell(j))/(TPPerWell(j)+TNPerWell(j));
            TPPerWell(j) = (TPPerWell(j))/(TPPerWell(j)+TNPerWell(j));
        end
    end

    TP = mean(TPPerWell);
    FP = mean(FPPerWell);