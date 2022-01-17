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

foldername = uigetdir;
allFiles = dir(foldername);
c3RegexpBig = strcat('(.*)Channel3Big.*');

for i=1:numel(allFiles)
    tokensc1 = regexpi(allFiles(i).name, c3RegexpBig, 'tokens');
    if(length(tokensc1) > 0)
        %Flipdim
        %Resize Image
        tokensc1 = tokensc1{1};
        selectedWell = tokensc1{1};
        currentImage=imread([foldername '\' allFiles(i).name]);
        currentImage=flipdim(currentImage,1);
        fileNameBig = [selectedWell 'OligoBig.png'];
        fileNameBig = [foldername '\' fileNameBig];
        imwrite(currentImage,fileNameBig);
        smallImage = imresize(currentImage,0.1);
        fileNameSmall = [foldername '\' selectedWell 'OligoSmall.png'];
        imwrite(smallImage,fileNameSmall);
        delete([foldername '/' allFiles(i).name]);
    end
end