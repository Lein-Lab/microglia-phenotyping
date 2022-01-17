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

classdef OnlineHelper
    
    methods (Access = public)
        function nucleusBinary = watershedNucleusImage(nucleusBinary,optionHandler)
            nucleusBinary(nucleusBinary<optionHandler.NucleusThreshold) = 0;
            nucleusPic(nucleusPic>=optionHandler.MigDistLowerNucleusThreshold) = 255;
            nucleusPic = uint8(nucleusPic);
            nucleusPicBinBigNuclei = xor(bwareaopen(nucleusPicBin,250),  bwareaopen(nucleusPicBin,15000));
            D = bwdist(~nucleusPicBinBigNuclei);
            D = -D;
            D(~nucleusPicBin) = -Inf;
            L = watershed(D);
            nucleusPicBin(find(~logical(L))) = 0;
            nucleusPicBin = xor(bwareaopen(nucleusPicBin,35),  bwareaopen(nucleusPicBin,15000));
            nucleusBinary = nucleusPicBin;            
        end    
        
        function imageMap=CutOutCircles(imageMap,selectedWell,CutInnerCircle, CutOuterCircle,refreshIM,optionHandler,foldername,sizeY,sizeX, NucleusM, NeuronM, imageHandler)
            %wellList = get(handles.lbWell, 'string');
            %selectedWell = wellList{selectedWellNumber};
            densityWidth = 50;
            densityHeight = 50;
            ringNumber = optionHandler.DensityDistributionRingNumber;
            
            %Get density distribution to exclude too dense areas

                 [filterDistance nonFilterDistance SphereArea markerPointCoordinates] = calculateMigrationDistance(selectedWell, densityWidth, densityHeight, optionHandler, NucleusM, NeuronM, foldername, sizeY, sizeX); 
           
            if(numel(markerPointCoordinates('10')./10) == 0)
                markerPointCoordinates=0;
            end


            if(markerPointCoordinates ~=0 && CutInnerCircle == 1 && CutOuterCircle == 1)    
                innerMask = ones(uint16(sizeY/10),uint16(sizeX/10));
                mk = markerPointCoordinates('10')./10;
                innerMask = roipoly(innerMask,mk(:,1),mk(:,2));
                outerMask = ones(uint16(sizeY/10),uint16(sizeX/10));
                mk = markerPointCoordinates(num2str(ringNumber * 10))./10;
                outerMask = roipoly(innerMask,mk(:,1),mk(:,2));
                %hInner = impoly(gca,double(markerPointCoordinates('10')./10));
                %innerMask = logical(createMask(hInner));
                %delete(hInner);


               % hOuter = impoly(gca,double(markerPointCoordinates(num2str(ringNumber * 10))./10));
               % outerMask = logical(createMask(hOuter));
               % delete(hOuter);
                currentRing = logical(outerMask-innerMask); 
            elseif(CutInnerCircle==1 && CutOuterCircle == 0 && markerPointCoordinates ~=0)
                innerMask = ones(uint16(sizeY/10),uint16(sizeX/10));
                mk = markerPointCoordinates('10')./10;
                innerMask = roipoly(innerMask,mk(:,1),mk(:,2));
                %hInner = impoly(gca,double(markerPointCoordinates('10')./10));
                %innerMask = logical(createMask(hInner));
                %delete(hInner);
                currentRing =logical(ones(int32(sizeY/10),int32(sizeX/10)) - innerMask);
            elseif(markerPointCoordinates ~=0);
                outerMask = ones(uint16(sizeY/10),uint16(sizeX/10));
                mk = markerPointCoordinates(num2str(ringNumber * 10))./10;
                outerMask = roipoly(innerMask,mk(:,1),mk(:,2));
                currentRing = outerMask;
            end
            if(markerPointCoordinates ~=0)
              currentRing = imresize(currentRing, [sizeY, sizeX]);
                if(size(imageMap,1) ~= sizeY)
                  for(i=1:size(imageMap,1))
                    image = imageMap(num2str(i));
                    imageMap(num2str(i)) = uint8(image) .* uint8(currentRing);
                  end
                else
                  imageMap = uint8(imageMap) .* uint8(currentRing);
                end
                else 
                    if(size(imageMap,1) ~= sizeY)
                  for(i=1:size(imageMap,1))
                    image = imageMap(num2str(i));
                    imageMap(num2str(i)) = uint8(image);
                  end
                else
                  imageMap = uint8(imageMap);
                end
            end           
        end            
    end
end