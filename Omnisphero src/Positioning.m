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

diffXCentroidl = x_centroid(1) - x_1;
               diffXCentroidr = x_2 - x_centroid(1);
               diffYCentroidl = y_centroid(1) - y_1;
               diffYCentroidr = y_2 - y_centroid(1);
               if(diffXCentroidl < diffXCentroidr)
               leftBoarder = 1;
               else 
               leftBoarder = 0;
               end
               
               if(diffYCentroidl <diffYCentroidr)
               lowerBoarder = 1;
               else 
               lowerBoarder = 0;
               end 
               
               
               if(leftBoarder == 0)
                   if(lowerBoarder == 0)
                        
                        yStart1 = 1;
                        yEnd1 = r1;
                        xStart1 = 1;
                        xEnd1 = c1;
                   else
                        
                        yStart1 = (7169-r1);
                        yEnd1 = 7168;
                        xStart1 = 1;
                        xEnd1 = c1;
                   end 
               else
                   if(lowerBoarder == 0)
                       
                        yStart1 = 1;
                        yEnd1 = r1;
                        xStart1 = (7169 - c1);
                        xEnd1 = 7168;
                   else
                       
                        yStart1 = (7169-r1);
                        yEnd1 = 7168;
                        xStart1 = (7169 - c1);
                        xEnd1 = 7168;
                   end
               end  
               
               
               imagesphereresize1(yStart1:yEnd1,xStart1:xEnd1)=imagesphere1(1:r1,1:c1);
               
               
               
               
               
               
               
               
               currentNM(yStart1:yEnd1,xStart1:xEnd1) = NucleusM(y_1:y_2,x_1:x_2);