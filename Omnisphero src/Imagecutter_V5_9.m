%Copyright (C) 2016-2021  Martin Schmuck

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

function[csvHandler wellList] = Imagecutter_V5_9(Neuron_Channel, Oligo_Channel, Astro_Channel, foldername, csvHandler, optionHandler)

%function csvHandler = Imagecutter_V5_8(Neuron_Channel, Oligo_Channel, Astro_Channel, foldername, csvHandler)
% This function is used to sepearte neurospheres from images of
% 8-chamber slide into seperate images, containing only one neurosphere.
% Therefore centroids of neurosphere cores were determined and distance
% between the boarders of the coverslip as well as distances among
% centroids are calculated. For both x and y coordinates the shortest
% distance among border and centroid centroid distances is either added or
% substracted from the centroid coordinate to generate a rectangular filter
% to cut out single spheres. In a last step those images are resized to a
% 96 well plate image.
n = 0;
wellList = cell(0);
% N gives the minimum required area for a neurosphere core
N_1 = 2000;
N_2 = 100000;
%Neuron_Channel = 0;
%Oligo_Channel = 0;
%Astro_Channel = 0;
% The MinimumDistance defines the minimum distance of two sphere cores for
% which migration areas should not intersect
MinimumDistance = 7168/2;
MinDistCentr = 7168/2;
MinDistanceCores = 1000;
%This is the maximum added or substracted from the center point. This
%should, if all spheres are at least 7168 pixels away from each other
%result in no overlapp of migration areas (3101 equals 2000�M)
maxaddedlength = 3101;
foldername = [foldername '/ConvertedCellomics'];
mkdir(foldername,'OT');
foldername1 = [foldername];
foldername2 = [foldername1 '/OT'];


%Well renamer:

% Sometimes the difference between x or y coordinates is very small. This
% is the result of two spheres with minimum distance smaller the defined
% positioned on the line (e.g. on a horizontal line or the upper two
% spheres of a five). Therefore a minimum ratio of x/y and y/x is defined
% as an additional criteria. To rule out the possibility that spheres are
% almost touching each other, further a mindist is defined at which the xy
% and yx cirtiria is applied:
minxy_ratio = 0.25;
mindist = 3000;
%Defines the source folder
%foldername = uigetdir;
%Define destination folder for saving processed images
%foldername1 = uigetdir;
%Iterate over all files in this folder
allfiles = dir(foldername);
% Iteration loop for all images in source folder
for(i=1:numel(allfiles))
    %Check if file ends with NucleusBig.png
    ind = strfind([foldername '/' allfiles(i).name],'NucleusBig.png');
        if(numel(ind) > 0)
            n=n+1;
            NucleusM = csvHandler.CellPosMatrix(strcat(allfiles(i).name(1),allfiles(i).name(3)));      
            if(Neuron_Channel>0);
                % Looks for cooresponding neurite image by replacing NucleusBig by
                % NeuriteBig in the name of image i
                newFileNeu = regexprep(allfiles(i).name,'NucleusBig','NeuriteBig','ignorecase'); 
                % New file is defined as imgneurite
                %imgneurite = imread([foldername '/' newFileNeu]);
            end
            
            %Same for OligoBig
            if(Oligo_Channel>0);
                newFileOli = regexprep(allfiles(i).name,'NucleusBig','OligoBig','ignorecase');        
                %imgoligo = imread([foldername '/' newFileOli]);
            end
                    
        %Same for AstroBig
            if(Astro_Channel >0)
                newFileAst = regexprep(allfiles(i).name,'NucleusBig','AstroBig','ignorecase');        
                %imgastro = imread([foldername '/' newFileAst]);
            end    
        
        % Loads NucleusBig  
        imgnucleus = imread([foldername '/' allfiles(i).name]);
        [lehisto x] = imhist(imgnucleus);
        level = triangle_th(lehisto,256);
        if(level>0.1)
            level = 0.099;
        end    
        imgnucleusbinary = im2bw(imgnucleus,level*2);
        %newFileWatershed = regexprep(allfiles(i).name,'NucleusBig','NucleusBigWatershed','ignorecase');
        %imgWatershed = imread([foldername '/' newFileWatershed]);
        %imgnucleusbinarysizedwater = bwareaopen(imgWatershed,N_1);
        %imgnucleusbinarysizedwater = bwmorph(imgnucleusbinarysizedwater, 'bridge');
        %imgnucleusbinarysizedwater = bwareaopen(imgnucleusbinarysizedwater,N_2);
        % Generates histogram of imgnucleus
        %[lehisto x] = imhist(imgnucleus);
        % Calculates background of imgnucleus
        %level = triangle_th(lehisto,256);
        % Converts imgnucleus to a binary 'imgnucleusbinary' image using
        % level as threshold
        %imgnucleusbinary = im2bw(imgnucleus,level);
        clear imgnucleus;
        % Convertes 'imgnucleusbinary' to an 8-bit image
        %imgnucleusbinary = uint8(imgnucleusbinary.*255);
        % Excludes all areas in 'imgnucleusbinary' smaller than 250000
        %imgnucleusbinarysizedwater = bwareaopen(imgWatershed,N);
        %clear imgnucleusbinary;
        % Excludes objects touching the rim of the cover slip to get rid of
        % bright edges of the cover slip
        imgnucleusbinarysizedwater = imclearborder(imgnucleusbinary);
        clear  imgnucleusbinary;
        imgnucleusbinarysizedwater = bwareaopen(imgnucleusbinarysizedwater,N_2);
        
        
        % Calculates all centroid coordinates of all areas in
        % 'imgnucleusbinarysized'
        sC = regionprops(imgnucleusbinarysizedwater,'centroid');
        sA = regionprops(imgnucleusbinarysizedwater,'area');
        
        %Function to delte big particles close to the sphere.
        s_LC = length(sC);
        nn = 0;
        s = 0;
        for(j=1:s_LC)
            x_centroid_1 = sC(j).Centroid(1);
            y_centroid_1 = sC(j).Centroid(2);
            
            
            for(o=1:s_LC)
                
                x_centroid_n = sC(o).Centroid(1);
                y_centroid_n = sC(o).Centroid(2);    
                diffCentroids = ((x_centroid_n - x_centroid_1)^2 + (y_centroid_n - y_centroid_1)^2)^0.5;
                if(diffCentroids < MinDistanceCores && diffCentroids ~= 0)
                    if(sA(o).Area(1) < sA(j).Area(1))
                    x_centroid_n = 0;
                    y_centroid_n = 0;
                    else
                    x_centroid_1 = 0;
                    y_centroid_1 = 0;
                    end
                end    
                
                   
            
            
                   
            end
            if(x_centroid_1>0)
                nn = nn + 1;
                s(nn,1) = x_centroid_1;
                s(nn,2) = y_centroid_1;
            end
        end
        
        % Determines the y-dimension of s
        s_L = length(s); 
        
        clear imgnucleusbinarysizedwater;
        
        % Case destinction is required, since through washing spheres can
        % be lost
        %Case 1: All spheres are lost, meaning no centroids and a y
        %dimension of s equals 0
            if(strfind(allfiles(i).name,  'D01'))
                wellind = '01';
            end   
            if(strfind(allfiles(i).name,  'C01'))
                wellind = '02';
            end 
            if(strfind(allfiles(i).name,  'B01'))
                wellind = '03';
            end 
            if(strfind(allfiles(i).name,  'A01'))
                wellind = '04';
            end 
            if(strfind(allfiles(i).name,  'D02'))
                wellind = '05';
            end 
            if(strfind(allfiles(i).name,  'C02'))
                wellind = '06';
            end 
            if(strfind(allfiles(i).name,  'B02'))
                wellind = '07';
            end 
            if(strfind(allfiles(i).name,  'A02'))
                wellind = '08';
            end
        if(s_L <=0);
            %Well renamer:
            
            % Creation of five epmty images with dimensions 7168x7168    
            imagesphereresize1 = zeros(7168,7168);
            newFile1 = regexprep(allfiles(i).name,'NucleusBig',['E' wellind 'NucleusBig'],'ignorecase');        
            imwrite(imagesphereresize1,[foldername1 '/' newFile1]);
            clear imagesphereresize1; 
            
            imagesphereresize2 = zeros(7168,7168);
            newFile2 = regexprep(allfiles(i).name,'NucleusBig',['F' wellind 'NucleusBig'],'ignorecase');        
            imwrite(imagesphereresize2,[foldername1 '/' newFile2]);
            clear imagesphereresize2; 
            
            imagesphereresize3 = zeros(7168,7168);
            newFile3 = regexprep(allfiles(i).name,'NucleusBig',['G' wellind 'NucleusBig'],'ignorecase');        
            imwrite(imagesphereresize3,[foldername1 '/' newFile3]);
            clear imagesphereresize3; 
            
            imagesphereresize4 = zeros(7168,7168);
            newFile4 = regexprep(allfiles(i).name,'NucleusBig',['H' wellind 'NucleusBig'],'ignorecase');        
            imwrite(imagesphereresize4,[foldername1 '/' newFile4]);
            clear imagesphereresize4; 
            
            imagesphereresize5 = zeros(7168,7168);
            newFile5 = regexprep(allfiles(i).name,'NucleusBig',['I' wellind 'NucleusBig'],'ignorecase');        
            imwrite(imagesphereresize5,[foldername1 '/' newFile5]);
            clear imagesphereresize5;
            
            if(Neuron_Channel > 0)
                   % Same for neurite channel


                   imagesphereresize6 = zeros(7168,7168);
            	   newFile6 = regexprep(allfiles(i).name,'NucleusBig',['E' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize6,[foldername1 '/' newFile6]);
                   clear imagesphereresize6;

                   imagesphereresize7 = zeros(7168,7168);
            	   newFile7 = regexprep(allfiles(i).name,'NucleusBig',['F' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize7,[foldername1 '/' newFile7]);
                   clear imagesphereresize7;

                   imagesphereresize8 = zeros(7168,7168);
            	   newFile8 = regexprep(allfiles(i).name,'NucleusBig',['G' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize8,[foldername1 '/' newFile8]);
                   clear imagesphereresize8;


                   imagesphereresize9 = zeros(7168,7168);
            	   newFile9 = regexprep(allfiles(i).name,'NucleusBig',['H' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize9,[foldername1 '/' newFile9]);
                   clear imagesphereresize9;

                   imagesphereresize10 = zeros(7168,7168);
            	   newFile10 = regexprep(allfiles(i).name,'NucleusBig',['I' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize10,[foldername1 '/' newFile10]);
                   clear imagesphereresize10;

               end
               
               if(Oligo_Channel > 0)
                   % Same for Oligos
                   imagesphereresize11 = zeros(7168,7168);
            	   newFile11 = regexprep(allfiles(i).name,'NucleusBig',['E' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize11,[foldername1 '/' newFile11]);
                   clear imagesphereresize11;

                   imagesphereresize12 = zeros(7168,7168);
            	   newFile12 = regexprep(allfiles(i).name,'NucleusBig',['F' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize12,[foldername1 '/' newFile12]);
                   clear imagesphereresize12;

                   imagesphereresize13 = zeros(7168,7168);
            	   newFile13 = regexprep(allfiles(i).name,'NucleusBig',['G' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize13,[foldername1 '/' newFile13]);
                   clear imagesphereresize13;

                   imagesphereresize14 = zeros(7168,7168);
            	   newFile14 = regexprep(allfiles(i).name,'NucleusBig',['H' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize14,[foldername1 '/' newFile14]);
                   clear imagesphereresize14;

                   imagesphereresize15 = zeros(7168,7168);
            	   newFile15 = regexprep(allfiles(i).name,'NucleusBig',['I' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize15,[foldername1 '/' newFile15]);
                   clear imagesphereresize15;
               end
               
               if(Astro_Channel >0)
                   % Same for Astros
                   imagesphereresize16 = zeros(7168,7168);
            	   newFile16 = regexprep(allfiles(i).name,'NucleusBig',['E' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize16,[foldername1 '/' newFile16]);
                   clear imagesphereresize16;

                   imagesphereresize17 = zeros(7168,7168);
            	   newFile17 = regexprep(allfiles(i).name,'NucleusBig',['F' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize17,[foldername1 '/' newFile17]);
                   clear imagesphereresize17;

                   imagesphereresize18 = zeros(7168,7168);
            	   newFile18 = regexprep(allfiles(i).name,'NucleusBig',['G' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize18,[foldername1 '/' newFile18]);
                   clear imagesphereresize18;

                   imagesphereresize19 = zeros(7168,7168);
            	   newFile19 = regexprep(allfiles(i).name,'NucleusBig',['H' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize19,[foldername1 '/' newFile19]);
                   clear imagesphereresize19;

                   imagesphereresize20 = zeros(7168,7168);
            	   newFile20 = regexprep(allfiles(i).name,'NucleusBig',['I' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize20,[foldername1 '/' newFile20]);
                   clear imagesphereresize20;
               end
               
            else
            end
        % Case 2: One sphere remained: s=1     
            if(0<s_L && s_L<2);
               % Extracts the x and y coordinates from s 
               x_centroid(1) = s(1,1);
               y_centroid(1) = s(1,2);
               
               % Calculates the distance between the centroid of the
               % neurosphere to all boarders of the cover slip
               
               % Distance between centroid and left x-boarder
               d_x0_1 = x_centroid(1)-1;
               % Distance between centroid and right x-boarder
               d_x10752_1 = 10752 - x_centroid(1);
               % Distance between centroid and upper y-boarder
               d_y0_1 = y_centroid(1)-1;
               % Distance between centroid and lower x-boarder
               d_y10752_1 = 13824 - y_centroid(1);
               
               %Generation of an vector containing all distances in x and y for  the sphere
               %core. 7168/2 is inserted as the maximal allowed value to
               %add or to substract
               % V_x1m contains all values left of the centroid
               V_x1m = [d_x0_1;(7168/2)];
               % V_x1p contains all values right of the centroid
               V_x1p = [d_x10752_1;(7168/2)];
               % V_y1m contains all values above the centroid
               V_y1m = [d_y0_1;(7168/2)];
               % V_y1p contains all values below the centroid
               V_y1p = [d_y10752_1;(7168/2)];
               
               
               % Determination of the smallest distance in all directions
               % of V
               minx_1m = min(V_x1m(1:2));
               minx_1p = min(V_x1p(1:2));
               miny_1m = min(V_y1m(1:2));
               miny_1p = min(V_y1p(1:2));
               
               % Calculates the x and y coordinates of the rectangular
               % around the sphere core. Always minimal distances were
               % substrated or added
               x_1 = x_centroid(1) - minx_1m;
               x_2 = x_centroid(1) + minx_1p - 1;
               y_1 = y_centroid(1) - miny_1m;
               y_2 = y_centroid(1) + miny_1p - 1;
               
               % Add decision on which side the black border will be added
               % later! Idea is that the rim will be added to the longer
               % side!
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
               
               imgnucleus = imread([foldername '/' allfiles(i).name]);
               %imwrite(imgnucleus,[foldername2 '/' allfiles(i).name]);
               %The rectangular image is cutted from the imgnucleus image containing the whole cover slip
               imagesphere1 = imgnucleus(y_1:y_2,x_1:x_2);
               % Dimensions of imagesphere1 are determined
               [r1 c1] = size(imagesphere1);
               % An empty matrix is generated
               imagesphereresize1 = zeros(7168,7168);
               % imagesphereresize1 is converted to 8-bit
               imagesphereresize1 = uint8(imagesphereresize1);
               %The empty indices in the imagesphereresize1 matrix are overwritten with imagesphere 1 entries
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
               % imagesphereresize1 is saved to final destination
               newFile1 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize1,[foldername1 '/' newFile1]);
               wellname = (newFile1(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(7168,7168);
               currentNM = csvHandler.CellPosMatrix(wellname);
                currentNM(yStart1:yEnd1,xStart1:xEnd1) = NucleusM(y_1:y_2,x_1:x_2);
               csvHandler.CellPosMatrix(wellname) = currentNM;
               imagesphereresize1_small = imresize(imagesphereresize1, optionHandler.ScalingFactor);
               newFile1s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize1_small,[foldername1 '/' newFile1s]);
               
               clear imagesphereresize1; 
               clear imagesphere1;
               clear imagesphereresize1_small
                
               imagesphereresize2 = zeros(7168,7168);
               newFile2 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize2,[foldername1 '/' newFile2]);
               imagesphereresize2_small = imresize(imagesphereresize2, optionHandler.ScalingFactor);
               newFile2s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize2_small,[foldername1 '/' newFile2s]);
               n = n + 1;
               wellname = (newFile2(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(zeros(7168,7168));
               clear imagesphereresize2; 
               clear imagesphereresize2_small
               
               
               imagesphereresize3 = zeros(7168,7168);
               newFile3 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize3,[foldername1 '/' newFile3]);
               imagesphereresize3_small = imresize(imagesphereresize3, optionHandler.ScalingFactor);
               newFile3s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize3_small,[foldername1 '/' newFile3s]);
               n = n + 1;
               wellname = (newFile3(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(zeros(7168,7168));
               clear imagesphereresize3; 
               clear imagesphereresize3_small
               imagesphereresize4 = zeros(7168,7168);
               newFile4 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize4,[foldername1 '/' newFile4]);
               imagesphereresize4_small = imresize(imagesphereresize4, optionHandler.ScalingFactor);
               newFile4s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize4_small,[foldername1 '/' newFile4s]);
               n = n + 1;
               wellname = (newFile4(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(zeros(7168,7168));
               clear imagesphereresize4; 
               clear imagesphereresize4_small
               imagesphereresize5 = zeros(7168,7168);
               newFile5 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize5,[foldername1 '/' newFile5]);
               imagesphereresize5_small = imresize(imagesphereresize5, optionHandler.ScalingFactor);
               newFile5s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize5_small,[foldername1 '/' newFile5s]);
               n = n + 1;
               wellname = (newFile5(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(zeros(7168,7168));
               
               clear imagesphereresize5;   
               clear imagesphereresize5_small;
               clear imgnucleus;
               
               
               if(Neuron_Channel > 0)
                   % Same for neurite channel
                        
                   imgneurite = imread([foldername '/' newFileNeu]);
                   %imwrite(imgneurite,[foldername2 '/' newFileNeu]);
                   imagesphere6 = imgneurite(y_1:y_2,x_1:x_2);
                   [r1 c1] = size(imagesphere6);
                   imagesphereresize6 = zeros(size(imagesphere6));
                   imagesphereresize6 = uint8(imagesphereresize6);
                   imagesphereresize6(yStart1:yEnd1,xStart1:xEnd1)=imagesphere6(1:r1,1:c1);
                   newFile6 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize6,[foldername1 '/' newFile6]);
                   imagesphereresize6_small = imresize(imagesphereresize6, optionHandler.ScalingFactor);
                   newFile6s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize6_small,[foldername1 '/' newFile6s]);
                   clear imagesphereresize6;
                   clear imagesphere6;
                   clear imagesphereresize6_small
                   imagesphereresize7 = zeros(7168,7168);
            	   newFile7 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize7,[foldername1 '/' newFile7]);
                   imagesphereresize7_small = imresize(imagesphereresize7, optionHandler.ScalingFactor);
                   newFile7s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize7_small,[foldername1 '/' newFile7s]);
                   clear imagesphereresize7;
                   clear imagesphereresize7_small
                   imagesphereresize8 = zeros(7168,7168);
            	   newFile8 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize8,[foldername1 '/' newFile8]);
                   imagesphereresize8_small = imresize(imagesphereresize8, optionHandler.ScalingFactor);
                   newFile8s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize8_small,[foldername1 '/' newFile8s]);
                   clear imagesphereresize8;
                   clear imagesphereresize8_small

                   imagesphereresize9 = zeros(7168,7168);
            	   newFile9 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize9,[foldername1 '/' newFile9]);
                   imagesphereresize9_small = imresize(imagesphereresize9, optionHandler.ScalingFactor);
                   newFile9s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize9_small,[foldername1 '/' newFile9s]);
                   clear imagesphereresize9;
                   clear imagesphereresize9_small
                   imagesphereresize10 = zeros(7168,7168);
            	   newFile10 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize10,[foldername1 '/' newFile10]);
                   imagesphereresize10_small = imresize(imagesphereresize10, optionHandler.ScalingFactor);
                   newFile10s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize10_small,[foldername1 '/' newFile10s]);
                   clear imagesphereresize10;
                   clear imagesphereresize10_small
                   clear imgneurite;
                  
               end
               if(Oligo_Channel > 0)
                   % Same for Oligo channel

                   imgoligo = imread([foldername '/' newFileOli]);
                   %imwrite(imgoligo,[foldername2 '/' newFileOli]);
                   imagesphere11 = imgoligo(y_1:y_2,x_1:x_2);
                   [r1 c1] = size(imagesphere11);
                   imagesphereresize11 = zeros(size(imagesphere11));
                   imagesphereresize11 = uint8(imagesphereresize11);
                   imagesphereresize11(yStart1:yEnd1,xStart1:xEnd1)=imagesphere11(1:r1,1:c1);
                   newFile11 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize11,[foldername1 '/' newFile11]);
                   imagesphereresize11_small = imresize(imagesphereresize11, optionHandler.ScalingFactor);
                   newFile11s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize11_small,[foldername1 '/' newFile11s]);
                   clear imagesphereresize11;
                   clear imagesphere11;
                   clear imagesphereresize11_small
                   imagesphereresize12 = zeros(7168,7168);
            	   newFile12 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize12,[foldername1 '/' newFile12]);
                   imagesphereresize12_small = imresize(imagesphereresize12, optionHandler.ScalingFactor);
                   newFile12s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize12_small,[foldername1 '/' newFile12s]);
                   clear imagesphereresize12;
                   clear imagesphereresize12_small
                   imagesphereresize13 = zeros(7168,7168);
            	   newFile13 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize13,[foldername1 '/' newFile13]);
                   imagesphereresize13_small = imresize(imagesphereresize13, optionHandler.ScalingFactor);
                   newFile13s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize13_small,[foldername1 '/' newFile13s]);
                   clear imagesphereresize13;
                   clear imagesphereresize13_small

                   imagesphereresize14 = zeros(7168,7168);
            	   newFile14 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize14,[foldername1 '/' newFile14]);
                   imagesphereresize14_small = imresize(imagesphereresize14, optionHandler.ScalingFactor);
                   newFile14s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize14_small,[foldername1 '/' newFile14s]);
                   clear imagesphereresize14;
                   clear imagesphereresize14_small
                   imagesphereresize15 = zeros(7168,7168);
            	   newFile15 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize15,[foldername1 '/' newFile15]);
                   imagesphereresize15_small = imresize(imagesphereresize15, optionHandler.ScalingFactor);
                   newFile15s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize15_small,[foldername1 '/' newFile15s]);
                   clear imagesphereresize15;
                   clear imagesphereresize15_small
                   clear imgoligo
               end
               if(Astro_Channel > 0)
                   % Same for Astro channel

                   imgastro = imread([foldername '/' newFileAst]);
                   %imwrite(imgastro,[foldername2 '/' newFileAst]);
                   imagesphere16 = imgastro(y_1:y_2,x_1:x_2);
                   [r1 c1] = size(imagesphere16);
                   imagesphereresize16 = zeros(size(imagesphere16));
                   imagesphereresize16 = uint8(imagesphereresize16);
                   imagesphereresize16(yStart1:yEnd1,xStart1:xEnd1)=imagesphere16(1:r1,1:c1);
                   newFile16 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize16,[foldername1 '/' newFile16]);
                   imagesphereresize16_small = imresize(imagesphereresize16, optionHandler.ScalingFactor);
                   newFile16s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize16_small,[foldername1 '/' newFile16s]);
                   clear imagesphereresize16;
                   clear imagesphere16;
                   clear imagesphereresize16_small
                   imagesphereresize17 = zeros(7168,7168);
            	   newFile17 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize17,[foldername1 '/' newFile17]);
                   imagesphereresize17_small = imresize(imagesphereresize17, optionHandler.ScalingFactor);
                   newFile17s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize17_small,[foldername1 '/' newFile17s]);
                   clear imagesphereresize17;
                   clear imagesphereresize17_small
                   imagesphereresize18 = zeros(7168,7168);
            	   newFile18 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize18,[foldername1 '/' newFile18]);
                   imagesphereresize18_small = imresize(imagesphereresize18, optionHandler.ScalingFactor);
                   newFile18s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize18_small,[foldername1 '/' newFile18s]);
                   clear imagesphereresize18;
                   clear imagesphereresize18_small

                   imagesphereresize19 = zeros(7168,7168);
            	   newFile19 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize19,[foldername1 '/' newFile19]);
                   imagesphereresize19_small = imresize(imagesphereresize19, optionHandler.ScalingFactor);
                   newFile19s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize19_small,[foldername1 '/' newFile19s]);
                   clear imagesphereresize19;
                   clear imagesphereresize19_small
                   imagesphereresize20 = zeros(7168,7168);
            	   newFile20 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize20,[foldername1 '/' newFile20]);
                   imagesphereresize20_small = imresize(imagesphereresize20, optionHandler.ScalingFactor);
                   newFile20s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize20_small,[foldername1 '/' newFile20s]);
                   clear imagesphereresize20;
                   clear imagesphereresize20_small
                   clear imgastro;
               end
               
                %delete(allfiles(i).name);
                %if(Neuron_Channel > 0)
                %delete(newFileNeu);
                %end
                %if(Oligo_Channel > 0)
                %delete(newFileOli);
                %end
                %if(Astro_Channel > 0)
                %delete(newFileAst);
                %end
            end    
        % Case 3: Two spheres remained: s=2    
            if(s_L>1 && s_L<3);
               % Extracts the x and y coordinates from s 
               x_centroid(1) = s(1,1);
               y_centroid(1) = s(1,2)
               x_centroid(2) = s(2,1);
               y_centroid(2) = s(2,2);
               % Calculates the distance between the centroid of the
               % neurosphere and all boarders of the cover slip
               
               % Distance between centroid(1) and centroid(2)
               d1_2 = 0.5* ((x_centroid(1)-x_centroid(2))^2+(y_centroid(1)-y_centroid(2))^2)^0.5;
               % Distance between centroid and left x-boarder
               d_x0_1 = x_centroid(1)-1;
               % Distance between centroid and right x-boarder
               d_x10752_1 = 10752 - x_centroid(1);
               % Distance between centroid and upper y-boarder
               d_y0_1 = y_centroid(1)-1;
               % Distance between centroid and lower y-boarder
               d_y10752_1 = 13824 - y_centroid(1);
               % Same for sphere core 2
               d_x0_2 = x_centroid(2)-1;
               d_x10752_2 = 10752 - x_centroid(2);
               d_y0_2 = y_centroid(2)-1;
               d_y10752_2 = 13824 - y_centroid(2);
               %Generation of an vector containing all distances in x and y for  the sphere
               %core. 7168/2 is inserted as the maximal allowed value to
               %add or to substract
               % V_x1m contains all values left of the centroid
               V_x1m = [d1_2;d_x0_1;(7168/2)];
               % V_x1p contains all values right of the centroid
               V_x1p = [d1_2;d_x10752_1;(7168/2)];
               % V_y1m contains all values above the centroid
               V_y1m = [d1_2;d_y0_1;(7168/2)];
               % V_y1p contains all values below the centroid
               V_y1p = [d1_2;d_y10752_1;(7168/2)];
               % Same for sphere core 2
               V_x2m = [d1_2;d_x0_2;(7168/2)];
               V_x2p = [d1_2;d_x10752_2;(7168/2)];
               V_y2m = [d1_2;d_y0_2;(7168/2)];
               V_y2p = [d1_2;d_y10752_2;(7168/2)];
               
                             
               % Determination of the smallest distance in all directions
               % of V
               
               minx_1m = min(V_x1m(1:3));
               minx_1p = min(V_x1p(1:3));
               miny_1m = min(V_y1m(1:3));
               miny_1p = min(V_y1p(1:3));
               minx_2m = min(V_x2m(1:3));
               minx_2p = min(V_x2p(1:3));
               miny_2m = min(V_y2m(1:3));
               miny_2p = min(V_y2p(1:3));
               
               
                imgnucleus = imread([foldername '/' allfiles(i).name]);
               %imwrite(imgnucleus,[foldername2 '/' allfiles(i).name]);
               % Calculates the x and y coordinates of the rectangular
               % around the sphere core. Always minimal distances were
               % substrated or added
               x_1 = x_centroid(1) - minx_1m;
               x_2 = x_centroid(1) + minx_1p;
               y_1 = y_centroid(1) - miny_1m;
               y_2 = y_centroid(1) + miny_1p;
               
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
               
               
               
               
               %The rectangular image is cutted from the imgnucleus image containing the whole cover slip
               imagesphere1 = imgnucleus(y_1:y_2,x_1:x_2);
               % Dimensions of imagesphere1 are determined
               [r1 c1] = size(imagesphere1);
               % An empty matrix is generated
               imagesphereresize1 = zeros(7168,7168);
               % imagesphereresize1 is converted to 8-bit
               imagesphereresize1 = uint8(imagesphereresize1);
               %The empty indices in the imagesphereresize1 matrix are overwritten with imagesphere 1 entries
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
               % imagesphereresize1 is saved to final destination
               newFile1 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize1,[foldername1 '/' newFile1]);
               wellname = (newFile1(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(7168,7168);
               currentNM = csvHandler.CellPosMatrix(wellname);
               currentNM(yStart1:yEnd1,xStart1:xEnd1) = NucleusM(y_1:y_2,x_1:x_2);
               csvHandler.CellPosMatrix(wellname) = currentNM;
               imagesphereresize1_small = imresize(imagesphereresize1, optionHandler.ScalingFactor);
               newFile1s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize1_small,[foldername1 '/' newFile1s]);
               clear imagesphereresize1_small
               clear imagesphereresize1; 
               clear imagesphere1;
               
                           
               % Same for second sphere
               
                
               x_3 = x_centroid(2) - minx_2m;
               x_4 = x_centroid(2) + minx_2p -1;
               y_3 = y_centroid(2) - miny_2m;
               y_4 = y_centroid(2) + miny_2p -1;
               
               
               diffXCentroidl = x_centroid(2) - x_3;
               diffXCentroidr = x_4 - x_centroid(2);
               diffYCentroidl = y_centroid(2) - y_3;
               diffYCentroidr = y_4 - y_centroid(2);
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
               
               imagesphere2 = imgnucleus(y_3:y_4,x_3:x_4);
               [r2 c2] = size(imagesphere2);
               imagesphereresize2 = zeros(7168,7168);
               imagesphereresize2 = uint8(imagesphereresize2);
               if(leftBoarder == 0)
                   if(lowerBoarder == 0)
                        
                        yStart2 = 1;
                        yEnd2 = r2;
                        xStart2 = 1;
                        xEnd2 = c2;
                   else
                        
                        yStart2 = (7169-r2);
                        yEnd2 = 7168;
                        xStart2 = 1;
                        xEnd2 = c2;
                   end 
               else
                   if(lowerBoarder == 0)
                       
                        yStart2 = 1;
                        yEnd2 = r2;
                        xStart2 = (7169 - c2);
                        xEnd2 = 7168;
                   else
                       
                        yStart2 = (7169-r2);
                        yEnd2 = 7168;
                        xStart2 = (7169 - c2);
                        xEnd2 = 7168;
                   end
               end  
               
               
               imagesphereresize2(yStart2:yEnd2,xStart2:xEnd2)=imagesphere2(1:r2,1:c2);
               newFile2 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize2,[foldername1 '/' newFile2]);
               n = n + 1;
               wellname = (newFile2(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(7168,7168);
               currentNM = csvHandler.CellPosMatrix(wellname);
               currentNM(yStart2:yEnd2,xStart2:xEnd2) = NucleusM(y_3:y_4,x_3:x_4);
               csvHandler.CellPosMatrix(wellname) = currentNM;
               imagesphereresize2_small = imresize(imagesphereresize2, optionHandler.ScalingFactor);
               newFile2s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize2_small,[foldername1 '/' newFile2s]);
               clear imagesphereresize2_small
               clear imagesphereresize2;
               clear imagesphere2;
               
               imagesphereresize3 = zeros(7168,7168);
               newFile3 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize3,[foldername1 '/' newFile3]);
               imagesphereresize3_small = imresize(imagesphereresize3, optionHandler.ScalingFactor);
               newFile3s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize3_small,[foldername1 '/' newFile3s]);
               n = n + 1;
               wellname = (newFile3(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(zeros(7168,7168));
               clear imagesphereresize3_small
               clear imagesphereresize3; 

               imagesphereresize4 = zeros(7168,7168);
               newFile4 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize4,[foldername1 '/' newFile4]);
                imagesphereresize4_small = imresize(imagesphereresize4, optionHandler.ScalingFactor);
               newFile4s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize4_small,[foldername1 '/' newFile4s]);
               n = n + 1;
               wellname = (newFile4(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(zeros(7168,7168));
               clear imagesphereresize4_small
               clear imagesphereresize4;
                

               imagesphereresize5 = zeros(7168,7168);
               newFile5 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NucleusBig'],'ignorecase');               
               imwrite(imagesphereresize5,[foldername1 '/' newFile5]);
               imagesphereresize5_small = imresize(imagesphereresize5, optionHandler.ScalingFactor);
               newFile5s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize5_small,[foldername1 '/' newFile5s]);
               n = n + 1;
               wellname = (newFile5(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(zeros(7168,7168));
               clear imagesphereresize5_small
               clear imagesphereresize5; 
               clear imgnucleus
               
               if(Neuron_Channel > 0)
                   % Same for neurite channel

                   imgneurite = imread([foldername '/' newFileNeu]);
                   %imwrite(imgneurite,[foldername2 '/' newFileNeu]);
                   imagesphere6 = imgneurite(y_1:y_2,x_1:x_2);
                   [r1 c1] = size(imagesphere6);
                   imagesphereresize6 = zeros(7168,7168);
                   imagesphereresize6 = uint8(imagesphereresize6);
                   imagesphereresize6(yStart1:yEnd1,xStart1:xEnd1)=imagesphere6(1:r1,1:c1);
                   newFile6 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize6,[foldername1 '/' newFile6]);
                   imagesphereresize6_small = imresize(imagesphereresize6, optionHandler.ScalingFactor);
                   newFile6s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize6_small,[foldername1 '/' newFile6s]);
                   clear imagesphere6;
                   clear imagesphereresize6;
                   clear imagesphereresize6_small

                   imagesphere7 = imgneurite(y_3:y_4,x_3:x_4);
                   [r2 c2] = size(imagesphere7);
                   imagesphereresize7 = zeros(7168,7168);
                   imagesphereresize7 = uint8(imagesphereresize7);
                   imagesphereresize7(yStart2:yEnd2,xStart2:xEnd2)=imagesphere7(1:r2,1:c2);
                   newFile7 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize7,[foldername1 '/' newFile7]);
                   imagesphereresize7_small = imresize(imagesphereresize7, optionHandler.ScalingFactor);
                   newFile7s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize7_small,[foldername1 '/' newFile7s]);
                   clear imagesphere7;
                   clear imagesphereresize7;
                   clear imagesphereresize7_small

                   imagesphereresize8 = zeros(7168,7168);
            	   newFile8 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize8,[foldername1 '/' newFile8]);
                   imagesphereresize8_small = imresize(imagesphereresize8, optionHandler.ScalingFactor);
                   newFile8s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize8_small,[foldername1 '/' newFile8s]);
                   clear imagesphereresize8;
                   clear imagesphereresize8_small


                   imagesphereresize9 = zeros(7168,7168);
            	   newFile9 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize9,[foldername1 '/' newFile9]);
                   imagesphereresize9_small = imresize(imagesphereresize9, optionHandler.ScalingFactor);
                   newFile9s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize9_small,[foldername1 '/' newFile9s]);
                   clear imagesphereresize9;
                   clear imagesphereresize9_small

                   imagesphereresize10 = zeros(7168,7168);
            	   newFile10 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize10,[foldername1 '/' newFile10]);
                   imagesphereresize10_small = imresize(imagesphereresize10, optionHandler.ScalingFactor);
                   newFile10s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize10_small,[foldername1 '/' newFile10s]);
                   clear imagesphereresize10;
                   clear imagesphereresize10_small
                   clear imgneurite;

               end
               
               if(Oligo_Channel > 0)
                   
                   imgoligo = imread([foldername '/' newFileOli]);
                   %imwrite(imgoligo,[foldername2 '/' newFileOli]);
                   % Same for Oligos
                   imagesphere11 = imgoligo(y_1:y_2,x_1:x_2);
                   [r1 c1] = size(imagesphere11);
                   imagesphereresize11 = zeros(7168,7168);
                   imagesphereresize11 = uint8(imagesphereresize11);
                   imagesphereresize11(yStart1:yEnd1,xStart1:xEnd1)=imagesphere11(1:r1,1:c1);
                   newFile11 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize11,[foldername1 '/' newFile11]);
                   imagesphereresize11_small = imresize(imagesphereresize11, optionHandler.ScalingFactor);
                   newFile11s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize11_small,[foldername1 '/' newFile11s]);
                   clear imagesphereresize11;
                   clear imagesphere11;
                   clear imagesphereresize11_small

                   imagesphere12 = imgoligo(y_3:y_4,x_3:x_4);
                   [r2 c2] = size(imagesphere12);
                   imagesphereresize12 = zeros(7168,7168);
                   imagesphereresize12 = uint8(imagesphereresize12);
                   imagesphereresize12(yStart2:yEnd2,xStart2:xEnd2)=imagesphere12(1:r2,1:c2);
                   newFile12 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize12,[foldername1 '/' newFile12]);
                   imagesphereresize12_small = imresize(imagesphereresize12, optionHandler.ScalingFactor);
                   newFile12s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize12_small,[foldername1 '/' newFile12s]);
                   clear imagesphereresize12;
                   clear imagesphere12;
                   clear imagesphereresize12_small

                   imagesphereresize13 = zeros(7168,7168);
            	   newFile13 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize13,[foldername1 '/' newFile13]);
                   imagesphereresize13_small = imresize(imagesphereresize13, optionHandler.ScalingFactor);
                   newFile13s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize13_small,[foldername1 '/' newFile13s]);
                   clear imagesphereresize13;
                   clear imagesphereresize13_small;

                   imagesphereresize14 = zeros(7168,7168);
            	   newFile14 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize14,[foldername1 '/' newFile14]);
                   imagesphereresize14_small = imresize(imagesphereresize14, optionHandler.ScalingFactor);
                   newFile14s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize14_small,[foldername1 '/' newFile14s]);
                   clear imagesphereresize14;
                   clear imagesphereresize14_small;

                   imagesphereresize15 = zeros(7168,7168);
            	   newFile15 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize15,[foldername1 '/' newFile15]);
                   imagesphereresize15_small = imresize(imagesphereresize15, optionHandler.ScalingFactor);
                   newFile15s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize15_small,[foldername1 '/' newFile15s]);
                   clear imagesphereresize15;
                   clear imagesphereresize15_small;
                   clear imgoligo;
               end
               
               if(Astro_Channel >0)
                   % Same for Astros
                   imgastro = imread([foldername '/' newFileAst]);
                   %imwrite(imgastro,[foldername2 '/' newFileAst]);
                   imagesphere16 = imgastro(y_1:y_2,x_1:x_2);
                   [r1 c1] = size(imagesphere16);
                   imagesphereresize16 = zeros(7168,7168);
                   imagesphereresize16 = uint8(imagesphereresize16);
                   imagesphereresize16(yStart1:yEnd1,xStart1:xEnd1)=imagesphere16(1:r1,1:c1);
                   newFile16 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize16,[foldername1 '/' newFile16]);
                   imagesphereresize16_small = imresize(imagesphereresize16, optionHandler.ScalingFactor);
                   newFile16s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize16_small,[foldername1 '/' newFile16s]);
                   clear imagesphereresize16;
                   clear imagesphere16;
                   clear imagesphereresize16_small

                   imagesphere17 = imgastro(y_3:y_4,x_3:x_4);
                   [r2 c2] = size(imagesphere17);
                   imagesphereresize17 = zeros(7168,7168);
                   imagesphereresize17 = uint8(imagesphereresize17);
                   imagesphereresize17(yStart2:yEnd2,xStart2:xEnd2)=imagesphere17(1:r2,1:c2);
                   newFile17 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize17,[foldername1 '/' newFile17]);
                   imagesphereresize17_small = imresize(imagesphereresize17, optionHandler.ScalingFactor);
                   newFile17s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize17_small,[foldername1 '/' newFile17s]);
                   clear imagesphereresize17;
                   clear imagesphere17;
                   clear imagesphereresize17_small

                   imagesphereresize18 = zeros(7168,7168);
            	   newFile18 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize18,[foldername1 '/' newFile18]);
                   imagesphereresize18_small = imresize(imagesphereresize18, optionHandler.ScalingFactor);
                   newFile18s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize18_small,[foldername1 '/' newFile18s]);
                   clear imagesphereresize18;
                   clear imagesphereresize18_small

                   imagesphereresize19 = zeros(7168,7168);
            	   newFile19 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize19,[foldername1 '/' newFile19]);
                   imagesphereresize19_small = imresize(imagesphereresize19, optionHandler.ScalingFactor);
                   newFile19s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize19_small,[foldername1 '/' newFile19s]);
                   clear imagesphereresize19;
                   clear imagesphereresize19_small

                   imagesphereresize20 = zeros(7168,7168);
            	   newFile20 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize20,[foldername1 '/' newFile20]);
                   imagesphereresize20_small = imresize(imagesphereresize20, optionHandler.ScalingFactor);
                   newFile20s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize20_small,[foldername1 '/' newFile20s]);
                   clear imagesphereresize20;
                   clear imagesphereresize20_small
                   clear imgastro

               end
               
               
                
            %delete(allfiles(i).name);
                %if(Neuron_Channel > 0)
                %delete(newFileNeu);
                %end
                %if(Oligo_Channel > 0)
                %delete(newFileOli);
                %end
                %if(Astro_Channel > 0)
                %delete(newFileAst);
                %end
                
            end
            
        % Case 4: Three spheres remained: s=3
        
            if(s_L >2 && s_L <4);
               % Extracts the x and y coordinates from s 
               x_centroid(1) = s(1,1);
               y_centroid(1) = s(1,2);
               x_centroid(2) = s(2,1);
               y_centroid(2) = s(2,2);
               x_centroid(3) = s(3,1);
               y_centroid(3) = s(3,2);
               % Definition of a centroid quadrant in which the centroid of a potential center neurosphere has to be located. 
               %All other centroids are considerd as rim coordinates 
               quadrant1_x1 = 4032;
               quadrant1_x2 = 6720;
               quadrant1_y1 = 5184;
               quadrant1_y2 = 8640;
               
               %Calculates the distance between the centroids of the
               % three neurosphere 
               d1_2 = 0.5* ((x_centroid(1)-x_centroid(2))^2+(y_centroid(1)-y_centroid(2))^2)^0.5;
               d1_3 = 0.5* ((x_centroid(1)-x_centroid(3))^2+(y_centroid(1)-y_centroid(3))^2)^0.5;
               d2_3 = 0.5* ((x_centroid(3)-x_centroid(2))^2+(y_centroid(3)-y_centroid(2))^2)^0.5;
               
               Vd1 = [d1_2;d1_3]
               Vd2 = [d1_2;d2_3]
               Vd3 = [d1_3;d2_3]
               % The shortest distance between neurosphere cores is the
               % distance between the rim spheres and the center sphere.
               % Therefore for all centroids it is checked wether the coordinates of 
               % their centroid is located within the center quadrant.
               if(x_centroid(1)>=quadrant1_x1 && x_centroid(1)<=quadrant1_x2 && y_centroid(1)>=quadrant1_y1 && y_centroid(1)<=quadrant1_y2);
               % Since more than one sphere might be located in the
               % centroid quadrant, further the distance between center
               % spheres and the middle of the cover slip (dc_x) is
               % calculated     
                    dc_1 = ((x_centroid(1)-5376)^2+(y_centroid(1)-6912)^2)^0.5;
                    
                else
                % If the centroid is not located in the center quadrant the distance is set to an abitrary high value of dc_1 = 1000000000    
                    dc_1 = 1000000000;
                end 
                % Same for sphere 2 and 3
                if(x_centroid(2)>=quadrant1_x1 && x_centroid(2)<=quadrant1_x2 && y_centroid(2)>=quadrant1_y1 && y_centroid(2)<=quadrant1_y2);
                    dc_2 = ((x_centroid(2)-5376)^2+(y_centroid(2)-6912)^2)^0.5;
                    
                else
                    
                    dc_2 = 1000000000;
                end
                
                if(x_centroid(3)>=quadrant1_x1 && x_centroid(3)<=quadrant1_x2 && y_centroid(3)>=quadrant1_y1 && y_centroid(3)<=quadrant1_y2);
                    dc_3 = ((x_centroid(3)-5376)^2+(y_centroid(3)-6912)^2)^0.5;
                    
                else
                    
                    dc_3 = 1000000000;
                end
               % It is checked wether dc_1 is smaller than 1000000000 (no centroid sphere according to the above checked criteria) and smaller
               % than the distances of other centroid spheres. If dc_1 is
               % the smallest distance the centroid of sphere 1 is defined
               % as the centroid of the center sphere
               if(dc_1 < 1000000000 && dc_1 < dc_2 && dc_1 < dc_3);
                   x_centroid_c = x_centroid(1);
                   y_centroid_c = y_centroid(1);
                   missingc1 = 0;
                 
                   
               else
               % If above criteria are not fullfilled the centroid sphere
               % is missing and the centroid of this sphere is treated as
               % a rim sphere
                   missingc1 = 1;
               end    
               %Same for sphere 2 and 3 
               if(dc_2 < 1000000000 && dc_2 < dc_1 && dc_2 < dc_3);
                   x_centroid_c = x_centroid(2);
                   y_centroid_c = y_centroid(2);
                   missingc2 = 0;
               else
                   missingc2 = 1;
               end              
                
               if(dc_3 < 1000000000 && dc_3 < dc_1 && dc_3 < dc_2);
                   x_centroid_c = x_centroid(3);
                   y_centroid_c = y_centroid(3);
                   missingc3 = 0;
           
               else
                   missingc3 = 1;
               end
               % If all 3 spheres were defined as rim spheres, an artifical
               % center point is defined as the center point of the cover
               % slip
               if((missingc1+missingc2+missingc3) > 2);
                   x_centroid_c = 5376;
                   y_centroid_c = 6912;
               else
               end 
               
               % If sphere 1 is no center sphere the shortest distance
               % between sphere1 and the center sphere is calculated
               if(x_centroid(1) ~= x_centroid_c);
                   d_1_c = 0.5* (((x_centroid(1)-x_centroid_c)^2+(y_centroid(1)-y_centroid_c)^2)^0.5);
               else
               end   
               % Centroid of sphere 2 is located between sphere 1 and 3
               if(x_centroid(2) ~= x_centroid_c);
               % If sphere 2 is no center sphere the shortest distance
               % between sphere2 and the center sphere is calculated    
                   d_2_c = 0.5* (((x_centroid(2)-x_centroid_c)^2+(y_centroid(2)-y_centroid_c)^2)^0.5);
               else
               end
               % Centroid of sphere 3 is located right of sphere 1 and 2
               if(x_centroid(3) ~= x_centroid_c);
                   d_3_c = 0.5* (((x_centroid(3)-x_centroid_c)^2+(y_centroid(3)-y_centroid_c)^2)^0.5);
               else
               end
               
            
               % Defines wether centroid of sphere 1 lays left or right of the center
               % sphere! If sphere 1 is not the center sphere
                if(x_centroid(1) ~= x_centroid_c);
                   % First case: Check whether x_centroid(1) is left of the center sphere 
                   if(x_centroid(1) < x_centroid_c);
                      % Check if y_centroid(1) is below y_centroid_c 
                      if(y_centroid(1)<y_centroid_c);
                       % If y_centroid(1)<y_centroid_c the sphere has to be located at the
                       % upper left corner of the cover slip 
                       x_centroid_l(1) = x_centroid(1);
                       y_centroid_l(1) = y_centroid(1);
                       % Sphere is located in the upper left of the cover
                       % slip
                       l_1_1 = 1;
                       % Sphere is located in the lower left of the cover
                       % slip
                       l_1_2 = 0;
                       % Sphere is located in the upper right of the cover
                       % slip
                       r_1_4 = 0;
                       % Sphere is located in the lower right of the cover
                       % slip
                       r_1_5 = 0;
                       
                      else
                       % If y_centroid(1)>y_centroid_c the sphere has to be located at the
                       % lower left corner of the cover slip  
                       x_centroid_l(2) = x_centroid(1);
                       y_centroid_l(2) = y_centroid(1);
                       l_1_1 = 0
                       l_1_2 = 1
                       r_1_4 = 0
                       r_1_5 = 0 
                      end
                   % Second case: Check whether x_centroid(1) is right of the center sphere    
                   else
                      % Check if y_centroid(1) is below y_centroid_c)                  
                      if(y_centroid(1)<y_centroid_c) 
                       x_centroid_r(4) = x_centroid(1);
                       y_centroid_r(4) = y_centroid(1);
                       % If y_centroid(1)<y_centroid_c the sphere has to be located at the
                       % upper right corner of the cover slip
                       l_1_1 = 0
                       l_1_2 = 0
                       r_1_4 = 1
                       r_1_5 = 0
                      else
                      % If y_centroid(1)>y_centroid_c the sphere has to be located at the
                       % lower right corner of the cover slip     
                       x_centroid_r(5) = x_centroid(1);
                       y_centroid_r(5) = x_centroid(1);
                       l_1_1 = 0
                       l_1_2 = 0
                       r_1_4 = 0
                       r_1_5 = 1 
                      end 
                   end    
                else
                    l_1_1 = 0
                    l_1_2 = 0
                    r_1_4 = 0
                    r_1_5 = 0
                end    
               % Same for sphere 2 and 3
               if(x_centroid(2) ~= x_centroid_c)
                   if(x_centroid(2) < x_centroid_c);
                      if(y_centroid(2)<y_centroid_c) 
                       x_centroid_l(1) = x_centroid(2);
                       y_centroid_l(1) = y_centroid(2);
                       l_2_1 = 1
                       l_2_2 = 0
                       r_2_4 = 0
                       r_2_5 = 0
                       
                      else
                       x_centroid_l(2) = x_centroid(2);
                       y_centroid_l(2) = y_centroid(2);
                       l_2_1 = 0
                       l_2_2 = 1
                       r_2_4 = 0
                       r_2_5 = 0 
                      end 
                   else
                      if(y_centroid(2)<y_centroid_c) 
                       x_centroid_r(4) = x_centroid(2);
                       y_centroid_r(4) = y_centroid(2);
                       l_2_1 = 0
                       l_2_2 = 0
                       r_2_4 = 1
                       r_2_5 = 0
                      else
                       x_centroid_r(5) = x_centroid(2);
                       y_centroid_r(5) = x_centroid(2);
                       l_2_1 = 0
                       l_2_2 = 0
                       r_2_4 = 0
                       r_2_5 = 1 
                      end 
                   end    
                else
                    l_2_1 = 0
                    l_2_2 = 0
                    r_2_4 = 0
                    r_2_5 = 0
                end 
               
               if(x_centroid(3) ~= x_centroid_c)
                   if(x_centroid(3) < x_centroid_c);
                      if(y_centroid(3)<y_centroid_c) 
                       x_centroid_l(1) = x_centroid(3);
                       y_centroid_l(1) = y_centroid(3);
                       l_3_1 = 1
                       l_3_2 = 0
                       r_3_4 = 0
                       r_3_5 = 0
                       
                      else
                       x_centroid_l(2) = x_centroid(3);
                       y_centroid_l(2) = y_centroid(3);
                       l_3_1 = 0
                       l_3_2 = 1
                       r_3_4 = 0
                       r_3_5 = 0 
                      end 
                   else
                      if(y_centroid(3)<y_centroid_c) 
                       x_centroid_r(4) = x_centroid(3);
                       y_centroid_r(4) = y_centroid(3);
                       l_3_1 = 0
                       l_3_2 = 0
                       r_3_4 = 1
                       r_3_5 = 0
                      else
                       x_centroid_r(5) = x_centroid(3);
                       y_centroid_r(5) = x_centroid(3);
                       l_3_1 = 0
                       l_3_2 = 0
                       r_3_4 = 0
                       r_3_5 = 1 
                      end 
                   end    
               else
                    l_3_1 = 0
                    l_3_2 = 0
                    r_3_4 = 0
                    r_3_5 = 0
               end
               
               
               % Check for double entry positions:
               
               %Check if sphere 1,2 and 3 have the same position 
               if(l_2_1 == l_1_1 && l_2_2 == l_1_2 && r_2_4 == r_1_4 && r_2_5 == r_1_5 && l_3_1 == l_1_1 && l_3_2 == l_1_2 && r_3_4 == r_1_4 && r_3_5 == r_1_5)
                   
                   l_1_1_t = l_1_1
                   l_1_2_t = l_1_2
                   r_1_4_t = r_1_4
                   r_1_5_t = r_1_5
                   l_2_1_t = l_2_1
                   l_2_2_t = l_2_2
                   r_2_4_t = r_2_4
                   r_2_5_t = r_2_5
                   l_3_1_t = l_3_1
                   l_3_2_t = l_3_2
                   r_3_4_t = r_3_4
                   r_3_5_t = r_3_5
                   l_1_1 = 0
                   l_1_2 = 0
                   r_1_4 = 0
                   r_1_5 = 0
                   l_2_1 = 0
                   l_2_2 = 0
                   r_2_4 = 0
                   r_2_5 = 0
                   l_3_1 = 0
                   l_3_2 = 0
                   r_3_4 = 0
                   r_3_5 = 0
                   
                   
               else
                   l_1_1_t = 0
                   l_1_2_t = 0
                   r_1_4_t = 0
                   r_1_5_t = 0
                   l_2_1_t = 0
                   l_2_2_t = 0
                   r_2_4_t = 0
                   r_2_5_t = 0
                   l_3_1_t = 0
                   l_3_2_t = 0
                   r_3_4_t = 0
                   r_3_5_t = 0 
                   
               end    
               
               % Check if sphere 1 and 2 share the same position
               if(l_2_1 == l_1_1 && l_2_2 == l_1_2 && r_2_4 == r_1_4 && r_2_5 == r_1_5)
                 
                   l_1_1_d12 = l_1_1
                   l_1_2_d12 = l_1_2
                   r_1_4_d12 = r_1_4
                   r_1_5_d12 = r_1_5
                   l_2_1_d12 = l_2_1
                   l_2_2_d12 = l_2_2
                   r_2_4_d12 = r_2_4
                   r_2_5_d12 = r_2_5
                   l_3_1_d12 = 0
                   l_3_2_d12 = 0
                   r_3_4_d12 = 0
                   r_3_5_d12 = 0
                   l_1_1 = 0
                   l_1_2 = 0
                   r_1_4 = 0
                   r_1_5 = 0
                   l_2_1 = 0
                   l_2_2 = 0
                   r_2_4 = 0
                   r_2_5 = 0
                                              
               else
                   l_1_1_d12 = 0
                   l_1_2_d12 = 0
                   r_1_4_d12 = 0
                   r_1_5_d12 = 0
                   l_2_1_d12 = 0
                   l_2_2_d12 = 0
                   r_2_4_d12 = 0
                   r_2_5_d12 = 0 
                   
                       
               end
                   
               % Check if sphere 1 and 3 share the same position
               if(l_3_1 == l_1_1 && l_3_2 == l_1_2 && r_3_4 == r_1_4 && r_3_5 == r_1_5)
                 
                   l_1_1_d13 = l_1_1
                   l_1_2_d13 = l_1_2
                   r_1_4_d13 = r_1_4
                   r_1_5_d13 = r_1_5
                   l_2_1_d13 = 0
                   l_2_2_d13 = 0
                   r_2_4_d13 = 0
                   r_2_5_d13 = 0
                   l_3_1_d13 = l_3_1
                   l_3_2_d13 = l_3_2
                   r_3_4_d13 = r_3_4
                   r_3_5_d13 = r_3_5
                   l_1_1 = 0
                   l_1_2 = 0
                   r_1_4 = 0
                   r_1_5 = 0
                   l_3_1 = 0
                   l_3_2 = 0
                   r_3_4 = 0
                   r_3_5 = 0
                                              
               else
                   l_1_1_d13 = 0
                   l_1_2_d13 = 0
                   r_1_4_d13 = 0
                   r_1_5_d13 = 0
                   l_3_1_d13 = 0
                   l_3_2_d13 = 0
                   r_3_4_d13 = 0
                   r_3_5_d13 = 0 
                   
                       
                   end
                       
                   
                 
               
               
               if(l_2_1 == l_3_1 && l_2_2 == l_3_2 && r_2_4 == r_3_4 && r_2_5 == r_3_5)
                   B = 1
                   l_2_1_d23 = l_2_1
                   l_2_2_d23 = l_2_2
                   r_2_4_d23 = r_2_4
                   r_2_5_d23 = r_2_5
                   l_3_1_d23 = l_3_1
                   l_3_2_d23 = l_3_2
                   r_3_4_d23 = r_3_4
                   r_3_5_d23 = r_3_5
                   l_2_1 = 0
                   l_2_2 = 0
                   r_2_4 = 0
                   r_2_5 = 0
                   l_3_1 = 0
                   l_3_2 = 0
                   r_3_4 = 0
                   r_3_5 = 0
                   
               else
                   l_2_1_d23 = 0
                   l_2_2_d23 = 0
                   r_2_4_d23 = 0
                   r_2_5_d23 = 0
                   l_3_1_d23 = 0
                   l_3_2_d23 = 0
                   r_3_4_d23 = 0
                   r_3_5_d23 = 0
                   
               end    
                  
               
               
                            
               
               
               % Sphere 1 with smallest value of x_centroid
               
               
               
               % Go through potential positions
               % Case 1: Sphere is no center sphere and is located at the
               % upper left
               if(l_1_1 >0 && l_1_1 < 2);
                       % Coordinates to cut out sphere 1
                       x_1 = 1;
                       x_2 = x_centroid(1) + d_1_c;
                       y_1 = 1;
                       y_2 = y_centroid(1) + d_1_c;
               else
               end    
               % Case 2: Sphere is no center sphere and is located at the
               % lower left
               if(l_1_2 >0 && l_1_2 < 2);
                       % Coordinates to cut out sphere 1
                       x_1 = 1;
                       x_2 = x_centroid(1) + d_1_c;
                       y_1 = y_centroid(1) - d_1_c;
                       y_2 = 13824;
               else
               end
               % Case 3: Sphere is no center sphere and is located at the
               % upper right
               if(r_1_4 >0 && r_1_4 < 2);
                       % Coordinates to cut out sphere 1
                       x_1 = x_centroid(1) - d_1_c;
                       x_2 = 10752;
                       y_1 = 1;
                       y_2 = y_centroid(1) + d_1_c;
               else
               end
               % Case 4: Sphere is no center sphere and is located at the
               % lower right
               if(r_1_5 >0 && r_1_5 < 2);
                       % Coordinates to cut out sphere 1 
                       x_1 = x_centroid(1) - d_1_c;
                       x_2 = 10752;
                       y_1 = y_centroid(1) - d_1_c;
                       y_2 = 13824;
               else
               end
               % Going through all iterations for the case sphere 1 =
               % center sphere
               if(x_centroid(1) ~= x_centroid_c);
               else
                   % Case 1: sphere 2 and 3 are located left of sphere 1
                   if(x_centroid(2) < x_centroid(1) && x_centroid(3) < x_centroid(1));
                       % Find shortest distance to sphere 1 to sphere 2 and
                       % 3 to a maximum of 7168/2. Since both sphere 2 and sphere 3
                       % are located left of sphere 1, distance is just
                       % limited to the left and not to the right
                       V_1m = [d1_2 ; d1_3 ; (7168/2)];
                       min1m = min(V_1m(1:3));
                       % Maximal added distance to the right to obtain
                       % maximal images of 7168x7168 pixel
                       min1p = (7168-min1m);
                   else
                   end
                   % Case 2: sphere 2 and 3 are located right of sphere 1
                   if(x_centroid(2) > x_centroid(1) && x_centroid(3) > x_centroid(1));
                       % Find shortest distance to sphere 1 to sphere 2 and
                       % 3 to a maximum of 7168/2. Since both sphere 2 and sphere 3
                       % are located right of sphere 1, distance is just
                       % limited to the right and not to the left
                       V_1p = [d1_2 ; d1_3 ; (7168/2)];
                       min1p = min(V_1p(1:3));
                       min1m = (7168-min1p);
                   else
                   end
                   % Case 3: sphere 2 is located right of sphere 1 and 3 is located left of sphere 1
                   if(x_centroid(2) > x_centroid(1) && x_centroid(3) < x_centroid(1));
                       % Find minimal distance to the left and to the right. To the left
                       % only the distance between sphere 1 and 2 and to
                       % the right only the distance of sphere 1 and 3
                       % matter.                     
                       V_1m = [d1_3 ; (7168/2)];
                       V_1p = [d1_2 ; (7168/2)];
                       min1p = min(V_1p(1:2));
                       min1m = min(V_1m(1:2));
                   else
                   end
                   % Case 3: sphere 2 is located left of sphere 1 and 3 is located right of sphere 1
                   if(x_centroid(2) < x_centroid(1) && x_centroid(3) > x_centroid(1));
                       % Find minimal distance to the left and to the right. To the left
                       % only the distance between sphere 1 and 3 and to
                       % the right only the distance of sphere 1 and 2
                       % matter. 
                       V_1m = [d1_2 ; (7168/2)];
                       V_1p = [d1_3 ; (7168/2)];
                       min1p = min(V_1p(1:2));
                       min1m = min(V_1m(1:2));
                   else
                   end
              % Coordinates to cut out sphere 1  
                   x_1 = x_centroid(1) - min1m 
                   x_2 = x_centroid(1) + min1p
                   y_1 = y_centroid(1) - min1m
                   y_2 = y_centroid(1) + min1p
               end    
               
               % Cases for having two or three spheres at one location
               
               % Sphere 1 and 2 are in the same location
               
               if(l_1_1_d12 == 1)
                  if(y_centroid(1) < y_centroid(2))
                  x_1 = 1;
                  x_2 = x_centroid(1) + d1_2;
                  y_1 = 1;
                  y_2 = y_centroid(1) + d1_2; 
                  else
                  x_1 = 1;
                  x_2 = x_centroid(1) + d1_2;
                  y_1 = y_centroid(1) - d1_2;
                  y_2 = y_centroid(1) + (7168 - d1_2);    
                  end    
               end   
               
               if(l_1_2_d12 == 1)
                  if(y_centroid(1) < y_centroid(2))
                  x_1 = 1;
                  x_2 = x_centroid(1) + d1_2;
                  y_1 = y_centroid(1) - (7168 - d1_2);
                  y_2 = y_centroid(1) + d1_2; 
                  else
                  x_1 = 1;
                  x_2 = x_centroid(1) + d1_2;
                  y_1 = y_centroid(1) - d1_2;
                  y_2 = 13824;    
                  end    
               end 
               
               if(r_1_4_d12 == 1)
                  if(y_centroid(1) < y_centroid(2))
                  x_1 = x_centroid(1) - (7168 - d1_2);
                  x_2 = x_centroid(1) + d1_2;
                  y_1 = 1;
                  y_2 = y_centroid(1) + d1_2;
                  else
                  x_1 = x_centroid(1) - (7168 - d1_2);
                  x_2 = x_centroid(1) + d1_2;
                  y_1 = y_centroid(1) - d1_2
                  y_2 = y_centroid(1) + (7168 - d1_2);    
                  end    
               end 
               
               if(r_1_5_d12 == 1)
                  if(y_centroid(1) < y_centroid(2))
                  x_1 = x_centroid(1) - (7168 - d1_2);
                  x_2 = x_centroid(1) + d1_2;
                  y_1 = y_centroid(1) + d1_2;
                  y_2 = 13824;
                  else
                  x_1 = x_centroid(1) - (7168 - d1_2);
                  x_2 = x_centroid(1) + d1_2;
                  y_1 = y_centroid(1) + (7168 - d1_2)
                  y_2 = y_centroid(1) + d1_2;    
                  end    
               end 
               
               % Sphere 1 and 3 are at the same quadrant
               
               if(l_1_1_d13 == 1)
                  if(y_centroid(1) < y_centroid(3))
                  x_1 = 1;
                  x_2 = x_centroid(1) + d1_3;
                  y_1 = 1;
                  y_2 = y_centroid(1) + d1_3; 
                  else
                  x_1 = 1;
                  x_2 = x_centroid(1) + d1_3;
                  y_1 = y_centroid(1) - d1_3;
                  y_2 = y_centroid(1) + (7168 - d1_3);    
                  end    
               end   
               
               if(l_1_2_d13 == 1)
                  if(y_centroid(1) < y_centroid(3))
                  x_1 = 1;
                  x_2 = x_centroid(1) + d1_3;
                  y_1 = y_centroid(1) - (7168 - d1_3);
                  y_2 = y_centroid(1) + d1_3; 
                  else
                  x_1 = 1;
                  x_2 = x_centroid(1) + d1_3;
                  y_1 = y_centroid(1) - d1_3;
                  y_2 = 13824;    
                  end    
               end 
               
               if(r_1_4_d13 == 1)
                  if(y_centroid(1) < y_centroid(3))
                  x_1 = x_centroid(1) - (7168 - d1_3);
                  x_2 = x_centroid(1) + d1_3;
                  y_1 = 1;
                  y_2 = y_centroid(1) + d1_3;
                  else
                  x_1 = x_centroid(1) - (7168 - d1_3);
                  x_2 = x_centroid(1) + d1_3;
                  y_1 = y_centroid(1) - d1_3
                  y_2 = y_centroid(1) + (7168 - d1_3);    
                  end    
               end 
               
               if(r_1_5_d13 == 1)
                  if(y_centroid(1) < y_centroid(3))
                  x_1 = x_centroid(1) - (7168 - d1_3);
                  x_2 = x_centroid(1) + d1_3;
                  y_1 = y_centroid(1) + d1_3;
                  y_2 = 13824;
                  else
                  x_1 = x_centroid(1) - (7168 - d1_3);
                  x_2 = x_centroid(1) + d1_3;
                  y_1 = y_centroid(1) + (7168 - d1_3)
                  y_2 = y_centroid(1) + d1_3;    
                  end    
               end 
               
               
               % Sphere 1 and 2 and 3 are in the same location
               
               if(l_1_1_t == 1)
                  if(y_centroid(1) < y_centroid(2) && y_centroid(1) < y_centroid(3))
                  minxy = min(Vd1(1:2))    
                  x_1 = 1;
                  x_2 = x_centroid(1) + minxy;
                  y_1 = 1;
                  y_2 = y_centroid(1) + minxy; 
                  end
                  if(y_centroid(1) > y_centroid(2) && y_centroid(1) > y_centroid(3))
                  minxy = min(Vd1(1:2))   
                  x_1 = 1;
                  x_2 = x_centroid(1) + minxy;
                  y_1 = y_centroid(1) - minxy;
                  y_2 = y_centroid(1) + (7168-minxy); 
                  end
                  if(y_centroid(1) < y_centroid(2) && y_centroid(1) > y_centroid(3))
                  minxy = min(Vd1(1:2))    
                  x_1 = 1;
                  x_2 = x_centroid(1) + minxy;
                  y_1 = y_centroid(1) - d1_3;
                  y_2 = y_centroid(1) + d1_2; 
                  end
                  if(y_centroid(1) > y_centroid(2) && y_centroid(1) < y_centroid(3))
                  minxy = min(Vd1(1:2))    
                  x_1 = 1;
                  x_2 = x_centroid(1) + minxy;
                  y_1 = y_centroid(1) - d1_2;
                  y_2 = y_centroid(1) + d1_3; 
                  end
               end 
               
               if(l_1_2_t == 1)
                  if(y_centroid(1) < y_centroid(2) && y_centroid(1) < y_centroid(3))
                  minxy = min(Vd1(1:2))    
                  x_1 = 1;
                  x_2 = x_centroid(1) + minxy;
                  y_1 = y_centroid(1) - (7168-minxy);
                  y_2 = y_centroid(1) + minxy; 
                  end
                  if(y_centroid(1) > y_centroid(2) && y_centroid(1) > y_centroid(3))
                  minxy = min(Vd1(1:2))    
                  x_1 = 1;
                  x_2 = x_centroid(1) + minxy;
                  y_1 = y_centroid(1) - minxy;
                  y_2 = 13824; 
                  end
                  if(y_centroid(1) < y_centroid(2) && y_centroid(1) > y_centroid(3))
                  minxy = min(Vd1(1:2))    
                  x_1 = 1;
                  x_2 = x_centroid(1) + minxy;
                  y_1 = y_centroid(1) - d1_3;
                  y_2 = y_centroid(1) + d1_2; 
                  end
                  if(y_centroid(1) > y_centroid(2) && y_centroid(1) < y_centroid(3))
                  minxy = min(Vd1(1:2))    
                  x_1 = 1;
                  x_2 = x_centroid(1) + minxy;
                  y_1 = y_centroid(1) - d1_2;
                  y_2 = y_centroid(1) + d1_3; 
                  end
               end
               
                if(r_1_4_t == 1)
                  if(y_centroid(1) < y_centroid(2) && y_centroid(1) < y_centroid(3))
                  minxy = min(Vd1(1:2))    
                  x_1 = x_centroid(1) - (7168 - minxy);
                  x_2 = x_centroid(1) + minxy;
                  y_1 = 1;
                  y_2 = y_centroid(1) + minxy; 
                  end
                  if(y_centroid(1) > y_centroid(2) && y_centroid(1) > y_centroid(3))
                  minxy = min(Vd1(1:2))    
                  x_1 = x_centroid(1) - (7168 - minxy);
                  x_2 = x_centroid(1) + minxy;
                  y_1 = y_centroid(1) - minxy;
                  y_2 = y_centroid(1) + (7168 - minxy); 
                  end
                  if(y_centroid(1) < y_centroid(2) && y_centroid(1) > y_centroid(3))
                  minxy = min(Vd1(1:2))   
                  x_1 = x_centroid(1) - (7168 - minxy);
                  x_2 = x_centroid(1) + minxy;
                  y_1 = y_centroid(1) - d1_3;
                  y_2 = y_centroid(1) + d1_2; 
                  end
                  if(y_centroid(1) > y_centroid(2) && y_centroid(1) < y_centroid(3))
                  minxy = min(Vd1(1:2))    
                  x_1 = x_centroid(1) - (7168 - minxy);
                  x_2 = x_centroid(1) + minxy;
                  y_1 = y_centroid(1) - d1_2;
                  y_2 = y_centroid(1) + d1_3; 
                  end
                end
                
                if(r_1_5_t == 1)
                  if(y_centroid(1) < y_centroid(2) && y_centroid(1) < y_centroid(3))
                  minxy = min(Vd1(1:2))    
                  x_1 = x_centroid(1) - (7168 - minxy);
                  x_2 = x_centroid(1) + minxy;
                  y_1 = y_centroid(1) - (7168-minxy);
                  y_2 = y_centroid(1) + minxy;
                  end
                  if(y_centroid(1) > y_centroid(2) && y_centroid(1) > y_centroid(3))
                  minxy = min(Vd1(1:2))   
                  x_1 = x_centroid(1) - (7168 - minxy);
                  x_2 = x_centroid(1) + minxy;
                  y_1 = y_centroid(1) - minxy;
                  y_2 = 13824;  
                  end
                  if(y_centroid(1) < y_centroid(2) && y_centroid(1) > y_centroid(3))
                  minxy = min(Vd1(1:2))    
                  x_1 = x_centroid(1) - (7168 - minxy);
                  x_2 = x_centroid(1) + minxy;
                  y_1 = y_centroid(1) - d1_3;
                  y_2 = y_centroid(1) + d1_2; 
                  end
                  if(y_centroid(1) > y_centroid(2) && y_centroid(1) < y_centroid(3))
                  minxy = min(Vd1(1:2))    
                  x_1 = x_centroid(1) - (7168 - minxy);
                  x_2 = x_centroid(1) + minxy;
                  y_1 = y_centroid(1) - d1_2;
                  y_2 = y_centroid(1) + d1_3; 
                  end
                end
                
                if(x_1 < 0)
                    x_1 =0;
                end    
                if(x_2 > 10752)
                    x_2 =10752;
                end
                if(y_1 < 0)
                    y_1 =0;
                end
                if(y_2 > 13824)
                    y_2 =13824;
                end
                
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
                
                
               imgnucleus = imread([foldername '/' allfiles(i).name]);
               %imwrite(imgnucleus,[foldername2 '/' allfiles(i).name]);
                %The rectangular image is cutted from the imgnucleus image containing the whole cover slip
               imagesphere1 = imgnucleus(y_1:y_2,x_1:x_2);
               % Dimensions of imagesphere1 are determined
               [r1 c1] = size(imagesphere1);
               % An empty matrix is generated
               imagesphereresize1 = zeros(7168,7168);
               % imagesphereresize1 is converted to 8-bit
               imagesphereresize1 = uint8(imagesphereresize1);
               %The empty indices in the imagesphereresize1 matrix are overwritten with imagesphere 1 entries
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
               % imagesphereresize1 is saved to final destination
               newFile1 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize1,[foldername1 '/' newFile1]);
               wellname = (newFile1(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(7168,7168);
               currentNM = csvHandler.CellPosMatrix(wellname);
               currentNM(yStart1:yEnd1,xStart1:xEnd1) = NucleusM(y_1:y_2,x_1:x_2);
               csvHandler.CellPosMatrix(wellname) = currentNM;
               imagesphereresize1_small = imresize(imagesphereresize1, optionHandler.ScalingFactor);
               newFile1s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize1_small,[foldername1 '/' newFile1s]);
               clear imagesphereresize1_small
               clear imagesphereresize1; 
               clear imagesphere1;
               
               % Sphere 2 with second smallest value of x_centroid
               
               
               
               
               if(l_2_1 >0 && l_2_1 < 2);
                       x_3 = 1;
                       x_4 = x_centroid(2) + d_2_c;
                       y_3 = 1;
                       y_4 = y_centroid(2) + d_2_c;
               else
               end    
               
               if(l_2_2 >0 && l_2_2 < 2);
                       x_3 = 1;
                       x_4 = x_centroid(2) + d_2_c;
                       y_3 = y_centroid(2) - d_2_c;
                       y_4 = 13824;
               else
               end
               
               if(r_2_4 >0 && r_2_4 < 2);
                       x_3 = x_centroid(2) - d_2_c;
                       x_4 = 10752;
                       y_3 = 1;
                       y_4 = y_centroid(2) + d_2_c;
               else
               end
               
               if(r_2_5 >0 && r_2_5 < 2);
                       x_3 = x_centroid(2) - d_2_c;
                       x_4 = 10752;
                       y_3 = y_centroid(2) - d_2_c;
                       y_4 = 13824;
               else
               end
               
               if(x_centroid(2) ~= x_centroid_c);
               else
                   if(x_centroid(1) < x_centroid(2) && x_centroid(3) < x_centroid(2));
                       V_1m = [d1_2 ; d2_3 ; (7168/2)];
                       min1m = min(V_1m(1:3));
                       min1p = (7168-min1m);
                   else
                   end
                   if(x_centroid(1) > x_centroid(2) && x_centroid(3) > x_centroid(2));
                       V_1p = [d1_2 ; d2_3 ; (7168/2)];
                       min1p = min(V_1p(1:3));
                       min1m = (7168-min1p);
                   else
                   end
                   if(x_centroid(1) > x_centroid(2) && x_centroid(3) < x_centroid(2));
                       V_1m = [d2_3 ; (7168/2)];
                       V_1p = [d1_2 ; (7168/2)];
                       min1p = min(V_1p(1:2));
                       min1m = min(V_1m(1:2));
                   else
                       
                   end
                   if(x_centroid(1) < x_centroid(2) && x_centroid(3) > x_centroid(2));
                       V_1m = [d1_2 ; (7168/2)];
                       V_1p = [d2_3 ; (7168/2)];
                       min1p = min(V_1p(1:2));
                       min1m = min(V_1m(1:2));
                        
                   else
                       
                   end
                   
                   x_3 = x_centroid(2) - min1m 
                   x_4 = x_centroid(2) + min1p
                   y_3 = y_centroid(2) - min1m
                   y_4 = y_centroid(2) + min1p
               end  
               
               % Cases for having two or three spheres at one location
               
               % Sphere 1 and 2 are in the same location
               
               if(l_2_1_d12 == 1)
                  if(y_centroid(2) < y_centroid(1))
                  x_3 = x_centroid(2) - d1_2;
                  x_4 = x_centroid(2) + (7168 - d1_2);
                  y_3 = 1;
                  y_4 = y_centroid(2) + d1_2; 
                  else
                  x_3 = x_centroid(2) - d1_2;
                  x_4 = x_centroid(2) + (7168 - d1_2);
                  y_3 = y_centroid(2) - d1_2;
                  y_4 = y_centroid(2) + (7168 - d1_2);    
                  end    
               end   
               
               if(l_2_2_d12 == 1)
                  if(y_centroid(2) < y_centroid(1))
                  x_3 = x_centroid(2) - d1_2;
                  x_4 = x_centroid(2) + (7168 - d1_2);
                  y_3 = y_centroid(2) - (7168 - d1_2);
                  y_4 = y_centroid(2) + d1_2; 
                  else
                  x_3 = x_centroid(2) - d1_2;
                  x_4 = x_centroid(2) + (7168 - d1_2);
                  y_3 = y_centroid(2) - d1_2;
                  y_4 = 13824;    
                  end    
               end 
               
               if(r_2_4_d12 == 1)
                  if(y_centroid(2) < y_centroid(1))
                  x_3 = x_centroid(2) - d1_2;
                  x_4 = 10752;
                  y_3 = 1;
                  y_4 = y_centroid(2) + d1_2;
                  else
                  x_3 = x_centroid(2) - d1_2;
                  x_4 = 10752;
                  y_3 = y_centroid(2) - d1_2
                  y_4 = y_centroid(2) + (7168 - d1_2);    
                  end    
               end 
               
               if(r_2_5_d12 == 1)
                  if(y_centroid(2) < y_centroid(1))
                  x_3 = x_centroid(2) - d1_2;
                  x_4 = 10752;
                  y_3 = y_centroid(2) + d1_2;
                  y_4 = 13824;
                  else
                  x_3 = x_centroid(2) - d1_2;
                  x_4 = 10752;
                  y_3 = y_centroid(2) + (7168 - d1_2)
                  y_4 = y_centroid(2) + d1_2;    
                  end    
               end 
               
               % Sphere 2 and 3 are in the same location
               
               if(l_2_1_d23 == 1)
                  if(y_centroid(2) < y_centroid(3))
                  x_3 = 1;
                  x_4 = x_centroid(2) + d2_3;
                  y_3 = 1;
                  y_4 = y_centroid(2) + d2_3; 
                  else
                  x_3 = 1;
                  x_4 = x_centroid(2) + d2_3;
                  y_3 = y_centroid(2) - d2_3;
                  y_4 = y_centroid(2) + (7168 - d2_3);    
                  end    
               end   
               
               if(l_2_2_d23 == 1)
                  if(y_centroid(2) < y_centroid(3))
                  x_3 = 1;
                  x_4 = x_centroid(2) + d2_3;
                  y_3 = y_centroid(2) - (7168 - d2_3);
                  y_4 = y_centroid(2) + d2_3; 
                  else
                  x_3 = 1;
                  x_4 = x_centroid(2) + d2_3;
                  y_3 = y_centroid(2) - d2_3;
                  y_4 = 13824;    
                  end    
               end 
               
               if(r_2_4_d23 == 1)
                  if(y_centroid(2) < y_centroid(3))
                  x_3 = x_centroid(2) - (7168 - d2_3);
                  x_4 = x_centroid(2) + d2_3;
                  y_3 = 1;
                  y_4 = y_centroid(2) + d2_3;
                  else
                  x_3 = x_centroid(2) - (7168 - d2_3);
                  x_4 = x_centroid(2) + d2_3;
                  y_3 = y_centroid(2) - d2_3
                  y_4 = y_centroid(2) + (7168 - d2_3);    
                  end    
               end 
               
               if(r_2_5_d23 == 1)
                  if(y_centroid(2) < y_centroid(3))
                  x_3 = x_centroid(2) - (7168 - d2_3);
                  x_4 = x_centroid(2) + d2_3;
                  y_3 = y_centroid(2) - (7168 - d2_3);
                  y_4 = y_centroid(2) + d2_3; 
                  else
                  x_3 = x_centroid(2) - (7168 - d2_3);
                  x_4 = x_centroid(2) + d2_3;
                  y_3 = y_centroid(2) - d2_3;
                  y_4 = 13824;     
                  end    
               end 
               
               % Sphere 1 and 2 and 3 are in the same location
               
               if(l_2_1_t == 1)
                  if(y_centroid(2) < y_centroid(1) && y_centroid(2) < y_centroid(3))
                  minxy = min(Vd2(1:2))    
                  x_3 = x_centroid(2) -d1_2;
                  x_4 = x_centroid(2) +d2_3;
                  y_3 = 1;
                  y_4 = y_centroid(2) + minxy; 
                  end
                  if(y_centroid(2) > y_centroid(1) && y_centroid(2) > y_centroid(3))
                  minxy = min(Vd2(1:2))     
                  x_3 = x_centroid(2) -d1_2;
                  x_4 = x_centroid(2) +d2_3;
                  y_3 = y_centroid(2) - minxy;
                  y_4 = y_centroid(2) + (7168 - minxy); 
                  end
                  if(y_centroid(2) < y_centroid(1) && y_centroid(2) > y_centroid(3))
                  minxy = min(Vd2(1:2))     
                  x_3 = x_centroid(2) -d1_2;
                  x_4 = x_centroid(2) +d2_3;
                  y_3 = y_centroid(2) - d2_3;
                  y_4 = y_centroid(2) + d1_2; 
                  end
                  if(y_centroid(2) > y_centroid(1) && y_centroid(2) < y_centroid(3))
                  minxy = min(Vd2(1:2))   
                  x_3 = x_centroid(2) -d1_2;
                  x_4 = x_centroid(2) +d2_3;
                  y_3 = y_centroid(2) - d1_2;
                  y_4 = y_centroid(2) + d2_3; 
                  end
               end 
               
               if(l_2_2_t == 1)
                  if(y_centroid(2) < y_centroid(1) && y_centroid(2) < y_centroid(3))
                  minxy = min(Vd2(1:2))     
                  x_3 = x_centroid(2) -d1_2;
                  x_4 = x_centroid(2) +d2_3;
                  y_3 = y_centroid(2) - (7168-minxy);
                  y_4 = y_centroid(2) + minxy; 
                  end
                  if(y_centroid(2) > y_centroid(1) && y_centroid(2) > y_centroid(3))
                  minxy = min(Vd2(1:2))     
                  x_3 = x_centroid(2) -d1_2;
                  x_4 = x_centroid(2) +d2_3;
                  y_3 = y_centroid(2) - minxy;
                  y_4 = 13824; 
                  end
                  if(y_centroid(2) < y_centroid(1) && y_centroid(2) > y_centroid(3))
                  minxy = min(Vd2(1:2))     
                  x_3 = x_centroid(2) -d1_2;
                  x_4 = x_centroid(2) +d2_3;
                  y_3 = y_centroid(2) - d2_3;
                  y_4 = y_centroid(2) + d1_2; 
                  end
                  if(y_centroid(2) > y_centroid(1) && y_centroid(2) < y_centroid(3))
                  minxy = min(Vd2(1:2))     
                  x_3 = x_centroid(2) -d1_2;
                  x_4 = x_centroid(2) +d2_3;
                  y_3 = y_centroid(2) - d1_2;
                  y_4 = y_centroid(2) + d2_3; 
                  end
               end
               
                if(r_2_4_t == 1)
                  if(y_centroid(2) < y_centroid(1) && y_centroid(2) < y_centroid(3))
                  minxy = min(Vd2(1:2))    
                  x_3 = x_centroid(2) -d1_2;
                  x_4 = x_centroid(2) +d2_3;
                  y_3 = 1;
                  y_4 = y_centroid(2) + minxy; 
                  end
                  if(y_centroid(2) > y_centroid(1) && y_centroid(2) > y_centroid(3))
                  minxy = min(Vd2(1:2))     
                  x_3 = x_centroid(2) -d1_2;
                  x_4 = x_centroid(2) +d2_3;
                  y_3 = y_centroid(2) - minxy;
                  y_4 = y_centroid(2) + (7168 - minxy); 
                  end
                  if(y_centroid(2) < y_centroid(1) && y_centroid(2) > y_centroid(3))
                  minxy = min(Vd2(1:2))     
                  x_3 = x_centroid(2) -d1_2;
                  x_4 = x_centroid(2) +d2_3;
                  y_3 = y_centroid(2) - d2_3;
                  y_4 = y_centroid(2) + d1_2; 
                  end
                  if(y_centroid(2) > y_centroid(1) && y_centroid(2) < y_centroid(3))
                  minxy = min(Vd2(1:2))     
                  x_3 = x_centroid(2) -d1_2;
                  x_4 = x_centroid(2) +d2_3;
                  y_3 = y_centroid(2) - d1_2;
                  y_4 = y_centroid(2) + d2_3; 
                  end
                end
                
                if(r_2_5_t == 1)
                  if(y_centroid(2) < y_centroid(1) && y_centroid(2) < y_centroid(3))
                  minxy = min(Vd2(1:2))    
                  x_3 = x_centroid(2) -d1_2;
                  x_4 = x_centroid(2) +d2_3;
                  y_3 = y_centroid(2) - (7168-minxy);
                  y_4 = y_centroid(2) + minxy;
                  end
                  if(y_centroid(2) > y_centroid(1) && y_centroid(2) > y_centroid(3))
                  minxy = min(Vd2(1:2))     
                  x_3 = x_centroid(2) -d1_2;
                  x_4 = x_centroid(2) +d2_3;
                  y_3 = y_centroid(2) - minxy;
                  y_4 = 13824;  
                  end
                  if(y_centroid(2) < y_centroid(1) && y_centroid(2) > y_centroid(3))
                  minxy = min(Vd2(1:2))     
                  x_3 = x_centroid(2) -d1_2;
                  x_4 = x_centroid(2) +d2_3;
                  y_3 = y_centroid(2) - d2_3;
                  y_4 = y_centroid(2) + d1_2; 
                  end
                  if(y_centroid(2) > y_centroid(1) && y_centroid(2) < y_centroid(3))
                  minxy = min(Vd2(1:2))     
                  x_3 = x_centroid(2) -d1_2;
                  x_4 = x_centroid(2) +d2_3;
                  y_3 = y_centroid(2) - d1_2;
                  y_4 = y_centroid(2) + d2_3; 
                  end
                end
               
                if(x_3 < 0)
                    x_3 =0;
                end    
                if(x_4 > 10752)
                    x_4 =10752;
                end
                if(y_3 < 0)
                    y_3 =0;
                end
                if(y_4 > 13824)
                    y_4 =13824;
                end
                
               diffXCentroidl = x_centroid(2) - x_3;
               diffXCentroidr = x_4 - x_centroid(2);
               diffYCentroidl = y_centroid(2) - y_3;
               diffYCentroidr = y_4 - y_centroid(2);
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
               %The rectangular image is cutted from the imgnucleus image containing the whole cover slip
               imagesphere2 = imgnucleus(y_3:y_4,x_3:x_4);
               % Dimensions of imagesphere1 are determined
               [r2 c2] = size(imagesphere2);
               % An empty matrix is generated
               imagesphereresize2 = zeros(7168,7168);
               % imagesphereresize1 is converted to 8-bit
               imagesphereresize2 = uint8(imagesphereresize2);
               %The empty indices in the imagesphereresize1 matrix are overwritten with imagesphere 1 entries
               if(leftBoarder == 0)
                   if(lowerBoarder == 0)
                        
                        yStart2 = 1;
                        yEnd2 = r2;
                        xStart2 = 1;
                        xEnd2 = c2;
                   else
                        
                        yStart2 = (7169-r2);
                        yEnd2 = 7168;
                        xStart2 = 1;
                        xEnd2 = c2;
                   end 
               else
                   if(lowerBoarder == 0)
                       
                        yStart2 = 1;
                        yEnd2 = r2;
                        xStart2 = (7169 - c2);
                        xEnd2 = 7168;
                   else
                       
                        yStart2 = (7169-r2);
                        yEnd2 = 7168;
                        xStart2 = (7169 - c2);
                        xEnd2 = 7168;
                   end
               end  
               
               
               imagesphereresize2(yStart2:yEnd2,xStart2:xEnd2)=imagesphere2(1:r2,1:c2);
               % imagesphereresize1 is saved to final destination
               newFile2 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize2,[foldername1 '/' newFile2]);
               n = n + 1;
               wellname = (newFile2(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(7168,7168);
               currentNM = csvHandler.CellPosMatrix(wellname);
               currentNM(yStart2:yEnd2,xStart2:xEnd2) = NucleusM(y_3:y_4,x_3:x_4);
               csvHandler.CellPosMatrix(wellname) = currentNM;
               imagesphereresize2_small = imresize(imagesphereresize2, optionHandler.ScalingFactor);
               newFile2s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize2_small,[foldername1 '/' newFile2s]);
               clear imagesphereresize2_small
               clear imagesphereresize2; 
               clear imagesphere2;
               
               
            % Sphere 3 with third smallest value of x_centroid
               
               
                             
               if(l_3_1 >0 && l_3_1 < 2);
                       x_5 = 1;
                       x_6 = x_centroid(3) + d_3_c;
                       y_5 = 1;
                       y_6 = y_centroid(3) + d_3_c;
               else
               end    
               
               if(l_3_2 >0 && l_3_2 < 2);
                       x_5 = 1;
                       x_6 = x_centroid(3) + d_3_c;
                       y_5 = y_centroid(3) - d_3_c;
                       y_6 = 13824;
               else
               end
               
               if(r_3_4 >0 && r_3_4 < 2);
                       x_5 = x_centroid(3) - d_3_c;
                       x_6 = 10752;
                       y_5 = 1;
                       y_6 = y_centroid(3) + d_3_c;
               else
               end
               
               if(r_3_5 >0 && r_3_5 < 2);
                       x_5 = x_centroid(3) - d_3_c;
                       x_6 = 10752;
                       y_5 = y_centroid(3) - d_3_c;
                       y_6 = 13824;
               else
               end
               
               if(x_centroid(3) ~= x_centroid_c);
               else
                   if(x_centroid(1) < x_centroid(3) && x_centroid(1) < x_centroid(3));
                       V_1m = [d1_2 ; d1_3 ; (7168/2)];
                       min1m = min(V_1m(1:3));
                       min1p = (7168-min1m);
                   else
                   end
                   if(x_centroid(1) > x_centroid(3) && x_centroid(2) > x_centroid(3));
                       V_1p = [d1_2 ; d1_3 ; (7168/2)];
                       min1p = min(V_1p(1:3));
                       min1m = (7168-min1p);
                   else
                   end
                   if(x_centroid(1) > x_centroid(3) && x_centroid(2) < x_centroid(3));
                       V_1m = [d1_3 ; (7168/2)];
                       V_1p = [d1_2 ; (7168/2)];
                       min1p = min(V_1p(1:2));
                       min1m = min(V_1m(1:2));
                   else
                       
                   end
                   if(x_centroid(1) < x_centroid(3) && x_centroid(2) > x_centroid(3));
                       V_1m = [d1_2 ; (7168/2)];
                       V_1p = [d1_3 ; (7168/2)];
                       min1p = min(V_1p(1:2));
                       min1m = min(V_1m(1:2));
                        
                   else
                       
                   end
                   
                   x_5 = x_centroid(3) - min1m 
                   x_6 = x_centroid(3) + min1p
                   y_5 = y_centroid(3) - min1m
                   y_6 = y_centroid(3) + min1p
               end  
               
               % Cases for having two or three spheres at one location
               
                % Sphere 1 and 3 are in the same location
               
               if(l_3_1_d13 == 1)
                  if(y_centroid(3) < y_centroid(1))
                  x_5 = x_centroid(3) - d1_3;
                  x_6 = x_centroid(3) + (7168 - d1_3);
                  y_5 = 1;
                  y_6 = y_centroid(3) + d1_3; 
                  else
                  x_5 = x_centroid(3) - d1_3;
                  x_6 = x_centroid(3) + (7168 - d1_3);
                  y_5 = y_centroid(3) - d1_3;
                  y_6 = y_centroid(3) + (7168 - d1_3);    
                  end    
               end   
               
               if(l_3_2_d13 == 1)
                  if(y_centroid(3) < y_centroid(1))
                  x_5 = x_centroid(3) - d1_3;
                  x_6 = x_centroid(3) + (7168 - d1_3);
                  y_5 = y_centroid(3) - (7168 - d1_3);
                  y_6 = y_centroid(3) + d1_3; 
                  else
                  x_5 = x_centroid(3) - d1_3;
                  x_6 = x_centroid(3) + (7168 - d1_3);
                  y_5 = y_centroid(3) - d1_3;
                  y_6 = 13824;    
                  end    
               end 
               
               if(r_3_4_d13 == 1)
                  if(y_centroid(3) < y_centroid(1))
                  x_5 = x_centroid(3)  - d1_3;
                  x_6 = 10752;
                  y_5 = 1;
                  y_6 = y_centroid(3) + d1_3;
                  else
                  x_5 = x_centroid(3)  - d1_3;
                  x_6 = 10752;
                  y_5 = y_centroid(3) - d1_3
                  y_6 = y_centroid(3) + (7168 - d1_3);    
                  end    
               end 
               
               % Sphere 2 and 3 are in the same location
               
               if(l_3_1_d23 == 1)
                  if(y_centroid(3) < y_centroid(2))
                  x_5 = x_centroid(3) - d2_3;
                  x_6 = x_centroid(3) + (7168 - d2_3);
                  y_5 = 1;
                  y_6 = y_centroid(3) + d2_3; 
                  else
                  x_5 = x_centroid(3) - d2_3;
                  x_6 = x_centroid(3) + (7168 - d2_3);
                  y_5 = y_centroid(3) - d1_3;
                  y_6 = y_centroid(3) + (7168 - d2_3);    
                  end    
               end   
               
               if(l_3_2_d23 == 1)
                  if(y_centroid(3) < y_centroid(2))
                  x_5 = x_centroid(3) - d2_3;
                  x_6 = x_centroid(3) + (7168 - d2_3);
                  y_5 = y_centroid(3) - (7168 - d2_3);
                  y_6 = y_centroid(3) + d2_3; 
                  else
                  x_5 = x_centroid(3) - d2_3;
                  x_6 = x_centroid(3) + (7168 - d2_3);
                  y_5 = y_centroid(3) - d2_3;
                  y_6 = 13824;    
                  end    
               end 
               
               if(r_3_4_d23 == 1)
                  if(y_centroid(3) < y_centroid(2))
                  x_5 = x_centroid(3)  - d2_3;
                  x_6 = 10752;
                  y_5 = 1;
                  y_6 = y_centroid(3) + d2_3;
                  else
                  x_5 = x_centroid(3)  - d2_3;
                  x_6 = 10752;
                  y_5 = y_centroid(3) - d2_3
                  y_6 = y_centroid(3) + (7168 - d2_3);    
                  end    
               end 
               
               if(r_3_5_d23 == 1)
                  if(y_centroid(3) < y_centroid(2))
                  x_5 = x_centroid(3)  - d2_3;
                  x_6 = 10752;
                  y_5 = y_centroid(3) - (7168 - d2_3);
                  y_6 = y_centroid(3) + d2_3; 
                  else
                  x_5 = x_centroid(3)  - d2_3;
                  x_6 = 10752;
                  y_5 = y_centroid(3) - d2_3;
                  y_6 = 13824;     
                  end    
               end 
               
               % Sphere 1 and 2 and 3 are in the same location
               
               if(l_3_1_t == 1)
                  if(y_centroid(3) < y_centroid(1) && y_centroid(3) < y_centroid(2))
                  minxy = min(Vd3(1:2))     
                  x_5 = x_centroid(3) - d2_3;
                  x_6 = x_centroid(3) + (7168 -d2_3);
                  y_5 = 1;
                  y_6 = y_centroid(3) + minxy; 
                  end
                  if(y_centroid(3) > y_centroid(1) && y_centroid(3) > y_centroid(2))
                  minxy = min(Vd3(1:2))    
                  x_5 = x_centroid(3) - d2_3;
                  x_6 = x_centroid(3) + (7168 -d2_3);
                  y_5 = y_centroid(3) - minxy;
                  y_6 = y_centroid(3) + (7168 - minxy); 
                  end
                  if(y_centroid(3) < y_centroid(1) && y_centroid(3) > y_centroid(2))
                  minxy = min(Vd3(1:2))    
                  x_5 = x_centroid(3) - d2_3;
                  x_6 = x_centroid(3) + (7168 -d2_3);
                  y_5 = y_centroid(3) - d1_3;
                  y_6 = y_centroid(3) + d1_2; 
                  end
                  if(y_centroid(3) > y_centroid(1) && y_centroid(3) < y_centroid(2))
                  minxy = min(Vd3(1:2))    
                  x_5 = x_centroid(3) - d2_3;
                  x_6 = x_centroid(3) + (7168 -d2_3);
                  y_5 = y_centroid(3) - d1_2;
                  y_6 = y_centroid(3) + d1_3; 
                  end
               end 
               
               if(l_3_2_t == 1)
                  if(y_centroid(3) < y_centroid(1) && y_centroid(3) < y_centroid(2))
                  minxy = min(Vd3(1:2))    
                  x_5 = x_centroid(3) - d2_3;
                  x_6 = x_centroid(3) + (7168 -d2_3);
                  y_5 = y_centroid(3) - (7168-minxy);
                  y_6 = y_centroid(3) + minxy; 
                  end
                  if(y_centroid(3) > y_centroid(1) && y_centroid(3) > y_centroid(2))
                  minxy = min(Vd3(1:2))    
                  x_5 = x_centroid(3) - d2_3;
                  x_6 = x_centroid(3) + (7168 -d2_3);
                  y_5 = y_centroid(3) - minxy;
                  y_6 = 13824; 
                  end
                  if(y_centroid(3) < y_centroid(1) && y_centroid(3) > y_centroid(2))
                  minxy = min(Vd3(1:2))    
                  x_5 = x_centroid(3) - d2_3;
                  x_6 = x_centroid(3) + (7168 -d2_3);
                  y_5 = y_centroid(3) - d1_3;
                  y_6 = y_centroid(3) + d1_2; 
                  end
                  if(y_centroid(3) > y_centroid(1) && y_centroid(3) < y_centroid(2))
                  minxy = min(Vd3(1:2))    
                  x_5 = x_centroid(3) - d2_3;
                  x_6 = x_centroid(3) + (7168 -d2_3);
                  y_5 = y_centroid(3) - d1_2;
                  y_6 = y_centroid(3) + d1_3; 
                  end
               end
               
                if(r_3_4_t == 1)
                  if(y_centroid(3) < y_centroid(1) && y_centroid(3) < y_centroid(2))
                  minxy = min(Vd3(1:2))    
                  x_5 = x_centroid(3) -d2_3;
                  x_6 = 10752;
                  y_5 = 1;
                  y_6 = y_centroid(3) + minxy; 
                  end
                  if(y_centroid(3) > y_centroid(1) && y_centroid(3) > y_centroid(2))
                  minxy = min(Vd3(1:2))    
                  x_5 = x_centroid(3) -d2_3;
                  x_6 = 10752;
                  y_5 = y_centroid(3) - minxy;
                  y_6 = y_centroid(3) + (7168 - minxy); 
                  end
                  if(y_centroid(3) < y_centroid(1) && y_centroid(3) > y_centroid(2))
                  minxy = min(Vd3(1:2))    
                  x_5 = x_centroid(3) -d2_3;
                  x_6 = 10752;
                  y_5 = y_centroid(3) - d1_3;
                  y_6 = y_centroid(3) + d1_2; 
                  end
                  if(y_centroid(3) > y_centroid(1) && y_centroid(3) < y_centroid(2))
                  minxy = min(Vd3(1:2))    
                  x_5 = x_centroid(3) -d2_3;
                  x_6 = 10752;
                  y_5 = y_centroid(3) - d1_2;
                  y_6 = y_centroid(3) + d1_3; 
                  end
                end
                
                if(r_3_5_t == 1)
                  if(y_centroid(3) < y_centroid(1) && y_centroid(3) < y_centroid(2))
                  minxy = min(Vd3(1:2))    
                  x_5 = x_centroid(3) -d2_3;
                  x_6 = 10752;
                  y_5 = y_centroid(3) - (7168-minxy);
                  y_6 = y_centroid(3) + minxy;
                  end
                  if(y_centroid(3) > y_centroid(1) && y_centroid(3) > y_centroid(2))
                  minxy = min(Vd3(1:2))    
                  x_5 = x_centroid(3) -d2_3;
                  x_6 = 10752;
                  y_5 = y_centroid(3) - minxy;
                  y_6 = 13824;  
                  end
                  if(y_centroid(3) < y_centroid(1) && y_centroid(3) > y_centroid(2))
                  minxy = min(Vd3(1:2))    
                  x_5 = x_centroid(3) -d2_3;
                  x_6 = 10752;
                  y_5 = y_centroid(3) - d1_3;
                  y_6 = y_centroid(3) + d1_2; 
                  end
                  if(y_centroid(3) > y_centroid(1) && y_centroid(3) < y_centroid(2))
                  minxy = min(Vd3(1:2))    
                  x_5 = x_centroid(3) -d2_3;
                  x_6 = 10752;
                  y_5 = y_centroid(3) - d1_2;
                  y_6 = y_centroid(3) + d1_3; 
                  end
                end
                
               if(x_5 < 0)
                    x_5 =0;
                end    
                if(x_6 > 10752)
                    x_6 =10752;
                end
                if(y_5 < 0)
                    y_5 =0;
                end
                if(y_6 > 13824)
                    y_6 =13824;
                end
               diffXCentroidl = x_centroid(3) - x_5;
               diffXCentroidr = x_6 - x_centroid(3);
               diffYCentroidl = y_centroid(3) - y_5;
               diffYCentroidr = y_6 - y_centroid(3);
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
                %The rectangular image is cutted from the imgnucleus image containing the whole cover slip
               imagesphere3 = imgnucleus(y_5:y_6,x_5:x_6);
               % Dimensions of imagesphere1 are determined
               [r3 c3] = size(imagesphere3);
               % An empty matrix is generated
               imagesphereresize3 = zeros(7168,7168);
               % imagesphereresize1 is converted to 8-bit
               imagesphereresize3 = uint8(imagesphereresize3);
               %The empty indices in the imagesphereresize1 matrix are overwritten with imagesphere 1 entries
               if(leftBoarder == 0)
                   if(lowerBoarder == 0)
                        
                        yStart3 = 1;
                        yEnd3 = r3;
                        xStart3 = 1;
                        xEnd3 = c3;
                   else
                        
                        yStart3 = (7169-r3);
                        yEnd3 = 7168;
                        xStart3 = 1;
                        xEnd3 = c3;
                   end 
               else
                   if(lowerBoarder == 0)
                       
                        yStart3 = 1;
                        yEnd3 = r3;
                        xStart3 = (7169 - c3);
                        xEnd3 = 7168;
                   else
                       
                        yStart3 = (7169-r3);
                        yEnd3 = 7168;
                        xStart3 = (7169 - c3);
                        xEnd3 = 7168;
                   end
               end  
               
               
               imagesphereresize3(yStart3:yEnd3,xStart3:xEnd3)=imagesphere3(1:r3,1:c3);
               % imagesphereresize1 is saved to final destination
               newFile3 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize3,[foldername1 '/' newFile3]);
               n = n + 1;
               wellname = (newFile3(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(7168,7168);
               currentNM = csvHandler.CellPosMatrix(wellname);
               currentNM(yStart3:yEnd3,xStart3:xEnd3) = NucleusM(y_5:y_6,x_5:x_6);
               csvHandler.CellPosMatrix(wellname) = currentNM;
               imagesphereresize3_small = imresize(imagesphereresize3, optionHandler.ScalingFactor);
               newFile3s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize3_small,[foldername1 '/' newFile3s]);
               clear imagesphereresize3_small
               clear imagesphereresize3; 
               clear imagesphere3;
               
               imagesphereresize4 = zeros(7168,7168);
               newFile4 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize4,[foldername1 '/' newFile4]);
               n = n + 1;
               wellname = (newFile4(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(logical(zeros(7168,7168)));
               imagesphereresize4_small = imresize(imagesphereresize4, optionHandler.ScalingFactor);
               newFile4s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize4_small,[foldername1 '/' newFile4s]);
               clear imagesphereresize4_small
               clear imagesphere4;
               clear imagesphereresize4;
               
               imagesphereresize5 = zeros(7168,7168);
               newFile5 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize5,[foldername1 '/' newFile5]);
               n = n + 1;
               wellname = (newFile5(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) =  sparse(logical(zeros(7168,7168)));
               imagesphereresize5_small = imresize(imagesphereresize5, optionHandler.ScalingFactor);
               newFile5s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize5_small,[foldername1 '/' newFile5s]);
               clear imagesphereresize5_small
               clear imagesphere5;
               clear imagesphereresize5;
               clear imgnucleus;
               
               if(Neuron_Channel > 0)
                   % Same for neurite channel

                   imgneurite = imread([foldername '/' newFileNeu]);
                   %imwrite(imgneurite,[foldername2 '/' newFileNeu]);
                   imagesphere6 = imgneurite(y_1:y_2,x_1:x_2);
                   [r1 c1] = size(imagesphere6);
                   imagesphereresize6 = zeros(7168,7168);
                   imagesphereresize6 = uint8(imagesphereresize6);
                   imagesphereresize6(yStart1:yEnd1,xStart1:xEnd1)=imagesphere6(1:r1,1:c1);
                   newFile6 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize6,[foldername1 '/' newFile6]);
                   imagesphereresize6_small = imresize(imagesphereresize6, optionHandler.ScalingFactor);
                   newFile6s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize6_small,[foldername1 '/' newFile6s]);
                   clear imagesphere6;
                   clear imagesphereresize6;
                   clear imagesphereresize6_small;

                   imagesphere7 = imgneurite(y_3:y_4,x_3:x_4);
                   [r2 c2] = size(imagesphere7);
                   imagesphereresize7 = zeros(7168,7168);
                   imagesphereresize7 = uint8(imagesphereresize7);
                   imagesphereresize7(yStart2:yEnd2,xStart2:xEnd2)=imagesphere7(1:r2,1:c2);
                   newFile7 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize7,[foldername1 '/' newFile7]);
                   imagesphereresize7_small = imresize(imagesphereresize7, optionHandler.ScalingFactor);
                   newFile7s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize7_small,[foldername1 '/' newFile7s]);
                   clear imagesphere7;
                   clear imagesphereresize7;
                   clear imagesphereresize7_small;

                   imagesphere8 = imgneurite(y_5:y_6,x_5:x_6);
                   [r3 c3] = size(imagesphere8);
                   imagesphereresize8 = zeros(7168,7168);
                   imagesphereresize8 = uint8(imagesphereresize8);
                   imagesphereresize8(yStart3:yEnd3,xStart3:xEnd3)=imagesphere8(1:r3,1:c3);
                   newFile8 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize8,[foldername1 '/' newFile8]);
                   imagesphereresize8_small = imresize(imagesphereresize8, optionHandler.ScalingFactor);
                   newFile8s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize8_small,[foldername1 '/' newFile8s]);
                   clear imagesphere8;
                   clear imagesphereresize8;
                   clear imagesphereresize8_small;

                   imagesphereresize9 = zeros(7168,7168);
            	   newFile9 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize9,[foldername1 '/' newFile9]);
                   imagesphereresize9_small = imresize(imagesphereresize9, optionHandler.ScalingFactor);
                   newFile9s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize9_small,[foldername1 '/' newFile9s]);
                   clear imagesphere9;
                   clear imagesphereresize9;
                   clear imagesphereresize9_small;

                   imagesphereresize10 = zeros(7168,7168);
            	   newFile10 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize10,[foldername1 '/' newFile10]);
                   imagesphereresize10_small = imresize(imagesphereresize10, optionHandler.ScalingFactor);
                   newFile10s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize10_small,[foldername1 '/' newFile10s]);
                   clear imagesphere10;
                   clear imagesphereresize10;
                   clear imagesphereresize10_small;
                   clear imgneurite;

               end
               
               if(Oligo_Channel > 0)
                   % Same for Oligos
                   imgoligo = imread([foldername '/' newFileOli]);
                   %imwrite(imgoligo,[foldername2 '/' newFileOli]);
                   imagesphere11 = imgoligo(y_1:y_2,x_1:x_2);
                   [r1 c1] = size(imagesphere11);
                   imagesphereresize11 = zeros(7168,7168);
                   imagesphereresize11 = uint8(imagesphereresize11);
                   imagesphereresize11(yStart1:yEnd1,xStart1:xEnd1)=imagesphere11(1:r1,1:c1);
                   newFile11 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize11,[foldername1 '/' newFile11]);
                   imagesphereresize11_small = imresize(imagesphereresize11, optionHandler.ScalingFactor);
                   newFile11s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize11_small,[foldername1 '/' newFile11s]);
                   clear imagesphere11;
                   clear imagesphereresize11;
                   clear imagesphereresize11_small;

                   imagesphere12 = imgoligo(y_3:y_4,x_3:x_4);
                   [r2 c2] = size(imagesphere12);
                   imagesphereresize12 =  zeros(7168,7168);
                   imagesphereresize12 = uint8(imagesphereresize12);
                   imagesphereresize12(yStart2:yEnd2,xStart2:xEnd2)=imagesphere12(1:r2,1:c2);
                   newFile12 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize12,[foldername1 '/' newFile12]);
                   imagesphereresize12_small = imresize(imagesphereresize12, optionHandler.ScalingFactor);
                   newFile12s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize12_small,[foldername1 '/' newFile12s]);
                   clear imagesphere12;
                   clear imagesphereresize12;
                   clear imagesphereresize12_small;

                   imagesphere13 = imgoligo(y_5:y_6,x_5:x_6);
                   [r3 c3] = size(imagesphere13);
                   imagesphereresize13 = zeros(7168,7168);
                   imagesphereresize13 = uint8(imagesphereresize13);
                   imagesphereresize13(yStart3:yEnd3,xStart3:xEnd3)=imagesphere13(1:r3,1:c3);
                   newFile13 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize13,[foldername1 '/' newFile13]);
                   imagesphereresize13_small = imresize(imagesphereresize13, optionHandler.ScalingFactor);
                   newFile13s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize13_small,[foldername1 '/' newFile13s]);
                   clear imagesphere13;
                   clear imagesphereresize13;
                   clear imagesphereresize13_small;


                   imagesphereresize14 = zeros(7168,7168);
            	   newFile14 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize14,[foldername1 '/' newFile14]);
                   imagesphereresize14_small = imresize(imagesphereresize14, optionHandler.ScalingFactor);
                   newFile14s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize14_small,[foldername1 '/' newFile14s]);
                   clear imagesphere14;
                   clear imagesphereresize14;
                   clear imagesphereresize14_small;

                   imagesphereresize15 = zeros(7168,7168);
            	   newFile15 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize15,[foldername1 '/' newFile15]);
                   imagesphereresize15_small = imresize(imagesphereresize15, optionHandler.ScalingFactor);
                   newFile15s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize15_small,[foldername1 '/' newFile15s]);
                   clear imagesphere15;
                   clear imagesphereresize15;
                   clear imagesphereresize15_small;
                   clear imgoligo
               end
               
               if(Astro_Channel >0)
                   % Same for Astros
                   
                   imgastro = imread([foldername '/' newFileAst]);
                   %imwrite(imgastro,[foldername2 '/' newFileAst]);
                   imagesphere16 = imgastro(y_1:y_2,x_1:x_2);
                   [r1 c1] = size(imagesphere16);
                   imagesphereresize16 = zeros(7168,7168);
                   imagesphereresize16 = uint8(imagesphereresize16);
                   imagesphereresize16(yStart1:yEnd1,xStart1:xEnd1)=imagesphere16(1:r1,1:c1);
                   newFile16 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize16,[foldername1 '/' newFile16]);
                   imagesphereresize16_small = imresize(imagesphereresize16, optionHandler.ScalingFactor);
                   newFile16s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize16_small,[foldername1 '/' newFile16s]);
                   clear imagesphere16;
                   clear imagesphereresize16;
                   clear imagesphereresize16_small;

                   imagesphere17 = imgastro(y_3:y_4,x_3:x_4);
                   [r2 c2] = size(imagesphere17);
                   imagesphereresize17 = zeros(7168,7168);
                   imagesphereresize17 = uint8(imagesphereresize17);
                   imagesphereresize17(yStart2:yEnd2,xStart2:xEnd2)=imagesphere17(1:r2,1:c2);
                   newFile17 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize17,[foldername1 '/' newFile17]);
                   imagesphereresize17_small = imresize(imagesphereresize17, optionHandler.ScalingFactor);
                   newFile17s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize17_small,[foldername1 '/' newFile17s]);
                   clear imagesphere17;
                   clear imagesphereresize17;
                   clear imagesphereresize17_small;

                   imagesphere18 = imgastro(y_5:y_6,x_5:x_6);
                   [r3 c3] = size(imagesphere18);
                   imagesphereresize18 = zeros(7168,7168);
                   imagesphereresize18 = uint8(imagesphereresize18);
                   imagesphereresize18(yStart3:yEnd3,xStart3:xEnd3)=imagesphere18(1:r3,1:c3);
                   newFile18 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize18,[foldername1 '/' newFile18]);
                   imagesphereresize18_small = imresize(imagesphereresize18, optionHandler.ScalingFactor);
                   newFile18s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize18_small,[foldername1 '/' newFile18s]);
                   clear imagesphere18;
                   clear imagesphereresize18;
                   clear imagesphereresize18_small;

                   imagesphereresize19 = zeros(7168,7168);
            	   newFile19 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize19,[foldername1 '/' newFile19]);
                   imagesphereresize19_small = imresize(imagesphereresize19, optionHandler.ScalingFactor);
                   newFile19s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize19_small,[foldername1 '/' newFile19s]);
                   clear imagesphere19;
                   clear imagesphereresize19;
                   clear imagesphereresize19_small;

                   imagesphereresize20 = zeros(7168,7168);
            	   newFile20 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize20,[foldername1 '/' newFile20]);
                   imagesphereresize20_small = imresize(imagesphereresize20, optionHandler.ScalingFactor);
                   newFile20s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize20_small,[foldername1 '/' newFile20s]);
                   clear imagesphere20;
                   clear imagesphereresize20;
                   clear imagesphereresize20_small;
                   clear imgastro;
               end
              
               %delete(allfiles(i).name);
                %if(Neuron_Channel > 0)
                %delete(newFileNeu);
                %end
                %if(Oligo_Channel > 0)
                %delete(newFileOli);
                %end
                %if(Astro_Channel > 0)
                %delete(newFileAst);
                %end
                
            end    
            
            % Case 5: Four spheres remained: s=4
            
            if(s_L>3 && s_L<5);
               % Extracts the x and y coordinates from s 
               x_centroid(1) = s(1,1);
               y_centroid(1) = s(1,2);
               x_centroid(2) = s(2,1);
               y_centroid(2) = s(2,2);
               x_centroid(3) = s(3,1);
               y_centroid(3) = s(3,2);
               x_centroid(4) = s(4,1);
               y_centroid(4) = s(4,2);
               % Definition of a centroid quadrant in which the centroid of a potential center neurosphere has to be located. 
               % All other centroids are considerd as rim coordinates
               quadrant1_x1 = 4032;
               quadrant1_x2 = 6720;
               quadrant1_y1 = 5184;
               quadrant1_y2 = 8640;
               % Calculates the distance between the centroids of the
               % four neurosphere
               d1_2 = 0.5* ((x_centroid(1)-x_centroid(2))^2+(y_centroid(1)-y_centroid(2))^2)^0.5;
               d1_3 = 0.5* ((x_centroid(1)-x_centroid(3))^2+(y_centroid(1)-y_centroid(3))^2)^0.5;
               d1_4 = 0.5* ((x_centroid(1)-x_centroid(4))^2+(y_centroid(1)-y_centroid(4))^2)^0.5;
               d2_3 = 0.5* ((x_centroid(3)-x_centroid(2))^2+(y_centroid(3)-y_centroid(2))^2)^0.5;
               d2_4 = 0.5* ((x_centroid(4)-x_centroid(2))^2+(y_centroid(4)-y_centroid(2))^2)^0.5;
               d3_4 = 0.5* ((x_centroid(3)-x_centroid(4))^2+(y_centroid(3)-y_centroid(4))^2)^0.5;
               % The shortest distance between neurosphere cores is the
               % distance between the rim spheres and the center sphere.
               % Therefore for all centroids it is checked wether the coordinates of 
               % their centroid is located within the center quadrant.
                if(x_centroid(1)>=quadrant1_x1 && x_centroid(1)<=quadrant1_x2 && y_centroid(1)>=quadrant1_y1 && y_centroid(1)<=quadrant1_y2);
                    % Since more than one sphere might be located in the
                    % centroid quadrant, further the distance between center
                    % spheres and the middle of the cover slip (dc_x) is
                    % calculated
                    dc_1 = ((x_centroid(1)-5376)^2+(y_centroid(1)-6912)^2)^0.5;
                else
                    % If the centroid is not located in the center quadrant the distance is set to an abitrary high value of dc_1 = 1000000000
                    dc_1 = 1000000000;
                end 
                % Same for sphere 2,3 and 4
                if(x_centroid(2)>=quadrant1_x1 && x_centroid(2)<=quadrant1_x2 && y_centroid(2)>=quadrant1_y1 && y_centroid(2)<=quadrant1_y2);
                    dc_2 = ((x_centroid(2)-5376)^2+(y_centroid(2)-6912)^2)^0.5;
                    
                else
                    
                    dc_2 = 1000000000;
                end
                
                if(x_centroid(3)>=quadrant1_x1 && x_centroid(3)<=quadrant1_x2 && y_centroid(3)>=quadrant1_y1 && y_centroid(3)<=quadrant1_y2);
                    dc_3 = ((x_centroid(3)-5376)^2+(y_centroid(3)-6912)^2)^0.5;
                    
                else
                    
                    dc_3 = 1000000000;
                end
                
                if(x_centroid(4)>=quadrant1_x1 && x_centroid(4)<=quadrant1_x2 && y_centroid(4)>=quadrant1_y1 && y_centroid(4)<=quadrant1_y2);
                    dc_4 = ((x_centroid(4)-5376)^2+(y_centroid(4)-6912)^2)^0.5;
                    
                else
                    
                    dc_4 = 1000000000;
                end
               % It is checked whether dc_1 is smaller than 1000000000 (no centroid sphere according to the above checked criteria) and smaller
               % than the distances of other centroid spheres. If dc_1 is
               % the smallest distance the centroid of sphere 1 is defined
               % as the centroid of the center sphere 
               if(dc_1 < 1000000000 && dc_1 < dc_2 && dc_1 < dc_3 && dc_1 < dc_4);
                   x_centroid_c = x_centroid(1);
                   y_centroid_c = y_centroid(1);
                   missingc1 = 0;
                 
                   
               else
                   % If above criteria are not fullfilled the centroid sphere
                   % is missing and the centroid of this sphere is treated as
                   % a rim sphere
                   missingc1 = 1;
               end    
               %Same for sphere 2, 3 and 4 
               if(dc_2 < 1000000000 && dc_2 < dc_1 && dc_2 < dc_3 && dc_2 < dc_4);
                   x_centroid_c = x_centroid(2);
                   y_centroid_c = y_centroid(2);
                   missingc2 = 0;
               else
                   missingc2 = 1;
               end              
                
               if(dc_3 < 1000000000 && dc_3 < dc_1 && dc_3 < dc_2 && dc_3 < dc_4);
                   x_centroid_c = x_centroid(3);
                   y_centroid_c = y_centroid(3);
                   missingc3 = 0;
           
               else
                   missingc3 = 1;
               end
               
               if(dc_4 < 1000000000 && dc_4 < dc_1 && dc_4 < dc_2 && dc_4 < dc_3);
                   x_centroid_c = x_centroid(4);
                   y_centroid_c = y_centroid(4);
                   missingc4 = 0;
           
               else
                   missingc4 = 1;
               end
               % If all 4 spheres were defined as rim spheres, an artifical
               % center point is defined as the center point of the cover
               % slip
               if((missingc1+missingc2+missingc3+missingc4) > 3);
                   x_centroid_c = 5376;
                   y_centroid_c = 6912;
               else
               end    
               % If sphere 1 is no center sphere the shortest distance
               % between sphere 1 and the center sphere is calculated
               
            
               % Defines wether centroid lays left or right of the center
               % sphere! If sphere 1 is not the center sphere
                if(x_centroid(1) ~= x_centroid_c);
                   % First case: Check whether x_centroid(1) is left of the center sphere
                   if(x_centroid(1) < x_centroid_c);
                      % Check if y_centroid(1) is below y_centroid_c 
                      if(y_centroid(1)<y_centroid_c); 
                       % If y_centroid(1)<y_centroid_c the sphere has to be located at the
                       % upper left corner of the cover slip
                       x_centroid_l(1) = x_centroid(1);
                       y_centroid_l(1) = y_centroid(1);
                       % Sphere is located in the upper left of the cover
                       % slip
                       l_1_1 = 1
                       % Sphere is located in the lower left of the cover
                       % slip
                       l_1_2 = 0
                       % Sphere is located in the upper right of the cover
                       % slip
                       r_1_4 = 0
                       % Sphere is located in the lower right of the cover
                       % slip
                       r_1_5 = 0
                       
                      else
                       % If y_centroid(1)>y_centroid_c the sphere has to be located at the
                       % lower left corner of the cover slip   
                       x_centroid_l(2) = x_centroid(1);
                       y_centroid_l(2) = y_centroid(1);
                       l_1_1 = 0
                       l_1_2 = 1
                       r_1_4 = 0
                       r_1_5 = 0 
                      end 
                   else
                   % Second case: Check whether x_centroid(1) is right of the center sphere    
                      % Check if y_centroid(1) is below y_centroid_c)
                      if(y_centroid(1)<y_centroid_c) 
                       % If y_centroid(1)<y_centroid_c the sphere has to be located at the
                       % upper right corner of the cover slip   
                       x_centroid_r(4) = x_centroid(1);
                       y_centroid_r(4) = y_centroid(1);
                       l_1_1 = 0
                       l_1_2 = 0
                       r_1_4 = 1
                       r_1_5 = 0
                      else
                       % If y_centroid(1)>y_centroid_c the sphere has to be located at the
                       % lower right corner of the cover slip   
                       x_centroid_r(5) = x_centroid(1);
                       y_centroid_r(5) = x_centroid(1);
                       l_1_1 = 0
                       l_1_2 = 0
                       r_1_4 = 0
                       r_1_5 = 1 
                      end 
                   end    
                else
                    l_1_1 = 0
                    l_1_2 = 0
                    r_1_4 = 0
                    r_1_5 = 0
                end    
               % Same for sphere 2, 3 and 4
               if(x_centroid(2) ~= x_centroid_c)
                   if(x_centroid(2) < x_centroid_c);
                      if(y_centroid(2)<y_centroid_c) 
                       x_centroid_l(1) = x_centroid(2);
                       y_centroid_l(1) = y_centroid(2);
                       l_2_1 = 1
                       l_2_2 = 0
                       r_2_4 = 0
                       r_2_5 = 0
                       
                      else
                       x_centroid_l(2) = x_centroid(2);
                       y_centroid_l(2) = y_centroid(2);
                       l_2_1 = 0
                       l_2_2 = 1
                       r_2_4 = 0
                       r_2_5 = 0 
                      end 
                   else
                      if(y_centroid(2)<y_centroid_c) 
                       x_centroid_r(4) = x_centroid(2);
                       y_centroid_r(4) = y_centroid(2);
                       l_2_1 = 0
                       l_2_2 = 0
                       r_2_4 = 1
                       r_2_5 = 0
                      else
                       x_centroid_r(5) = x_centroid(2);
                       y_centroid_r(5) = x_centroid(2);
                       l_2_1 = 0
                       l_2_2 = 0
                       r_2_4 = 0
                       r_2_5 = 1 
                      end 
                   end    
                else
                    l_2_1 = 0
                    l_2_2 = 0
                    r_2_4 = 0
                    r_2_5 = 0
                end 
               
               if(x_centroid(3) ~= x_centroid_c)
                   if(x_centroid(3) < x_centroid_c);
                      if(y_centroid(3)<y_centroid_c) 
                       x_centroid_l(1) = x_centroid(3);
                       y_centroid_l(1) = y_centroid(3);
                       l_3_1 = 1
                       l_3_2 = 0
                       r_3_4 = 0
                       r_3_5 = 0
                       
                      else
                       x_centroid_l(2) = x_centroid(3);
                       y_centroid_l(2) = y_centroid(3);
                       l_3_1 = 0
                       l_3_2 = 1
                       r_3_4 = 0
                       r_3_5 = 0 
                      end 
                   else
                      if(y_centroid(3)<y_centroid_c) 
                       x_centroid_r(4) = x_centroid(3);
                       y_centroid_r(4) = y_centroid(3);
                       l_3_1 = 0
                       l_3_2 = 0
                       r_3_4 = 1
                       r_3_5 = 0
                      else
                       x_centroid_r(5) = x_centroid(3);
                       y_centroid_r(5) = x_centroid(3);
                       l_3_1 = 0
                       l_3_2 = 0
                       r_3_4 = 0
                       r_3_5 = 1 
                      end 
                   end    
               else
                    l_3_1 = 0
                    l_3_2 = 0
                    r_3_4 = 0
                    r_3_5 = 0
               end
               
               
               if(x_centroid(4) ~= x_centroid_c);
                   if(x_centroid(4) < x_centroid_c);
                      if(y_centroid(4)<y_centroid_c); 
                       x_centroid_l(1) = x_centroid(4);
                       y_centroid_l(1) = y_centroid(4);
                       l_4_1 = 1
                       l_4_2 = 0
                       r_4_4 = 0
                       r_4_5 = 0
                       
                      else
                       x_centroid_l(2) = x_centroid(4);
                       y_centroid_l(2) = y_centroid(4);
                       l_4_1 = 0
                       l_4_2 = 1
                       r_4_4 = 0
                       r_4_5 = 0 
                      end 
                   else
                      if(y_centroid(4)<y_centroid_c) 
                       x_centroid_r(4) = x_centroid(4);
                       y_centroid_r(4) = y_centroid(4);
                       l_4_1 = 0
                       l_4_2 = 0
                       r_4_4 = 1
                       r_4_5 = 0
                      else
                       x_centroid_r(5) = x_centroid(4);
                       y_centroid_r(5) = x_centroid(4);
                       l_4_1 = 0
                       l_4_2 = 0
                       r_4_4 = 0
                       r_4_5 = 1 
                      end 
                   end    
                else
                    l_4_1 = 0
                    l_4_2 = 0
                    r_4_4 = 0
                    r_4_5 = 0
               end    
               
               
               % Check for multiple entry positions:
               
               % Check if all spheres at the same position
               if(l_2_1 == l_1_1 && l_2_2 == l_1_2 && r_2_4 == r_1_4 && r_2_5 == r_1_5 && l_3_1 == l_1_1 && l_3_2 == l_1_2 && r_3_4 == r_1_4 && r_3_5 == r_1_5 && l_4_1 == l_1_1 && l_4_2 == l_1_2 && r_4_4 == r_1_4 && r_4_5 == r_1_5)
                   
                   l_1_1_q = l_1_1
                   l_1_2_q = l_1_2
                   r_1_4_q = r_1_4
                   r_1_5_q = r_1_5
                   l_2_1_q = l_2_1
                   l_2_2_q = l_2_2
                   r_2_4_q = r_2_4
                   r_2_5_q = r_2_5
                   l_3_1_q = l_3_1
                   l_3_2_q = l_3_2
                   r_3_4_q = r_3_4
                   r_3_5_q = r_3_5
                   l_4_1_q = l_4_1
                   l_4_2_q = l_4_2
                   r_4_4_q = r_4_4
                   r_4_5_q = r_4_5
                   l_1_1 = 0
                   l_1_2 = 0
                   r_1_4 = 0
                   r_1_5 = 0
                   l_2_1 = 0
                   l_2_2 = 0
                   r_2_4 = 0
                   r_2_5 = 0
                   l_3_1 = 0
                   l_3_2 = 0
                   r_3_4 = 0
                   r_3_5 = 0
                   l_4_1 = 0
                   l_4_2 = 0
                   r_4_4 = 0
                   r_4_5 = 0
                   
                else
                   l_1_1_q = 0
                   l_1_2_q = 0
                   r_1_4_q = 0
                   r_1_5_q = 0
                   l_2_1_q = 0
                   l_2_2_q = 0
                   r_2_4_q = 0
                   r_2_5_q = 0
                   l_3_1_q = 0
                   l_3_2_q = 0
                   r_3_4_q = 0
                   r_3_5_q = 0 
                   l_4_1_q = 0
                   l_4_2_q = 0
                   r_4_4_q = 0
                   r_4_5_q = 0
                   
               end    
               
               % Check for Triplets
               % 1,2,3 Triplet
               
               if(l_1_1 == l_2_1 && l_1_2 == l_2_2 && r_1_4 == r_2_4 && r_1_5 == r_2_5 && l_1_1 == l_3_1 && l_1_2 == l_3_2 && r_1_4 == r_3_4 && r_1_5 == r_3_5)
                   l_1_1_t = l_1_1
                   l_1_2_t = l_1_2
                   r_1_4_t = r_1_4
                   r_1_5_t = r_1_5
                   l_2_1_t = l_2_1
                   l_2_2_t = l_2_2
                   r_2_4_t = r_2_4
                   r_2_5_t = r_2_5
                   l_3_1_t = l_3_1
                   l_3_2_t = l_3_2
                   r_3_4_t = r_3_4
                   r_3_5_t = r_3_5
                   l_4_1_t = 0
                   l_4_2_t = 0
                   r_4_4_t = 0
                   r_4_5_t = 0
                   l_1_1 = 0
                   l_1_2 = 0
                   r_1_4 = 0
                   r_1_5 = 0
                   l_2_1 = 0
                   l_2_2 = 0
                   r_2_4 = 0
                   r_2_5 = 0
                   l_3_1 = 0
                   l_3_2 = 0
                   r_3_4 = 0
                   r_3_5 = 0
                   
                   
                   
               else
                   l_1_1_t = 0
                   l_1_2_t = 0
                   r_1_4_t = 0
                   r_1_5_t = 0
                   l_2_1_t = 0
                   l_2_2_t = 0
                   r_2_4_t = 0
                   r_2_5_t = 0
                   l_3_1_t = 0
                   l_3_2_t = 0
                   r_3_4_t = 0
                   r_3_5_t = 0 
                 
                   
               end  
               
               % 1,2,4 Triplet
               if(l_1_1 == l_2_1 && l_1_2 == l_2_2 && r_1_4 == r_2_4 && r_1_5 == r_2_5 && l_1_1 == l_3_1 && l_1_2 == l_3_2 && r_1_4 == r_3_4 && r_1_5 == r_3_5)
                   t = 1
               else    
               if(l_1_1 == l_2_1 && l_1_2 == l_2_2 && r_1_4 == r_2_4 && r_1_5 == r_2_5 && l_1_1 == l_4_1 && l_1_2 == l_4_2 && r_1_4 == r_4_4 && r_1_5 == r_4_5)
                   l_1_1_t = l_1_1
                   l_1_2_t = l_1_2
                   r_1_4_t = r_1_4
                   r_1_5_t = r_1_5
                   l_2_1_t = l_2_1
                   l_2_2_t = l_2_2
                   r_2_4_t = r_2_4
                   r_2_5_t = r_2_5
                   l_3_1_t = 0
                   l_3_2_t = 0
                   r_3_4_t = 0
                   r_3_5_t = 0
                   l_4_1_t = l_4_1
                   l_4_2_t = l_4_2
                   r_4_4_t = r_4_4
                   r_4_5_t = r_4_5
                   l_1_1 = 0
                   l_1_2 = 0
                   r_1_4 = 0
                   r_1_5 = 0
                   l_2_1 = 0
                   l_2_2 = 0
                   r_2_4 = 0
                   r_2_5 = 0
                   l_4_1 = 0
                   l_4_2 = 0
                   r_4_4 = 0
                   r_4_5 = 0
                   
                   t = 1
                   
               else
                   l_1_1_t = 0
                   l_1_2_t = 0
                   r_1_4_t = 0
                   r_1_5_t = 0
                   l_2_1_t = 0
                   l_2_2_t = 0
                   r_2_4_t = 0
                   r_2_5_t = 0
                   l_4_1_t = 0
                   l_4_2_t = 0
                   r_4_4_t = 0
                   r_4_5_t = 0 
                   t = 0
                   
               end 
               end
               
               
               % 1,3,4 Triplet
               if(t == 1)
                   t = 1
               else    
               if(l_1_1 == l_3_1 && l_1_2 == l_3_2 && r_1_4 == r_3_4 && r_1_5 == r_3_5 && l_1_1 == l_4_1 && l_1_2 == l_4_2 && r_1_4 == r_4_4 && r_1_5 == r_4_5)
                   l_1_1_t = l_1_1
                   l_1_2_t = l_1_2
                   r_1_4_t = r_1_4
                   r_1_5_t = r_1_5
                   l_2_1_t = 0
                   l_2_2_t = 0
                   r_2_4_t = 0
                   r_2_5_t = 0
                   l_3_1_t = l_3_1
                   l_3_2_t = l_3_2
                   r_3_4_t = r_3_4
                   r_3_5_t = r_3_5
                   l_4_1_t = l_4_1
                   l_4_2_t = l_4_2
                   r_4_4_t = r_4_4
                   r_4_5_t = r_4_5
                   l_1_1 = 0
                   l_1_2 = 0
                   r_1_4 = 0
                   r_1_5 = 0
                   l_3_1 = 0
                   l_3_2 = 0
                   r_3_4 = 0
                   r_3_5 = 0
                   l_4_1 = 0
                   l_4_2 = 0
                   r_4_4 = 0
                   r_4_5 = 0
                   t=1
                   
                   
               else
                   l_1_1_t = 0
                   l_1_2_t = 0
                   r_1_4_t = 0
                   r_1_5_t = 0
                   l_3_1_t = 0
                   l_3_2_t = 0
                   r_3_4_t = 0
                   r_3_5_t = 0
                   l_4_1_t = 0
                   l_4_2_t = 0
                   r_4_4_t = 0
                   r_4_5_t = 0 
                   t = 0
                   
                   
               end 
               end
               
               % 2,3,4 Triplet
               if(t==1)
                  
               else    
               if(l_3_1 == l_2_1 && l_3_2 == l_2_2 && r_3_4 == r_2_4 && r_3_5 == r_2_5 && l_4_1 == l_2_1 && l_4_2 == l_2_2 && r_4_4 == r_2_4 && r_4_5 == r_2_5)
                   
                   l_1_1_t = 0
                   l_1_2_t = 0
                   r_1_4_t = 0
                   r_1_5_t = 0
                   l_2_1_t = l_2_1
                   l_2_2_t = l_2_2
                   r_2_4_t = r_2_4
                   r_2_5_t = r_2_5
                   l_3_1_t = l_3_1
                   l_3_2_t = l_3_2
                   r_3_4_t = r_3_4
                   r_3_5_t = r_3_5
                   l_4_1_t = l_4_1
                   l_4_2_t = l_4_2
                   r_4_4_t = r_4_4
                   r_4_5_t = r_4_5
                   l_2_1 = 0
                   l_2_2 = 0
                   r_2_4 = 0
                   r_2_5 = 0
                   l_3_1 = 0
                   l_3_2 = 0
                   r_3_4 = 0
                   r_3_5 = 0
                   l_4_1 = 0
                   l_4_2 = 0
                   r_4_4 = 0
                   r_4_5 = 0
                   
                   
               else
                   l_2_1_t = 0
                   l_2_2_t = 0
                   r_2_4_t = 0
                   r_2_5_t = 0
                   l_3_1_t = 0
                   l_3_2_t = 0
                   r_3_4_t = 0
                   r_3_5_t = 0 
                   l_4_1_t = 0
                   l_4_2_t = 0
                   r_4_4_t = 0
                   r_4_5_t = 0
                   
               end  
               end
               % Check for doublets
               
               % Same position of 1,2
               
               if(l_1_1 == l_2_1 && l_1_2 == l_2_2 && r_1_4 == r_2_4 && r_1_5 == r_2_5)
                 
                   l_1_1_d12 = l_1_1
                   l_1_2_d12 = l_1_2
                   r_1_4_d12 = r_1_4
                   r_1_5_d12 = r_1_5
                   l_2_1_d12 = l_2_1
                   l_2_2_d12 = l_2_2
                   r_2_4_d12 = r_2_4
                   r_2_5_d12 = r_2_5
                   l_1_1 = 0
                   l_1_2 = 0
                   r_1_4 = 0
                   r_1_5 = 0
                   l_2_1 = 0
                   l_2_2 = 0
                   r_2_4 = 0
                   r_2_5 = 0
                                              
               else
                   l_1_1_d12 = 0
                   l_1_2_d12 = 0
                   r_1_4_d12 = 0
                   r_1_5_d12 = 0
                   l_2_1_d12 = 0
                   l_2_2_d12 = 0
                   r_2_4_d12 = 0
                   r_2_5_d12 = 0 
                   
               end
              
               % Same position of 1,3
               
               if(l_1_1 == l_3_1 && l_1_2 == l_3_2 && r_1_4 == r_3_4 && r_1_5 == r_3_5)
                 
                   l_1_1_d13 = l_1_1
                   l_1_2_d13 = l_1_2
                   r_1_4_d13 = r_1_4
                   r_1_5_d13 = r_1_5
                   l_3_1_d13 = l_2_1
                   l_3_2_d13 = l_2_2
                   r_3_4_d13 = r_2_4
                   r_3_5_d13 = r_2_5
                   l_1_1 = 0
                   l_1_2 = 0
                   r_1_4 = 0
                   r_1_5 = 0
                   l_3_1 = 0
                   l_3_2 = 0
                   r_3_4 = 0
                   r_3_5 = 0
                                              
               else
                   l_1_1_d13 = 0
                   l_1_2_d13 = 0
                   r_1_4_d13 = 0
                   r_1_5_d13 = 0
                   l_3_1_d13 = 0
                   l_3_2_d13 = 0
                   r_3_4_d13 = 0
                   r_3_5_d13 = 0 
                   
               end
              
               % Same position of 1,4
               
               if(l_1_1 == l_4_1 && l_1_2 == l_4_2 && r_1_4 == r_4_4 && r_1_5 == r_4_5)
                 
                   l_1_1_d14 = l_1_1
                   l_1_2_d14 = l_1_2
                   r_1_4_d14 = r_1_4
                   r_1_5_d14 = r_1_5
                   l_4_1_d14 = l_2_1
                   l_4_2_d14 = l_2_2
                   r_4_4_d14 = r_2_4
                   r_4_5_d14 = r_2_5
                   l_1_1 = 0
                   l_1_2 = 0
                   r_1_4 = 0
                   r_1_5 = 0
                   l_4_1 = 0
                   l_4_2 = 0
                   r_4_4 = 0
                   r_4_5 = 0
                                              
               else
                   l_1_1_d14 = 0
                   l_1_2_d14 = 0
                   r_1_4_d14 = 0
                   r_1_5_d14 = 0
                   l_4_1_d14 = 0
                   l_4_2_d14 = 0
                   r_4_4_d14 = 0
                   r_4_5_d14 = 0 
                   
              end
                       
                   
              % Same position of 2,3   
               
               
               if(l_2_1 == l_3_1 && l_2_2 == l_3_2 && r_2_4 == r_3_4 && r_2_5 == r_3_5)
                   B = 1
                   l_2_1_d23 = l_2_1
                   l_2_2_d23 = l_2_2
                   r_2_4_d23 = r_2_4
                   r_2_5_d23 = r_2_5
                   l_3_1_d23 = l_3_1
                   l_3_2_d23 = l_3_2
                   r_3_4_d23 = r_3_4
                   r_3_5_d23 = r_3_5
                   l_2_1 = 0
                   l_2_2 = 0
                   r_2_4 = 0
                   r_2_5 = 0
                   l_3_1 = 0
                   l_3_2 = 0
                   r_3_4 = 0
                   r_3_5 = 0
                   
               else
                   l_2_1_d23 = 0
                   l_2_2_d23 = 0
                   r_2_4_d23 = 0
                   r_2_5_d23 = 0
                   l_3_1_d23 = 0
                   l_3_2_d23 = 0
                   r_3_4_d23 = 0
                   r_3_5_d23 = 0
                   
               end    
               
               % Same position of 2,4   
               
               
               if(l_2_1 == l_4_1 && l_2_2 == l_4_2 && r_2_4 == r_4_4 && r_2_5 == r_4_5)
                   B = 1
                   l_2_1_d24 = l_2_1
                   l_2_2_d24 = l_2_2
                   r_2_4_d24 = r_2_4
                   r_2_5_d24 = r_2_5
                   l_4_1_d24 = l_3_1
                   l_4_2_d24 = l_3_2
                   r_4_4_d24 = r_3_4
                   r_4_5_d24 = r_3_5
                   l_2_1 = 0
                   l_2_2 = 0
                   r_2_4 = 0
                   r_2_5 = 0
                   l_4_1 = 0
                   l_4_2 = 0
                   r_4_4 = 0
                   r_4_5 = 0
                   
               else
                   l_2_1_d24 = 0
                   l_2_2_d24 = 0
                   r_2_4_d24 = 0
                   r_2_5_d24 = 0
                   l_4_1_d24 = 0
                   l_4_2_d24 = 0
                   r_4_4_d24 = 0
                   r_4_5_d24 = 0
                   
               end  
               % Same position of 3,4   
               
               
               if(l_4_1 == l_3_1 && l_4_2 == l_3_2 && r_4_4 == r_3_4 && r_4_5 == r_3_5)
                   B = 1
                   l_3_1_d34 = l_3_1
                   l_3_2_d34 = l_3_2
                   r_3_4_d34 = r_3_4
                   r_3_5_d34 = r_3_5
                   l_4_1_d34 = l_4_1
                   l_4_2_d34 = l_4_2
                   r_4_4_d34 = r_4_4
                   r_4_5_d34 = r_4_5
                   l_3_1 = 0
                   l_3_2 = 0
                   r_3_4 = 0
                   r_3_5 = 0
                   l_4_1 = 0
                   l_4_2 = 0
                   r_4_4 = 0
                   r_4_5 = 0
                   
               else
                   l_3_1_d34 = 0
                   l_3_2_d34 = 0
                   r_3_4_d34 = 0
                   r_3_5_d34 = 0
                   l_4_1_d34 = 0
                   l_4_2_d34 = 0
                   r_4_4_d34 = 0
                   r_4_5_d34 = 0
                   
               end 
               
               
               %Calculate max distance in each direction
               if(x_centroid(1) ~= x_centroid_c);
                  if((missingc1+missingc2+missingc3+missingc4) > 3);
                      % If there is no center sphere identified, the
                      % distance to the center of the image can be adjusted
                      % to help not to cut far migrated cells at day 5!
             
                      d_1_c = MinDistCentr
                  else
                      d_1_c = 0.5* (((x_centroid(1)-x_centroid_c)^2+(y_centroid(1)-y_centroid_c)^2)^0.5);
                  end    
               else
               end   
               % Same for sphere 2, 3 and 4
               if(x_centroid(2) ~= x_centroid_c);
                   if((missingc1+missingc2+missingc3+missingc4) > 3);
                       d_2_c = MinDistCentr
                   else
                       d_2_c = 0.5* (((x_centroid(2)-x_centroid_c)^2+(y_centroid(2)-y_centroid_c)^2)^0.5);
                   end    
               else
               end
               if(x_centroid(3) ~= x_centroid_c);
                   if((missingc1+missingc2+missingc3+missingc4) > 3);
                       d_3_c = MinDistCentr
                   else    
                       d_3_c = 0.5* (((x_centroid(3)-x_centroid_c)^2+(y_centroid(3)-y_centroid_c)^2)^0.5);
                   end    
               else
               end
               
               if(x_centroid(4) ~= x_centroid_c);
                   if((missingc1+missingc2+missingc3+missingc4) > 3);
                    d_4_c = MinDistCentr
                   else
                    d_4_c = 0.5* (((x_centroid(4)-x_centroid_c)^2+(y_centroid(4)-y_centroid_c)^2)^0.5);
                   end    
               else
                   V_4 = [d1_4 ; d2_4 ; d3_4 ; (7168/2)];
                   min4 = min(V_4(1:4));
                   d_4_c = min4;
               end
               % Sphere 1 with smallest value of x_centroid
               
               % Go through potential positions
               % Case 1: Sphere is no center sphere and is located at the
               % upper left
               if(l_1_1 >0 && l_1_1 < 2);
                       % Coordinates to cut out sphere 1
                       x_1 = 1;
                       x_2 = x_centroid(1) + d_1_c;
                       y_1 = 1;
                       y_2 = y_centroid(1) + d_1_c;
               else
               end    
               % Case 2: Sphere is no center sphere and is located at the
               % lower left
               if(l_1_2 >0 && l_1_2 < 2);
                       % Coordinates to cut out sphere 1
                       x_1 = 1;
                       x_2 = x_centroid(1) + d_1_c;
                       y_1 = y_centroid(1) - d_1_c;
                       y_2 = 13824;
               else
               end
               % Case 3: Sphere is no center sphere and is located at the
               % upper right
               if(r_1_4 >0 && r_1_4 < 2);
                       % Coordinates to cut out sphere 1    
                       x_1 = x_centroid(1) - d_1_c;
                       x_2 = 10752;
                       y_1 = 1;
                       y_2 = y_centroid(1) + d_1_c;
               else
               end
               % Case 4: Sphere is no center sphere and is located at the
               % lower right
               if(r_1_5 >0 && r_1_5 < 2);
                       % Coordinates to cut out sphere 1
                       x_1 = x_centroid(1) - d_1_c;
                       x_2 = 10752;
                       y_1 = y_centroid(1) - d_1_c;
                       y_2 = 13824;
               else
               end
               
               % Going through all iterations for the case sphere 1 =
               % center sphere
               if(x_centroid(1) ~= x_centroid_c);
               else
                   % Case 1: sphere 2, 3, 4 are located left of sphere 1
                   if(x_centroid(2) < x_centroid(1) && x_centroid(3) < x_centroid(1) && x_centroid(4) < x_centroid(1));
                       % Find shortest distance to sphere 1 to sphere 2,3
                       % and 3 to a maximum of 7168/2. Since sphere 2, 3
                       % and 4 are located left of sphere 1, distance is just
                       % limited to the left and not to the right
                       V_1m = [d1_2 ; d1_3 ; d1_4 ; (7168/2)];
                       min1m = min(V_1m);
                       % Maximal added distance to the right to obtain
                       % maximal images of 7168x7168 pixel
                       min1p = (7168-min1m);
                   else
                   end
                   % Case 2: sphere 2, 3 and 4 are located right of sphere 1
                   if(x_centroid(2) > x_centroid(1) && x_centroid(3) > x_centroid(1) && x_centroid(4) > x_centroid(1));
                       % Find shortest distance to sphere 1 to sphere 2, 3 and
                       % 4 to a maximum of 7168/2. Since sphere 2, 3
                       % and 4 are located right of sphere 1, distance is just
                       % limited to the right and not to the left
                       V_1p = [d1_2 ; d1_3 ; d1_4 (7168/2)];
                       min1p = min(V_1p);
                       % Maximal added distance to the left to obtain
                       % maximal images of 7168x7168 pixel
                       min1m = (7168-min1p);
                   else
                   end
                   % Case 3: sphere 2 and 3 are located left of sphere 1 and 4 is located right of sphere 1
                   if(x_centroid(2) < x_centroid(1) && x_centroid(3) < x_centroid(1) && x_centroid(4) > x_centroid(1));
                       % Find minimal distance to the left and to the right. To the left
                       % only the distance between sphere 1 and 2 and 1 and 3 and to
                       % the right only the distance of sphere 1 and 4
                       % matter. 
                       V_1m = [d1_2 ; d1_3 ; (7168/2)];
                       V_1p = [d1_4 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);  
                   else  
                   end
                   % Case 4: sphere 2 and 4 are located left of sphere 1 and 3 is located right of sphere 1
                   if(x_centroid(2) < x_centroid(1) && x_centroid(3) > x_centroid(1) && x_centroid(4) < x_centroid(1));
                       % Find minimal distance to the left and to the right. To the left
                       % only the distance between sphere 1 and 2 and 1 and 4 and to
                       % the right only the distance of sphere 1 and 4
                       % matter.
                       V_1m = [d1_2 ; d1_4 ; (7168/2)];
                       V_1p = [d1_3 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   % Case 5: sphere 3 and 4 are located left of sphere 1 and 2 is located right of sphere 1
                   if(x_centroid(2) > x_centroid(1) && x_centroid(3) < x_centroid(1) && x_centroid(4) < x_centroid(1));
                       % Find minimal distance to the left and to the right. To the left
                       % only the distance between sphere 1 and 3 and 1 and 4 and to
                       % the right only the distance of sphere 1 and 2
                       % matter.
                       V_1m = [d1_3 ; d1_4 ; (7168/2)];
                       V_1p = [d1_2 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   % Case 6: sphere 2 is located left of sphere 1 and 3 and 4 are located right of sphere 1
                   if(x_centroid(2) < x_centroid(1) && x_centroid(3) > x_centroid(1) && x_centroid(4) > x_centroid(1));
                       % Find minimal distance to the left and to the right. To the left
                       % only the distance between sphere 1 and 2 and to
                       % the right only the distance of sphere 1 and 3 and
                       % 1 and 4 matter.
                       V_1m = [d1_2 ; (7168/2)];
                       V_1p = [d1_3 ; d1_4 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   % Case 7: sphere 4 is located left of sphere 1 and 2 and 3 are located right of sphere 1
                   if(x_centroid(2) > x_centroid(1) && x_centroid(3) > x_centroid(1) && x_centroid(4) < x_centroid(1));
                       % Find minimal distance to the left and to the right. To the left
                       % only the distance between sphere 1 and 4 and to
                       % the right only the distance of sphere 1 and 2 and
                       % 1 and 3 matter.
                       V_1m = [d1_4 ; (7168/2)];
                       V_1p = [d1_2 ; d1_3 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   % Case 8: sphere 3 is located left of sphere 1 and 2 and 4 are located right of sphere 1
                   if(x_centroid(2) > x_centroid(1) && x_centroid(3) < x_centroid(1) && x_centroid(4) > x_centroid(1));
                       % Find minimal distance to the left and to the right. To the left
                       % only the distance between sphere 1 and 3 and to
                       % the right only the distance of sphere 1 and 2 and
                       % 1 and 4 matter.
                       V_1m = [d1_3 ; (7168/2)];
                       V_1p = [d1_2 ; d1_4 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   
                   % Coordinates to cut out sphere 1
                   x_1 = x_centroid(1) - min1m 
                   x_2 = x_centroid(1) + min1p
                   y_1 = y_centroid(1) - min1m
                   y_2 = y_centroid(1) + min1p
                   
                    
               end    
               % Case for all four spheres in the same spot
                   
                   if(l_1_1_q == 1)
                       Vdq_x = [d1_2 ; d1_3 ; d1_4 ; 7168/2]
                       Vdq_y_m = [d1_2 ; d1_3 ; d1_4 ; 7168/2 ; (y_centroid(1)-1)]
                       Vdq_y_p = [d1_2 ; d1_3 ; d1_4 ; 7168/2 ]
                       miny_m = min(Vdq_y_m(1:5))
                       miny_p = min(Vdq_y_m(1:4))
                       minx = min(Vdq_x(1:4))
                       x_1 = 1
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end   
                       
                   if(l_1_2_q == 1)
                       Vdq_x = [d1_2 ; d1_3 ; d1_4 ; 7168/2]
                       Vdq_y_m = [d1_2 ; d1_3 ; d1_4 ; 7168/2]
                       Vdq_y_p = [d1_2 ; d1_3 ; d1_4 ; 7168/2 ; (13824 - y_centroid(1))]
                       miny_m = min(Vdq_y_m(1:4))
                       miny_p = min(Vdq_y_p(1:5))
                       minx = min(Vdq_x(1:4))
                       x_1 = 1
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end 
                   
                   if(r_1_4_q == 1)
                       Vdq_x = [d1_2 ; d1_3 ; d1_4 ; 7168/2]
                       Vdq_y_m = [d1_2 ; d1_3 ; d1_4 ; 7168/2 ; (y_centroid(1)-1)]
                       Vdq_y_p = [d1_2 ; d1_3 ; d1_4 ; 7168/2 ]
                       miny_m = min(Vdq_y_m(1:5))
                       miny_p = min(Vdq_y_p(1:4))
                       minx = min(Vdq_x(1:4))
                       x_1 = x_centroid(1) - (7168 - minx)
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end 
                   
                   if(r_1_5_q == 1)
                       Vdq_x = [d1_2 ; d1_3 ; d1_4 ; 7168/2]
                        Vdq_y_m = [d1_2 ; d1_3 ; d1_4 ; 7168/2]
                       Vdq_y_p = [d1_2 ; d1_3 ; d1_4 ; 7168/2 ; (13824 - y_centroid(1))]
                       miny_m = min(Vdq_y_m(1:4))
                       miny_p = min(Vdq_y_p(1:5))
                       minx = min(Vdq_x(1:4))
                       x_1 = x_centroid(1) - (7168 - minx)
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end
               
                   % Case for 1,2,3 at the same position
                   if(l_1_1_t == 1)
                       Vdt_x = [d1_2 ; d1_3 ; 7168/2]
                       Vdt_y_m = [d1_2 ; d1_3  ; 7168/2 ; (y_centroid(1)-1)]
                       Vdt_y_p = [d1_2 ; d1_3  ; 7168/2 ]
                       miny_m = min(Vdt_y_m(1:4))
                       miny_p = min(Vdt_y_m(1:3))
                       minx = min(Vdt_x(1:3))
                       x_1 = 1
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end   
                       
                   if(l_1_2_t == 1)
                       Vdt_x = [d1_2 ; d1_3 ;  7168/2]
                       Vdt_y_m = [d1_2 ; d1_3 ;  7168/2]
                       Vdt_y_p = [d1_2 ; d1_3 ;  7168/2 ; (13824 - y_centroid(1))]
                       miny_m = min(Vdt_y_m(1:3))
                       miny_p = min(Vdt_y_p(1:4))
                       minx = min(Vdt_x(1:3))
                       x_1 = 1
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end 
                   
                   if(r_1_4_t == 1)
                       Vdt_x = [d1_2 ; d1_3 ;  7168/2]
                       Vdt_y_m = [d1_2 ; d1_3 ;  7168/2 ; (y_centroid(1)-1)]
                       Vdt_y_p = [d1_2 ; d1_3 ;  7168/2 ]
                       miny_m = min(Vdt_y_m(1:4))
                       miny_p = min(Vdt_y_p(1:3))
                       minx = min(Vdt_x(1:3))
                       x_1 = x_centroid(1) - (7168 - minx)
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end 
                   
                   if(r_1_5_t == 1)
                       Vdt_x = [d1_2 ; d1_3 ;  7168/2]
                        Vdt_y_m = [d1_2 ; d1_3 ;  7168/2]
                       Vdt_y_p = [d1_2 ; d1_3 ;  7168/2 ; (13824 - y_centroid(1))]
                       miny_m = min(Vdt_y_m(1:3))
                       miny_p = min(Vdt_y_p(1:4))
                       minx = min(Vdt_x(1:3))
                       x_1 = x_centroid(1) - (7168 - minx)
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end
                   
                   
                   % Case fox 1,2 are the same quadrant
                   
                   if(l_1_1_d12 == 1)
                       Vdd_x = [d1_2 ; 7168/2]
                       Vdd_y_m = [d1_2; 7168/2 ; (y_centroid(1)-1)]
                       Vdd_y_p = [d1_2 7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_m(1:2))
                       minx = min(Vdd_x(1:2))
                       x_1 = 1
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end   
                       
                   if(l_1_2_d12 == 1)
                       Vdd_x = [d1_2 ; 7168/2]
                       Vdd_y_m = [d1_2  ;  7168/2]
                       Vdd_y_p = [d1_2;  7168/2 ; (13824 - y_centroid(1))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_1 = 1
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end 
                   
                   if(r_1_4_d12 == 1)
                       Vdd_x = [d1_2;   7168/2]
                       Vdd_y_m = [d1_2;   7168/2 ; (y_centroid(1)-1)]
                       Vdd_y_p = [d1_2;   7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_p(1:2))
                       minx = min(Vdd_x(1:2))
                       x_1 = x_centroid(1) - (7168 - minx)
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end 
                   
                   if(r_1_5_d12 == 1)
                       Vdd_x = [d1_2;  7168/2]
                        Vdd_y_m = [d1_2;  7168/2]
                       Vdd_y_p = [d1_2;   7168/2 ; (13824 - y_centroid(1))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdt_x(1:2))
                       x_1 = x_centroid(1) - (7168 - minx)
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end
                   
                % Case fox 1,3 are the same quadrant
                   
                   if(l_1_1_d13 == 1)
                       Vdd_x = [d1_3 ; 7168/2]
                       Vdd_y_m = [d1_3; 7168/2 ; (y_centroid(1)-1)]
                       Vdd_y_p = [d1_3 7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_m(1:2))
                       minx = min(Vdd_x(1:2))
                       x_1 = 1
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end   
                       
                   if(l_1_2_d13 == 1)
                       Vdd_x = [d1_3 ; 7168/2]
                       Vdd_y_m = [d1_3  ;  7168/2]
                       Vdd_y_p = [d1_3;  7168/2 ; (13824 - y_centroid(1))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_1 = 1
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end 
                   
                   if(r_1_4_d13 == 1)
                       Vdd_x = [d1_3;   7168/2]
                       Vdd_y_m = [d1_3;   7168/2 ; (y_centroid(1)-1)]
                       Vdd_y_p = [d1_3;   7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_p(1:2))
                       minx = min(Vdd_x(1:2))
                       x_1 = x_centroid(1) - (7168 - minx)
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end 
                   
                   if(r_1_5_d13 == 1)
                       Vdd_x = [d1_3;  7168/2]
                        Vdd_y_m = [d1_3;  7168/2]
                       Vdd_y_p = [d1_3;   7168/2 ; (13824 - y_centroid(1))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdt_x(1:2))
                       x_1 = x_centroid(1) - (7168 - minx)
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end   
                  % Case fox 1,4 are the same quadrant
                   
                   if(l_1_1_d14 == 1)
                       Vdd_x = [d1_4 ; 7168/2]
                       Vdd_y_m = [d1_4; 7168/2 ; (y_centroid(1)-1)]
                       Vdd_y_p = [d1_4 7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_m(1:2))
                       minx = min(Vdd_x(1:2))
                       x_1 = 1
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end   
                       
                   if(l_1_2_d14 == 1)
                       Vdd_x = [d1_4 ; 7168/2]
                       Vdd_y_m = [d1_4  ;  7168/2]
                       Vdd_y_p = [d1_4;  7168/2 ; (13824 - y_centroid(1))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_1 = 1
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end 
                   
                   if(r_1_4_d14 == 1)
                       Vdd_x = [d1_4;   7168/2]
                       Vdd_y_m = [d1_4;   7168/2 ; (y_centroid(1)-1)]
                       Vdd_y_p = [d1_4;   7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_p(1:2))
                       minx = min(Vdd_x(1:2))
                       x_1 = x_centroid(1) - (7168 - minx)
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end 
                   
                   if(r_1_5_d14 == 1)
                       Vdd_x = [d1_4;  7168/2]
                        Vdd_y_m = [d1_4;  7168/2]
                       Vdd_y_p = [d1_4;   7168/2 ; (13824 - y_centroid(1))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdt_x(1:2))
                       x_1 = x_centroid(1) - (7168 - minx)
                       x_2 = x_centroid(1) + minx
                       y_1 = y_centroid(1) - miny_m
                       y_2 = y_centroid(1) + miny_p
                   end 
                   
                    
                if(x_1 < 1)
                    x_1 = 1
                end    
                
                if(x_2 > 10752)
                    x_2 = 10752
                end
                
                if(y_1 < 1)
                    y_1 = 1
                end
                if(y_2 < 1)
                    y_2 = 1
                end   
                if(x_2-x_1>7167)
                    x_1 = x_2-7167;
                end
                if(y_2-y_1>7167)
                   y_1 = y_2-7167; 
                end   
                if(x_2-x_1>7167)
                   x_1 = x_2-7167; 
                end
               
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
                
               imgnucleus = imread([foldername '/' allfiles(i).name]);
               %imwrite(imgnucleus,[foldername2 '/' allfiles(i).name]);
                %The rectangular image is cutted from the imgnucleus image containing the whole cover slip
               imagesphere1 = imgnucleus(y_1:y_2,x_1:x_2);
               % Dimensions of imagesphere1 are determined
               [r1 c1] = size(imagesphere1);
               % An empty matrix is generated
               imagesphereresize1 = zeros(7168,7168);
               % imagesphereresize1 is converted to 8-bit
               imagesphereresize1 = uint8(imagesphereresize1);
               %The empty indices in the imagesphereresize1 matrix are overwritten with imagesphere 1 entries
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
               % imagesphereresize1 is saved to final destination
               newFile1 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize1,[foldername1 '/' newFile1]);
               wellname = (newFile1(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(size(imagesphere1));
               currentNM = csvHandler.CellPosMatrix(wellname);
               currentNM(yStart1:yEnd1,xStart1:xEnd1) = NucleusM(y_1:y_2,x_1:x_2);
               csvHandler.CellPosMatrix(wellname) = currentNM;
               imagesphereresize1_small = imresize(imagesphereresize1, optionHandler.ScalingFactor);
               newFile1s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize1_small,[foldername1 '/' newFile1s]);
               clear imagesphereresize1_small
               clear imagesphereresize1; 
               clear imagesphere1;
            
               
               % Sphere 2 with smallest value of x_centroid
               
               if(l_2_1 >0 && l_2_1 < 2);
                       x_3 = 1;
                       x_4 = x_centroid(2) + d_2_c;
                       y_3 = 1;
                       y_4 = y_centroid(2) + d_2_c;
               else
               end    
               
               if(l_2_2 >0 && l_2_2 < 2);
                       x_3 = 1;
                       x_4 = x_centroid(2) + d_2_c;
                       y_3 = y_centroid(2) - d_2_c;
                       y_4 = 13824;
               else
               end
               
               if(r_2_4 >0 && r_2_4 < 2);
                       x_3 = x_centroid(2) - d_2_c;
                       x_4 = 10752;
                       y_3 = 1;
                       y_4 = y_centroid(2) + d_2_c;
               else
               end
               
               if(r_2_5 >0 && r_2_5 < 2);
                       x_3 = x_centroid(2) - d_2_c;
                       x_4 = 10752;
                       y_3 = y_centroid(2) - d_2_c;
                       y_4 = 13824;
               else
               end
               
               % Going over all possible distributions of sphere cores
               if(x_centroid(2) ~= x_centroid_c);
               else
                   % 1st permutation
                   if(x_centroid(1) < x_centroid(2) && x_centroid(3) < x_centroid(2) && x_centroid(4) < x_centroid(2));
                       V_1m = [d1_2 ; d2_3 ; d2_4 ; (7168/2)];
                       min1m = min(V_1m);
                       min1p = (7168-min1m);
                   else
                   end
                   % 2nd permutation
                   if(x_centroid(1) > x_centroid(2) && x_centroid(3) > x_centroid(2) && x_centroid(4) > x_centroid(2));
                       V_1p = [d1_2 ; d2_3 ; d2_4 (7168/2)];
                       min1p = min(V_1p);
                       min1m = (7168-min1p);
                   else
                   end
                   %3rd permutation
                   if(x_centroid(1) < x_centroid(2) && x_centroid(3) < x_centroid(2) && x_centroid(4) > x_centroid(2));
                       V_1m = [d1_2 ; d2_3 ; (7168/2)];
                       V_1p = [d2_4 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                       
                   else
                       
                   end
                   %4th permutation
                   if(x_centroid(1) < x_centroid(2) && x_centroid(3) > x_centroid(2) && x_centroid(4) < x_centroid(2));
                       V_1m = [d1_2 ; d2_4 ; (7168/2)];
                       V_1p = [d2_3 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   %5th permutation
                   if(x_centroid(1) > x_centroid(2) && x_centroid(3) < x_centroid(2) && x_centroid(4) < x_centroid(2));
                       V_1m = [d2_3 ; d2_4 ; (7168/2)];
                       V_1p = [d1_2 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   %6th permutation
                   if(x_centroid(1) < x_centroid(2) && x_centroid(3) > x_centroid(2) && x_centroid(4) > x_centroid(2));
                       V_1m = [d1_2 ; (7168/2)];
                       V_1p = [d2_3 ; d2_4 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   %7th permutation
                   if(x_centroid(1) > x_centroid(2) && x_centroid(3) > x_centroid(2) && x_centroid(4) < x_centroid(2));
                       V_1m = [d2_4 ; (7168/2)];
                       V_1p = [d1_2 ; d2_3 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   %8th permutation
                   if(x_centroid(1) > x_centroid(2) && x_centroid(3) < x_centroid(2) && x_centroid(4) > x_centroid(2));
                       V_1m = [d2_3 ; (7168/2)];
                       V_1p = [d1_2 ; d2_4 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   x_3 = x_centroid(2) - min1m 
                   x_4 = x_centroid(2) + min1p
                   y_3 = y_centroid(2) - min1m
                   y_4 = y_centroid(2) + min1p
                   
                   
                   
                   
                   
                   
               end    
               
               % Case for all four spheres in the same spot
                   
                   if(l_2_1_q == 1)
                       Vdq_x = [d1_2 ; d2_3 ; d2_4 ; 7168/2]
                       Vdq_y_m = [d1_2 ; d2_3 ; d2_4 ; 7168/2 ; (y_centroid(2)-1)]
                       Vdq_y_p = [d1_2 ; d2_3 ; d2_4 ; 7168/2 ]
                       miny_m = min(Vdq_y_m(1:5))
                       miny_p = min(Vdq_y_p(1:4))
                       minx = min(Vdq_x(1:4))
                       x_3 = x_centroid(2) - minx 
                       x_4 = x_centroid(2) + minx
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end   
                       
                   if(l_2_2_q == 1)
                       Vdq_x = [d1_2 ; d2_3 ; d2_4 ; 7168/2]
                       Vdq_y_m = [d1_2 ; d2_3 ; d2_4 ; 7168/2]
                       Vdq_y_p = [d1_2 ; d2_3 ; d2_4 ; 7168/2 ; (13824 - y_centroid(2))]
                       miny_m = min(Vdq_y_m(1:4))
                       miny_p = min(Vdq_y_p(1:5))
                       minx = min(Vdq_x(1:4))
                       x_3 = x_centroid(2) - minx
                       x_4 = x_centroid(2) + minx
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end 
                   
                   if(r_2_4_q == 1)
                       Vdq_x = [d1_2 ; d2_3 ; d2_4 ; 7168/2]
                       Vdq_y_m = [d1_2 ; d2_3 ; d2_4 ; 7168/2 ; (y_centroid(2)-1)]
                       Vdq_y_p = [d1_2 ; d2_3 ; d2_4 ; 7168/2 ]
                       miny_m = min(Vdq_y_m(1:5))
                       miny_p = min(Vdq_y_p(1:4))
                       minx = min(Vdq_x(1:4))
                       x_3 = x_centroid(2) - minx
                       x_4 = x_centroid(2) + minx
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end 
                   
                   if(r_2_5_q == 1)
                       Vdq_x = [d1_2 ; d2_3 ; d2_4 ; 7168/2]
                        Vdq_y_m = [d1_2 ; d2_3 ; d2_4 ; 7168/2]
                       Vdq_y_p = [d1_2 ; d2_3 ; d2_4 ; 7168/2 ; (13824 - y_centroid(2))]
                       miny_m = min(Vdq_y_m(1:4))
                       miny_p = min(Vdq_y_p(1:5))
                       minx = min(Vdq_x(1:4))
                       x_3 = x_centroid(2) - minx
                       x_4 = x_centroid(2) + minx
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end 
                   
                   
               % Case for 1,2,3 at the same position
                   if(l_2_1_t == 1)
                       Vdt_x = [d1_2 ; d2_3 ; 7168/2]
                       Vdt_y_m = [d1_2 ; d2_3  ; 7168/2 ; (y_centroid(2)-1)]
                       Vdt_y_p = [d1_2 ; d2_3  ; 7168/2 ]
                       miny_m = min(Vdt_y_m(1:4))
                       miny_p = min(Vdt_y_m(1:3))
                       minx = min(Vdt_x(1:3))
                       x_3 = x_centroid(2) - minx
                       x_4 = x_centroid(2) + minx
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end   
                       
                   if(l_2_2_t == 1)
                       Vdt_x = [d1_2 ; d2_3 ;  7168/2]
                       Vdt_y_m = [d1_2 ; d2_3 ;  7168/2]
                       Vdt_y_p = [d1_2 ; d2_3 ;  7168/2 ; (13824 - y_centroid(2))]
                       miny_m = min(Vdt_y_m(1:3))
                       miny_p = min(Vdt_y_p(1:4))
                       minx = min(Vdt_x(1:3))
                       x_3 = x_centroid(2) - minx
                       x_4 = x_centroid(2) + minx
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end 
                   
                   if(r_2_4_t == 1)
                       Vdt_x = [d1_2 ; d2_3 ;  7168/2]
                       Vdt_y_m = [d1_2 ; d2_3 ;  7168/2 ; (y_centroid(2)-1)]
                       Vdt_y_p = [d1_2 ; d2_3 ;  7168/2 ]
                       miny_m = min(Vdt_y_m(1:4))
                       miny_p = min(Vdt_y_p(1:3))
                       minx = min(Vdt_x(1:3))
                       x_3 = x_centroid(2) - minx
                       x_4 = x_centroid(2) + minx
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end 
                   
                   if(r_2_5_t == 1)
                       Vdt_x = [d1_2 ; d2_3 ;  7168/2]
                        Vdt_y_m = [d1_2 ; d2_3 ;  7168/2]
                       Vdt_y_p = [d1_2 ; d2_3 ;  7168/2 ; (13824 - y_centroid(2))]
                       miny_m = min(Vdt_y_m(1:3))
                       miny_p = min(Vdt_y_p(1:4))
                       minx = min(Vdt_x(1:3))
                       x_3 = x_centroid(2) - minx
                       x_4 = x_centroid(2) + minx
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end   
                   
                    % Case for 2,3,4 at the same position
                   if(l_2_1_t == 1)
                       Vdt_x = [d2_4 ; d2_3 ; 7168/2]
                       Vdt_y_m = [d2_4 ; d2_3  ; 7168/2 ; (y_centroid(2)-1)]
                       Vdt_y_p = [d2_4 ; d2_3  ; 7168/2 ]
                       miny_m = min(Vdt_y_m(1:4))
                       miny_p = min(Vdt_y_m(1:3))
                       minx = min(Vdt_x(1:3))
                       x_3 = x_centroid(2) - minx
                       x_4 = x_centroid(2) + minx
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end   
                       
                   if(l_2_2_t == 1)
                       Vdt_x = [d2_4 ; d2_3 ;  7168/2]
                       Vdt_y_m = [d2_4 ; d2_3 ;  7168/2]
                       Vdt_y_p = [d2_4 ; d2_3 ;  7168/2 ; (13824 - y_centroid(2))]
                       miny_m = min(Vdt_y_m(1:3))
                       miny_p = min(Vdt_y_p(1:4))
                       minx = min(Vdt_x(1:3))
                       x_3 = x_centroid(2) - minx
                       x_4 = x_centroid(2) + minx
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end 
                   
                   if(r_2_4_t == 1)
                       Vdt_x = [d2_4 ; d2_3 ;  7168/2]
                       Vdt_y_m = [d2_4 ; d2_3 ;  7168/2 ; (y_centroid(2)-1)]
                       Vdt_y_p = [d2_4 ; d2_3 ;  7168/2 ]
                       miny_m = min(Vdt_y_m(1:4))
                       miny_p = min(Vdt_y_p(1:3))
                       minx = min(Vdt_x(1:3))
                       x_3 = x_centroid(2) - minx
                       x_4 = x_centroid(2) + minx
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end 
                   
                   if(r_2_5_t == 1)
                       Vdt_x = [d2_4 ; d2_3 ;  7168/2]
                        Vdt_y_m = [d2_4 ; d2_3 ;  7168/2]
                       Vdt_y_p = [d2_4 ; d2_3 ;  7168/2 ; (13824 - y_centroid(2))]
                       miny_m = min(Vdt_y_m(1:3))
                       miny_p = min(Vdt_y_p(1:4))
                       minx = min(Vdt_x(1:3))
                       x_3 = x_centroid(2) - minx
                       x_4 = x_centroid(2) + minx
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end  
                   
                % Case fox 1,2 are the same quadrant
                   
                   if(l_2_1_d12 == 1)
                       Vdd_x = [d1_2 ; 7168/2]
                       Vdd_y_m = [d1_2 ; 7168/2 ; (y_centroid(2)-1)]
                       Vdd_y_p = [d1_2 ;  7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_m(1:2))
                       minx = min(Vdd_x(1:2))
                       x_3 = x_centroid(2) - minx
                       x_4 = x_centroid(2) + (7168-minx)
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end   
                       
                   if(l_2_2_d12 == 1)
                       Vdd_x = [d1_2 ;  7168/2]
                       Vdd_y_m = [d1_2  ;  7168/2]
                       Vdd_y_p = [d1_2 ;  7168/2 ; (13824 - y_centroid(2))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_3 = x_centroid(2) - minx
                       x_4 = x_centroid(2) + (7168-minx)
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end 
                   
                   if(r_2_4_d12 == 1)
                       Vdd_x = [d1_2 ;   7168/2]
                       Vdd_y_m = [d1_2 ;   7168/2 ; (y_centroid(2)-1)]
                       Vdd_y_p = [d1_2 ;   7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_p(1:2))
                       minx = min(Vdd_x(1:2))
                       x_3 = x_centroid(2) - minx
                       x_4 = 10752
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end 
                   
                   if(r_2_5_d12 == 1)
                       Vdd_x = [d1_2 ;  7168/2]
                        Vdd_y_m = [d1_2 ;  7168/2]
                       Vdd_y_p = [d1_2 ;   7168/2 ; (13824 - y_centroid(2))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdt_x(1:2))
                       x_3 = x_centroid(2) - minx
                       x_4 = 10752
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end  
                   
                   % Case fox 2,3 are the same quadrant
                   
                               
                   if(l_2_1_d23 == 1)
                       Vdd_x = [d2_3 ; 7168/2]
                       Vdd_y_m = [d2_3 ; 7168/2 ; (y_centroid(2)-1)]
                       Vdd_y_p = [d2_3 ;  7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_m(1:2))
                       minx = min(Vdd_x(1:2))
                       x_3 = x_centroid(2) - minx
                       x_4 = x_centroid(2) + (7168-minx)
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end   
                       
                   if(l_2_2_d23 == 1)
                       Vdd_x = [d2_3 ;  7168/2]
                       Vdd_y_m = [d2_3  ;  7168/2]
                       Vdd_y_p = [d2_3 ;  7168/2 ; (13824 - y_centroid(2))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_3 = x_centroid(2) - minx
                       x_4 = x_centroid(2) + (7168-minx)
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end 
                   
                   if(r_2_4_d23 == 1)
                       Vdd_x = [d2_3 ;   7168/2]
                       Vdd_y_m = [d2_3 ;   7168/2 ; (y_centroid(2)-1)]
                       Vdd_y_p = [d2_3 ;   7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_p(1:2))
                       minx = min(Vdd_x(1:2))
                       x_3 = x_centroid(2) - minx
                       x_4 = 10752
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end 
                   
                   if(r_2_5_d23 == 1)
                       Vdd_x = [d2_3 ;  7168/2]
                        Vdd_y_m = [d2_3 ;  7168/2]
                       Vdd_y_p = [d2_3 ;   7168/2 ; (13824 - y_centroid(2))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_3 = x_centroid(2) - minx
                       x_4 = 10752
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end  
                   
                   % Case fox 2,4 are the same quadrant
                   
                               
                   if(l_2_1_d24 == 1)
                       Vdd_x = [d2_4 ; 7168/2]
                       Vdd_y_m = [d2_4 ; 7168/2 ; (y_centroid(2)-1)]
                       Vdd_y_p = [d2_4 ;  7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_m(1:2))
                       minx = min(Vdd_x(1:2))
                       x_3 = x_centroid(2) - minx
                       x_4 = x_centroid(2) + (7168-minx)
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end   
                       
                   if(l_2_2_d24 == 1)
                       Vdd_x = [d2_4 ;  7168/2]
                       Vdd_y_m = [d2_4  ;  7168/2]
                       Vdd_y_p = [d2_4 ;  7168/2 ; (13824 - y_centroid(2))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_3 = x_centroid(2) - minx
                       x_4 = x_centroid(2) + (7168-minx)
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end 
                   
                   if(r_2_4_d24 == 1)
                       Vdd_x = [d2_4 ;   7168/2]
                       Vdd_y_m = [d2_4 ;   7168/2 ; (y_centroid(2)-1)]
                       Vdd_y_p = [d2_4 ;   7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_p(1:2))
                       minx = min(Vdd_x(1:2))
                       x_3 = x_centroid(2) - minx
                       x_4 = 10752
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end 
                   
                   if(r_2_5_d24 == 1)
                       Vdd_x = [d2_4 ;  7168/2]
                        Vdd_y_m = [d2_4 ;  7168/2]
                       Vdd_y_p = [d2_4 ;   7168/2 ; (13824 - y_centroid(2))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_3 = x_centroid(2) - minx
                       x_4 = 10752
                       y_3 = y_centroid(2) - miny_m
                       y_4 = y_centroid(2) + miny_p
                   end
                if(x_3 < 1)
                    x_3 = 1
                end    
                
                if(x_4 > 10752)
                    x_4 = 10752
                end
                
                if(y_3 < 1)
                    y_3 = 1
                end
                if(y_4 < 1)
                    y_4 = 1
                end 
                if(x_4-x_3>7167)
                    x_3 = x_4-7167;
                end
                if(y_4-y_3>7167)
                   y_3 = y_4-7167; 
                end 
                if(x_4-x_3>7167)
                   x_3 = x_4-7167; 
                end
                
               diffXCentroidl = x_centroid(2) - x_3;
               diffXCentroidr = x_4 - x_centroid(2);
               diffYCentroidl = y_centroid(2) - y_3;
               diffYCentroidr = y_4 - y_centroid(2);
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
                
                
                
                
                
                
                %The rectangular image is cutted from the imgnucleus image containing the whole cover slip
               imagesphere2 = imgnucleus(y_3:y_4,x_3:x_4);
               % Dimensions of imagesphere1 are determined
               [r2 c2] = size(imagesphere2);
               % An empty matrix is generated
               imagesphereresize2 = zeros(7168,7168);
               % imagesphereresize1 is converted to 8-bit
               imagesphereresize2 = uint8(imagesphereresize2);
               %The empty indices in the imagesphereresize1 matrix are overwritten with imagesphere 1 entries
               if(leftBoarder == 0)
                   if(lowerBoarder == 0)
                        
                        yStart2 = 1;
                        yEnd2 = r2;
                        xStart2 = 1;
                        xEnd2 = c2;
                   else
                        
                        yStart2 = (7169-r2);
                        yEnd2 = 7168;
                        xStart2 = 1;
                        xEnd2 = c2;
                   end 
               else
                   if(lowerBoarder == 0)
                       
                        yStart2 = 1;
                        yEnd2 = r2;
                        xStart2 = (7169 - c2);
                        xEnd2 = 7168;
                   else
                       
                        yStart2 = (7169-r2);
                        yEnd2 = 7168;
                        xStart2 = (7169 - c2);
                        xEnd2 = 7168;
                   end
               end  
               
               
               imagesphereresize2(yStart2:yEnd2,xStart2:xEnd2)=imagesphere2(1:r2,1:c2);
               % imagesphereresize1 is saved to final destination
               newFile2 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize2,[foldername1 '/' newFile2]);
               n = n + 1;
               wellname = (newFile2(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(size(imagesphere2));
               currentNM = csvHandler.CellPosMatrix(wellname);
               currentNM(yStart2:yEnd2,xStart2:xEnd2) = NucleusM(y_3:y_4,x_3:x_4);
               csvHandler.CellPosMatrix(wellname) = currentNM;
               imagesphereresize2_small = imresize(imagesphereresize2, optionHandler.ScalingFactor);
               newFile2s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize2_small,[foldername1 '/' newFile2s]);
               clear imagesphereresize2_small
               clear imagesphereresize2; 
               clear imagesphere2;
              
               % Sphere 3 with smallest value of x_centroid
              
               if(l_3_1 >0 && l_3_1 < 2);
                       x_5 = 1;
                       x_6 = x_centroid(3) + d_3_c;
                       y_5 = 1;
                       y_6 = y_centroid(3) + d_3_c;
               else
               end    
               
               if(l_3_2 >0 && l_3_2 < 2);
                       x_5 = 1;
                       x_6 = x_centroid(3) + d_3_c;
                       y_5 = y_centroid(3) - d_3_c;
                       y_6 = 13824;
               else
               end
               
               if(r_3_4 >0 && r_3_4 < 2);
                       x_5 = x_centroid(3) - d_3_c;
                       x_6 = 10752;
                       y_5 = 1;
                       y_6 = y_centroid(3) + d_3_c;
               else
               end
               
               if(r_3_5 >0 && r_3_5 < 2);
                       x_5 = x_centroid(3) - d_3_c;
                       x_6 = 10752;
                       y_5 = y_centroid(3) - d_3_c;
                       y_6 = 13824;
               else
               end
               
               % Going over all possible distributions of sphere cores
               if(x_centroid(3) ~= x_centroid_c);
               else
                   % 1st permutation
                   if(x_centroid(1) < x_centroid(3) && x_centroid(2) < x_centroid(3) && x_centroid(4) < x_centroid(3));
                       V_1m = [d1_3 ; d2_3 ; d3_4 ; (7168/2)];
                       min1m = min(V_1m);
                       min1p = (7168-min1m);
                   else
                   end
                   % 2nd permutation
                   if(x_centroid(1) > x_centroid(3) && x_centroid(2) > x_centroid(3) && x_centroid(4) > x_centroid(3));
                       V_1p = [d1_3 ; d2_3 ; d3_4 (7168/2)];
                       min1p = min(V_1p);
                       min1m = (7168-min1p);
                   else
                   end
                   %3rd permutation
                   if(x_centroid(1) < x_centroid(3) && x_centroid(2) < x_centroid(3) && x_centroid(4) > x_centroid(3));
                       V_1m = [d1_3 ; d2_3 ; (7168/2)];
                       V_1p = [d3_4 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                       
                   else
                       
                   end
                   %4th permutation
                   if(x_centroid(1) < x_centroid(3) && x_centroid(2) > x_centroid(3) && x_centroid(4) < x_centroid(3));
                       V_1m = [d1_3 ; d3_4 ; (7168/2)];
                       V_1p = [d2_3 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   %5th permutation
                   if(x_centroid(1) > x_centroid(3) && x_centroid(2) < x_centroid(3) && x_centroid(4) < x_centroid(3));
                       V_1m = [d2_3 ; d3_4 ; (7168/2)];
                       V_1p = [d1_3 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   %6th permutation
                   if(x_centroid(1) < x_centroid(3) && x_centroid(2) > x_centroid(3) && x_centroid(4) > x_centroid(3));
                       V_1m = [d1_3 ; (7168/2)];
                       V_1p = [d2_3 ; d3_4 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   %7th permutation
                   if(x_centroid(1) > x_centroid(3) && x_centroid(2) > x_centroid(3) && x_centroid(4) < x_centroid(3));
                       V_1m = [d3_4 ; (7168/2)];
                       V_1p = [d1_3 ; d2_3 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   %8th permutation
                   if(x_centroid(1) > x_centroid(3) && x_centroid(2) < x_centroid(3) && x_centroid(4) > x_centroid(3));
                       V_1m = [d2_3 ; (7168/2)];
                       V_1p = [d1_3 ; d3_4 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   x_5 = x_centroid(3) - min1m 
                   x_6 = x_centroid(3) + min1p
                   y_5 = y_centroid(3) - min1m
                   y_6 = y_centroid(3) + min1p
                   
                   
                   
               end    
               
               % Case for all four spheres in the same spot
                   
                   if(l_3_1_q == 1)
                       Vdq_x = [d1_3 ; d2_3 ; d3_4 ; 7168/2]
                       Vdq_y_m = [d1_3 ; d2_3 ; d3_4 ; 7168/2 ; (y_centroid(3)-1)]
                       Vdq_y_p = [d1_3 ; d2_3 ; d3_4 ; 7168/2 ]
                       miny_m = min(Vdq_y_m(1:5))
                       miny_p = min(Vdq_y_p(1:4))
                       minx = min(Vdq_x(1:4))
                       x_5 = x_centroid(3) - minx 
                       x_6 = x_centroid(3) + minx
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end   
                       
                   if(l_3_2_q == 1)
                       Vdq_x = [d1_3 ; d2_3 ; d3_4 ; 7168/2]
                       Vdq_y_m = [d1_3 ; d2_3 ; d3_4 ; 7168/2]
                       Vdq_y_p = [d1_3 ; d2_3 ; d3_4 ; 7168/2 ; (13824 - y_centroid(3))]
                       miny_m = min(Vdq_y_m(1:4))
                       miny_p = min(Vdq_p_m(1:5))
                       minx = min(Vdq_x(1:4))
                       x_5 = x_centroid(3) - minx
                       x_6 = x_centroid(3) + minx
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end 
                   
                   if(r_3_4_q == 1)
                       Vdq_x = [d1_3 ; d2_3 ; d3_4 ; 7168/2]
                       Vdq_y_m = [d1_3 ; d2_3 ; d3_4 ; 7168/2 ; (y_centroid(3)-1)]
                       Vdq_y_p = [d1_3 ; d2_3 ; d3_4 ; 7168/2 ]
                       miny_m = min(Vdq_y_m(1:5))
                       miny_p = min(Vdq_y_p(1:4))
                       minx = min(Vdq_x(1:4))
                       x_5 = x_centroid(3) - minx
                       x_6 = x_centroid(3) + minx
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end 
                   
                   if(r_3_5_q == 1)
                       Vdq_x = [d1_3 ; d2_3 ; d3_4 ; 7168/2]
                        Vdq_y_m = [d1_3 ; d2_3 ; d3_4 ; 7168/2]
                       Vdq_y_p = [d1_3 ; d2_3 ; d3_4 ; 7168/2 ; (13824 - y_centroid(3))]
                       miny_m = min(Vdq_y_m(1:4))
                       miny_p = min(Vdq_y_p(1:5))
                       minx = min(Vdq_x(1:4))
                       x_5 = x_centroid(3) - minx
                       x_6 = x_centroid(3) + minx
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end 
                 
                 % Case for 1,2,3 at the same position
                 
                   if(l_3_1_t == 1)
                       Vdt_x = [d1_3 ; d2_3 ; 7168/2]
                       Vdt_y_m = [d1_3 ; d2_3  ; 7168/2 ; (y_centroid(3)-1)]
                       Vdt_y_p = [d1_3 ; d2_3  ; 7168/2 ]
                       miny_m = min(Vdt_y_m(1:4))
                       miny_p = min(Vdt_y_m(1:3))
                       minx = min(Vdt_x(1:3))
                       x_5 = x_centroid(3) - minx
                       x_6 = x_centroid(3) + minx
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end   
                       
                   if(l_3_2_t == 1)
                       Vdt_x = [d1_3 ; d2_3 ;  7168/2]
                       Vdt_y_m = [d1_3 ; d2_3 ;  7168/2]
                       Vdt_y_p = [d1_3 ; d2_3 ;  7168/2 ; (13824 - y_centroid(3))]
                       miny_m = min(Vdt_y_m(1:3))
                       miny_p = min(Vdt_y_p(1:4))
                       minx = min(Vdt_x(1:3))
                       x_5 = x_centroid(3) - minx
                       x_6 = x_centroid(3) + minx
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end 
                   
                   if(r_3_4_t == 1)
                       Vdt_x = [d1_3 ; d2_3 ;  7168/2]
                       Vdt_y_m = [d1_3 ; d2_3 ;  7168/2 ; (y_centroid(3)-1)]
                       Vdt_y_p = [d1_3 ; d2_3 ;  7168/2 ]
                       miny_m = min(Vdt_y_m(1:4))
                       miny_p = min(Vdt_y_p(1:3))
                       minx = min(Vdt_x(1:3))
                       x_5 = x_centroid(3) - minx
                       x_6 = x_centroid(3) + minx
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end 
                   
                   if(r_3_5_t == 1)
                       Vdt_x = [d1_3 ; d2_3 ;  7168/2]
                        Vdt_y_m = [d1_3 ; d2_3 ;  7168/2]
                       Vdt_y_p = [d1_3 ; d2_3 ;  7168/2 ; (13824 - y_centroid(3))]
                       miny_m = min(Vdt_y_m(1:3))
                       miny_p = min(Vdt_y_p(1:4))
                       minx = min(Vdt_x(1:3))
                       x_5 = x_centroid(3) - minx
                       x_6 = x_centroid(3) + minx
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end  
                 
                 % Case for 2,3,4 at the same position
                   if(l_3_1_t == 1)
                       Vdt_x = [d3_4 ; d2_3 ; 7168/2]
                       Vdt_y_m = [d3_4 ; d2_3  ; 7168/2 ; (y_centroid(3)-1)]
                       Vdt_y_p = [d3_4 ; d2_3  ; 7168/2 ]
                       miny_m = min(Vdt_y_m(1:4))
                       miny_p = min(Vdt_y_m(1:3))
                       minx = min(Vdt_x(1:3))
                       x_5 = x_centroid(3) - minx
                       x_6 = x_centroid(3) + minx
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end   
                       
                   if(l_3_2_t == 1)
                       Vdt_x = [d3_4 ; d2_3 ;  7168/2]
                       Vdt_y_m = [d3_4 ; d2_3 ;  7168/2]
                       Vdt_y_p = [d3_4 ; d2_3 ;  7168/2 ; (13824 - y_centroid(3))]
                       miny_m = min(Vdt_y_m(1:3))
                       miny_p = min(Vdt_y_p(1:4))
                       minx = min(Vdt_x(1:3))
                       x_5 = x_centroid(3) - minx
                       x_6 = x_centroid(3) + minx
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end 
                   
                   if(r_3_4_t == 1)
                       Vdt_x = [d3_4 ; d2_3 ;  7168/2]
                       Vdt_y_m = [d3_4 ; d2_3 ;  7168/2 ; (y_centroid(3)-1)]
                       Vdt_y_p = [d3_4 ; d2_3 ;  7168/2 ]
                       miny_m = min(Vdt_y_m(1:4))
                       miny_p = min(Vdt_y_p(1:3))
                       minx = min(Vdt_x(1:3))
                       x_5 = x_centroid(3) - minx
                       x_6 = x_centroid(3) + minx
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end 
                   
                   if(r_3_5_t == 1)
                       Vdt_x = [d3_4 ; d2_3 ;  7168/2]
                        Vdt_y_m = [d3_4 ; d2_3 ;  7168/2]
                       Vdt_y_p = [d3_4 ; d2_3 ;  7168/2 ; (13824 - y_centroid(3))]
                       miny_m = min(Vdt_y_m(1:3))
                       miny_p = min(Vdt_y_p(1:4))
                       minx = min(Vdt_x(1:3))
                       x_5 = x_centroid(3) - minx
                       x_6 = x_centroid(3) + minx
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end    
               % Case fox 1,3 are the same quadrant
                   
                   if(l_3_1_d13 == 1)
                       Vdd_x = [d1_3 ; 7168/2]
                       Vdd_y_m = [d1_3 ; 7168/2 ; (y_centroid(3)-1)]
                       Vdd_y_p = [d1_3 ;  7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_m(1:2))
                       minx = min(Vdd_x(1:2))
                       x_5 = x_centroid(3) - minx
                       x_6 = x_centroid(3) + (7168 - minx)
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end   
                       
                   if(l_3_2_d13 == 1)
                       Vdd_x = [d1_3 ;  7168/2]
                       Vdd_y_m = [d1_3  ;  7168/2]
                       Vdd_y_p = [d1_3 ;  7168/2 ; (13824 - y_centroid(3))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_5 = x_centroid(3) - minx
                       x_6 = x_centroid(3) + (7168 - minx)
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end 
                   
                   if(r_3_4_d13 == 1)
                       Vdd_x = [d1_3 ;   7168/2]
                       Vdd_y_m = [d1_3 ;   7168/2 ; (y_centroid(3)-1)]
                       Vdd_y_p = [d1_3 ;   7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_p(1:2))
                       minx = min(Vdd_x(1:2))
                       x_5 = x_centroid(3) - minx
                       x_6 = 10752
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end 
                   
                   if(r_3_5_d13 == 1)
                       Vdd_x = [d1_3 ;  7168/2]
                        Vdd_y_m = [d1_3 ;  7168/2]
                       Vdd_y_p = [d1_3 ;   7168/2 ; (13824 - y_centroid(3))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_5 = x_centroid(3) - minx
                       x_6 = 10752
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end
                   
                % Case fox 2,3 are the same quadrant
                   
                   if(l_3_1_d23 == 1)
                       Vdd_x = [d2_3 ; 7168/2]
                       Vdd_y_m = [d2_3;7168/2 ; (y_centroid(3)-1)]
                       Vdd_y_p = [d2_3 ;  7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_m(1:2))
                       minx = min(Vdd_x(1:2))
                       x_5 = x_centroid(3) - minx
                       x_6 = x_centroid(3) + (7168 - minx)
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end   
                       
                   if(l_3_2_d23 == 1)
                       Vdd_x = [d2_3 ;  7168/2]
                       Vdd_y_m = [d2_3  ;  7168/2]
                       Vdd_y_p = [d2_3 ;  7168/2 ; (13824 - y_centroid(3))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_5 = x_centroid(3) - minx
                       x_6 = x_centroid(3) + (7168 - minx)
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end 
                   
                   if(r_3_4_d23 == 1)
                       Vdd_x = [d2_3 ;   7168/2]
                       Vdd_y_m = [d2_3 ;   7168/2 ; (y_centroid(3)-1)]
                       Vdd_y_p = [d2_3 ;   7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_p(1:2))
                       minx = min(Vdd_x(1:2))
                       x_5 = x_centroid(3) - minx
                       x_6 = 10752
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end 
                   
                   if(r_3_5_d23 == 1)
                       Vdd_x = [d2_3 ;  7168/2]
                        Vdd_y_m = [d2_3 ;  7168/2]
                       Vdd_y_p = [d2_3 ;   7168/2 ; (13824 - y_centroid(3))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_5 = x_centroid(3) - minx
                       x_6 = 10752
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end  
                   
                   % Case fox 3,4 are the same quadrant
                   
                   if(l_3_1_d34 == 1)
                       Vdd_x = [d3_4 ; 7168/2]
                       Vdd_y_m = [d3_4 ; 7168/2 ; (y_centroid(3)-1)]
                       Vdd_y_p = [d3_4 ;  7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_m(1:2))
                       minx = min(Vdd_x(1:2))
                       x_5 = x_centroid(3) - minx
                       x_6 = x_centroid(3) + (7168 - minx)
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end   
                       
                   if(l_3_2_d34 == 1)
                       Vdd_x = [d3_4 ;  7168/2]
                       Vdd_y_m = [d3_4  ;  7168/2]
                       Vdd_y_p = [d3_4 ;  7168/2 ; (13824 - y_centroid(3))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_5 = x_centroid(3) - minx
                       x_6 = x_centroid(3) + (7168 - minx)
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end 
                   
                   if(r_3_4_d34 == 1)
                       Vdd_x = [d3_4 ;   7168/2]
                       Vdd_y_m = [d3_4 ;   7168/2 ; (y_centroid(3)-1)]
                       Vdd_y_p = [d3_4 ;   7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_p(1:2))
                       minx = min(Vdd_x(1:2))
                       x_5 = x_centroid(3) - minx
                       x_6 = 10752
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end 
                   
                   if(r_3_5_d34 == 1)
                       Vdd_x = [d3_4 ;  7168/2]
                        Vdd_y_m = [d3_4 ;  7168/2]
                       Vdd_y_p = [d3_4 ;   7168/2 ; (13824 - y_centroid(3))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_5 = x_centroid(3) - minx
                       x_6 = 10752
                       y_5 = y_centroid(3) - miny_m
                       y_6 = y_centroid(3) + miny_p
                   end 
                if(x_5 < 1)
                    x_5 = 1
                end    
                
                if(x_6 > 10752)
                    x_6 = 10752
                end
                
                if(y_5 < 1)
                    y_5 = 1
                end
                if(y_6 < 1)
                    y_6 = 1
                end
                if(x_6-x_5>7167)
                    x_5 = x_6-7167;
                end
                if(y_6-y_5>7167)
                   y_5 = y_6-7167; 
                end 
                if(x_6-x_5>7167)
                   x_5 = x_6-7167; 
                end
                
               diffXCentroidl = x_centroid(3) - x_5;
               diffXCentroidr = x_6 - x_centroid(3);
               diffYCentroidl = y_centroid(3) - y_5;
               diffYCentroidr = y_6 - y_centroid(3);
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
                
                
                
                %The rectangular image is cutted from the imgnucleus image containing the whole cover slip
               imagesphere3 = imgnucleus(y_5:y_6,x_5:x_6);
               % Dimensions of imagesphere1 are determined
               [r3 c3] = size(imagesphere3);
               % An empty matrix is generated
               imagesphereresize3 = zeros(7168,7168);
               % imagesphereresize1 is converted to 8-bit
               imagesphereresize3 = uint8(imagesphereresize3);
               %The empty indices in the imagesphereresize1 matrix are overwritten with imagesphere 1 entries
               if(leftBoarder == 0)
                   if(lowerBoarder == 0)
                        
                        yStart3 = 1;
                        yEnd3 = r3;
                        xStart3 = 1;
                        xEnd3 = c3;
                   else
                        
                        yStart3 = (7169-r3);
                        yEnd3 = 7168;
                        xStart3 = 1;
                        xEnd3 = c3;
                   end 
               else
                   if(lowerBoarder == 0)
                       
                        yStart3 = 1;
                        yEnd3 = r3;
                        xStart3 = (7169 - c3);
                        xEnd3 = 7168;
                   else
                       
                        yStart3 = (7169-r3);
                        yEnd3 = 7168;
                        xStart3 = (7169 - c3);
                        xEnd3 = 7168;
                   end
               end  
               imagesphereresize3(yStart3:yEnd3,xStart3:xEnd3)=imagesphere3(1:r3,1:c3);
               % imagesphereresize1 is saved to final destination
               newFile3 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize3,[foldername1 '/' newFile3]);
               n = n + 1;
               wellname = (newFile3(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(size(imagesphere3));
               currentNM = csvHandler.CellPosMatrix(wellname);
               currentNM(yStart3:yEnd3,xStart3:xEnd3) = NucleusM(y_5:y_6,x_5:x_6);
               csvHandler.CellPosMatrix(wellname) = currentNM;
               imagesphereresize3_small = imresize(imagesphereresize3, optionHandler.ScalingFactor);
               newFile3s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize3_small,[foldername1 '/' newFile3s]);
               clear imagesphereresize3_small
               clear imagesphereresize3; 
               clear imagesphere3;
               
               
               % Sphere 4 with smallest value of x_centroid
              
               if(l_4_1 >0 && l_4_1 < 2);
                       x_7 = 1;
                       x_8 = x_centroid(4) + d_4_c;
                       y_7 = 1;
                       y_8 = y_centroid(4) + d_4_c;
               else
               end    
               
               if(l_4_2 >0 && l_4_2 < 2);
                       x_7 = 1;
                       x_8 = x_centroid(4) + d_4_c;
                       y_7 = y_centroid(4) - d_4_c;
                       y_8 = 13824;
               else
               end
               
               if(r_4_4 >0 && r_4_4 < 2);
                       x_7 = x_centroid(4) - d_4_c;
                       x_8 = 10752;
                       y_7 = 1;
                       y_8 = y_centroid(4) + d_4_c;
               else
               end
               
               if(r_4_5 >0 && r_4_5 < 2);
                       x_7 = x_centroid(4) - d_4_c;
                       x_8 = 10752;
                       y_7 = y_centroid(4) - d_4_c;
                       y_8 = 13824;
               else
               end
               
               % Going over all possible distributions of sphere cores
               if(x_centroid(4) ~= x_centroid_c);
               else
                   % 1st permutation
                   if(x_centroid(1) < x_centroid(4) && x_centroid(2) < x_centroid(4) && x_centroid(3) < x_centroid(4));
                       V_1m = [d1_4 ; d2_4 ; d3_4 ; (7168/2)];
                       min1m = min(V_1m);
                       min1p = (7168-min1m);
                   else
                   end
                   % 2nd permutation
                   if(x_centroid(1) > x_centroid(4) && x_centroid(2) > x_centroid(4) && x_centroid(3) > x_centroid(4));
                       V_1p = [d1_4 ; d2_4 ; d3_4 (7168/2)];
                       min1p = min(V_1p);
                       min1m = (7168-min1p);
                   else
                   end
                   %3rd permutation
                   if(x_centroid(1) < x_centroid(4) && x_centroid(2) < x_centroid(4) && x_centroid(3) > x_centroid(4));
                       V_1m = [d1_4 ; d2_4 ; (7168/2)];
                       V_1p = [d3_4 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                       
                   else
                       
                   end
                   %4th permutation
                   if(x_centroid(1) < x_centroid(4) && x_centroid(2) > x_centroid(4) && x_centroid(3) < x_centroid(4));
                       V_1m = [d1_4 ; d3_4 ; (7168/2)];
                       V_1p = [d2_4 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   %5th permutation
                   if(x_centroid(1) > x_centroid(4) && x_centroid(2) < x_centroid(4) && x_centroid(3) < x_centroid(4));
                       V_1m = [d2_4 ; d3_4 ; (7168/2)];
                       V_1p = [d1_4 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   %6th permutation
                   if(x_centroid(1) < x_centroid(4) && x_centroid(2) > x_centroid(4) && x_centroid(3) > x_centroid(4));
                       V_1m = [d1_4 ; (7168/2)];
                       V_1p = [d2_4 ; d3_4 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   %7th permutation
                   if(x_centroid(1) > x_centroid(4) && x_centroid(2) > x_centroid(4) && x_centroid(3) < x_centroid(4));
                       V_1m = [d3_4 ; (7168/2)];
                       V_1p = [d1_4 ; d2_4 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   %8th permutation
                   if(x_centroid(1) > x_centroid(4) && x_centroid(2) < x_centroid(4) && x_centroid(3) > x_centroid(4));
                       V_1m = [d2_4 ; (7168/2)];
                       V_1p = [d1_4 ; d3_4 ; (7168/2)];
                       min1p = min(V_1p);
                       min1m = min(V_1m);
                   else
                   end
                   x_7 = x_centroid(4) - min1m 
                   x_8 = x_centroid(4) + min1p
                   y_7 = y_centroid(4) - min1m
                   y_8 = y_centroid(4) + min1p
                   
                   
                  
                   
                   
               end    
               
                % Case for all four spheres in the same spot
                   
                   if(l_4_1_q == 1)
                       Vdq_x = [d1_4 ; d2_4 ; d3_4 ; 7168/2]
                       Vdq_y_m = [d1_4 ; d2_4 ; d3_4 ; 7168/2 ; (y_centroid(4)-1)]
                       Vdq_y_p = [d1_4 ; d2_4 ; d3_4 ; 7168/2 ]
                       miny_m = min(Vdq_y_m(1:5))
                       miny_p = min(Vdq_y_p(1:4))
                       minx = min(Vdq_x(1:4))
                       x_7 = x_centroid(4) - minx 
                       x_8 = x_centroid(4) + minx
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end   
                       
                   if(l_4_2_q == 1)
                       Vdq_x = [d1_4 ; d2_4 ; d3_4 ; 7168/2]
                       Vdq_y_m = [d1_4 ; d2_4 ; d3_4 ; 7168/2]
                       Vdq_y_p = [d1_4 ; d2_4 ; d3_4 ; 7168/2 ; (13824 - y_centroid(4))]
                       miny_m = min(Vdq_y_m(1:4))
                       miny_p = min(Vdq_y_p(1:5))
                       minx = min(Vdq_x(1:4))
                       x_7 = x_centroid(4) - minx
                       x_8 = x_centroid(4) + minx
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end 
                   
                   if(r_4_4_q == 1)
                       Vdq_x = [d1_4 ; d2_4 ; d3_4 ; 7168/2]
                       Vdq_y_m = [d1_4 ; d2_4 ; d3_4 ; 7168/2 ; (y_centroid(4)-1)]
                       Vdq_y_p = [d1_4 ; d2_4 ; d3_4 ; 7168/2 ]
                       miny_m = min(Vdq_y_m(1:5))
                       miny_p = min(Vdq_y_p(1:4))
                       minx = min(Vdq_x(1:4))
                       x_7 = x_centroid(4) - minx
                       x_8 = x_centroid(4) + minx
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end 
                   
                   if(r_4_5_q == 1)
                       Vdq_x = [d1_4 ; d2_4 ; d3_4 ; 7168/2]
                        Vdq_y_m = [d1_4 ; d2_4 ; d3_4 ; 7168/2]
                       Vdq_y_p = [d1_4 ; d2_4 ; d3_4 ; 7168/2 ; (13824 - y_centroid(4))]
                       miny_m = min(Vdq_y_m(1:4))
                       miny_p = min(Vdq_y_p(1:5))
                       minx = min(Vdq_x(1:4))
                       x_7 = x_centroid(4) - minx
                       x_8 = x_centroid(4) + minx
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end 
                   
                   
                 % Case for 2,3,4 at the same position
                   if(l_4_1_t == 1)
                       Vdt_x = [d3_4 ; d2_4 ; 7168/2]
                       Vdt_y_m = [d3_4 ; d2_4  ; 7168/2 ; (y_centroid(4)-1)]
                       Vdt_y_p = [d3_4 ; d2_4  ; 7168/2 ]
                       miny_m = min(Vdt_y_m(1:4))
                       miny_p = min(Vdt_y_m(1:3))
                       minx = min(Vdt_x(1:3))
                       x_7 = x_centroid(4) - minx
                       x_8 = x_centroid(4) + minx
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end   
                       
                   if(l_4_2_t == 1)
                       Vdt_x = [d3_4 ; d2_4 ;  7168/2]
                       Vdt_y_m = [d3_4 ; d2_4 ;  7168/2]
                       Vdt_y_p = [d3_4 ; d2_4 ;  7168/2 ; (13824 - y_centroid(4))]
                       miny_m = min(Vdt_y_m(1:3))
                       miny_p = min(Vdt_y_p(1:4))
                       minx = min(Vdt_x(1:3))
                       x_7 = x_centroid(4) - minx
                       x_8 = x_centroid(4) + minx
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end 
                   
                   if(r_4_4_t == 1)
                       Vdt_x = [d3_4 ; d2_4 ;  7168/2]
                       Vdt_y_m = [d3_4 ; d2_4 ;  7168/2 ; (y_centroid(4)-1)]
                       Vdt_y_p = [d3_4 ; d2_4 ;  7168/2 ]
                       miny_m = min(Vdt_y_m(1:4))
                       miny_p = min(Vdt_y_p(1:3))
                       minx = min(Vdt_x(1:3))
                       x_7 = x_centroid(4) - minx
                       x_8 = x_centroid(4) + minx
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end 
                   
                   if(r_4_5_t == 1)
                       Vdt_x = [d3_4 ; d2_4 ;  7168/2]
                        Vdt_y_m = [d3_4 ; d2_4 ;  7168/2]
                       Vdt_y_p = [d3_4 ; d2_4 ;  7168/2 ; (13824 - y_centroid(4))]
                       miny_m = min(Vdt_y_m(1:3))
                       miny_p = min(Vdt_y_p(1:4))
                       minx = min(Vdt_x(1:3))
                       x_7 = x_centroid(4) - minx
                       x_8 = x_centroid(4) + minx
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end   
                    % Case fox 1,4 are the same quadrant
                   
                   if(l_4_1_d14 == 1)
                       Vdd_x = [d1_4 ; 7168/2]
                       Vdd_y_m = [d1_4 ; 7168/2 ; (y_centroid(4)-1)]
                       Vdd_y_p = [d1_4 ;  7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_m(1:2))
                       minx = min(Vdd_x(1:2))
                       x_7 = x_centroid(4) - minx
                       x_8 = x_centroid(4) + (7168 - minx)
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end   
                       
                   if(l_4_2_d14 == 1)
                       Vdd_x = [d1_4 ;  7168/2]
                       Vdd_y_m = [d1_4  ;  7168/2]
                       Vdd_y_p = [d1_4 ;  7168/2 ; (13824 - y_centroid(4))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_7 = x_centroid(4) - minx
                       x_8 = x_centroid(4) + (7168 - minx)
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end 
                   
                   if(r_4_4_d14 == 1)
                       Vdd_x = [d1_4 ;   7168/2]
                       Vdd_y_m = [d1_4 ;   7168/2 ; (y_centroid(4)-1)]
                       Vdd_y_p = [d1_4 ;   7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_p(1:2))
                       minx = min(Vdd_x(1:2))
                       x_7 = x_centroid(4) - minx
                       x_8 = 10752
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end 
                   
                   if(r_4_5_d14 == 1)
                       Vdd_x = [d1_4 ;  7168/2]
                        Vdd_y_m = [d1_4 ;  7168/2]
                       Vdd_y_p = [d1_4 ;   7168/2 ; (13824 - y_centroid(4))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_7 = x_centroid(4) - minx
                       x_8 = 10752
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end 
                   
                    % Case fox 2,4 are the same quadrant
                   
                   if(l_4_1_d24 == 1)
                       Vdd_x = [d2_4 ; 7168/2]
                       Vdd_y_m = [d2_4 ; 7168/2 ; (y_centroid(4)-1)]
                       Vdd_y_p = [d2_4 ;  7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_m(1:2))
                       minx = min(Vdd_x(1:2))
                       x_7 = x_centroid(4) - minx
                       x_8 = x_centroid(4) + (7168 - minx)
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end   
                       
                   if(l_4_2_d24 == 1)
                       Vdd_x = [d2_4 ;  7168/2]
                       Vdd_y_m = [d2_4  ;  7168/2]
                       Vdd_y_p = [d2_4 ;  7168/2 ; (13824 - y_centroid(4))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_7 = x_centroid(4) - minx
                       x_8 = x_centroid(4) + (7168 - minx)
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end 
                   
                   if(r_4_4_d24 == 1)
                       Vdd_x = [d2_4 ;   7168/2]
                       Vdd_y_m = [d2_4 ;   7168/2 ; (y_centroid(4)-1)]
                       Vdd_y_p = [d2_4 ;   7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_p(1:2))
                       minx = min(Vdd_x(1:2))
                       x_7 = x_centroid(4) - minx
                       x_8 = 10752
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end 
                   
                   if(r_4_5_d24 == 1)
                       Vdd_x = [d2_4 ;  7168/2]
                        Vdd_y_m = [d2_4 ;  7168/2]
                       Vdd_y_p = [d2_4 ;   7168/2 ; (13824 - y_centroid(4))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_7 = x_centroid(4) - minx
                       x_8 = 10752
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end 
                   
                   % Case fox 3,4 are the same quadrant
                   
                   if(l_4_1_d34 == 1)
                       Vdd_x = [d3_4 ; 7168/2]
                       Vdd_y_m = [d3_4 ; 7168/2 ; (y_centroid(4)-1)]
                       Vdd_y_p = [d3_4 ;  7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_m(1:2))
                       minx = min(Vdd_x(1:2))
                       x_7 = x_centroid(4) - minx
                       x_8 = x_centroid(4) + (7168 - minx)
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end   
                       
                   if(l_4_2_d34 == 1)
                       Vdd_x = [d3_4 ;  7168/2]
                       Vdd_y_m = [d3_4  ;  7168/2]
                       Vdd_y_p = [d3_4 ;  7168/2 ; (13824 - y_centroid(4))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_7 = x_centroid(4) - minx
                       x_8 = x_centroid(4) + (7168 - minx)
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end 
                   
                   if(r_4_4_d34 == 1)
                       Vdd_x = [d3_4 ;   7168/2]
                       Vdd_y_m = [d3_4 ;   7168/2 ; (y_centroid(4)-1)]
                       Vdd_y_p = [d3_4 ;   7168/2 ]
                       miny_m = min(Vdd_y_m(1:3))
                       miny_p = min(Vdd_y_p(1:2))
                       minx = min(Vdd_x(1:2))
                       x_7 = x_centroid(4) - minx
                       x_8 = 10752
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end 
                   
                   if(r_4_5_d34 == 1)
                       Vdd_x = [d3_4 ;  7168/2]
                        Vdd_y_m = [d3_4 ;  7168/2]
                       Vdd_y_p = [d3_4 ;   7168/2 ; (13824 - y_centroid(4))]
                       miny_m = min(Vdd_y_m(1:2))
                       miny_p = min(Vdd_y_p(1:3))
                       minx = min(Vdd_x(1:2))
                       x_7 = x_centroid(4) - minx
                       x_8 = 10752
                       y_7 = y_centroid(4) - miny_m
                       y_8 = y_centroid(4) + miny_p
                   end 
                   
                if(x_7 < 1)
                    x_7 = 1
                end    
                
                if(x_8 > 10752)
                    x_8 = 10752
                end
                
                if(y_7 < 1)
                    y_7 = 1
                end
                if(y_8 < 1)
                    y_8 = 1
                end
                if(x_8-x_7>7167)
                    x_7 = x_8-7167;
                end
                if(y_8-y_7>7167)
                   y_7 = y_8-7167; 
                end 
                if(x_8-x_7>7167)
                   x_7 = x_8-7167; 
                end
                
               diffXCentroidl = x_centroid(4) - x_7;
               diffXCentroidr = x_8 - x_centroid(4);
               diffYCentroidl = y_centroid(4) - y_7;
               diffYCentroidr = y_8 - y_centroid(4);
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
                
                
                %The rectangular image is cutted from the imgnucleus image containing the whole cover slip
               imagesphere4 = imgnucleus(y_7:y_8,x_7:x_8);
               % Dimensions of imagesphere1 are determined
               [r4 c4] = size(imagesphere4);
               % An empty matrix is generated
               imagesphereresize4 = zeros(7168,7168);
               % imagesphereresize1 is converted to 8-bit
               imagesphereresize4 = uint8(imagesphereresize4);
               %The empty indices in the imagesphereresize1 matrix are overwritten with imagesphere 1 entries
               if(leftBoarder == 0)
                   if(lowerBoarder == 0)
                        
                        yStart4 = 1;
                        yEnd4 = r4;
                        xStart4 = 1;
                        xEnd4 = c4;
                   else
                        
                        yStart4 = (7169-r4);
                        yEnd4 = 7168;
                        xStart4 = 1;
                        xEnd4 = c4;
                   end 
               else
                   if(lowerBoarder == 0)
                       
                        yStart4 = 1;
                        yEnd4 = r4;
                        xStart4 = (7169 - c4);
                        xEnd4 = 7168;
                   else
                       
                        yStart4 = (7169-r4);
                        yEnd4 = 7168;
                        xStart4 = (7169 - c4);
                        xEnd4 = 7168;
                   end
               end  
               
               
               imagesphereresize4(yStart4:yEnd4,xStart4:xEnd4)=imagesphere4(1:r4,1:c4);
               % imagesphereresize1 is saved to final destination
               newFile4 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize4,[foldername1 '/' newFile4]);
               n = n + 1;
               wellname = (newFile4(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(size(imagesphere4));
               currentNM = csvHandler.CellPosMatrix(wellname);
               currentNM(yStart4:yEnd4,xStart4:xEnd4) = NucleusM(y_7:y_8,x_7:x_8);
               csvHandler.CellPosMatrix(wellname) = currentNM;
               imagesphereresize4_small = imresize(imagesphereresize4, optionHandler.ScalingFactor);
               newFile4s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize4_small,[foldername1 '/' newFile4s]);
               clear imagesphereresize4_small
               clear imagesphereresize4; 
               clear imagesphere4;
               
               imagesphereresize5 = zeros(7168,7168);
               newFile5 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize5,[foldername1 '/' newFile5]);
               newFile5s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NucleusSmall'],'ignorecase');  
               imagesphereresize5_small = imresize(imagesphereresize5, optionHandler.ScalingFactor);
               clear imagesphereresize5; 
               imwrite(imagesphereresize5_small,[foldername1 '/' newFile5s]);
               n = n + 1;
               wellname = (newFile5(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(zeros(7168,7168));
               clear imagesphereresize5_small
               clear imgnucleus;
               
               if(Neuron_Channel > 0)
                   % Same for neurite channel


                   imgneurite = imread([foldername '/' newFileNeu]);
                   %imwrite(imgneurite,[foldername2 '/' newFileNeu]);
                   imagesphere6 = imgneurite(y_1:y_2,x_1:x_2);
                   [r1 c1] = size(imagesphere6);
                   imagesphereresize6 = zeros(7168,7168);
                   imagesphereresize6 = uint8(imagesphereresize6);
                   imagesphereresize6(yStart1:yEnd1,xStart1:xEnd1)=imagesphere6(1:r1,1:c1);
                   newFile6 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize6,[foldername1 '/' newFile6]);
                   imagesphereresize6_small = imresize(imagesphereresize6, optionHandler.ScalingFactor);
                   newFile6s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize6_small,[foldername1 '/' newFile6s]);
                   clear imagespher6;
                   clear imagesphereresize6; 
                   clear imagesphereresize6_small

                   imagesphere7 = imgneurite(y_3:y_4,x_3:x_4);
                   [r2 c2] = size(imagesphere7);
                   imagesphereresize7 = zeros(7168,7168);
                   imagesphereresize7 = uint8(imagesphereresize7);
                   imagesphereresize7(yStart2:yEnd2,xStart2:xEnd2)=imagesphere7(1:r2,1:c2);
                   newFile7 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize7,[foldername1 '/' newFile7]);
                   imagesphereresize7_small = imresize(imagesphereresize7, optionHandler.ScalingFactor);
                   newFile7s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize7_small,[foldername1 '/' newFile7s]);
                   clear imagesphere7;
                   clear imagesphereresize7; 
                   clear imagesphereresize7_small

                   imagesphere8 = imgneurite(y_5:y_6,x_5:x_6);
                   [r3 c3] = size(imagesphere8);
                   imagesphereresize8 = zeros(7168,7168);
                   imagesphereresize8 = uint8(imagesphereresize8);
                   imagesphereresize8(yStart3:yEnd3,xStart3:xEnd3)=imagesphere8(1:r3,1:c3);
                   newFile8 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize8,[foldername1 '/' newFile8]);
                   imagesphereresize8_small = imresize(imagesphereresize8, optionHandler.ScalingFactor);
                   newFile8s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize8_small,[foldername1 '/' newFile8s]);
                   clear imagesphere8;
                   clear imagesphereresize8; 
                   clear imagesphereresize8_small

                   imagesphere9 = imgneurite(y_7:y_8,x_7:x_8);
                   [r4 c4] = size(imagesphere9);
                   imagesphereresize9 = zeros(7168,7168);
                   imagesphereresize9 = uint8(imagesphereresize9);
                   imagesphereresize9(yStart4:yEnd4,xStart4:xEnd4)=imagesphere9(1:r4,1:c4);
                   newFile9 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize9,[foldername1 '/' newFile9]);
                   imagesphereresize9_small = imresize(imagesphereresize9, optionHandler.ScalingFactor);
                   newFile9s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize9_small,[foldername1 '/' newFile9s]);
                   clear imagesphere9;
                   clear imagesphereresize9; 
                   clear imagesphereresize9_small

                   imagesphereresize10 = zeros(7168,7168);
            	   newFile10 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize10,[foldername1 '/' newFile10]);
                   imagesphereresize10_small = imresize(imagesphereresize10, optionHandler.ScalingFactor);
                   clear imagesphereresize10; 
                   newFile10s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize10_small,[foldername1 '/' newFile10s]);
                   clear imagesphereresize10_small
                   clear imgneurite
               end
               
               if(Oligo_Channel > 0)
                   % Same for Oligos
                   
                   imgoligo = imread([foldername '/' newFileOli]);
                   %imwrite(imgoligo,[foldername2 '/' newFileOli]);
                   imagesphere11 = imgoligo(y_1:y_2,x_1:x_2);
                   [r1 c1] = size(imagesphere11);
                   imagesphereresize11 = zeros(7168,7168);
                   imagesphereresize11 = uint8(imagesphereresize11);
                   imagesphereresize11(yStart1:yEnd1,xStart1:xEnd1)=imagesphere11(1:r1,1:c1);
                   newFile11 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize11,[foldername1 '/' newFile11]);
                   imagesphereresize11_small = imresize(imagesphereresize11, optionHandler.ScalingFactor);
                   newFile11s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize11_small,[foldername1 '/' newFile11s]);
                   clear imagesphere11;
                   clear imagesphereresize11; 
                   clear imagesphereresize11_small

                   imagesphere12 = imgoligo(y_3:y_4,x_3:x_4);
                   [r2 c2] = size(imagesphere12);
                   imagesphereresize12 = zeros(7168,7168);
                   imagesphereresize12 = uint8(imagesphereresize12);
                   imagesphereresize12(yStart2:yEnd2,xStart2:xEnd2)=imagesphere12(1:r2,1:c2);
                   newFile12 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize12,[foldername1 '/' newFile12]);
                   imagesphereresize12_small = imresize(imagesphereresize12, optionHandler.ScalingFactor);
                   newFile12s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize12_small,[foldername1 '/' newFile12s]);
                   clear imagesphere12;
                   clear imagesphereresize12; 
                   clear imagesphereresize12_small

                   imagesphere13 = imgoligo(y_5:y_6,x_5:x_6);
                   [r3 c3] = size(imagesphere13);
                   imagesphereresize13 = zeros(7168,7168);
                   imagesphereresize13 = uint8(imagesphereresize13);
                   imagesphereresize13(yStart3:yEnd3,xStart3:xEnd3)=imagesphere13(1:r3,1:c3);
                   newFile13 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize13,[foldername1 '/' newFile13]);
                   imagesphereresize13_small = imresize(imagesphereresize13, optionHandler.ScalingFactor);
                   newFile13s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize13_small,[foldername1 '/' newFile13s]);
                   clear imagesphere13;
                   clear imagesphereresize13; 
                   clear imagesphereresize13_small


                   imagesphere14 = imgoligo(y_7:y_8,x_7:x_8);
                   [r4 c4] = size(imagesphere14);
                   imagesphereresize14 = zeros(7168,7168);
                   imagesphereresize14 = uint8(imagesphereresize14);
                   imagesphereresize14(yStart4:yEnd4,xStart4:xEnd4)=imagesphere14(1:r4,1:c4);
                   newFile14 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize14,[foldername1 '/' newFile14]);
                   imagesphereresize14_small = imresize(imagesphereresize14, optionHandler.ScalingFactor);
                   newFile14s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize14_small,[foldername1 '/' newFile14s]);
                   clear imagesphere14;
                   clear imagesphereresize14; 
                   clear imagesphereresize14_small

                   imagesphereresize15 = zeros(7168,7168);
            	   newFile15 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize15,[foldername1 '/' newFile15]);
                   imagesphereresize15_small = imresize(imagesphereresize15, optionHandler.ScalingFactor);
                   clear imagesphereresize15;
                   newFile15s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize15_small,[foldername1 '/' newFile15s]);
                   clear imagesphereresize15_small 
                   clear imgoligo;
               end
               
               if(Astro_Channel >0)
                   % Same for Astros
                   imgastro = imread([foldername '/' newFileAst]);
                   %imwrite(imgastro,[foldername2 '/' newFileAst]);
                   imagesphere16 = imgastro(y_1:y_2,x_1:x_2);
                   [r1 c1] = size(imagesphere16);
                   imagesphereresize16 = zeros(7168,7168);
                   imagesphereresize16 = uint8(imagesphereresize16);
                   imagesphereresize16(yStart1:yEnd1,xStart1:xEnd1)=imagesphere16(1:r1,1:c1);
                   newFile16 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize16,[foldername1 '/' newFile16]);
                   imagesphereresize16_small = imresize(imagesphereresize16, optionHandler.ScalingFactor);
                   newFile16s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize16_small,[foldername1 '/' newFile16s]);
                   clear imagesphere16;
                   clear imagesphereresize16; 
                   clear imagesphereresize16_small

                   imagesphere17 = imgastro(y_3:y_4,x_3:x_4);
                   [r2 c2] = size(imagesphere17);
                   imagesphereresize17 = zeros(7168,7168);
                   imagesphereresize17 = uint8(imagesphereresize17);
                   imagesphereresize17(yStart2:yEnd2,xStart2:xEnd2)=imagesphere17(1:r2,1:c2);
                   newFile17 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize17,[foldername1 '/' newFile17]);
                   imagesphereresize17_small = imresize(imagesphereresize17, optionHandler.ScalingFactor);
                   newFile17s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize17_small,[foldername1 '/' newFile17s]);
                   clear imagesphere17;
                   clear imagesphereresize17; 
                   clear imagesphereresize17_small

                   imagesphere18 = imgastro(y_5:y_6,x_5:x_6);
                   [r3 c3] = size(imagesphere18);
                   imagesphereresize18 = zeros(7168,7168);
                   imagesphereresize18 = uint8(imagesphereresize18);
                   imagesphereresize18(yStart3:yEnd3,xStart3:xEnd3)=imagesphere18(1:r3,1:c3);
                   newFile18 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize18,[foldername1 '/' newFile18]);
                   imagesphereresize18_small = imresize(imagesphereresize18, optionHandler.ScalingFactor);
                   newFile18s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize18_small,[foldername1 '/' newFile18s]);
                   clear imagesphere18;
                   clear imagesphereresize18; 
                   clear imagesphereresize18_small

                   imagesphere19 = imgastro(y_7:y_8,x_7:x_8);
                   [r4 c4] = size(imagesphere19);
                   imagesphereresize19 = zeros(7168,7168);
                   imagesphereresize19 = uint8(imagesphereresize19);
                   imagesphereresize19(yStart4:yEnd4,xStart4:xEnd4)=imagesphere19(1:r4,1:c4);
                   newFile19 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize19,[foldername1 '/' newFile19]);
                   imagesphereresize19_small = imresize(imagesphereresize19, optionHandler.ScalingFactor);
                   newFile19s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize19_small,[foldername1 '/' newFile19s]);
                   clear imagesphere19;
                   clear imagesphereresize19; 
                   clear imagesphereresize19_small
                   
                   imagesphereresize20 = zeros(7168,7168);
            	   newFile20 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize20,[foldername1 '/' newFile20]);
                   imagesphereresize20_small = imresize(imagesphereresize20, optionHandler.ScalingFactor);
                   clear imagesphereresize20;
                   newFile20s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize20_small,[foldername1 '/' newFile20s]);
                   clear imagesphereresize20_small;
                   clear imgastro
               end
                %delete(allfiles(i).name);
                %if(Neuron_Channel > 0)
                %delete(newFileNeu);
                %end
                %if(Oligo_Channel > 0)
                %delete(newFileOli);
                %end
                %if(Astro_Channel > 0)
                %delete(newFileAst);
                %end
                
            end  
            
            if(s_L >= 5);    
               % Extracts the x and y coordinates from s  
               x_centroid(1) = s(1,1);
               y_centroid(1) = s(1,2);
               x_centroid(2) = s(2,1);
               y_centroid(2) = s(2,2);
               x_centroid(3) = s(3,1);
               y_centroid(3) = s(3,2);
               x_centroid(4) = s(4,1);
               y_centroid(4) = s(4,2);
               x_centroid(5) = s(5,1);
               y_centroid(5) = s(5,2);
               % Calculates the distance between the centroids of the
               % neurosphere an all boarders of the cover slip      
               
               
               
               
               
               % Distance between centroid(1) and centroid(2) 
               d1_2 = ((x_centroid(1)-x_centroid(2))^2+(y_centroid(1)-y_centroid(2))^2)^0.5; 
               % Distance between centroid(1) and centroid(3)
               d1_3 = ((x_centroid(1)-x_centroid(3))^2+(y_centroid(1)-y_centroid(3))^2)^0.5;
               % Distance between centroid(1) and centroid(4)
               d1_4 = ((x_centroid(1)-x_centroid(4))^2+(y_centroid(1)-y_centroid(4))^2)^0.5;
               % Distance between centroid(1) and centroid(5)
               d1_5 = ((x_centroid(1)-x_centroid(5))^2+(y_centroid(1)-y_centroid(5))^2)^0.5;
               % Distance between centroid(2) and centroid(3)
               d2_3 = ((x_centroid(3)-x_centroid(2))^2+(y_centroid(3)-y_centroid(2))^2)^0.5;
               % Distance between centroid(2) and centroid(4)
               d2_4 = ((x_centroid(4)-x_centroid(2))^2+(y_centroid(4)-y_centroid(2))^2)^0.5;
               % Distance between centroid(2) and centroid(5)
               d2_5 = ((x_centroid(5)-x_centroid(2))^2+(y_centroid(5)-y_centroid(2))^2)^0.5;
               % Distance between centroid(3) and centroid(4)
               d3_4 = ((x_centroid(4)-x_centroid(3))^2+(y_centroid(4)-y_centroid(3))^2)^0.5;
               % Distance between centroid(3) and centroid(5)
               d3_5 = ((x_centroid(5)-x_centroid(3))^2+(y_centroid(5)-y_centroid(3))^2)^0.5;
               % Distance between centroid(4) and centroid(5)
               d4_5 = ((x_centroid(4)-x_centroid(5))^2+(y_centroid(4)-y_centroid(5))^2)^0.5;
               %Distances to rims for sphere 1
               
               d1_x_1 = ((x_centroid(1)-1)^2+(y_centroid(1)-y_centroid(1))^2)^0.5
               d1_x10752 = ((10752 - x_centroid(1)-1)+(y_centroid(1)-y_centroid(1))^2)^0.5
               d1_y0 = ((x_centroid(1)-x_centroid(1))^2+(y_centroid(1)-1)^2)^0.5
               d1_y13824 = ((x_centroid(1)-x_centroid(1))^2+(13824-y_centroid(1))^2)^0.5
               
               
               
         
               
               
                              
               %Sphere 1 and 2:
               if(d1_2<MinimumDistance);
                   S12 = 1;
                   % Calculate line equation for d_12 
                   m_12 = (y_centroid(2)-y_centroid(1))/(x_centroid(2)-x_centroid(1));
                   b_12 = y_centroid(1)-(m_12*x_centroid(1));
                   % 2) Calculate Coordinate of d_12/2
                   x_c_12 = x_centroid(1) + (x_centroid(2)-x_centroid(1))/2;
                   y_c_12 = m_12*x_c_12 + b_12
                   % Calculate max distance in x from sphere 2 to 1:
                   x_12_1 = x_centroid(1) - (maxaddedlength);
                   x_12_2 = x_c_12
                   %Introduce ratio of x/y (to rule out cases where two
                   %spheres have a very small differences in x or y but not
                   %in the other coordinate
                   xy_1 = abs(x_centroid(2)-x_centroid(1))/abs(y_centroid(2)-y_centroid(1))
                   yx_1 = abs(y_centroid(2)-y_centroid(1))/abs(x_centroid(2)-x_centroid(1))
                   if(xy_1 < minxy_ratio && d1_2 > mindist)
                       x_c_12 = 1000000;
                       y_c_12 = 1000000;
                       x_12_1 = x_centroid(1) - (maxaddedlength);
                       x_12_2 = x_centroid(1) + (maxaddedlength);
                       y_12_1 = y_centroid(1)-(maxaddedlength);
                       y_12_2 = y_centroid(1)+(maxaddedlength);
                   end 
                   if(yx_1 < minxy_ratio && d1_2 > mindist)
                       x_c_12 = 1000000;
                       y_c_12 = 1000000;
                       x_12_1 = x_centroid(1) - (maxaddedlength);
                       x_12_2 = x_centroid(1) + (maxaddedlength);
                       y_12_1 = y_centroid(1)-(maxaddedlength);
                       y_12_2 = y_centroid(1)+(maxaddedlength);
                   end
                   if(y_c_12 > y_centroid(1))
                       y_12_1 = y_centroid(1)-(maxaddedlength)
                       y_12_2 = y_c_12
                   else 
                       y_12_1 = y_c_12
                       y_12_2 = y_centroid(1)+(maxaddedlength)
                   end    
               else
                   S12 = 0;
                   x_c_12 = 1000000;
                   y_c_12 = 1000000;
                   x_12_1 = x_centroid(1) - (maxaddedlength);
                   x_12_2 = x_centroid(1) + (maxaddedlength);
                   y_12_1 = y_centroid(1)-(maxaddedlength);
                   y_12_2 = y_centroid(1)+(maxaddedlength);
               end
               
               %Sphere 1 and 3:
               if(d1_3<MinimumDistance);
                   S13 = 1;
                   % Calculate line equation for d_13 
                   m_13 = (y_centroid(3)-y_centroid(1))/(x_centroid(3)-x_centroid(1));
                   b_13 = y_centroid(1)-(m_13*x_centroid(1));
                   % 2) Calculate Coordinate of d_13/2
                   x_c_13 = x_centroid(1) + (x_centroid(3)-x_centroid(1))/2;
                   y_c_13 = m_13*x_c_13 + b_13
                   % Calculate max distance in x from sphere 2 to 1:
                   x_13_1 = x_centroid(1) - (maxaddedlength);
                   x_13_2 = x_c_13
                   %Introduce ratio of x/y (to rule out cases where two
                   %spheres have a very small differences in x or y but not
                   %in the other coordinate
                   xy_2 = abs(x_centroid(3)-x_centroid(1))/abs(y_centroid(3)-y_centroid(1))
                   yx_2 = abs(y_centroid(3)-y_centroid(1))/abs(x_centroid(3)-x_centroid(1))
                   if(xy_2 < minxy_ratio && d1_3 > mindist)
                       x_c_13 = 1000000;
                       y_c_13 = 1000000;
                       x_13_1 = x_centroid(1) - (maxaddedlength);
                       x_13_2 = x_centroid(1) + (maxaddedlength);
                       y_13_1 = y_centroid(1)-(maxaddedlength);
                       y_13_2 = y_centroid(1)+(maxaddedlength);
                   end 
                   if(yx_2 < minxy_ratio && d1_3 > mindist)
                       x_c_13 = 1000000;
                       y_c_13 = 1000000;
                       x_13_1 = x_centroid(1) - (maxaddedlength);
                       x_13_2 = x_centroid(1) + (maxaddedlength);
                       y_13_1 = y_centroid(1)-(maxaddedlength);
                       y_13_2 = y_centroid(1)+(maxaddedlength);
                   end
                   if(y_c_13 > y_centroid(1))
                       y_13_1 = y_centroid(1)-(maxaddedlength)
                       y_13_2 = y_c_13
                   else 
                       y_13_1 = y_c_13
                       y_13_2 = y_centroid(1)+(maxaddedlength)
                   end    
               else
                   S13 = 0;
                   x_c_13 = 1000000;
                   y_c_13 = 1000000;
                   x_13_1 = x_centroid(1) - (maxaddedlength);
                   x_13_2 = x_centroid(1) + (maxaddedlength);
                   y_13_1 = y_centroid(1)-(maxaddedlength);
                   y_13_2 = y_centroid(1)+(maxaddedlength);
               end
               
               %Sphere 1 and 4:
               if(d1_4<MinimumDistance);
                   S14 = 1;
                   % Calculate line equation for d_14 
                   m_14 = (y_centroid(4)-y_centroid(1))/(x_centroid(4)-x_centroid(1));
                   b_14 = y_centroid(1)-(m_14*x_centroid(1));
                   % 2) Calculate Coordinate of d_14/2
                   x_c_14 = x_centroid(1) + (x_centroid(4)-x_centroid(1))/2;
                   y_c_14 = m_14*x_c_14 + b_14
                   % Calculate max distance in x from sphere 2 to 1:
                   x_14_1 = x_centroid(1) - (maxaddedlength);
                   x_14_2 = x_c_14
                   %Introduce ratio of x/y (to rule out cases where two
                   %spheres have a very small differences in x or y but not
                   %in the other coordinate
                   xy_3 = abs(x_centroid(4)-x_centroid(1))/abs(y_centroid(4)-y_centroid(1))
                   yx_3 = abs(y_centroid(4)-y_centroid(1))/abs(x_centroid(4)-x_centroid(1))
                   if(xy_3 < minxy_ratio && d1_4 > mindist)
                       x_c_14 = 1000000;
                       y_c_14 = 1000000;
                       x_14_1 = x_centroid(1) - (maxaddedlength);
                       x_14_2 = x_centroid(1) + (maxaddedlength);
                       y_14_1 = y_centroid(1)-(maxaddedlength);
                       y_14_2 = y_centroid(1)+(maxaddedlength);
                   end   
                   if(yx_3 < minxy_ratio && d1_4 > mindist)
                       x_c_14 = 1000000;
                       y_c_14 = 1000000;
                       x_14_1 = x_centroid(1) - (maxaddedlength);
                       x_14_2 = x_centroid(1) + (maxaddedlength);
                       y_14_1 = y_centroid(1)-(maxaddedlength);
                       y_14_2 = y_centroid(1)+(maxaddedlength);
                   end
                   if(y_c_14 > y_centroid(1))
                       y_14_1 = y_centroid(1)-(maxaddedlength)
                       y_14_2 = y_c_14
                   else 
                       y_14_1 = y_c_14
                       y_14_2 = y_centroid(1)+(maxaddedlength)
                   end    
               else
                   S14 = 0;
                   x_c_14 = 1000000;
                   y_c_14 = 1000000;
                   x_14_1 = x_centroid(1) - (maxaddedlength);
                   x_14_2 = x_centroid(1) + (maxaddedlength);
                   y_14_1 = y_centroid(1)-(maxaddedlength);
                   y_14_2 = y_centroid(1)+(maxaddedlength);
               end
               
               %Sphere 1 and 5:
               if(d1_5<MinimumDistance);
                   S15 = 1;
                   % Calculate line equation for d_15 
                   m_15 = (y_centroid(5)-y_centroid(1))/(x_centroid(5)-x_centroid(1));
                   b_15 = y_centroid(1)-(m_15*x_centroid(1));
                   % 2) Calculate Coordinate of d_15/2
                   x_c_15 = x_centroid(1) + (x_centroid(5)-x_centroid(1))/2;
                   y_c_15 = m_15*x_c_15 + b_15
                   % Calculate max distance in x from sphere 2 to 1:
                   x_15_1 = x_centroid(1) - (maxaddedlength);
                   x_15_2 = x_c_15
                   %Introduce ratio of x/y (to rule out cases where two
                   %spheres have a very small differences in x or y but not
                   %in the other coordinate
                   xy_4 = abs(x_centroid(5)-x_centroid(1))/abs(y_centroid(5)-y_centroid(1))
                   yx_4 = abs(y_centroid(5)-y_centroid(1))/abs(x_centroid(5)-x_centroid(1))
                   if(xy_4 < minxy_ratio && d1_5 > mindist)
                       x_c_15 = 1000000;
                       y_c_15 = 1000000;
                       x_15_1 = x_centroid(1) - (maxaddedlength);
                       x_15_2 = x_centroid(1) + (maxaddedlength);
                       y_15_1 = y_centroid(1)-(maxaddedlength);
                       y_15_2 = y_centroid(1)+(maxaddedlength);
                   end  
                   if(yx_4 < minxy_ratio && d1_5 > mindist)
                       x_c_15 = 1000000;
                       y_c_15 = 1000000;
                       x_15_1 = x_centroid(1) - (maxaddedlength);
                       x_15_2 = x_centroid(1) + (maxaddedlength);
                       y_15_1 = y_centroid(1)-(maxaddedlength);
                       y_15_2 = y_centroid(1)+(maxaddedlength);
                   end  
                   if(y_c_15 > y_centroid(1))
                       y_15_1 = y_centroid(1)-(maxaddedlength)
                       y_15_2 = y_c_15
                   else 
                       y_15_1 = y_c_15
                       y_15_2 = y_centroid(1)+(maxaddedlength)
                   end    
               else
                   S15 = 0;
                   x_c_15 = 1000000;
                   y_c_15 = 1000000;
                   x_15_1 = x_centroid(1) - (maxaddedlength);
                   x_15_2 = x_centroid(1) + (maxaddedlength);
                   y_15_1 = y_centroid(1)-(maxaddedlength);
                   y_15_2 = y_centroid(1)+(maxaddedlength);
               end
               
             
               
               min_x_1 = x_centroid(1) - (maxaddedlength);
               
               V_x6 = [x_c_12;x_c_13;x_c_14;x_c_15];
               min_x_2 = min(V_x6(1:4));
               
                            
               if(min_x_2 == 1000000)
                   min_x_2 = x_centroid(1) + (maxaddedlength);
               end    
                                             
               if(y_c_12 < 1000000 || y_c_13 < 1000000 || y_c_14 < 1000000 || y_c_15 < 1000000)
                   
                   if(y_c_12 - y_centroid(1)>0 && y_c_12 - y_centroid(1)<500000);
                       y_c_12_s = 1000000;
                   else
                       y_c_12_s = y_c_12; 
                   end 
                   if(y_c_13 - y_centroid(1)>0 && y_c_13 - y_centroid(1)<500000);
                       y_c_13_s = 1000000;
                   else
                       y_c_13_s = y_c_13;     
                   end
                   if(y_c_14 - y_centroid(1)>0 && y_c_14 - y_centroid(1)<500000);
                       y_c_14_s = 1000000;
                   else
                       y_c_14_s = y_c_14;     
                   end
                   if(y_c_15 - y_centroid(1)>0 && y_c_15 - y_centroid(1)<500000);
                       y_c_15_s = 1000000;
                   else
                       y_c_15_s = y_c_15;    
                   end
                   V_y5=[y_c_12_s;y_c_13_s;y_c_14_s;y_c_15_s];
                   min_y_1 = min(V_y5(1:4));
                   if(min_y_1 == 1000000)
                       min_y_1 = y_centroid(1)- (maxaddedlength);
                   end
                   
               else
                   min_y_1 = y_centroid(1)- (maxaddedlength);
               end
               
               if(y_c_12 < 1000000 || y_c_13 < 1000000 || y_c_14 < 1000000 || y_c_15 < 1000000)
                   if(y_c_12 - y_centroid(1)<0 && y_c_12 - y_centroid(1)<500000);
                       y_c_12_s = 1000000;
                   else
                       y_c_12_s = y_c_12;
                   end    
                   if(y_c_13 - y_centroid(1)<0 && y_c_13 - y_centroid(1)<500000);
                       y_c_13_s = 1000000;
                   else
                       y_c_13_s = y_c_13;
                   end
                   if(y_c_14 - y_centroid(1)<0 && y_c_14 - y_centroid(1)<500000);
                       y_c_14_s = 1000000;
                   else
                       y_c_14_s = y_c_14;
                   end
                   if(y_c_15 - y_centroid(1)<0 && y_c_15 - y_centroid(1)<500000);
                       y_c_15_s = 1000000;
                   else
                       y_c_15_s = y_c_15;
                   end
                   V_y6=[y_c_12_s;y_c_13_s;y_c_14_s;y_c_15_s];
                   min_y_2 = min(V_y6(1:4));
                   if(min_y_2 == 1000000)
                       min_y_2 = y_centroid(1)+ (maxaddedlength);
                   end    
                   
               else
                   min_y_2 = y_centroid(1)+ (maxaddedlength);
               end
               
                                  
               if(min_x_1<1)
                   x_1 = 1
               else
                   x_1 = min_x_1;
               end
               
               if(min_x_2>10752)
                   x_2 = 10752
               else
                   x_2 = min_x_2;
               end
               
               if(min_y_1<1)
                   y_1 = 1
               else
                   y_1 = min_y_1;
               end
               
               if(min_y_2>13824)
                   y_2 = 13824
               else
                   y_2 = min_y_2;
               end
                     
               diffXCentroidl = x_centroid(1) - x_1;
               diffXCentroidr = x_2 - x_centroid(1);
               diffYCentroidl = y_centroid(1) - y_1;
               diffYCentroidr = y_2 - y_centroid(1);
               if(diffXCentroidl < diffXCentroidr)
               leftBoarder = 1;
               else 
               leftBoarder = 0;
               end
               
               imgnucleus = imread([foldername '/' allfiles(i).name]);
               %imwrite(imgnucleus,[foldername2 '/' allfiles(i).name]);
               %The rectangular image is cutted from the imgnucleus image containing the whole cover slip
               imagesphere1 = imgnucleus(y_1:y_2,x_1:x_2);
               % Dimensions of imagesphere1 are determined
               [r1 c1] = size(imagesphere1);
               % An empty matrix is generated
               imagesphereresize1 = zeros(7168,7168);
               % imagesphereresize1 is converted to 8-bit
               imagesphereresize1 = uint8(imagesphereresize1);
               %The empty indices in the imagesphereresize1 matrix are overwritten with imagesphere 1 entries
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
               % imagesphereresize1 is saved to final destination
               newFile1 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize1,[foldername1 '/' newFile1]);
               wellname = (newFile1(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(7168,7168);
               currentNM = csvHandler.CellPosMatrix(wellname);
               currentNM(yStart1:yEnd1,xStart1:xEnd1) = NucleusM(y_1:y_2,x_1:x_2);
               csvHandler.CellPosMatrix(wellname) = currentNM;
               imagesphereresize1_small = imresize(imagesphereresize1, optionHandler.ScalingFactor);
               newFile1s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize1_small,[foldername1 '/' newFile1s]);
               clear imagesphere1;
               clear imagesphereresize1;
               clear imagesphereresize1_small;
               
               % Same for sphere core 2
              if(d1_2<MinimumDistance);
                   S21 = 1;
                   % Calculate line equation for d_12 
                   m_12 = (y_centroid(2)-y_centroid(1))/(x_centroid(2)-x_centroid(1));
                   b_12 = y_centroid(1)-(m_12*x_centroid(1));
                   % 2) Calculate Coordinate of d_12/2
                   x_c_12 = x_centroid(1) + (x_centroid(2)-x_centroid(1))/2;
                   y_c_12 = m_12*x_c_12 + b_12
                   % Calculate max distance in x from sphere 2 to 1:
                   x_12_1 = x_c_12;
                   x_12_2 = x_centroid(2) + (maxaddedlength);
                   %Introduce ratio of x/y (to rule out cases where two
                   %spheres have a very small differences in x or y but not
                   %in the other coordinate
                   xy_1 = abs(x_centroid(2)-x_centroid(1))/abs(y_centroid(2)-y_centroid(1))
                   yx_1 = abs(y_centroid(2)-y_centroid(1))/abs(x_centroid(2)-x_centroid(1))
                   if(xy_1 < minxy_ratio && d1_2 > mindist)
                       x_c_12 = 1000000;
                       y_c_12 = 1000000;
                       x_12_1 = x_centroid(2) - (maxaddedlength);
                       x_12_2 = x_centroid(2) + (maxaddedlength);
                       y_12_1 = y_centroid(2)-(maxaddedlength);
                       y_12_2 = y_centroid(2)+(maxaddedlength);
                   end  
                   if(yx_1 < minxy_ratio && d1_2 > mindist)
                       x_c_12 = 1000000;
                       y_c_12 = 1000000;
                       x_12_1 = x_centroid(2) - (maxaddedlength);
                       x_12_2 = x_centroid(2) + (maxaddedlength);
                       y_12_1 = y_centroid(2)-(maxaddedlength);
                       y_12_2 = y_centroid(2)+(maxaddedlength);
                   end 
                   if(y_c_12 > y_centroid(2))
                       y_12_1 = y_centroid(2)-(maxaddedlength)
                       y_12_2 = y_c_12
                   else 
                       y_12_1 = y_c_12
                       y_12_2 = y_centroid(2)+(maxaddedlength)
                   end    
               else
                   S12 = 0;
                   x_c_12 = 1000000;
                   y_c_12 = 1000000;
                   x_12_1 = x_centroid(2) - (maxaddedlength);
                   x_12_2 = x_centroid(2) + (maxaddedlength);
                   y_12_1 = y_centroid(2)-(maxaddedlength);
                   y_12_2 = y_centroid(2)+(maxaddedlength);
              end
               
              if(d2_3<MinimumDistance);
                   S23 = 1;
                   % Calculate line equation for d_12 
                   m_23 = (y_centroid(3)-y_centroid(2))/(x_centroid(3)-x_centroid(2));
                   b_23 = y_centroid(2)-(m_23*x_centroid(2));
                   % 2) Calculate Coordinate of d_12/2
                   x_c_23 = x_centroid(2) + (x_centroid(3)-x_centroid(2))/2;
                   y_c_23 = m_23*x_c_23 + b_23
                   % Calculate max distance in x from sphere 2 to 1:
                   x_23_1 = x_centroid(2) - (maxaddedlength);
                   x_23_2 = x_c_23;
                   xy_2 = abs(x_centroid(2)-x_centroid(3))/abs(y_centroid(2)-y_centroid(3))
                   yx_2 = abs(y_centroid(2)-y_centroid(3))/abs(x_centroid(2)-x_centroid(3))
                   if(xy_2 < minxy_ratio && d2_3 > mindist)
                       x_c_23 = 1000000;
                       y_c_23 = 1000000;
                       x_23_1 = x_centroid(2) - (maxaddedlength);
                       x_23_2 = x_centroid(2) + (maxaddedlength);
                       y_23_1 = y_centroid(2)-(maxaddedlength);
                       y_23_2 = y_centroid(2)+(maxaddedlength);
                   end
                   if(yx_2 < minxy_ratio && d2_3 > mindist)
                       x_c_23 = 1000000;
                       y_c_23 = 1000000;
                       x_23_1 = x_centroid(2) - (maxaddedlength);
                       x_23_2 = x_centroid(2) + (maxaddedlength);
                       y_23_1 = y_centroid(2)-(maxaddedlength);
                       y_23_2 = y_centroid(2)+(maxaddedlength);
                   end
                   if(y_c_23 > y_centroid(2))
                       y_23_1 = y_centroid(2)-(maxaddedlength)
                       y_23_2 = y_c_23
                   else 
                       y_23_1 = y_c_23
                       y_23_2 = y_centroid(2)+(maxaddedlength)
                   end    
               else
                   S12 = 0;
                   x_c_23 = 1000000;
                   y_c_23 = 1000000;
                   x_23_1 = x_centroid(2) - (maxaddedlength);
                   x_23_2 = x_centroid(2) + (maxaddedlength);
                   y_23_1 = y_centroid(2)-(maxaddedlength);
                   y_23_2 = y_centroid(2)+(maxaddedlength);
              end
              
              if(d2_4<MinimumDistance);
                   S24 = 1;
                   % Calculate line equation for d_12 
                   m_24 = (y_centroid(4)-y_centroid(2))/(x_centroid(4)-x_centroid(2));
                   b_24 = y_centroid(2)-(m_24*x_centroid(2));
                   % 2) Calculate Coordinate of d_12/2
                   x_c_24 = x_centroid(2) + (x_centroid(4)-x_centroid(2))/2;
                   y_c_24 = m_24*x_c_24 + b_24
                   % Calculate max distance in x from sphere 2 to 1:
                   x_24_1 = x_centroid(2) - (maxaddedlength);
                   x_24_2 = x_c_24;
                   xy_3 = abs(x_centroid(2)-x_centroid(4))/abs(y_centroid(2)-y_centroid(4))
                   yx_3 = abs(y_centroid(2)-y_centroid(4))/abs(x_centroid(2)-x_centroid(4))
                   if(xy_3 < minxy_ratio && d2_4 > mindist)
                       x_c_24 = 1000000;
                       y_c_24 = 1000000;
                       x_24_1 = x_centroid(2) - (maxaddedlength);
                       x_24_2 = x_centroid(2) + (maxaddedlength);
                       y_24_1 = y_centroid(2)-(maxaddedlength);
                       y_24_2 = y_centroid(2)+(maxaddedlength);
                   end   
                   if(yx_3 < minxy_ratio && d2_4 > mindist)
                       x_c_24 = 1000000;
                       y_c_24 = 1000000;
                       x_24_1 = x_centroid(2) - (maxaddedlength);
                       x_24_2 = x_centroid(2) + (maxaddedlength);
                       y_24_1 = y_centroid(2)-(maxaddedlength);
                       y_24_2 = y_centroid(2)+(maxaddedlength);
                   end  
                   if(y_c_24 > y_centroid(2))
                       y_24_1 = y_centroid(2)-(maxaddedlength)
                       y_24_2 = y_c_24
                   else 
                       y_24_1 = y_c_24
                       y_24_2 = y_centroid(2)+(maxaddedlength)
                   end    
               else
                   S24 = 0;
                   x_c_24 = 1000000;
                   y_c_24 = 1000000;
                   x_24_1 = x_centroid(2) - (maxaddedlength);
                   x_24_2 = x_centroid(2) + (maxaddedlength);
                   y_24_1 = y_centroid(2)-(maxaddedlength);
                   y_24_2 = y_centroid(2)+(maxaddedlength);
              end
               if(d2_5<MinimumDistance);
                   S25 = 1;
                   % Calculate line equation for d_12 
                   m_25 = (y_centroid(5)-y_centroid(2))/(x_centroid(5)-x_centroid(2));
                   b_25 = y_centroid(2)-(m_25*x_centroid(2));
                   % 2) Calculate Coordinate of d_12/2
                   x_c_25 = x_centroid(2) + (x_centroid(5)-x_centroid(2))/2;
                   y_c_25 = m_25*x_c_25 + b_25
                   % Calculate max distance in x from sphere 2 to 1:
                   x_25_1 = x_centroid(2) - (maxaddedlength);
                   x_25_2 = x_c_25;
                   xy_4 = abs(x_centroid(2)-x_centroid(5))/abs(y_centroid(2)-y_centroid(5))
                   yx_4 = abs(y_centroid(2)-y_centroid(5))/abs(x_centroid(2)-x_centroid(5))
                   if(xy_4 < minxy_ratio && d2_5 > mindist)
                       x_c_25 = 1000000;
                       y_c_25 = 1000000;
                       x_25_1 = x_centroid(2) - (maxaddedlength);
                       x_25_2 = x_centroid(2) + (maxaddedlength);
                       y_25_1 = y_centroid(2)-(maxaddedlength);
                       y_25_2 = y_centroid(2)+(maxaddedlength);
                   end  
                   if(yx_4 < minxy_ratio && d2_5 > mindist)
                       x_c_25 = 1000000;
                       y_c_25 = 1000000;
                       x_25_1 = x_centroid(2) - (maxaddedlength);
                       x_25_2 = x_centroid(2) + (maxaddedlength);
                       y_25_1 = y_centroid(2)-(maxaddedlength);
                       y_25_2 = y_centroid(2)+(maxaddedlength);
                   end
                   if(y_c_25 > y_centroid(2))
                       y_25_1 = y_centroid(2)-(maxaddedlength)
                       y_25_2 = y_c_25
                   else 
                       y_25_1 = y_c_25
                       y_25_2 = y_centroid(2)+(maxaddedlength)
                   end    
               else
                   S25 = 0;
                   x_c_25 = 1000000;
                   y_c_25 = 1000000;
                   x_25_1 = x_centroid(2) - (maxaddedlength);
                   x_25_2 = x_centroid(2) + (maxaddedlength);
                   y_25_1 = y_centroid(2)-(maxaddedlength);
                   y_25_2 = y_centroid(2)+(maxaddedlength);
               end
               
               if(x_c_12 < 1000000)
                   min_x_3 = x_c_12;
               else
                   min_x_3 = x_centroid(2) - (maxaddedlength);
               end    
               
               V_x6 = [x_c_23;x_c_24;x_c_25];
               min_x_4 = min(V_x6(1:3));
               
               if(min_x_4 == 1000000)
                   min_x_4 = x_centroid(2) + (maxaddedlength);
               end    
                                             
               if(y_c_12 < 1000000 || y_c_23 < 1000000 || y_c_24 < 1000000 || y_c_25 < 1000000)
                   
                   if(y_c_12 - y_centroid(2)>0 && y_c_12 - y_centroid(2)<500000);
                       y_c_12_s = 1000000;
                   else
                       y_c_12_s = y_c_12; 
                   end 
                   if(y_c_23 - y_centroid(2)>0 && y_c_23 - y_centroid(2)<500000);
                       y_c_23_s = 1000000;
                   else
                       y_c_23_s = y_c_23;     
                   end
                   if(y_c_24 - y_centroid(2)>0 && y_c_24 - y_centroid(2)<500000);
                       y_c_24_s = 1000000;
                   else
                       y_c_24_s = y_c_24;     
                   end
                   if(y_c_25 - y_centroid(2)>0 && y_c_25 - y_centroid(2)<500000);
                       y_c_25_s = 1000000;
                   else
                       y_c_25_s = y_c_25;    
                   end
                   V_y5=[y_c_12_s;y_c_23_s;y_c_24_s;y_c_25_s];
                   min_y_3 = min(V_y5(1:4));
                   if(min_y_3 == 1000000)
                       min_y_3 = y_centroid(2)- (maxaddedlength);
                   end
                   
               else
                   min_y_3 = y_centroid(2)- (maxaddedlength);
               end
               
               if(y_c_12 < 1000000 || y_c_23 < 1000000 || y_c_24 < 1000000 || y_c_25 < 1000000)
                   if(y_c_12 - y_centroid(2)<0 && y_c_12 - y_centroid(2)<500000);
                       y_c_12_s = 1000000;
                   else
                       y_c_12_s = y_c_12;
                   end    
                   if(y_c_23 - y_centroid(2)<0 && y_c_23 - y_centroid(2)<500000);
                       y_c_23_s = 1000000;
                   else
                       y_c_23_s = y_c_23;
                   end
                   if(y_c_24 - y_centroid(2)<0 && y_c_24 - y_centroid(2)<500000);
                       y_c_24_s = 1000000;
                   else
                       y_c_24_s = y_c_24;
                   end
                   if(y_c_25 - y_centroid(2)<0 && y_c_25 - y_centroid(2)<500000);
                       y_c_25_s = 1000000;
                   else
                       y_c_25_s = y_c_25;
                   end
                   V_y6=[y_c_12_s;y_c_23_s;y_c_24_s;y_c_25_s];
                   min_y_4 = min(V_y6(1:4));
                   if(min_y_4 == 1000000)
                       min_y_4 = y_centroid(2)+ (maxaddedlength);
                   end    
                   
               else
                   min_y_4 = y_centroid(2)+ (maxaddedlength);
               end
                   
               if(min_x_3<1)
                   x_3 = 1
               else
                   x_3 = min_x_3;
               end
               
               if(min_x_4>10752)
                   x_4 = 10752
               else
                   x_4 = min_x_4;
               end
               
               if(min_y_3<1)
                   y_3 = 1
               else
                   y_3 = min_y_3;
               end
               
               if(min_y_4>13824)
                   y_4 = 13824
               else
                   y_4 = min_y_4;
               end
               
               diffXCentroidl = x_centroid(2) - x_3;
               diffXCentroidr = x_4 - x_centroid(2);
               diffYCentroidl = y_centroid(2) - y_3;
               diffYCentroidr = y_4 - y_centroid(2);
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
               
               imagesphere2 = imgnucleus(y_3:y_4,x_3:x_4);
               [r2 c2] = size(imagesphere2);
               imagesphereresize2 = zeros(7168,7168);
               imagesphereresize2 = uint8(imagesphereresize2);
               if(leftBoarder == 0)
                   if(lowerBoarder == 0)
                        
                        yStart2 = 1;
                        yEnd2 = r2;
                        xStart2 = 1;
                        xEnd2 = c2;
                   else
                        
                        yStart2 = (7169-r2);
                        yEnd2 = 7168;
                        xStart2 = 1;
                        xEnd2 = c2;
                   end 
               else
                   if(lowerBoarder == 0)
                       
                        yStart2 = 1;
                        yEnd2 = r2;
                        xStart2 = (7169 - c2);
                        xEnd2 = 7168;
                   else
                       
                        yStart2 = (7169-r2);
                        yEnd2 = 7168;
                        xStart2 = (7169 - c2);
                        xEnd2 = 7168;
                   end
               end  
               
               
               imagesphereresize2(yStart2:yEnd2,xStart2:xEnd2)=imagesphere2(1:r2,1:c2);
               newFile2 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize2,[foldername1 '/' newFile2]);
               n=n+1;
               wellname = (newFile2(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(7168,7168);
               currentNM = csvHandler.CellPosMatrix(wellname);
               currentNM(yStart2:yEnd2,xStart2:xEnd2) = NucleusM(y_3:y_4,x_3:x_4);
               csvHandler.CellPosMatrix(wellname) = currentNM;
               imagesphereresize2_small = imresize(imagesphereresize2, optionHandler.ScalingFactor);
               newFile2s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize2_small,[foldername1 '/' newFile2s]);
               clear imagesphere2;
               clear imagesphereresize2;
               clear imagesphereresize2_small;
               
               % Same for sphere core 3
               
               if(d1_3<MinimumDistance);
                   S31 = 1;
                   % Calculate line equation for d_12 
                   m_13 = (y_centroid(3)-y_centroid(1))/(x_centroid(3)-x_centroid(1));
                   b_13 = y_centroid(1)-(m_13*x_centroid(1));
                   % 2) Calculate Coordinate of d_12/2
                   x_c_13 = x_centroid(1) + (x_centroid(3)-x_centroid(1))/2;
                   y_c_13 = m_13*x_c_13 + b_13
                   % Calculate max distance in x from sphere 2 to 1:
                   x_13_1 = x_c_13;
                   x_13_2 = x_centroid(3) + (maxaddedlength);
                   %Introduce ratio of x/y (to rule out cases where two
                   %spheres have a very small differences in x or y but not
                   %in the other coordinate
                   xy_1 = abs(x_centroid(3)-x_centroid(1))/abs(y_centroid(3)-y_centroid(1))
                   yx_1 = abs(y_centroid(3)-y_centroid(1))/abs(x_centroid(3)-x_centroid(1))
                   if(xy_1 < minxy_ratio && d1_3 > mindist)
                       x_c_13 = 1000000;
                       y_c_13 = 1000000;
                       x_13_1 = x_centroid(3) - (maxaddedlength);
                       x_13_2 = x_centroid(3) + (maxaddedlength);
                       y_13_1 = y_centroid(3)-(maxaddedlength);
                       y_13_2 = y_centroid(3)+(maxaddedlength);
                   end  
                   if(yx_1 < minxy_ratio && d1_3 > mindist)
                       x_c_13 = 1000000;
                       y_c_13 = 1000000;
                       x_13_1 = x_centroid(3) - (maxaddedlength);
                       x_13_2 = x_centroid(3) + (maxaddedlength);
                       y_13_1 = y_centroid(3)-(maxaddedlength);
                       y_13_2 = y_centroid(3)+(maxaddedlength);
                   end  
                   if(y_c_13 > y_centroid(3))
                       y_13_1 = y_centroid(3)-(maxaddedlength)
                       y_13_2 = y_c_13
                   else 
                       y_13_1 = y_c_13
                       y_13_2 = y_centroid(3)+(maxaddedlength)
                   end    
               else
                   S31 = 0;
                   x_c_13 = 1000000;
                   y_c_13 = 1000000;
                   x_13_1 = x_centroid(3) - (maxaddedlength);
                   x_13_2 = x_centroid(3) + (maxaddedlength);
                   y_13_1 = y_centroid(3)-(maxaddedlength);
                   y_13_2 = y_centroid(3)+(maxaddedlength);
              end
               
              if(d2_3<MinimumDistance);
                   S32 = 1;
                   % Calculate line equation for d_12 
                   m_23 = (y_centroid(3)-y_centroid(2))/(x_centroid(3)-x_centroid(2));
                   b_23 = y_centroid(2)-(m_23*x_centroid(2));
                   % 2) Calculate Coordinate of d_12/2
                   x_c_23 = x_centroid(2) + (x_centroid(3)-x_centroid(2))/2;
                   y_c_23 = m_23*x_c_23 + b_23
                   % Calculate max distance in x from sphere 2 to 1:
                   x_23_1 = x_c_23;
                   x_23_2 = x_centroid(3) + (maxaddedlength);
                   xy_2 = abs(x_centroid(3)-x_centroid(2))/abs(y_centroid(3)-y_centroid(2))
                   yx_2 = abs(y_centroid(3)-y_centroid(2))/abs(x_centroid(3)-x_centroid(2))
                   if(xy_2 < minxy_ratio && d2_3 > mindist)
                       x_c_23 = 1000000;
                       y_c_23 = 1000000;
                       x_23_1 = x_centroid(3) - (maxaddedlength);
                       x_23_2 = x_centroid(3) + (maxaddedlength);
                       y_23_1 = y_centroid(3)-(maxaddedlength);
                       y_23_2 = y_centroid(3)+(maxaddedlength);
                   end  
                   if(yx_2 < minxy_ratio && d2_3 > mindist)
                       x_c_23 = 1000000;
                       y_c_23 = 1000000;
                       x_23_1 = x_centroid(3) - (maxaddedlength);
                       x_23_2 = x_centroid(3) + (maxaddedlength);
                       y_23_1 = y_centroid(3)-(maxaddedlength);
                       y_23_2 = y_centroid(3)+(maxaddedlength);
                   end
                   if(y_c_23 > y_centroid(3))
                       y_23_1 = y_centroid(3)-(maxaddedlength)
                       y_23_2 = y_c_23
                   else 
                       y_23_1 = y_c_23
                       y_23_2 = y_centroid(3)+(maxaddedlength)
                   end    
               else
                   S32 = 0;
                   x_c_23 = 1000000;
                   y_c_23 = 1000000;
                   x_23_1 = x_centroid(3) - (maxaddedlength);
                   x_23_2 = x_centroid(3) + (maxaddedlength);
                   y_23_1 = y_centroid(3)-(maxaddedlength);
                   y_23_2 = y_centroid(3)+(maxaddedlength);
              end
              
              if(d3_4<MinimumDistance);
                   S34 = 1;
                   % Calculate line equation for d_12 
                   m_34 = (y_centroid(4)-y_centroid(3))/(x_centroid(4)-x_centroid(3));
                   b_34 = y_centroid(3)-(m_34*x_centroid(3));
                   % 2) Calculate Coordinate of d_12/2
                   x_c_34 = x_centroid(3) + (x_centroid(4)-x_centroid(3))/2;
                   y_c_34 = m_34*x_c_34 + b_34
                   % Calculate max distance in x from sphere 2 to 1:
                   x_34_1 = x_centroid(3) - (maxaddedlength);
                   x_34_2 = x_c_34;
                   xy_3 = abs(x_centroid(3)-x_centroid(4))/abs(y_centroid(3)-y_centroid(4))
                   yx_3 = abs(y_centroid(3)-y_centroid(4))/abs(x_centroid(3)-x_centroid(4))
                   if(xy_3 < minxy_ratio && d3_4 > mindist)
                       x_c_34 = 1000000;
                       y_c_34 = 1000000;
                       x_34_1 = x_centroid(3) - (maxaddedlength);
                       x_34_2 = x_centroid(3) + (maxaddedlength);
                       y_34_1 = y_centroid(3)-(maxaddedlength);
                       y_34_2 = y_centroid(3)+(maxaddedlength);
                   end  
                   if(yx_3 < minxy_ratio && d3_4 > mindist)
                       x_c_34 = 1000000;
                       y_c_34 = 1000000;
                       x_34_1 = x_centroid(3) - (maxaddedlength);
                       x_34_2 = x_centroid(3) + (maxaddedlength);
                       y_34_1 = y_centroid(3)-(maxaddedlength);
                       y_34_2 = y_centroid(3)+(maxaddedlength);
                   end
                   if(y_c_34 > y_centroid(3))
                       y_34_1 = y_centroid(3)-(maxaddedlength)
                       y_34_2 = y_c_34
                   else 
                       y_34_1 = y_c_34
                       y_34_2 = y_centroid(3)+(maxaddedlength)
                   end    
               else
                   S34 = 0;
                   x_c_34 = 1000000;
                   y_c_34 = 1000000;
                   x_34_1 = x_centroid(3) - (maxaddedlength);
                   x_34_2 = x_centroid(3) + (maxaddedlength);
                   y_34_1 = y_centroid(3)-(maxaddedlength);
                   y_34_2 = y_centroid(3)+(maxaddedlength);
              end
               if(d3_5<MinimumDistance);
                   S35 = 1;
                   % Calculate line equation for d_12 
                   m_35 = (y_centroid(5)-y_centroid(3))/(x_centroid(5)-x_centroid(3));
                   b_35 = y_centroid(3)-(m_35*x_centroid(3));
                   % 2) Calculate Coordinate of d_12/2
                   x_c_35 = x_centroid(3) + (x_centroid(5)-x_centroid(3))/2;
                   y_c_35 = m_35*x_c_35 + b_35
                   % Calculate max distance in x from sphere 2 to 1:
                   x_35_1 = x_centroid(3) - (maxaddedlength);
                   x_35_2 = x_c_35;
                   xy_4 = abs(x_centroid(3)-x_centroid(5))/abs(y_centroid(3)-y_centroid(5))
                   yx_4 = abs(y_centroid(3)-y_centroid(5))/abs(x_centroid(3)-x_centroid(5))
                   if(xy_4 < minxy_ratio && d3_5 > mindist)
                       x_c_35 = 1000000;
                       y_c_35 = 1000000;
                       x_35_1 = x_centroid(3) - (maxaddedlength);
                       x_35_2 = x_centroid(3) + (maxaddedlength);
                       y_35_1 = y_centroid(3)-(maxaddedlength);
                       y_35_2 = y_centroid(3)+(maxaddedlength);
                   end  
                   if(yx_4 < minxy_ratio && d3_5 > mindist)
                       x_c_35 = 1000000;
                       y_c_35 = 1000000;
                       x_35_1 = x_centroid(3) - (maxaddedlength);
                       x_35_2 = x_centroid(3) + (maxaddedlength);
                       y_35_1 = y_centroid(3)-(maxaddedlength);
                       y_35_2 = y_centroid(3)+(maxaddedlength);
                   end
                   if(y_c_35 > y_centroid(3))
                       y_35_1 = y_centroid(3)-(maxaddedlength)
                       y_35_2 = y_c_35
                   else 
                       y_35_1 = y_c_35
                       y_35_2 = y_centroid(3)+(maxaddedlength)
                   end    
               else
                   S35 = 0;
                   x_c_35 = 1000000;
                   y_c_35 = 1000000;
                   x_35_1 = x_centroid(3) - (maxaddedlength);
                   x_35_2 = x_centroid(3) + (maxaddedlength);
                   y_35_1 = y_centroid(3)-(maxaddedlength);
                   y_35_2 = y_centroid(3)+(maxaddedlength);
               end
                              
               if(x_c_13 >0 && x_c_23>0 && x_c_13 <1000000 && x_c_23<1000000 )
                   V_x5=[x_c_13;x_c_23];
                   min_x_5 = min(V_x5(1:2));
               else
                   if(x_c_13>0 && x_c_13 <1000000)
                       min_x_5 = x_c_13
                   else
                       if(x_c_23>0 && x_c_23<1000000)
                       min_x_5 = x_c_23
                       else
                           min_x_5 = x_centroid(3)-(maxaddedlength);
                       end
                   end
               end 
               
               if(x_c_34 >0 && x_c_35>0 && x_c_34 <1000000 && x_c_35<1000000)
                   V_x6=[x_c_34;x_c_35];
                   min_x_6 = min(V_x6(1:2));
               else
                   if(x_c_34>0 && x_c_34 <1000000)
                       min_x_6 = x_c_34
                   else
                       if(x_c_35>0 && x_c_35<1000000)
                       min_x_6 = x_c_35
                       else
                           min_x_6 = x_centroid(3)+(maxaddedlength);
                       end
                   end
               end 
               
               
                                             
               if(y_c_13 < 1000000 || y_c_23 < 1000000 || y_c_34 < 1000000 || y_c_35 < 1000000)
                   
                   if(y_c_13 - y_centroid(3)>0 && y_c_13 - y_centroid(3)<500000);
                       y_c_13_s = 1000000;
                   else
                       y_c_13_s = y_c_13; 
                   end 
                   if(y_c_23 - y_centroid(3)>0 && y_c_23 - y_centroid(3)<500000);
                       y_c_23_s = 1000000;
                   else
                       y_c_23_s = y_c_23;     
                   end
                   if(y_c_34 - y_centroid(3)>0 && y_c_34 - y_centroid(3)<500000);
                       y_c_34_s = 1000000;
                   else
                       y_c_34_s = y_c_34;     
                   end
                   if(y_c_35 - y_centroid(3)>0 && y_c_35 - y_centroid(3)<500000);
                       y_c_35_s = 1000000;
                   else
                       y_c_35_s = y_c_35;    
                   end
                   V_y5=[y_c_13_s;y_c_23_s;y_c_34_s;y_c_35_s];
                   min_y_5 = min(V_y5(1:4));
                   if(min_y_5 == 1000000)
                       min_y_5 = y_centroid(3)- (maxaddedlength);
                   end
                   
               else
                   min_y_5 = y_centroid(3)- (maxaddedlength);
               end
               
               if(y_c_13 < 1000000 || y_c_23 < 1000000 || y_c_34 < 1000000 || y_c_35 < 1000000)
                   if(y_c_13 - y_centroid(3)<0 && y_c_13 - y_centroid(3)<500000);
                       y_c_13_s = 1000000;
                   else
                       y_c_13_s = y_c_13;
                   end    
                   if(y_c_23 - y_centroid(3)<0 && y_c_23 - y_centroid(3)<500000);
                       y_c_23_s = 1000000;
                   else
                       y_c_23_s = y_c_23;
                   end
                   if(y_c_34 - y_centroid(3)<0 && y_c_34 - y_centroid(3)<500000);
                       y_c_34_s = 1000000;
                   else
                       y_c_34_s = y_c_34;
                   end
                   if(y_c_35 - y_centroid(3)<0 && y_c_35 - y_centroid(3)<500000);
                       y_c_35_s = 1000000;
                   else
                       y_c_35_s = y_c_35;
                   end
                   V_y6=[y_c_13_s;y_c_23_s;y_c_34_s;y_c_35_s];
                   min_y_6 = min(V_y6(1:4));
                   if(min_y_6 == 1000000)
                       min_y_6 = y_centroid(3)+ (maxaddedlength);
                   end    
                   
               else
                   min_y_6 = y_centroid(3)+ (maxaddedlength);
               end
                   
               if(min_x_5<1)
                   x_5 = 1
               else
                   x_5 = min_x_5;
               end
               
               if(min_x_6>10752)
                   x_6 = 10752
               else
                   x_6 = min_x_6;
               end
               
               if(min_y_5<1)
                   y_5 = 1
               else
                   y_5 = min_y_5;
               end
               
               if(min_y_6>13824)
                   y_6 = 13824
               else
                   y_6 = min_y_6;
               end
               
               
               diffXCentroidl = x_centroid(3) - x_5;
               diffXCentroidr = x_6 - x_centroid(3);
               diffYCentroidl = y_centroid(3) - y_5;
               diffYCentroidr = y_6 - y_centroid(3);
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
               
               
               
               imagesphere3 = imgnucleus(y_5:y_6,x_5:x_6);
               [r3 c3] = size(imagesphere3);
               imagesphereresize3 = zeros(7168,7168);
               imagesphereresize3 = uint8(imagesphereresize3);
               if(leftBoarder == 0)
                   if(lowerBoarder == 0)
                        
                        yStart3 = 1;
                        yEnd3 = r3;
                        xStart3 = 1;
                        xEnd3 = c3;
                   else
                        
                        yStart3 = (7169-r3);
                        yEnd3 = 7168;
                        xStart3 = 1;
                        xEnd3 = c3;
                   end 
               else
                   if(lowerBoarder == 0)
                       
                        yStart3 = 1;
                        yEnd3 = r3;
                        xStart3 = (7169 - c3);
                        xEnd3 = 7168;
                   else
                       
                        yStart3 = (7169-r3);
                        yEnd3 = 7168;
                        xStart3 = (7169 - c3);
                        xEnd3 = 7168;
                   end
               end  
               
               
               imagesphereresize3(yStart3:yEnd3,xStart3:xEnd3)=imagesphere3(1:r3,1:c3);
               newFile3 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize3,[foldername1 '/' newFile3]);
               n=n+1;
               wellname = (newFile3(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(7168,7168);
               currentNM = csvHandler.CellPosMatrix(wellname);
               currentNM(yStart3:yEnd3,xStart3:xEnd3) = NucleusM(y_5:y_6,x_5:x_6);
               csvHandler.CellPosMatrix(wellname) = currentNM;
               imagesphereresize3_small = imresize(imagesphereresize3, optionHandler.ScalingFactor);
               newFile3s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize3_small,[foldername1 '/' newFile3s]);
               clear imagesphere3;
               clear imagesphereresize3;
               clear imagesphereresize3_small;
               
               % Same for sphere core 4
               
               if(d1_4<MinimumDistance);
                   S41 = 1;
                   % Calculate line equation for d_12 
                   m_14 = (y_centroid(4)-y_centroid(1))/(x_centroid(4)-x_centroid(1));
                   b_14 = y_centroid(1)-(m_14*x_centroid(1));
                   % 2) Calculate Coordinate of d_12/2
                   x_c_14 = x_centroid(1) + (x_centroid(4)-x_centroid(1))/2;
                   y_c_14 = m_14*x_c_14 + b_14
                   % Calculate max distance in x from sphere 2 to 1:
                   x_14_1 = x_c_14;
                   x_14_2 = x_centroid(4) + (maxaddedlength);
                   %Introduce ratio of x/y (to rule out cases where two
                   %spheres have a very small differences in x or y but not
                   %in the other coordinate
                   xy_1 = abs(x_centroid(4)-x_centroid(1))/abs(y_centroid(4)-y_centroid(1))
                   yx_1 = abs(y_centroid(4)-y_centroid(1))/abs(x_centroid(4)-x_centroid(1))
                   if(xy_1 < minxy_ratio && d1_4 > mindist)
                       x_c_14 = 1000000;
                       y_c_14 = 1000000;
                       x_14_1 = x_centroid(4) - (maxaddedlength);
                       x_14_2 = x_centroid(4) + (maxaddedlength);
                       y_14_1 = y_centroid(4)-(maxaddedlength);
                       y_14_2 = y_centroid(4)+(maxaddedlength);
                   end  
                   if(yx_1 < minxy_ratio && d1_4 > mindist)
                       x_c_14 = 1000000;
                       y_c_14 = 1000000;
                       x_14_1 = x_centroid(4) - (maxaddedlength);
                       x_14_2 = x_centroid(4) + (maxaddedlength);
                       y_14_1 = y_centroid(4)-(maxaddedlength);
                       y_14_2 = y_centroid(4)+(maxaddedlength);
                   end  
                   if(y_c_14 > y_centroid(4))
                       y_14_1 = y_centroid(4)-(maxaddedlength)
                       y_14_2 = y_c_14
                   else 
                       y_14_1 = y_c_14
                       y_14_2 = y_centroid(4)+(maxaddedlength)
                   end    
               else
                   S41 = 0;
                   x_c_14 = 1000000;
                   y_c_14 = 1000000;
                   x_14_1 = x_centroid(4) - (maxaddedlength);
                   x_14_2 = x_centroid(4) + (maxaddedlength);
                   y_14_1 = y_centroid(4)-(maxaddedlength);
                   y_14_2 = y_centroid(4)+(maxaddedlength);
              end
               
              if(d2_4<MinimumDistance);
                   S42 = 1;
                   % Calculate line equation for d_12 
                   m_24 = (y_centroid(4)-y_centroid(2))/(x_centroid(4)-x_centroid(2));
                   b_24 = y_centroid(2)-(m_24*x_centroid(2));
                   % 2) Calculate Coordinate of d_12/2
                   x_c_24 = x_centroid(2) + (x_centroid(4)-x_centroid(2))/2;
                   y_c_24 = m_24*x_c_24 + b_24
                   % Calculate max distance in x from sphere 2 to 1:
                   x_24_1 = x_c_24;
                   x_24_2 = x_centroid(4) + (maxaddedlength);
                   xy_2 = abs(x_centroid(4)-x_centroid(2))/abs(y_centroid(4)-y_centroid(2))
                   yx_2 = abs(y_centroid(4)-y_centroid(2))/abs(x_centroid(4)-x_centroid(2))
                   if(xy_2 < minxy_ratio && d2_4 > mindist)
                       x_c_24 = 1000000;
                       y_c_24 = 1000000;
                       x_24_1 = x_centroid(4) - (maxaddedlength);
                       x_24_2 = x_centroid(4) + (maxaddedlength);
                       y_24_1 = y_centroid(4)-(maxaddedlength);
                       y_24_2 = y_centroid(4)+(maxaddedlength);
                   end  
                   if(yx_2 < minxy_ratio && d2_4 > mindist)
                       x_c_24 = 1000000;
                       y_c_24 = 1000000;
                       x_24_1 = x_centroid(4) - (maxaddedlength);
                       x_24_2 = x_centroid(4) + (maxaddedlength);
                       y_24_1 = y_centroid(4)-(maxaddedlength);
                       y_24_2 = y_centroid(4)+(maxaddedlength);
                   end
                   if(y_c_24 > y_centroid(4))
                       y_24_1 = y_centroid(4)-(maxaddedlength)
                       y_24_2 = y_c_24
                   else 
                       y_24_1 = y_c_24
                       y_24_2 = y_centroid(4)+(maxaddedlength)
                   end    
               else
                   S42 = 0;
                   x_c_24 = 1000000;
                   y_c_24 = 1000000;
                   x_24_1 = x_centroid(4) - (maxaddedlength);
                   x_24_2 = x_centroid(4) + (maxaddedlength);
                   y_24_1 = y_centroid(4)-(maxaddedlength);
                   y_24_2 = y_centroid(4)+(maxaddedlength);
              end
              
              if(d3_4<MinimumDistance);
                   S43 = 1;
                   % Calculate line equation for d_12 
                   m_34 = (y_centroid(4)-y_centroid(3))/(x_centroid(4)-x_centroid(3));
                   b_34 = y_centroid(3)-(m_34*x_centroid(3));
                   % 2) Calculate Coordinate of d_12/2
                   x_c_34 = x_centroid(3) + (x_centroid(4)-x_centroid(3))/2;
                   y_c_34 = m_34*x_c_34 + b_34
                   % Calculate max distance in x from sphere 2 to 1:
                   x_34_1 = x_c_34;
                   x_34_2 = x_centroid(4) - (maxaddedlength);
                   xy_3 = abs(x_centroid(3)-x_centroid(4))/abs(y_centroid(3)-y_centroid(4))
                   yx_3 = abs(y_centroid(3)-y_centroid(4))/abs(x_centroid(3)-x_centroid(4))
                   if(xy_3 < minxy_ratio && d3_4 > mindist)
                       x_c_34 = 1000000;
                       y_c_34 = 1000000;
                       x_34_1 = x_centroid(4) - (maxaddedlength);
                       x_34_2 = x_centroid(4) + (maxaddedlength);
                       y_34_1 = y_centroid(4)-(maxaddedlength);
                       y_34_2 = y_centroid(4)+(maxaddedlength);
                   end  
                   if(yx_3 < minxy_ratio && d3_4 > mindist)
                       x_c_34 = 1000000;
                       y_c_34 = 1000000;
                       x_34_1 = x_centroid(4) - (maxaddedlength);
                       x_34_2 = x_centroid(4) + (maxaddedlength);
                       y_34_1 = y_centroid(4)-(maxaddedlength);
                       y_34_2 = y_centroid(4)+(maxaddedlength);
                   end
                   if(y_c_34 > y_centroid(4))
                       y_34_1 = y_centroid(4)-(maxaddedlength)
                       y_34_2 = y_c_34
                   else 
                       y_34_1 = y_c_34
                       y_34_2 = y_centroid(4)+(maxaddedlength)
                   end    
               else
                   S34 = 0;
                   x_c_34 = 1000000;
                   y_c_34 = 1000000;
                   x_34_1 = x_centroid(4) - (maxaddedlength);
                   x_34_2 = x_centroid(4) + (maxaddedlength);
                   y_34_1 = y_centroid(4)-(maxaddedlength);
                   y_34_2 = y_centroid(4)+(maxaddedlength);
              end
               if(d4_5<MinimumDistance);
                   S54 = 1;
                   % Calculate line equation for d_12 
                   m_45 = (y_centroid(5)-y_centroid(4))/(x_centroid(5)-x_centroid(4));
                   b_45 = y_centroid(4)-(m_45*x_centroid(4));
                   % 2) Calculate Coordinate of d_12/2
                   x_c_45 = x_centroid(4) + (x_centroid(5)-x_centroid(4))/2;
                   y_c_45 = m_45*x_c_45 + b_45
                   % Calculate max distance in x from sphere 2 to 1:
                   x_45_1 = x_centroid(4) - (maxaddedlength);
                   x_45_2 = x_c_45;
                   xy_4 = abs(x_centroid(4)-x_centroid(5))/abs(y_centroid(4)-y_centroid(5))
                   yx_4 = abs(y_centroid(4)-y_centroid(5))/abs(x_centroid(4)-x_centroid(5))
                   if(xy_4 < minxy_ratio && d4_5 > mindist)
                       x_c_45 = 1000000;
                       y_c_45 = 1000000;
                       x_45_1 = x_centroid(4) - (maxaddedlength);
                       x_45_2 = x_centroid(4) + (maxaddedlength);
                       y_45_1 = y_centroid(4)-(maxaddedlength);
                       y_45_2 = y_centroid(4)+(maxaddedlength);
                   end  
                   if(yx_4 < minxy_ratio && d4_5 > mindist)
                       x_c_45 = 1000000;
                       y_c_45 = 1000000;
                       x_45_1 = x_centroid(4) - (maxaddedlength);
                       x_45_2 = x_centroid(4) + (maxaddedlength);
                       y_45_1 = y_centroid(4)-(maxaddedlength);
                       y_45_2 = y_centroid(4)+(maxaddedlength);
                   end
                   if(y_c_45 > y_centroid(4))
                       y_45_1 = y_centroid(4)-(maxaddedlength)
                       y_45_2 = y_c_45
                   else 
                       y_45_1 = y_c_45
                       y_45_2 = y_centroid(4)+(maxaddedlength)
                   end    
               else
                   S45 = 0;
                   x_c_45 = 1000000;
                   y_c_45 = 1000000;
                   x_45_1 = x_centroid(4) - (maxaddedlength);
                   x_45_2 = x_centroid(4) + (maxaddedlength);
                   y_45_1 = y_centroid(4)-(maxaddedlength);
                   y_45_2 = y_centroid(4)+(maxaddedlength);
               end
                              
               V_x7 = [x_c_14;x_c_24;x_c_35];
               min_x_7 = min(V_x7(1:3));
               if(min_x_7 ==1000000)
                   min_x_7 = x_centroid(4)-(maxaddedlength);
               end    
               min_x_8 = x_c_45;
               if(min_x_8 == 1000000)
                   min_x_8 = x_centroid(4) + (maxaddedlength);
               end    
               
                                             
               if(y_c_14 < 1000000 || y_c_24 < 1000000 || y_c_34 < 1000000 || y_c_45 < 1000000)
                   
                   if(y_c_14 - y_centroid(4)>0 && y_c_14 - y_centroid(4)<500000);
                       y_c_14_s = 1000000;
                   else
                       y_c_14_s = y_c_14; 
                   end 
                   if(y_c_24 - y_centroid(4)>0 && y_c_24 - y_centroid(4)<500000);
                       y_c_24_s = 1000000;
                   else
                       y_c_24_s = y_c_24;     
                   end
                   if(y_c_34 - y_centroid(4)>0 && y_c_34 - y_centroid(4)<500000);
                       y_c_34_s = 1000000;
                   else
                       y_c_34_s = y_c_34;     
                   end
                   if(y_c_45 - y_centroid(4)>0 && y_c_45 - y_centroid(4)<500000);
                       y_c_45_s = 1000000;
                   else
                       y_c_45_s = y_c_45;    
                   end
                   V_y5=[y_c_14_s;y_c_24_s;y_c_34_s;y_c_45_s];
                   min_y_7 = min(V_y5(1:4));
                   if(min_y_7 == 1000000)
                       min_y_7 = y_centroid(4)- (maxaddedlength);
                   end
                   
               else
                   min_y_7 = y_centroid(4)- (maxaddedlength);
               end
               
               if(y_c_14 < 1000000 || y_c_24 < 1000000 || y_c_34 < 1000000 || y_c_45 < 1000000)
                   if(y_c_14 - y_centroid(4)<0 && y_c_14 - y_centroid(4)<500000);
                       y_c_14_s = 1000000;
                   else
                       y_c_14_s = y_c_14;
                   end    
                   if(y_c_24 - y_centroid(4)<0 && y_c_24 - y_centroid(4)<500000);
                       y_c_24_s = 1000000;
                   else
                       y_c_24_s = y_c_24;
                   end
                   if(y_c_34 - y_centroid(4)<0 && y_c_34 - y_centroid(4)<500000);
                       y_c_34_s = 1000000;
                   else
                       y_c_34_s = y_c_34;
                   end
                   if(y_c_45 - y_centroid(4)<0 && y_c_45 - y_centroid(4)<500000);
                       y_c_45_s = 1000000;
                   else
                       y_c_45_s = y_c_45;
                   end
                   V_y6=[y_c_14_s;y_c_24_s;y_c_34_s;y_c_45_s];
                   min_y_8 = min(V_y6(1:4));
                   if(min_y_8 == 1000000)
                       min_y_8 = y_centroid(4)+ (maxaddedlength);
                   end    
                   
               else
                   min_y_8 = y_centroid(4)+ (maxaddedlength);
               end
                   
               if(min_x_7<1)
                   x_7 = 1
               else
                   x_7 = min_x_7;
               end
               
               if(min_x_8>10752)
                   x_8 = 10752
               else
                   x_8 = min_x_8;
               end
               
               if(min_y_7<1)
                   y_7 = 1
               else
                   y_7 = min_y_7;
               end
               
               if(min_y_8>13824)
                   y_8 = 13824
               else
                   y_8 = min_y_8;
               end
               
               diffXCentroidl = x_centroid(4) - x_7;
               diffXCentroidr = x_8 - x_centroid(4);
               diffYCentroidl = y_centroid(4) - y_7;
               diffYCentroidr = y_8 - y_centroid(4);
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
               
               
               imagesphere4 = imgnucleus(y_7:y_8,x_7:x_8);
               [r4 c4] = size(imagesphere4);
               imagesphereresize4 = zeros(7168,7168);
               imagesphereresize4 = uint8(imagesphereresize4);
               if(leftBoarder == 0)
                   if(lowerBoarder == 0)
                        
                        yStart4 = 1;
                        yEnd4 = r4;
                        xStart4 = 1;
                        xEnd4 = c4;
                   else
                        
                        yStart4 = (7169-r4);
                        yEnd4 = 7168;
                        xStart4 = 1;
                        xEnd4 = c4;
                   end 
               else
                   if(lowerBoarder == 0)
                       
                        yStart4 = 1;
                        yEnd4 = r4;
                        xStart4 = (7169 - c4);
                        xEnd4 = 7168;
                   else
                       
                        yStart4 = (7169-r4);
                        yEnd4 = 7168;
                        xStart4 = (7169 - c4);
                        xEnd4 = 7168;
                   end
               end  
               
               
               imagesphereresize4(yStart4:yEnd4,xStart4:xEnd4)=imagesphere4(1:r4,1:c4);
               newFile4 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize4,[foldername1 '/' newFile4]);
                n=n+1;
               wellname = (newFile4(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(7168,7168);
               currentNM = csvHandler.CellPosMatrix(wellname);
               currentNM(yStart4:yEnd4,xStart4:xEnd4) = NucleusM(y_7:y_8,x_7:x_8);
               csvHandler.CellPosMatrix(wellname) = currentNM;
               imagesphereresize4_small = imresize(imagesphereresize4, optionHandler.ScalingFactor);
               newFile4s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize4_small,[foldername1 '/' newFile4s]);
               clear imagesphere4;
               clear imagesphereresize4;
               clear imagesphereresize4_small;

               if(d1_5<MinimumDistance);
                   S51 = 1;
                   % Calculate line equation for d_12 
                   m_15 = (y_centroid(5)-y_centroid(1))/(x_centroid(5)-x_centroid(1));
                   b_15 = y_centroid(1)-(m_15*x_centroid(1));
                   % 2) Calculate Coordinate of d_12/2
                   x_c_15 = x_centroid(1) + (x_centroid(5)-x_centroid(1))/2;
                   y_c_15 = m_15*x_c_15 + b_15
                   % Calculate max distance in x from sphere 2 to 1:
                   x_15_1 = x_c_15;
                   x_15_2 = x_centroid(5) + (maxaddedlength);
                   %Introduce ratio of x/y (to rule out cases where two
                   %spheres have a very small differences in x or y but not
                   %in the other coordinate
                   xy_1 = abs(x_centroid(5)-x_centroid(1))/abs(y_centroid(5)-y_centroid(1))
                   yx_1 = abs(y_centroid(5)-y_centroid(1))/abs(x_centroid(5)-x_centroid(1))
                   if(xy_1 < minxy_ratio && d1_5 > mindist)
                       x_c_15 = 1000000;
                       y_c_15 = 1000000;
                       x_15_1 = x_centroid(5) - (maxaddedlength);
                       x_15_2 = x_centroid(5) + (maxaddedlength);
                       y_15_1 = y_centroid(5)-(maxaddedlength);
                       y_15_2 = y_centroid(5)+(maxaddedlength);
                   end  
                   if(yx_1 < minxy_ratio && d1_5 > mindist)
                       x_c_15 = 1000000;
                       y_c_15 = 1000000;
                       x_15_1 = x_centroid(5) - (maxaddedlength);
                       x_15_2 = x_centroid(5) + (maxaddedlength);
                       y_15_1 = y_centroid(5)-(maxaddedlength);
                       y_15_2 = y_centroid(5)+(maxaddedlength);
                   end  
                   if(y_c_15 > y_centroid(5))
                       y_15_1 = y_centroid(5)-(maxaddedlength)
                       y_15_2 = y_c_15
                   else 
                       y_15_1 = y_c_15
                       y_15_2 = y_centroid(5)+(maxaddedlength)
                   end    
               else
                   S51 = 0;
                   x_c_15 = 1000000;
                   y_c_15 = 1000000;
                   x_15_1 = x_centroid(5) - (maxaddedlength);
                   x_15_2 = x_centroid(5) + (maxaddedlength);
                   y_15_1 = y_centroid(5)-(maxaddedlength);
                   y_15_2 = y_centroid(5)+(maxaddedlength);
              end
               
              if(d2_5<MinimumDistance);
                   S52 = 1;
                   % Calculate line equation for d_12 
                   m_25 = (y_centroid(5)-y_centroid(2))/(x_centroid(5)-x_centroid(2));
                   b_25 = y_centroid(2)-(m_25*x_centroid(2));
                   % 2) Calculate Coordinate of d_12/2
                   x_c_25 = x_centroid(2) + (x_centroid(5)-x_centroid(2))/2;
                   y_c_25 = m_25*x_c_25 + b_25
                   % Calculate max distance in x from sphere 2 to 1:
                   x_25_1 = x_c_25;
                   x_25_2 = x_centroid(5) + (maxaddedlength);
                   xy_2 = abs(x_centroid(5)-x_centroid(2))/abs(y_centroid(5)-y_centroid(2))
                   yx_2 = abs(y_centroid(5)-y_centroid(2))/abs(x_centroid(5)-x_centroid(2))
                   if(xy_2 < minxy_ratio && d2_5 > mindist)
                       x_c_25 = 1000000;
                       y_c_25 = 1000000;
                       x_25_1 = x_centroid(5) - (maxaddedlength);
                       x_25_2 = x_centroid(5) + (maxaddedlength);
                       y_25_1 = y_centroid(5)-(maxaddedlength);
                       y_25_2 = y_centroid(5)+(maxaddedlength);
                   end  
                   if(yx_2 < minxy_ratio && d2_5 > mindist)
                       x_c_25 = 1000000;
                       y_c_25 = 1000000;
                       x_25_1 = x_centroid(5) - (maxaddedlength);
                       x_25_2 = x_centroid(5) + (maxaddedlength);
                       y_25_1 = y_centroid(5)-(maxaddedlength);
                       y_25_2 = y_centroid(5)+(maxaddedlength);
                   end
                   if(y_c_25 > y_centroid(5))
                       y_25_1 = y_centroid(5)-(maxaddedlength)
                       y_25_2 = y_c_25
                   else 
                       y_25_1 = y_c_25
                       y_25_2 = y_centroid(5)+(maxaddedlength)
                   end    
               else
                   S52 = 0;
                   x_c_25 = 1000000;
                   y_c_25 = 1000000;
                   x_25_1 = x_centroid(5) - (maxaddedlength);
                   x_25_2 = x_centroid(5) + (maxaddedlength);
                   y_25_1 = y_centroid(5)-(maxaddedlength);
                   y_25_2 = y_centroid(5)+(maxaddedlength);
              end
              
              if(d3_5<MinimumDistance);
                   S53 = 1;
                   % Calculate line equation for d_12 
                   m_35 = (y_centroid(5)-y_centroid(3))/(x_centroid(5)-x_centroid(3));
                   b_35 = y_centroid(3)-(m_35*x_centroid(3));
                   % 2) Calculate Coordinate of d_12/2
                   x_c_35 = x_centroid(3) + (x_centroid(5)-x_centroid(3))/2;
                   y_c_35 = m_35*x_c_35 + b_35
                   % Calculate max distance in x from sphere 2 to 1:
                   x_35_1 = x_c_35;
                   x_35_2 = x_centroid(5) - (maxaddedlength);
                   xy_3 = abs(x_centroid(3)-x_centroid(5))/abs(y_centroid(3)-y_centroid(5))
                   yx_3 = abs(y_centroid(3)-y_centroid(5))/abs(x_centroid(3)-x_centroid(5))
                   if(xy_3 < minxy_ratio && d3_5 > mindist)
                       x_c_35 = 1000000;
                       y_c_35 = 1000000;
                       x_35_1 = x_centroid(5) - (maxaddedlength);
                       x_35_2 = x_centroid(5) + (maxaddedlength);
                       y_35_1 = y_centroid(5)-(maxaddedlength);
                       y_35_2 = y_centroid(5)+(maxaddedlength);
                   end  
                   if(yx_3 < minxy_ratio && d3_5 > mindist)
                       x_c_35 = 1000000;
                       y_c_35 = 1000000;
                       x_35_1 = x_centroid(5) - (maxaddedlength);
                       x_35_2 = x_centroid(5) + (maxaddedlength);
                       y_35_1 = y_centroid(5)-(maxaddedlength);
                       y_35_2 = y_centroid(5)+(maxaddedlength);
                   end
                   if(y_c_35 > y_centroid(5))
                       y_35_1 = y_centroid(5)-(maxaddedlength)
                       y_35_2 = y_c_35
                   else 
                       y_35_1 = y_c_35
                       y_35_2 = y_centroid(5)+(maxaddedlength)
                   end    
               else
                   S53 = 0;
                   x_c_35 = 1000000;
                   y_c_35 = 1000000;
                   x_35_1 = x_centroid(5) - (maxaddedlength);
                   x_35_2 = x_centroid(5) + (maxaddedlength);
                   y_35_1 = y_centroid(5)-(maxaddedlength);
                   y_35_2 = y_centroid(5)+(maxaddedlength);
              end
               if(d4_5<MinimumDistance);
                   S54 = 1;
                   % Calculate line equation for d_12 
                   m_45 = (y_centroid(5)-y_centroid(4))/(x_centroid(5)-x_centroid(4));
                   b_45 = y_centroid(4)-(m_45*x_centroid(4));
                   % 2) Calculate Coordinate of d_12/2
                   x_c_45 = x_centroid(4) + (x_centroid(5)-x_centroid(4))/2;
                   y_c_45 = m_45*x_c_45 + b_45
                   % Calculate max distance in x from sphere 2 to 1:
                   x_45_1 = x_c_45;
                   x_45_2 = x_centroid(4) - (maxaddedlength);
                   xy_4 = abs(x_centroid(4)-x_centroid(5))/abs(y_centroid(4)-y_centroid(5))
                   yx_4 = abs(y_centroid(4)-y_centroid(5))/abs(x_centroid(4)-x_centroid(5))
                   if(xy_4 < minxy_ratio && d4_5 > mindist)
                       x_c_45 = 1000000;
                       y_c_45 = 1000000;
                       x_45_1 = x_centroid(5) - (maxaddedlength);
                       x_45_2 = x_centroid(5) + (maxaddedlength);
                       y_45_1 = y_centroid(5)-(maxaddedlength);
                       y_45_2 = y_centroid(5)+(maxaddedlength);
                   end  
                   if(yx_4 < minxy_ratio && d4_5 > mindist)
                       x_c_45 = 1000000;
                       y_c_45 = 1000000;
                       x_45_1 = x_centroid(5) - (maxaddedlength);
                       x_45_2 = x_centroid(5) + (maxaddedlength);
                       y_45_1 = y_centroid(5)-(maxaddedlength);
                       y_45_2 = y_centroid(5)+(maxaddedlength);
                   end
                   if(y_c_45 > y_centroid(5))
                       y_45_1 = y_centroid(5)-(maxaddedlength)
                       y_45_2 = y_c_45
                   else 
                       y_45_1 = y_c_45
                       y_45_2 = y_centroid(5)+(maxaddedlength)
                   end    
               else
                   S45 = 0;
                   x_c_45 = 1000000;
                   y_c_45 = 1000000;
                   x_45_1 = x_centroid(5) - (maxaddedlength);
                   x_45_2 = x_centroid(5) + (maxaddedlength);
                   y_45_1 = y_centroid(5)-(maxaddedlength);
                   y_45_2 = y_centroid(5)+(maxaddedlength);
               end
                              
               V_x9 = [x_c_15;x_c_25;x_c_35;x_c_45];
               min_x_9 = min(V_x9(1:4));
               if(min_x_9 ==1000000)
                   min_x_9 = x_centroid(5)-(maxaddedlength);
               end    
               min_x_10 = x_centroid(5) + (maxaddedlength);
               
               
                                             
               if(y_c_15 < 1000000 || y_c_25 < 1000000 || y_c_35 < 1000000 || y_c_45 < 1000000)
                   
                   if(y_c_15 - y_centroid(5)>0 && y_c_15 - y_centroid(5)<500000);
                       y_c_15_s = 1000000;
                   else
                       y_c_15_s = y_c_15; 
                   end 
                   if(y_c_25 - y_centroid(5)>0 && y_c_25 - y_centroid(5)<500000);
                       y_c_25_s = 1000000;
                   else
                       y_c_25_s = y_c_25;     
                   end
                   if(y_c_35 - y_centroid(5)>0 && y_c_35 - y_centroid(5)<500000);
                       y_c_35_s = 1000000;
                   else
                       y_c_35_s = y_c_35;     
                   end
                   if(y_c_45 - y_centroid(5)>0 && y_c_45 - y_centroid(5)<500000);
                       y_c_45_s = 1000000;
                   else
                       y_c_45_s = y_c_45;    
                   end
                   V_y9=[y_c_15_s;y_c_25_s;y_c_35_s;y_c_45_s];
                   min_y_9 = min(V_y9(1:4));
                   if(min_y_9 == 1000000)
                       min_y_9 = y_centroid(5)- (maxaddedlength);
                   end
                   
               else
                   min_y_9 = y_centroid(5)- (maxaddedlength);
               end
               
               if(y_c_15 < 1000000 || y_c_25 < 1000000 || y_c_35 < 1000000 || y_c_45 < 1000000)
                   if(y_c_15 - y_centroid(5)<0 && y_c_15 - y_centroid(5)<500000);
                       y_c_15_s = 1000000;
                   else
                       y_c_15_s = y_c_15;
                   end    
                   if(y_c_25 - y_centroid(5)<0 && y_c_25 - y_centroid(5)<500000);
                       y_c_25_s = 1000000;
                   else
                       y_c_25_s = y_c_25;
                   end
                   if(y_c_35 - y_centroid(5)<0 && y_c_35 - y_centroid(5)<500000);
                       y_c_35_s = 1000000;
                   else
                       y_c_35_s = y_c_35;
                   end
                   if(y_c_45 - y_centroid(5)<0 && y_c_45 - y_centroid(5)<500000);
                       y_c_45_s = 1000000;
                   else
                       y_c_45_s = y_c_45;
                   end
                   V_y10=[y_c_15_s;y_c_25_s;y_c_35_s;y_c_45_s];
                   min_y_10 = min(V_y10(1:4));
                   if(min_y_10 == 1000000)
                       min_y_10 = y_centroid(5)+ (maxaddedlength);
                   end    
                   
               else
                   min_y_10 = y_centroid(5)+ (maxaddedlength);
               end
                   
               if(min_x_9<1)
                   x_9 = 1
               else
                   x_9 = min_x_9;
               end
               
               if(min_x_10>10752)
                   x_10 = 10752
               else
                   x_10 = min_x_10;
               end
               
               if(min_y_9<1)
                   y_9 = 1
               else
                   y_9 = min_y_9;
               end
               
               if(min_y_10>13824)
                   y_10 = 13824
               else
                   y_10 = min_y_10;
               end
               
               diffXCentroidl = x_centroid(5) - x_9;
               diffXCentroidr = x_10 - x_centroid(5);
               diffYCentroidl = y_centroid(5) - y_9;
               diffYCentroidr = y_10 - y_centroid(5);
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
               
               imagesphere5 = imgnucleus(y_9:y_10,x_9:x_10);
               [r5 c5] = size(imagesphere5);
               imagesphereresize5 = zeros(7168,7168);
               imagesphereresize5 = uint8(imagesphereresize5);
               if(leftBoarder == 0)
                   if(lowerBoarder == 0)
                        
                        yStart5 = 1;
                        yEnd5 = r5;
                        xStart5 = 1;
                        xEnd5 = c5;
                   else
                        
                        yStart5 = (7169-r5);
                        yEnd5 = 7168;
                        xStart5 = 1;
                        xEnd5 = c5;
                   end 
               else
                   if(lowerBoarder == 0)
                       
                        yStart5 = 1;
                        yEnd5 = r5;
                        xStart5 = (7169 - c5);
                        xEnd5 = 7168;
                   else
                       
                        yStart5 = (7169-r5);
                        yEnd5 = 7168;
                        xStart5 = (7169 - c5);
                        xEnd5 = 7168;
                   end
               end  
               
               
               imagesphereresize5(yStart5:yEnd5,xStart5:xEnd5)=imagesphere5(1:r5,1:c5);
               newFile5 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NucleusBig'],'ignorecase');        
               imwrite(imagesphereresize5,[foldername1 '/' newFile5]);
                n=n+1;
               wellname = (newFile5(1:3));
               wellList{n} = wellname;
               csvHandler.CellPosMatrix(wellname) = sparse(7168,7168);
               currentNM = csvHandler.CellPosMatrix(wellname);
               currentNM(yStart5:yEnd5,xStart5:xEnd5) = NucleusM(y_9:y_10,x_9:x_10);
               csvHandler.CellPosMatrix(wellname) = currentNM;
               imagesphereresize5_small = imresize(imagesphereresize5, optionHandler.ScalingFactor);
               newFile5s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NucleusSmall'],'ignorecase');        
               imwrite(imagesphereresize5_small,[foldername1 '/' newFile5s]);
               clear imagesphere5;
               clear imagesphereresize5;
               clear imagesphereresize5_small;
               clear imgnucleus;
               
               
               if(Neuron_Channel > 0)
                   % Same for neurite channel
                    
                   imgneurite = imread([foldername '/' newFileNeu]);
                   %imwrite(imgneurite,[foldername2 '/' newFileNeu]);
                   imagesphere6 = imgneurite(y_1:y_2,x_1:x_2);
                   [r1 c1] = size(imagesphere6);
                   imagesphereresize6 = zeros(7168,7168);
                   imagesphereresize6 = uint8(imagesphereresize6);
                   imagesphereresize6(yStart1:yEnd1,xStart1:xEnd1)=imagesphere6(1:r1,1:c1);
                   newFile6 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize6,[foldername1 '/' newFile6]);
                   imagesphereresize6_small = imresize(imagesphereresize6, optionHandler.ScalingFactor);
                   newFile6s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize6_small,[foldername1 '/' newFile6s]);
                   clear imagesphere6;
                   clear imagesphereresize6;
                   clear imagesphereresize6_small;

                   imagesphere7 = imgneurite(y_3:y_4,x_3:x_4);
                   [r2 c2] = size(imagesphere7);
                   imagesphereresize7 = zeros(7168,7168);
                   imagesphereresize7 = uint8(imagesphereresize7);
                   imagesphereresize7(yStart2:yEnd2,xStart2:xEnd2)=imagesphere7(1:r2,1:c2);
                   newFile7 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize7,[foldername1 '/' newFile7]);
                   imagesphereresize7_small = imresize(imagesphereresize7, optionHandler.ScalingFactor);
                   newFile7s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize7_small,[foldername1 '/' newFile7s]);
                   clear imagesphere7;
                   clear imagesphereresize7;
                   clear imagesphereresize7_small;

                   imagesphere8 = imgneurite(y_5:y_6,x_5:x_6);
                   [r3 c3] = size(imagesphere8);
                   imagesphereresize8 = zeros(7168,7168);
                   imagesphereresize8 = uint8(imagesphereresize8);
                   imagesphereresize8(yStart3:yEnd3,xStart3:xEnd3)=imagesphere8(1:r3,1:c3);
                   newFile8 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize8,[foldername1 '/' newFile8]);
                   imagesphereresize8_small = imresize(imagesphereresize8, optionHandler.ScalingFactor);
                   newFile8s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize8_small,[foldername1 '/' newFile8s]);
                   clear imagesphere8;
                   clear imagesphereresize8;
                   clear imagesphereresize8_small;

                   imagesphere9 = imgneurite(y_7:y_8,x_7:x_8);
                   [r4 c4] = size(imagesphere9);
                   imagesphereresize9 = zeros(7168,7168);
                   imagesphereresize9 = uint8(imagesphereresize9);
                   imagesphereresize9(yStart4:yEnd4,xStart4:xEnd4)=imagesphere9(1:r4,1:c4);
                   newFile9 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize9,[foldername1 '/' newFile9]);
                   imagesphereresize9_small = imresize(imagesphereresize9, optionHandler.ScalingFactor);
                   newFile9s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize9_small,[foldername1 '/' newFile9s]);
                   clear imagesphere9;
                   clear imagesphereresize9;
                   clear imagesphereresize9_small;

                   imagesphere10 = imgneurite(y_9:y_10,x_9:x_10);
                   [r5 c5] = size(imagesphere10);
                   imagesphereresize10 = zeros(7168,7168);
                   imagesphereresize10 = uint8(imagesphereresize10);
                   imagesphereresize10(yStart5:yEnd5,xStart5:xEnd5)=imagesphere10(1:r5,1:c5);
                   newFile5 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NeuriteBig'],'ignorecase');        
                   imwrite(imagesphereresize10,[foldername1 '/' newFile5]);
                   imagesphereresize10_small = imresize(imagesphereresize10, optionHandler.ScalingFactor);
                   newFile10s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'NeuriteSmall'],'ignorecase');        
                   imwrite(imagesphereresize10_small,[foldername1 '/' newFile10s]);
                   clear imagesphere10;
                   clear imagesphereresize10;
                   clear imagesphereresize10_small;
                   clear imgneurite
               end
               
               if(Oligo_Channel > 0)
                   % Same for Oligos
                   imgoligo = imread([foldername '/' newFileOli]);
                   %imwrite(imgoligo,[foldername2 '/' newFileOli]);
                   imagesphere11 = imgoligo(y_1:y_2,x_1:x_2);
                   [r1 c1] = size(imagesphere11);
                   imagesphereresize11 = zeros(7168,7168);
                   imagesphereresize11 = uint8(imagesphereresize11);
                   imagesphereresize11(yStart1:yEnd1,xStart1:xEnd1)=imagesphere11(1:r1,1:c1);
                   newFile11 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize11,[foldername1 '/' newFile11]);
                   imagesphereresize11_small = imresize(imagesphereresize11, optionHandler.ScalingFactor);
                   newFile11s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize11_small,[foldername1 '/' newFile11s]);
                   clear imagesphere11;
                   clear imagesphereresize11;
                   clear imagesphereresize11_small;
                  
                   
                   
                   
                   imagesphere12 = imgoligo(y_3:y_4,x_3:x_4);
                   [r2 c2] = size(imagesphere12);
                   imagesphereresize12 = zeros(7168,7168);
                   imagesphereresize12 = uint8(imagesphereresize12);
                   imagesphereresize12(yStart2:yEnd2,xStart2:xEnd2)=imagesphere12(1:r2,1:c2);
                   newFile12 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize12,[foldername1 '/' newFile12]);
                   imagesphereresize12_small = imresize(imagesphereresize12, optionHandler.ScalingFactor);
                   newFile12s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize12_small,[foldername1 '/' newFile12s]);
                   clear imagesphere12;
                   clear imagesphereresize12;
                   clear imagesphereresize12_small;

                   imagesphere13 = imgoligo(y_5:y_6,x_5:x_6);
                   [r3 c3] = size(imagesphere13);
                   imagesphereresize13 = zeros(7168,7168);
                   imagesphereresize13 = uint8(imagesphereresize13);
                   imagesphereresize13(yStart3:yEnd3,xStart3:xEnd3)=imagesphere13(1:r3,1:c3);
                   newFile13 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize13,[foldername1 '/' newFile13]);
                   imagesphereresize13_small = imresize(imagesphereresize13, optionHandler.ScalingFactor);
                   newFile13s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize13_small,[foldername1 '/' newFile13s]);
                   clear imagesphere13;
                   clear imagesphereresize13;
                   clear imagesphereresize13_small;

                   imagesphere14 = imgoligo(y_7:y_8,x_7:x_8);
                   [r4 c4] = size(imagesphere14);
                   imagesphereresize14 = zeros(7168,7168);
                   imagesphereresize14 = uint8(imagesphereresize14);
                   imagesphereresize14(yStart4:yEnd4,xStart4:xEnd4)=imagesphere14(1:r4,1:c4);
                   newFile14 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize14,[foldername1 '/' newFile14]);
                   imagesphereresize14_small = imresize(imagesphereresize14, optionHandler.ScalingFactor);
                   newFile14s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize14_small,[foldername1 '/' newFile14s]);
                   clear imagesphere14;
                   clear imagesphereresize14;
                   clear imagesphereresize14_small;

                   imagesphere15 = imgoligo(y_9:y_10,x_9:x_10);
                   [r5 c5] = size(imagesphere15);
                   imagesphereresize15 = zeros(7168,7168);
                   imagesphereresize15 = uint8(imagesphereresize15);
                   imagesphereresize15(yStart5:yEnd5,xStart5:xEnd5)=imagesphere15(1:r5,1:c5);
                   newFile15 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'OligoBig'],'ignorecase');        
                   imwrite(imagesphereresize15,[foldername1 '/' newFile15]);
                   imagesphereresize15_small = imresize(imagesphereresize15, optionHandler.ScalingFactor);
                   newFile15s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'OligoSmall'],'ignorecase');        
                   imwrite(imagesphereresize15_small,[foldername1 '/' newFile15s]);
                   clear imagesphere15;
                   clear imagesphereresize15;
                   clear imagesphereresize15_small;
                   clear imgoligo;
               end
               
               if(Astro_Channel >0)
                   % Same for Astros
                   
                   
                   imgastro = imread([foldername '/' newFileAst]);
                   %imwrite(imgastro,[foldername2 '/' newFileAst]);
                   imagesphere16 = imgastro(y_1:y_2,x_1:x_2);
                   [r1 c1] = size(imagesphere16);
                   imagesphereresize16 = zeros(7168,7168);
                   imagesphereresize16 = uint8(imagesphereresize16);
                   imagesphereresize16(yStart1:yEnd1,xStart1:xEnd1)=imagesphere16(1:r1,1:c1);
                   newFile16 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize16,[foldername1 '/' newFile16]);
                   csvHandler.CellPosMatrix(wellname) = sparse(logical(zeros(7168,7168)));
                   imagesphereresize16_small = imresize(imagesphereresize16, optionHandler.ScalingFactor);
                   newFile16s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['E' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize16_small,[foldername1 '/' newFile16s]);
                   clear imagesphere16;
                   clear imagesphereresize16;
                   clear imagesphereresize16_small;

                   imagesphere17 = imgastro(y_3:y_4,x_3:x_4);
                   [r2 c2] = size(imagesphere17);
                   imagesphereresize17 = zeros(7168,7168);
                   imagesphereresize17 = uint8(imagesphereresize17);
                   imagesphereresize17(yStart2:yEnd2,xStart2:xEnd2)=imagesphere17(1:r2,1:c2);
                   newFile17 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize17,[foldername1 '/' newFile17]);
                   csvHandler.CellPosMatrix(wellname) = sparse(logical(zeros(7168,7168)));
                   imagesphereresize17_small = imresize(imagesphereresize17, optionHandler.ScalingFactor);
                   newFile17s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['F' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize17_small,[foldername1 '/' newFile17s]);
                   clear imagesphere17;
                   clear imagesphereresize17;
                   clear imagesphereresize17_small;


                   imagesphere18 = imgastro(y_5:y_6,x_5:x_6);
                   [r3 c3] = size(imagesphere18);
                   imagesphereresize18 = zeros(7168,7168);
                   imagesphereresize18 = uint8(imagesphereresize18);
                   imagesphereresize18(yStart3:yEnd3,xStart3:xEnd3)=imagesphere18(1:r3,1:c3);
                   newFile18 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize18,[foldername1 '/' newFile18]);
                   csvHandler.CellPosMatrix(wellname) = sparse(logical(zeros(7168,7168)));
                   imagesphereresize18_small = imresize(imagesphereresize18, optionHandler.ScalingFactor);
                   newFile18s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['G' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize18_small,[foldername1 '/' newFile18s]);
                   clear imagesphere18;
                   clear imagesphereresize18;
                   clear imagesphereresize18_small;

                   imagesphere19 = imgastro(y_7:y_8,x_7:x_8);
                   [r4 c4] = size(imagesphere19);
                   imagesphereresize19 = zeros(7168,7168);
                   imagesphereresize19 = uint8(imagesphereresize19);
                   imagesphereresize19(yStart4:yEnd4,xStart4:xEnd4)=imagesphere19(1:r4,1:c4);
                   newFile19 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize19,[foldername1 '/' newFile19]);
                   csvHandler.CellPosMatrix(wellname) = sparse(logical(zeros(7168,7168)));
                   imagesphereresize19_small = imresize(imagesphereresize19, optionHandler.ScalingFactor);
                   newFile19s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['H' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize19_small,[foldername1 '/' newFile19s]);
                   clear imagesphere19;
                   clear imagesphereresize19;
                   clear imagesphereresize19_small;

                   imagesphere20 = imgastro(y_9:y_10,x_9:x_10);
                   [r5 c5] = size(imagesphere20);
                   imagesphereresize20 = zeros(7168,7168);
                   imagesphereresize20 = uint8(imagesphereresize20);
                   imagesphereresize20(yStart5:yEnd5,xStart5:xEnd5)=imagesphere20(1:r5,1:c5);
                   newFile20 = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'AstroBig'],'ignorecase');        
                   imwrite(imagesphereresize20,[foldername1 '/' newFile20]);
                   csvHandler.CellPosMatrix(wellname) = sparse(logical(zeros(7168,7168)));
                   imagesphereresize20_small = imresize(imagesphereresize20, optionHandler.ScalingFactor);
                   newFile20s = regexprep(allfiles(i).name,allfiles(i).name(1:13),['I' wellind 'AstroSmall'],'ignorecase');        
                   imwrite(imagesphereresize20_small,[foldername1 '/' newFile20s]);
                   clear imagesphere20;
                   clear imagesphereresize20;
                   clear imagesphereresize20_small;
                   clear imgastro
               end
                
        
               %delete(allfiles(i).name);
                %if(Neuron_Channel > 0)
                %delete(newFileNeu);
                %end
                %if(Oligo_Channel > 0)
                %delete(newFileOli);
                %end
                %if(Astro_Channel > 0)
                %delete(newFileAst);
                %end
                          
            end    
            
               
        
    
        
   %csvHandler.CellPosMatrix(allfiles(i).name) = 0;   
   
end
        
     %CopyImage = imread([foldername '/' allfiles(i).name]);
     %imwrite(CopyImage,[foldername2 '/' allfiles(i).name]);
     %delete(allfiles(i).name);  
     %Delete original Matrix
     
end

  
        


end

