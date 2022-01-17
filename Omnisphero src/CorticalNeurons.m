%Copyright (C) 2017-2021 Martin Schmuck

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

% This function is designed to pre process images from stacks to oabtain
% binary images of neurons and to isolate cell soma and to identify the centroids of the cell soma

% Afer function are the variables obtained from Omnisphero in [], after =
% CorticalNeurons is variable which is returned to Omnisphero ()
function [wellList, foldername, ScalingFactor, csvHandler,   Th] = CorticalNeurons (csvHandler) 
% This is a quick script to convert 16-bit images of cortical neurons into
% binary image, to assess their center point and to generate a
% 'NucleusImage'.
n = 0;
% This variable describes the minimal overlapp the cell soma has to have in
% a quadratic box in order to further expand the box and therby the area of
% the cell soma
RatioFill = 0.7;
% This is the maximum allowed size of the box sourronding the cell soma
Maxlengthsquare = 100;
% Minimus required size of a structure to be considered as a neuron:
FinalNeuronSize = 4000;
% This factor is used to identify large bubbles, which will create huge
% areas with lower background compared to the rest of the image
MedianFactorMask = 0.3;
% This is the maximal area which is allowed to be filled by the imfill
% function used to fill small holes within neuronal structures
HoleArea = 50;
% This area describes the minimal area of foreground pixles required to not
% delete the image
MinReqArea = 700;
% This variable is used to remove neuronal structures, which should be
% bigger than this values, from the opened image for background removal 
HardBackgroundArea = 1000;
% Areas filtered from the edge processed image using canny method
AreaCanny = 100;
% Areas filtered from the edge processed image using prewitt method
AreaPrewitt = 10;
% Areas filtered from the edge processed image using sobel method
AreaSobel = 10;
% Areas filtered from the edge processed image using Roberts method
AreaRoberts = 5;
% Areas filtered from the edge processed image using Log method
AreaLogHigh = 100;
AreaLogLow = 50;
% Minmal ration of added feature size to area of surrounding
% box
RatioAreaValue = 3;
%Minal length per endpoint of added features connecting mutiple sites
LengthPerEndpointInner = 12;
%Minal length per endpoint of added features connecting one
%site
LengthPerEndpointOuter = 20;
% This variable gives the aspect ratio of the sourroundingbox of the cell
% soma!
ShapeCellSoma= 1.2;
%Get Image folder
foldername = uigetdir;
% Creates directory within the parent folder
allfiles = dir(foldername);
% Creates new folder within parent folder
mkdir(foldername,'ConvertedCellomics');
% 
%Create function to seperate image stack
%1) Load Stack
 nI = 0;
% Gets file name in folder (Be carefull, this will pick up the first
% available image, so do not place multipe stacks in one folder!
 fname = [foldername '/' allfiles(3).name];
% Get info, like number of images in the stack
 info = imfinfo(fname);
% Extracts number of images from stack!
 num_images = numel(info);
 for(w=1:num_images)
    % Numeration of images in stack 
    nI = nI + 1;
    % Reads first image from stack
    subImage = imread(fname, w, 'Info', info);
    % Binning to 8 bit
    subImage = double(subImage)./4095;
    subImage = uint8(subImage.*255);
    % Naming images in a onward fashion (Max images now 99!)
   if(w>9)
   a{1} = 'A';
   else
   a{1} = 'A0';
   end
   a{2} = int2str(nI);
   a{3} = '-';
   a{4} = allfiles(3).name;
   newfileI = [a{1} a{2} a{3} a{4}];
   newFileI = regexprep(newfileI, '.tif','.png','ignorecase');
   imwrite(subImage,[foldername '/' newFileI]);
end    


% Create folder path: Subfolder of parent directory
% 'foldername'
foldername1 = [foldername '/ConvertedCellomics'];
% Creeate directory in parent foldername
allfiles = dir(foldername);
% Create an empty welllist
wellList = cell(0);
% Iterate over all images in parent folder
for(i=1:numel(allfiles))
    %Check if file ends with .png. This will ignore the original stack,
    %which is in tiff.
    ind = strfind([foldername '/' allfiles(i).name],'.png');
    % Check of image i is a png or not
    if(numel(ind) > 0)%
        % Count number of images
         n = n + 1;
        %Load image n from parent folder
        imgNeurite = imread([foldername '/' allfiles(i).name]);    
        % Create a nonflat structuring element with ball shape for dilation
        % of a gray scale image
        SE = strel('ball',10,10);
        % Performs a morphological opening on the gray scale image. This
        % averages pixel. This function is used for morphological noise removal 
        backgroundball = imopen(imgNeurite,SE);
        % Convert image into 8 - bit
        backgroundball = uint8(backgroundball);
        
        
        
        % This function is exculde neuronal structures form the opened
        % image, to circumvent deletion of neuron pixels
        level = isodata(backgroundball);
        % Check if isodate can be computed. If not use triangle method
        if(isnan(level) == 1)
            % Create histogram of image
            [lehisto x] = imhist(backgroundball);
            % Calculate level using the triangel method
            level = triangle_th(lehisto,256);
        end 
        % Creates a binary image of opened image, with a high value for the
        % fixed background. If the threshold would exceed 1, we set it 0.9
        if(level*1.5 >1)
            level = 0.9;
            HardBackgroundBall = im2bw(backgroundball,level);
        else
            HardBackgroundBall = im2bw(backgroundball,level*1.5);
        end    
        % Creates a binary image of opened image, with a much higher value for the
        % fixed background
        if(level*5 >1)
            level = 0.9;
            EvenHarderBackgroundBall = im2bw(backgroundball,level)
        else
            EvenHarderBackgroundBall = im2bw(backgroundball,level*5);
        end
        %Sort out too big particles which would correspond to neuronal
        %structures
        % This gets all areas larger than 1000
        Blobs = bwareaopen(HardBackgroundBall,HardBackgroundArea);
        % Repeat the same for the EvenHarderBackgroundBall
        BlobsHard = bwareaopen(EvenHarderBackgroundBall,HardBackgroundArea);
        % First substract small areas from less stringent thresholded
        % binary image, resulting in only areas smaller than 100. Then
        % adding the areas of the binary image with the higher threshold.
        % This areas normally contribute to cell debris or 
        HardBackgroundBall = HardBackgroundBall - Blobs + BlobsHard;
        % Multiply image with 255 (saturate all pixels)
        HardBackgroundBall = 255*HardBackgroundBall;
        % Convert image to 8-bit
        HardBackgroundBall = uint8(HardBackgroundBall);
        % Substract artifacts from original opened image
        backgroundball = backgroundball - HardBackgroundBall;
        % Now substract noise from original grayscale image
        imgNeuriteBall = imgNeurite - backgroundball;
        % Now create binary image from corrected gray scale image using the
        % calculated threshold *1
        binaryNeuriteFill = im2bw(imgNeuriteBall,level);
        
        % Function to delete Bright artificat spots (cell debris or antibody participation) from binary image 
        %Remove small particles
        Artifact = bwareaopen(binaryNeuriteFill,100);
        %Everything bigger than 1000 should be a neuron structure!
        NonArtifact = bwareaopen(binaryNeuriteFill,1000);
        % Retain only particles smaller than 1000. This could be neurite
        % parts or artifacts
        Artifact = Artifact - NonArtifact;
        % COnvert image to 8 bit
        Artifact = uint8(Artifact);
        % Convert image to logical
        Artifact = logical(Artifact);
        %Calculate all major and minor axis of objects (Neurite structures should be longer than high, thus giving a high axis ratio from minor to major)
        Amj = regionprops(Artifact,'majorAxis');
        Ami = regionprops(Artifact,'minorAxis');
        % Returns all minaor and major axis of remianing structures
        An = size(Amj);
        An = An(1);
        % Label remianing structures structure
        ArtefactLabel = bwlabel(Artifact);
        % Create empty image later filled with removable artifacts
        ArtefactRemove = zeros(size(binaryNeuriteFill));
        % This loop should identify all round structures and fill them into
        % the removal image
        for(t=1:An)
            % Create image containing one particle from ArtefactLabel
            ArtefactOne = ArtefactLabel;
            % Index all pixels within the image ArtefactOne with a value
            % equal to t
            ii=ArtefactOne==t;
            % Set all these pixel values to 255
            ArtefactOne(ii)=255;
            % Substract 250, resulting in only pixel of index to be
            % positive
            ArtefactOne = ArtefactOne -250;
            % Convert image to 8-bit to get rid of negative pixle values
            ArtefactOne = uint8(ArtefactOne);
            % Convert image back to logical to have max pixle value to one
            ArtefactOne = logical(ArtefactOne);
            % Get minor and major axis of single object
            AmjN = Amj(t).MajorAxisLength(1);
            AmiN = Ami(t).MinorAxisLength(1);
            % Calculate the ratio of major and minor axis
            ArtifactAx = AmjN/AmiN;
            % Create matrix conating all ratios (helpful to identify
            % usefull vlaue for filtering)
            VAAX(t,1) = ArtifactAx;
            % Search for particles with a ratio of axis bigger than 2
            if(ArtifactAx>2)
                ArtefactOne = zeros(size(ArtefactRemove));
            end
            % Add particle or empty image to final image used for
            % substraction of aritacts
            ArtefactRemove = ArtefactRemove + ArtefactOne;
        end   
        % Remove artifacts from original primary image
        binaryNeuriteFill = binaryNeuriteFill -     ArtefactRemove;
        
        
%         This function is designed to decide, whether an image contains a
%         neuron or not. If there are less than a user defined number of
%         pixles presnet in the binary image, the image is excluded from
%         further analysis and treated as an empty image
        if(sum(sum(bwareaopen(binaryNeuriteFill,MinReqArea)))<MinReqArea)
            % Create empty image
            binaryNeuriteFill = zeros(size(imgNeurite));
            % Create empty opened image
            imgNeuriteBall = zeros(size(imgNeurite));
        end  
        
        
        
        % This function now performs the actual thresholding of the
        % corrected gray scale image. The basic idea is to use different edge detector algorithms, namely, canny, prewitt, sobel, Roberts and log, and filter the obtained imgaes for lines exceeding user defined thresholds, in order to remove non neuronal structures. 
        % Perform edge Canny  
        
        E1 = edge(imgNeuriteBall,'canny');
        % Filter image E1 for particles bigger AreaCanny
        E1old = bwareaopen(E1,AreaCanny);
        % Perform edge Prewitt
        E2 = edge(imgNeuriteBall,'prewitt');
        % Filter image E2 for particles bigger AreaPrewitt
        E2old = bwareaopen(E2,AreaPrewitt);
        % Perform edge Sobel
        E3 = edge(imgNeuriteBall,'sobel');
        % Filter image E3 for particles bigger AreaSobel
        E3old = bwareaopen(E3,AreaSobel);
        % Perform edge Roberts
        E4 = edge(imgNeuriteBall,'Roberts');
        % Filter image E4 for particles bigger AreaRoberts
        E4old = bwareaopen(E4,AreaRoberts);
        % Perform edge log
        E5 = edge(imgNeuriteBall,'log');
        % Filter image E5 for particles bigger AreaLogHigh
        E5old = bwareaopen(E5,AreaLogHigh);
        % Filter image E5 for particles bigger AreaLogLow
        E5 = bwareaopen(E5,AreaLogLow);
        
        % We use will again reomove identified artifacts from above, but
        % this time growing them before
        ArtefactRemove = bwmorph(ArtefactRemove, 'thicken', 5);
        % Growing leads to artifical gaps, which are connected with
        % 'bridge'
        ArtefactRemove = bwmorph( ArtefactRemove,'bridge');
        % Now we add all filtered edge detector images
         EdgeBinaryold = E1old + E2old + E3old + E4old + E5old; 
         % Convert image to logical
         EdgeBinaryold = logical(EdgeBinaryold);
         % Add binary image to the binary Edge image
         EdgeBinaryold = EdgeBinaryold + binaryNeuriteFill;
         % Convert image to logical
         EdgeBinaryold = logical(EdgeBinaryold);
         % Substract Artifcats from image
         EdgeBinaryold = EdgeBinaryold - ArtefactRemove;
         % Convert image to 8-bit
         EdgeBinaryold = uint8(EdgeBinaryold);
         % Convert image to logical
         EdgeBinaryold = logical(EdgeBinaryold);
         % Perform bridge to connect areas in close proximity
         EdgeBinaryold = bwmorph(EdgeBinaryold,'bridge');
         % Sort for particles bigger than FinalNeuronSize
         EdgeBinaryold = bwareaopen(EdgeBinaryold,FinalNeuronSize);
         % Create image in which all holes are filled
         BW = imfill(EdgeBinaryold,'holes');
         % Substract original EdgeBinary from filled image to obtain only
         % areas of holes
         BW2 = BW - EdgeBinaryold;
         % Perform image opening
         BW2 = bwmorph(BW2,'open');
         % Only filter Holes bigger than HoleArea
         Holes = bwareaopen(BW2,HoleArea);
         % Substract Holes from BW to obtain only the small holes
         EdgeBinaryold = BW - Holes;
         % Perform brdige again to connect objects
         EdgeBinaryold = bwmorph(EdgeBinaryold,'bridge');
         % Repeat a second time to fill new holes
         BW = imfill(EdgeBinaryold,'holes');
         BW2 = BW - EdgeBinaryold;
         BW2 = bwmorph(BW2,'open');
         Holes = bwareaopen(BW2,HoleArea);
         EdgeBinaryold = BW - Holes;
         % Perform majority to fill small one pixel missing holes
         EdgeBinaryold = bwmorph(EdgeBinaryold,'majority');
         % Sort again for big particles
         EdgeBinaryold = bwareaopen(EdgeBinaryold,2000);
        
        % Now create less stringent edge detector (No area filtering)
        
        % Sum all edge dector images
        EdgeBinary = E1+E2+E3+E4+E5;
        % Convert image to logical
        EdgeBinary = logical(EdgeBinary);
        % Add binary neurite image
        EdgeBinary = EdgeBinary + binaryNeuriteFill;
        % Convert image to logical
        EdgeBinary = logical(EdgeBinary);
        % Remove artifacts
        EdgeBinary = EdgeBinary - ArtefactRemove;
        % Convert image to 8-bit
        EdgeBinary = uint8(EdgeBinary);
        % Convert image to logical
        EdgeBinary = logical(EdgeBinary);
        % Perform bridge to connect areas
        EdgeBinary = bwmorph(EdgeBinary,'bridge');
        % Create image in which all holes are filled 
        BW = imfill(EdgeBinary,'holes');
        % Substract original EdgeBinary from filled image to obtain only
        % areas of holes
        BW2 = BW - EdgeBinary;
        % Perform image opening
        BW2 = bwmorph(BW2,'open');
        % Only filter Holes bigger than HoleArea
        Holes = bwareaopen(BW2,HoleArea);
        % Substract Holes from BW to obtain only the small holes
        EdgeBinary = BW - Holes;
        % Perform bridge
        EdgeBinary = bwmorph(EdgeBinary,'bridge');
        % Repeat a second time to fill new holes
        BW = imfill(EdgeBinary,'holes');
        BW2 = BW - EdgeBinary;
        BW2 = bwmorph(BW2,'open');
        Holes = bwareaopen(BW2,HoleArea);
        EdgeBinary = BW - Holes;
        % Perform majority to fill small one pixel missing holes
        EdgeBinary = bwmorph(EdgeBinary,'majority');
        % Only filter very small particles
        EdgeBinary = bwareaopen(EdgeBinary,10);
        % Create a copy of the original EdgeBinary before further
        % processing
        EdgeOriginal = EdgeBinary;
        
        
        % Gap filling algorithm:
        %1) Find all gaps!
        %1) dilatate image, so that gaps are closed
        SE = strel('disk',5);
        Gap = imdilate(EdgeBinary,SE);
        %2) Then erode again:
        Gap = imerode(Gap,SE);
        % Now get difference resulting in the actual 'gaps'
        Gap = Gap - EdgeBinary;
        % Convert image to 8-bit
        Gap = uint8(Gap);
        % Convert image to logical
        Gap = logical(Gap);
        % Dilatate Gaps
        SE = strel('disk',4);
        GapBig = imdilate(Gap,SE);
        % Identify big gaps, which are usually not real connections
        GapBig = bwareaopen(GapBig,200);
        % Substract big gaps to only obtain small true gaps
        Gap = Gap - GapBig;
        % Convert image to 8-bit
        Gap = uint8(Gap);
        % Convert image to logical
        Gap = logical(Gap);
        % Remove gaps attached to the image border
        Gap = imclearborder(Gap);
        % We should now further reduce reduntant gaps, by adding the
        % dialted image as a mask and only keep gaps within this mask
        % Imdialte image again, but this time even more to really close all
        % gaps
        SE2 = strel('disk',6);
        GapCorrected = imdilate(EdgeBinary,SE2);
        % Filter for big areas 
        GapCorrected = bwareaopen(GapCorrected,10000);
        % Add gaps to this mask. Gaps overlapping with the mask will have a
        % value of 2
        GapCorrected = GapCorrected + Gap;
        % Substraction of 1 results in gaps solely located on the mask
        GapCorrected = GapCorrected -1;
        % Convert image to 8-bit
        GapCorrected = uint8(GapCorrected);
        % Convert image to logical
        GapCorrected = logical(GapCorrected);
        % Rename GapCorrected to Gap 
        Gap = GapCorrected;
        %Now check expand each area for one pixel:
        GapFat = bwmorph(Gap,'thicken',1);
        % Now test for each area, wheather it is overlapping with two areas
        % from EdgeBinary which is only the case for a gap ar a pit:
        % Label all Gaps original
        GapLabel = bwlabel(Gap);
        % Label dilated Gaps
        GapLabelFat = bwlabel(GapFat);
        % Determine area with highest label
        sG = max(max(GapLabel));
        % Create image for GapAdd
        GapAdd = zeros(size(EdgeBinary));
        for(p=1:sG)
            %Extract one of the gap areas
            GapKeep = GapLabel;
            % Index all pixles in GapKeep equal to p
            ii = GapKeep == p;
            % Since we could have more than 255 areas, we create a value
            % depending on sG
            ValueAdd = (100*sG);
            % Set all indexed pixles to ValueAdd
            GapKeep(ii) = ValueAdd;
            ValueSubstract = ValueAdd -5; 
            % Substract ValueSubstract form Gap Keep to only have positve
            % pixel values for indexed pixels
            GapKeep = GapKeep -ValueSubstract;
            % Convert image to 8-bit
            GapKeep = uint8(GapKeep);
             % Convert image to logical
            GapKeep = logical(GapKeep);
            % Repeat for the dilated Gap imahe
            GapOne = GapLabelFat;
            ii = GapOne == p;
            GapOne(ii) = ValueAdd;
            GapOne = GapOne -ValueSubstract;
            GapOne = uint8(GapOne);
            GapOne = logical(GapOne);
            %Now add labeled EdgeBinary and check if there are three
            %or more different pixel values (means it connects two areas):
            % Add GapOne*150 to have a high value which can be substracted
            ConnectOne = (GapOne*150) + EdgeBinary;
            % Substract 149 to only keep pixels from GapOne, with values of
            % 1 (no overlapp) and 2 (overlapp with EdgeBinary)
            ConnectOne = ConnectOne - 149;
            % Convert image to 8-bit to get rid of negative values
            ConnectOne = uint8(ConnectOne);
            % Substract one to only keep pixell of overlapp
            ConnectOne = ConnectOne - 1;
            % Convert image to 8-bit to get rid of negative values
            ConnectOne = uint8(ConnectOne);
            % Perform bridge to connecte very near overlapp areas
            ConnectOne = bwmorph(ConnectOne,'bridge');
            % Look for number of areas within the Connected image
            AreaOne = regionprops(ConnectOne,'Area');
            AreaOne = size(AreaOne);
            AreaOne = AreaOne(1);
            % If the number of Areas is bigger than 1, the original gap is
            % kept
            if(AreaOne>1)
                if(sum(sum(ConnectOne))<21)
                GapAdd = GapAdd + GapKeep;
                end
            end
        end    
        % Converti image to 8-bit
        EdgeBinary = uint8(EdgeBinary);
        %Convert image to logical
        EdgeBinary = logical(EdgeBinary);
        % Add GapAdd to EdgeBinary to close Gaps
        EdgeBinary = EdgeBinary + GapAdd;
        % Converti image to 8-bit
        EdgeBinary = uint8(EdgeBinary);
        %Convert image to logical
        EdgeBinary = logical(EdgeBinary);
        % In order to remove artifical gap fills within loops of neurites,
        % we substract now all areas within neuronal loops from the original EdgeOriginal.
        % We perform imfill on EdgeOriginal
        HolesB = imfill(EdgeOriginal,'holes');
        %Seperate Areas by substracting EdgeOriginal
        HolesB = HolesB - EdgeOriginal;
        % Substract HolesB from EdgeBinary to eleimnate artifical gap fills
        % within loops
        EdgeBinary = EdgeBinary - HolesB;
        % Converti image to 8-bit
        EdgeBinary = uint8(EdgeBinary);
        % Converti image to logical
        EdgeBinary = logical(EdgeBinary);
        EdgeBinary = bwareaopen(EdgeBinary,FinalNeuronSize);
        
        % Now check whether added part to strict EdgeBinary is artifact or
        % real neurite extension. Idea here would be to skeletonize
        % additional parts and check for number of branches or endpoints:
        % Convert image to logical
        EdgeBinaryold = logical(EdgeBinaryold);
        % Substract EdgeBinary from EdgeBinaryold, to only obtain added
        % features
        AddedFeature = EdgeBinary - EdgeBinaryold;
        % Converti image to 8-bit
        AddedFeature = uint8(AddedFeature);
        % Convert image to logical
        AddedFeature = logical(AddedFeature);
        % Sort for features bigger than 50 (Everything else ist propably
        % not really enhacing anything
        AddedFeature = bwareaopen(AddedFeature,50);
        % Create Skeleton of added features
        AddedFeatureSkel  = bwmorph(AddedFeature,'thin',inf);
        % Label all added feature areas
        AddedFeatureLabel = bwlabel(AddedFeature);
        % Label all skeleton of added features
        AddedFeatureLabelSkel = bwlabel(AddedFeatureSkel);
        % Assess maximal label within added features image
        sAdd = max(max(AddedFeatureLabel));
        % Create empty image for features added which should be kept
        AddFeatureKeep = zeros(size(EdgeBinary));
        % This variable will be used down to create a matrix with different
        % morphological features of the skeletons
        VLP = 0;
        % Sort out features which ......
        for(u=1:sAdd)
            AddFeatureOne = AddedFeatureLabel;
            % index all pixels with value u in AddFeatureOne
            ii = AddFeatureOne == u;
            AddN = (10*sAdd);
            % Set all indexed pixel values in AddFeatureOne to AddN
            AddFeatureOne(ii) = AddN;
            SubN = AddN - 5;
            % Substract SubN from AddFeatureOne, to maintain only positive
            % values for pixel indexed
            AddFeatureOne = AddFeatureOne - SubN;
            % Convert image to 8-bit
            AddFeatureOne = uint8(AddFeatureOne);
            % Convert image to logical
            AddFeatureOne = logical(AddFeatureOne);
            % Dilation of AddFeatureOne by one
            AddFeatureOneFat = bwmorph(AddFeatureOne,'thicken',1);
            % Repeat the same for AddFeatureSkelOne
            AddFeatureSkelOne = AddedFeatureLabelSkel;
            AddN = (10*sAdd);
            ii = AddFeatureSkelOne == u;
            AddFeatureSkelOne(ii) = AddN;
            SubN = AddN - 5;
            AddFeatureSkelOne = AddFeatureSkelOne - SubN;
            AddFeatureSkelOne = uint8(AddFeatureSkelOne);
            AddFeatureSkelOne = logical(AddFeatureSkelOne);
            % Identify endpoints in skeletons
            EndOne = bwmorph(AddFeatureSkelOne,'endpoints');
            % Get number of endpoints
            EndOne = sum(sum(EndOne));
            % Get length of skeleton
            LengthOne = sum(sum(AddFeatureSkelOne));
            % Caluclate average length per endpoint (Usually cell debris
            % will result in high number of endpoints, resulting in a low
            % length per endpoint
            LengthPerEndpoint = LengthOne/EndOne;
            % Get minimal dimensions of a box sourrounding the added
            % feature
            OneBox = regionprops(AddFeatureOne,'BoundingBox');
            % Extract width of the box
            Width = OneBox(1).BoundingBox(3);
            % Extract length of the box
            Length = OneBox(1).BoundingBox(4);
            % Calulate area of the box
            AreaBox = Width*Length;
            % Calculate area of added feature
            AreaOne = sum(sum(AddFeatureOne));
            % Calculate ratio of area box to box sourrounding area
            RatioArea = AreaBox/AreaOne;
            % Write obtained values in matrix VLP
            % This variable is taking into account the raio of branches to
            % length, the ratio of box fill and the total area of the
            % object:
            WeightFactor = (RatioArea*LengthPerEndpoint)/AreaOne;
            %Check if AddFeatureOne is an inner part of EdgeBinary by
            %overlapping the dilated added feature to EdgeBinaryold
            Overlapp = AddFeatureOneFat + EdgeBinaryold;
            % Substract 1 to only keep overlay images
            Overlapp = Overlapp -1;
            % Converti image to 8bit
            Overlapp = uint8(Overlapp);
            % Convert image to logical
            Overlapp = logical(Overlapp);
            % Get areas within Overlapp
            sO = regionprops(Overlapp,'Area');
            % Get number of areas
            sO = size(sO);
            sO = sO(1);
            VLP (u,1) = RatioArea;
            VLP (u,2) = AreaOne;
            VLP (u,3) = LengthPerEndpoint;
            VLP (u,4) = WeightFactor;
            VLP (u,5) = sO;
            % If the added feature is within the inner part of EdgeBinaryold (sO>1), more endpoints area allowed, since the structure might combine different lose ends. If the added feature is non an interal addition, less endpoints are allowed, since added neurite parts should be mostly linear            
            if(sO >1)
                if(WeightFactor>0.07)
                  % we now need to figure, whhether the added Feature is closing a hole.
                  TestFill = imfill(EdgeBinaryold,'holes');
                  TestFill = TestFill - EdgeBinaryold;
                  TestFill = logical(TestFill);
                  TestFill = TestFill + AddFeatureOne;
                  if(max(max(TestFill))==1)
                  AddFeatureKeep = AddFeatureKeep +   AddFeatureOne;
                  end
                end 
            else
                if(WeightFactor>0.1)
                  TestFill = imfill(EdgeBinaryold,'holes');
                  TestFill = TestFill - EdgeBinaryold;
                  TestFill = logical(TestFill);
                  TestFill = TestFill + AddFeatureOne;
                  if(max(max(TestFill))==1)
                  AddFeatureKeep = AddFeatureKeep +   AddFeatureOne;
                  end
                end   
            end    
        end   
        % Add Only features fullfilling above criteria to EdgeBinaryold
        EdgeBinary = EdgeBinaryold + AddFeatureKeep;
        % Now we will bridge to connect big areas! But onlt if this
        % combines two areas!
        TestBridge = bwmorph(EdgeBinary,'bridge');
        %This will will give us all the bridge points!
        EdgeBinary = logical(EdgeBinary);
        BridgePoints = TestBridge - EdgeBinary;
        BridgePointsLabel = bwlabel(BridgePoints);
        NumberofBridgePoints = max(max(BridgePointsLabel));
        % Now we check for each single bridge Objects whether it reduces the number of binary particles     
        for(d=1:NumberofBridgePoints)
            BinaryRepair = BridgePointsLabel;
            ii = BinaryRepair == d;
            BinaryRepair(ii) = 255;
            BinaryRepair = BinaryRepair - 250;
            BinaryRepair = uint8(BinaryRepair);
            BinaryRepair = logical(BinaryRepair);
            EdgeBinaryLabel = bwlabel(EdgeBinary);
            NumberofParticels = max(max(EdgeBinaryLabel ));
            %Now we test if addition of the bridge object will actually
            %reduce the number oif binary particles of Neurons
            BinaryTest = EdgeBinary + BinaryRepair;
            BinaryTest = bwlabel(BinaryTest);
            BinaryTest = max(max(BinaryTest));
            if(BinaryTest<NumberofParticels)
                EdgeBinary = EdgeBinary + BinaryRepair ;
            end
        end    
        
        
        
        
        %Bubble Detector: Basic idea would be to identify dark areas,
        %define their rim and substract it from EdgeBinary, and thereby
        %delete bright rims between bubble and neuron
        % Converts image to single
        imgNeurite = single(imgNeurite);
        % Calulates median
        Median = median(median(imgNeurite));
        % Substract very low fixed threshold based on median (all non dark
        % areas of bubbles should be picked up!). The resulting binary
        % image will be used as a mask in which only neurons can be
        % searched for
        Mask = imgNeurite - Median*MedianFactorMask;
        % Convert image to 8-bit
        Mask = uint8(Mask);
        % Convert image to 8-bit
        imgNeurite = uint8(imgNeurite);
        % Index all pixel from Mask which are equal to 0. This are pixels
        % located within the bubble
        ii = Mask == 0;
        % Set indexed pixels to 255
        Mask(ii) = 255;
        % Subtract 254 to only have positve pixel values for indexed pixles
        Mask = Mask - 254;
        % Convert image to 8-bit
        Mask = uint8(Mask);
        % Convert image to logical
        Mask = logical(Mask);
        % Fill holes in Mask
        Mask = imfill(Mask,'holes');
        % Sort only for very big regions (bubbles have to be big)
        Mask = bwareaopen(Mask,10000);
        % Get edge of Mask
        Mask = edge(Mask,'canny');
        % Thicken Edge by 10
        Mask = bwmorph(Mask,'thicken',10);
        % Close generated gaps with bridge
        Mask = bwmorph(Mask,'bridge');
        % Correct EdgeBinary with Mask
        EdgeBinary = EdgeBinary - Mask;
        % Convert image to 8-bit
        EdgeBinary = uint8(EdgeBinary);
        
        
          
        % This fuction aims to get rid of well rim parts and bubbles and huge areas of
        % different intensities
        % Filter EdgeBinary again for only sizes bigger than FinalNeuronSize and rename it to imgbinary
        imgbinary = bwareaopen(EdgeBinary,FinalNeuronSize);
        % Convert image to 8-bit
        imgbinary = uint8(imgbinary);
        % Label particles in imgbinary
        imgbinaryLabel = bwlabel(imgbinary);
        % Number of particles in imgbinary
        NN = max(max(imgbinaryLabel));
        % Size of image
        [r c] = size(imgbinary);
        % Create empty imaga which will be filled by neuronal structures
        imgBinarySort = zeros(size(imgbinary));
        imgBinarySort = logical(imgBinarySort);
        % Check for each particle whithin imgbinary, whether it is an
        % artifact or a neuronal structure
        for(l=1:NN)
            % Extract inly pixle with index l from labeld image
            NeuronCheck =  imgbinaryLabel;
            NeuronCheck = uint8(NeuronCheck);
            ii = NeuronCheck == l;
            % Set indexed pixels to 255
            NeuronCheck(ii) = 255;
            % Substract 25-, so that only indexed pixel remain positive
            NeuronCheck = NeuronCheck -250;
            % Convert to 8 bit to delete all negative values
            NeuronCheck = uint8(NeuronCheck);
            % Convert image to logical
            NeuronCheck = logical(NeuronCheck);
            % Get area of NeuronCheck
            sOneArea = regionprops(NeuronCheck,'area');
            AreaOne = sOneArea(1).Area;
            %Get area of smallest surrounding box
            sOneBox = regionprops(NeuronCheck,'BoundingBox');
            Width = sOneBox(1).BoundingBox(3);
            Length = sOneBox(1).BoundingBox(4);
            AreaBox = Width*Length;
            % Calculate ratios of area box to area structure
            RatioAreas = AreaOne/AreaBox;
            % Create skeleton of structure
            NeuronCheckSkel = bwmorph(NeuronCheck,'thin',Inf);
            % Get number of branches
            BranchCheck = bwmorph(NeuronCheckSkel,'branch');
            BranchCheck = sum(sum(BranchCheck));
            % Get skeleton length
            LengthCheck = sum(sum(NeuronCheckSkel));
            % Get ratio of branch to length
            LBRatio = BranchCheck/LengthCheck;
            % 1) Check for Ratio areas
            if (RatioAreas > 0.3)
            else
            %2) Then check for ratio of branch to length    
            if(LBRatio < 0.04) 
            imgBinarySort = imgBinarySort +     NeuronCheck;
            end
            end
        end
        imgbinary = imgBinarySort;
        % All neurons touching the edge of the image are removed
        imgbinary = imclearborder(imgbinary);
        
              
        % This function will be used to identify the centroids of the
        % neuron soma
        %Index neurons in the image
        imgbinary = bwlabel(imgbinary);
        % Get maximal index
        NN = max(max(imgbinary));
        % Create empty image later filled with cell somas
        CutNeuriteImg = zeros(size(imgbinary));
        %Convert image to 8-bit
        CutNeuriteImg = uint8(CutNeuriteImg);
        % Get size of image
        [r c] = size(imgbinary);
        % Create a spinning disk for imdilation
        SE = strel('disk',4);
        for(d=1:NN)
            % Create empty image for neuron
            OneNeuron =  zeros(size(imgbinary));
            %Convert image to 8-bit
            OneNeuron = uint8(OneNeuron);
            % Use imgbinary with pixle index d as mask to only get pixle
            % values from imgNeurite within this mask
            for(h=1:r)
            for(o=1:c)
                if(imgbinary(h,o)==d)
                 
                    OneNeuron(h,o) = imgNeurite(h,o);
                end
            end
            end 
            %Convert image to singel
            OneNeuron = single(OneNeuron);
            %Get median from OneNeuron
            Median =  median(OneNeuron((OneNeuron>0)));
            % Subtract 1.5 Median to have only the birght part of the cell body 
            OneNeuron = OneNeuron - (1.5*Median);
            % Get binary image with super low threshold
            OneNeuron = im2bw(OneNeuron,0.00001);
            % Fill holes by generating a filled image
            Holes = imfill(OneNeuron,'holes');
            % Subtract original image from filled image to oly obtain holes
            Holes = Holes - OneNeuron;
            % Get only holes bigger than HoleArea
            HolesBig = bwareaopen(Holes, HoleArea);
            % Subtract HolesBig from Holes
            Holes = Holes - HolesBig;
            % Fill small holes
            OneNeuron = OneNeuron + Holes;
            %Convert image to 8-bit
            OneNeuron = uint8(OneNeuron);
            %Cnvert image to logical
            OneNeuron = logical(OneNeuron);
            % Get number of areas in One Neuron
            Os = regionprops(OneNeuron,'Area');
            sA = size(Os);
            AreaSum = zeros(sA(1),1);
            for(q=1:sA)
                Area = Os(q).Area(1);
                AreaSum(q,1) = Area;
            end
            % Identify largest area and subtract one to keep this area
            imF = max(AreaSum) - 1;
            % Special case if imF is zero it is set to 1
            LimF = length(imF);
            if(LimF < 1);
                imF = 1;
            end    
            % Only largest particele is chosen
            OneNeuron = bwareaopen(OneNeuron, imF);
            % If the particle is bigger than 500
            if(imF>500)
            % Image is eroded
            TestErode = imerode(OneNeuron,SE);
            % If erosion completly deletes the area, a smaller disk is used
            % for erosion
            if(sum(sum(TestErode)) == 0);
                SE = strel('disk',2); 
                OneNeuron = imerode(OneNeuron,SE);
                SE = strel('disk',4); 
            else
                OneNeuron = imerode(OneNeuron,SE);
            end    
            % Again number of areas are checked
            OsE = regionprops(OneNeuron,'Area');
            sAE = size(OsE);
            sAE = sAE(1);
            AreaSumE = zeros(sAE,1);
            % Biggest area is searched
            for(f=1:sAE)
                AreaE = OsE(f).Area(1);
                AreaSumE(f,1) = AreaE;
            end
            % Maximal area -1 to maintain particle after bwareaopen
            imFE = max(AreaSumE) - 1;
            LimFE = length(imFE);
            %Special case if imFE is 0
            if(LimFE < 1);
                imFE = 1;
            end
            % Filter One neuron for largest particle
            OneNeuron = bwareaopen(OneNeuron, imFE);
            end
            %Add particle to CutNeuriteImg
            OneNeuron = uint8(OneNeuron);
            CutNeuriteImg = CutNeuriteImg + OneNeuron;
        end    
        % CutNeuriteImg becomses binaryimgerode       
        binaryimgerode = CutNeuriteImg;
        %Create an empty image for thickend nuclei centroids
        NucleusDummy = zeros(size(imgNeurite));
        % Convert image to logical
        binaryimgerode = logical(binaryimgerode);
        % Test if image contains particles
        if(sum(sum(imgbinary)) >0)
        % If yes determine number of particles
        s = regionprops(binaryimgerode, 'Centroid');
        sdimensions = size(s);
        else
        % If no, set number to 0
        sdimensions = 0;
        end
        % For each particle assess the centroid (x,y coordinates and fill
        % the NucleusM
        for(k=1:sdimensions)
            x_centroid = s(k).Centroid(1);
            y_centroid = s(k).Centroid(2);
            x_centroid = round(x_centroid);
            y_centroid = round(y_centroid);
            NucleusDummy(y_centroid,x_centroid)= 1;
        end
       
        
        %This function is used to cut out the cell soma from the original
        %grayscale image
        % Number of neuron cell bodies
        number_particles = sdimensions;
        %Create empty image for nucleus gray scale image
        imgNucleus = zeros(size(imgNeurite));
        %Create empty image for nucleus binary image
        imgNucleusWatershed = zeros(size(imgNeurite));
        % Convert image to logical
        imgbinary = logical(imgbinary);
        % Now we thicken the NucleusMatrix to serve as an identificator for
        % right area for the cell some
        NucleusDummy = bwmorph(NucleusDummy,'thicken',15);
        %Convert to logical:
        NucleusDummy = logical(NucleusDummy);
        %Convert image to logical
        imgbinary = logical(imgbinary);
        % Now we have to label the particles
        imgbinaryLabel = bwlabel(imgbinary);
        % Now we iterate over all binary components and isolate the cell
        % soma. The basic idea is.....
        NeuronNucOneProcess11 = zeros(size(imgNeurite));
        NeuronNucOneProcess11Fat = zeros(size(imgNeurite));
        NeuronNucOneProcess11FatFat = zeros(size(imgNeurite));
        if(number_particles > 0)
        for(y=1:number_particles)
            %Isolate one neuron binary component
            NeuronNucOne = imgbinaryLabel;
            ii = NeuronNucOne == y;
            NeuronNucOne(ii) = 255;
            NeuronNucOne = NeuronNucOne -250;
            NeuronNucOne = uint8(NeuronNucOne);
            NeuronNucOne = logical(NeuronNucOne);
            %Now we find the centroid lying on the isolated neuronal
            %structure. Since we do not need the entire area, we will first
            %add the NucleusDummy and the subtract 1 to only get the close
            %connected component:
            CCOne = NucleusDummy + NeuronNucOne;
            CCOne = CCOne -1;
            CCOne = uint8(CCOne);
            CCOne = logical(CCOne);
            % Get nucleus from NucleusDummy:
            NucleusDummyLabel = bwlabel(NucleusDummy);
            CentroidThick = zeros(size(imgbinary));
            nD = max(max(NucleusDummyLabel));
              for(z=1:nD)
                  CCOneCom = NucleusDummyLabel;
                  ii = CCOneCom == z;
                  CCOneCom(ii) = 255;
                  CCOneCom = CCOneCom -250;
                  CCOneCom = uint8(CCOneCom);
                  CCOneCom = logical(CCOneCom);
                  if(max(max(CCOneCom+ CCOne))==2)
                      CentroidThick = CCOneCom;
                            % Here we can now also fill holes in imgbinary, which
                            % will be here NeuronNucOne!
                            % The idea would be to chech if CentroidThick is overlapping the holes in th imifilled image. Of course we need to limit the size of the hole!
                            HolesCellBody = imfill(NeuronNucOne,'holes');
                            HolesCellBody = HolesCellBody - NeuronNucOne;
                            HolesCellBodyBig = bwareaopen(HolesCellBody, 200);
                            HolesCellBody = HolesCellBody - HolesCellBodyBig;
                            if(max(max(HolesCellBody+CentroidThick))==2)
                                NeuronNucOne = NeuronNucOne+HolesCellBody;
                            end    
                  end
              end    
            % Now the actual processing takes place
            %1) Thinning of the objects, so that the neurite structures
            %will become lines
            NeuronNucOneProcess = bwmorph(NeuronNucOne,'thin',3);
            % Now we will call the majority function to delete the lines:
            NeuronNucOneProcess1 = bwmorph(NeuronNucOneProcess,'majority');
            % Now we check which area is overlapping with CCOne. In case we
            % have multiple we take the largest:
            %1) We filter for small areas to optimize the sorting!
            NeuronNucOneProcess2 = bwareaopen(NeuronNucOneProcess1,20);
            %Now we label and look for each area whether it overlapped with
            %CCone
            NeuronNucOneProcessLabel = bwlabel(NeuronNucOneProcess2);
            nNucOne = max(max(NeuronNucOneProcessLabel));
            NeuronNucOneProcess3 = zeros(size(imgNeurite));
            for(j=1:nNucOne)
                NeuronNucOnePart = NeuronNucOneProcessLabel;
                ii = NeuronNucOnePart ==j;
                NeuronNucOnePart(ii) = 255;
                NeuronNucOnePart = NeuronNucOnePart -250;
                NeuronNucOnePart = uint8(NeuronNucOnePart);
                NeuronNucOnePart = logical(NeuronNucOnePart);
                if(max(max(NeuronNucOnePart+CCOne)) == 2)
                    % We only include particle if bigger than
                    % NeuronNucOnePorcess3
                    if(sum(sum(NeuronNucOnePart))> sum(sum(NeuronNucOneProcess3)))
                    NeuronNucOneProcess3 = NeuronNucOnePart;
                    end
                end
            end
            %if no particle is left after this procedure, we jsut take the
            %largest One!
            if(sum(sum(NeuronNucOneProcess3))==0)
                SAreaNuc = regionprops(NeuronNucOneProcess2, 'Area');
                SAreaNuc = struct2cell(SAreaNuc);
                SAreaNuc = cell2mat(SAreaNuc);
                SAreaNuc = max(max(SAreaNuc));
                if(SAreaNuc>0)
                NeuronNucOneProcess3 = bwareaopen(NeuronNucOneProcess2,SAreaNuc);
                end
            end    
            %Now we proceed with processing
            %Since our thresholding might include holes, we need to close
            %these now.
            % We will perform image dilataion and use the result as a mask
            % for the original binary image!
            SE = strel('disk',10);
            NeuronNucOneProcess4 = imdilate(NeuronNucOneProcess3,SE);
            NeuronNucOneProcess5 = NeuronNucOneProcess4 + imgbinary;
            NeuronNucOneProcess5 = NeuronNucOneProcess5 -1;
            NeuronNucOneProcess5 = uint8(NeuronNucOneProcess5);
            NeuronNucOneProcess5 = logical(NeuronNucOneProcess5);
            % Now we close potential holes:
            NeuronNucOneProcess6 = imfill(NeuronNucOneProcess5,'holes');
            % Now we will thin again, to further shring the area. However,
            % since we will have differences in thicknes of structures, we
            % will need to check for n iterations of thinning. Basic idea
            % would be to thin, then majority, then check for areas, If
            % there are more than one area repeat until only one area
            % remains!
            Continue = 1;
            %This variable is used to dertermine whther the area of the
            %particle is assessed or not!
            AreaThin = 1;
            for(u=3:50)
                if(Continue ==1)
                Thin = u;
                NeuronNucOneProcess8 = bwmorph(NeuronNucOneProcess6,'thin',u);
                NeuronNucOneProcess9 = bwmorph(NeuronNucOneProcess8,'majority');
                NeuronNucOneProcess10 = bwareaopen(NeuronNucOneProcess9,50);
                % In some instances this will lead to two areas, which area
                % both bigger than 50. In this case, we have to keep the
                % area overlapping with CCOneCom!
                TestParticleThin = bwlabel(NeuronNucOneProcess10);
                NTestParticleThin = max(max(TestParticleThin));
                if(NTestParticleThin>1)
                NeuronNucOneProcess10 = zeros(size(imgbinary));
                    for(r=1:NTestParticleThin)
                    TestParticleThinOne =   TestParticleThin;
                    ii = TestParticleThinOne == r;
                    TestParticleThinOne(ii) = 255;
                    TestParticleThinOne = TestParticleThinOne - 250;
                    TestParticleThinOne = uint8(TestParticleThinOne);
                    TestParticleThinOne = logical(TestParticleThinOne);
                    if(max(max(TestParticleThinOne+CentroidThick))==2);
                        NeuronNucOneProcess10 = TestParticleThinOne;
                    end
                    end
                end 
                %In Case we would have the case that shrinking is leading
                %to now overlapp withany particle, we have to use the
                %particle from the previous iteration!
                if(sum(sum(NeuronNucOneProcess10))==0);
                    Thin = Thin -1;
                    NeuronNucOneProcess8 = bwmorph(NeuronNucOneProcess6,'thin',Thin);
                    NeuronNucOneProcess9 = bwmorph(NeuronNucOneProcess8,'majority');
                    NeuronNucOneProcess10 = bwareaopen(NeuronNucOneProcess9,50);
                    % In some instances this will lead to two areas, which area
                % both bigger than 50. In this case, we have to keep the
                % area overlapping with CCOneCom!
                TestParticleThin = bwlabel(NeuronNucOneProcess10);
                NTestParticleThin = max(max(TestParticleThin));
                if(NTestParticleThin>1)
                NeuronNucOneProcess10 = zeros(size(imgbinary));
                    for(r=1:NTestParticleThin)
                    TestParticleThinOne =   TestParticleThin;
                    ii = TestParticleThinOne == r;
                    TestParticleThinOne(ii) = 255;
                    TestParticleThinOne = TestParticleThinOne - 250;
                    TestParticleThinOne = uint8(TestParticleThinOne);
                    TestParticleThinOne = logical(TestParticleThinOne);
                    if(max(max(TestParticleThinOne+CentroidThick))==2);
                        NeuronNucOneProcess10 = TestParticleThinOne;
                    end
                    end
                end 
                % Now we also need to cancel the loop and need to skip the
                % the Box; As an easy fix we simply set the Area to below
                % the if and will add a variable which will determine
                % whether area is determined or not
                AreaThin = 0;
                BoxSomaArea = 1;
                end
                % After the the first thinning we look for the dimensions
                % of the object in a box. If the box becomes more a square,
                % it is more likely to be at the soma, in case it is more
                % longitutanal, we have still processess!
                if(AreaThin==1)
                BoxSoma = regionprops(NeuronNucOneProcess10,'BoundingBox');
                BoxSomaWidth = BoxSoma(1).BoundingBox(3);
                BoxSomaLength = BoxSoma(1).BoundingBox(4);
                BoxSomaArea = BoxSomaWidth*BoxSomaLength;
                end
                if(BoxSomaArea <1000)
                   ThinFinal = Thin;
                   % We will also thicken the preliminary soma to use it
                   % later to add it the imgbinary to fill holes in the
                   % cell soma!
                   NeuronNucOneProcess10Fat = bwmorph(NeuronNucOneProcess10,'thicken',ThinFinal);
                   NeuronNucOneProcess10Fat = bwmorph(NeuronNucOneProcess10Fat,'bridge');
                   NeuronNucOneProcess10FatFat = bwmorph(NeuronNucOneProcess10,'thicken',ThinFinal+1);
                   NeuronNucOneProcess10FatFat = bwmorph(NeuronNucOneProcess10FatFat,'bridge');
                   %Here we have to cancel the loop!
                   Continue = 0;
                end
                end 
            end
            NeuronNucOneProcess11 = NeuronNucOneProcess11 + NeuronNucOneProcess10;
            NeuronNucOneProcess11Fat = NeuronNucOneProcess11Fat + NeuronNucOneProcess10Fat;
            NeuronNucOneProcess11FatFat = NeuronNucOneProcess11FatFat + NeuronNucOneProcess10FatFat;
        end
        end
        
       % This function will simply correct imgbinary, by adding the
       % NeuronNucOneProcess11Fat to the imgbinary:
       imgbinary = imgbinary + NeuronNucOneProcess11Fat;
       
       % Now we need to reassess the centroid coordinates from
       % NeuronNucOneProcess11:
        NeuronNucOneProcess11 = logical(NeuronNucOneProcess11);
        % Create empty matrix for nucleus matrix
        NucleusM = sparse(zeros(size(imgNeurite)));
        if(sum(sum(imgbinary)) >0)
        % If yes determine number of particles
        sR = regionprops(NeuronNucOneProcess11, 'Centroid');
        sRdimensions = size(sR);
        else
        % If no, set number to 0
        sRdimensions = 0;
        end
        
        % For each particle assess the centroid (x,y coordinates and fill
        % the NucleusM
        for(k=1:sRdimensions)
            x_centroid = sR(k).Centroid(1);
            y_centroid = sR(k).Centroid(2);
            NucleusM(uint16(y_centroid),uint16(x_centroid)) = 1;
        end
       
        % This function will generate the binary and gray scale image from
        % img Nucleus. If we wanna apply filter, the image must be intact.
        % Therefore we will just increase the intensity in the grayscale
        % image for the cell soma
        imgNucleus = zeros(size(imgNeurite));
        imgNucleus = uint8(imgNucleus);
        [r c] = size(imgbinary);
        for(e=1:r)
            for(f=1:c)
                if(NeuronNucOneProcess11FatFat(e,f) == 0)
                    imgNucleusWatershed(e,f) = 0;
                    imgNucleus(e,f) = imgNeurite(e,f);
                else
                    imgNucleusWatershed(e,f) = imgbinary(e,f);
                    imgNucleus(e,f) = (imgNeurite(e,f))*5;
                end
            end
        end    
      
        
        
        
        %This part is for saving the images
        %Wells are named in a linear order (max 99)
        if(n>9)
        a{1} = 'A';
        else
        a{1} = 'A0';
        end
        a{2} = int2str(n);
        a{3} = 'space';
        newfile = [a{1} a{2} a{3}];
        % Names of images
        newFile1 = regexprep(newfile, 'space','NucleusSmall.png','ignorecase');
        newFile2 = regexprep(newfile, 'space','NeuriteBig.png','ignorecase'); 
        newFile3 = regexprep(newfile, 'space','NeuriteSmall.png','ignorecase');
        newFile4 = regexprep(newfile, 'space','Binary.png','ignorecase');
        newFile5 = regexprep(newfile, 'space','NucleusBigWatershed.png','ignorecase');
        newFile7 = regexprep(newfile, 'space','OligoBig.png','ignorecase');
        newFile8 = regexprep(newfile, 'space','OligoSmall.png','ignorecase');
        newFile9 = regexprep(newfile, 'space','NucleusBig.png','ignorecase');
        newFile10 = regexprep(newfile, 'space','Template.png','ignorecase');
        % Saving of Images
        imwrite(imgNucleus,[foldername1 '/' newFile9]);
        imwrite(imgNucleus,[foldername1 '/' newFile1]);
        imwrite(imgNeurite,[foldername1 '/' newFile2]);
        imwrite(imgNeurite,[foldername1 '/' newFile3]);
        imwrite(imgbinary,[foldername1 '/' newFile4]);
        imwrite(imgNucleusWatershed,[foldername1 '/' newFile5]);
        % Well name gets name from filename
        wellname = (newFile1(1:3));
        % Welllist is filled with well names
        wellList{n} = wellname;
        % Cell centroids are send to csvHandler for Omnisphero GUI
        csvHandler.CellPosMatrix(wellname) = NucleusM;
        % ScalingFactor is set to 1
        ScalingFactor = 1;
        clear newfile newfile1 newfile2 newfile3 newfile4 newfile5 newfile6 newfile7 newfile8 
        
    end
    
end

% Comments: It would be good to include a suitable filter for bubble rims