%Copyright (C) 2017-2021  Martin Schmuck

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
function [wellList, foldername, ScalingFactor, csvHandler,   Th] = CorticalNeuronsRat (csvHandler,Auto) 
% This is a quick script to convert 16-bit images of cortical neurons into
% binary image, to assess their center point and to generate a
% 'NucleusImage'.
n = 0;
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
%Get Image folder
foldername = uigetdir;
% Creates directory within the parent folder
allfiles = dir(foldername);
% Creates new folder within parent folder
mkdir(foldername,'ConvertedCellomics');

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
          
        % We use will again reomove identified artifacts from above, but
        % this time growing them before
        ArtefactRemove = bwmorph(ArtefactRemove, 'thicken', 5);
        % Growing leads to artifical gaps, which are connected with
        % 'bridge'
        ArtefactRemove = bwmorph( ArtefactRemove,'bridge');
        % Now we add all filtered edge detector images
         EdgeBinary = E1old + E2old + E3old + E4old + E5old; 
         % Convert image to logical
         EdgeBinary = logical(EdgeBinary);
         Keep = EdgeBinary;
         % Add binary image to the binary Edge image
         EdgeBinary = EdgeBinary + binaryNeuriteFill;
         % Convert image to logical
         EdgeBinary = logical(EdgeBinary);
         % Substract Artifcats from image
         EdgeBinary = EdgeBinary - ArtefactRemove;
         % Convert image to 8-bit
         EdgeBinary = uint8(EdgeBinary);
         % Convert image to logical
         EdgeBinary = logical(EdgeBinary);
         % Perform bridge to connect areas in close proximity
         EdgeBinary = bwmorph(EdgeBinary,'bridge');
         % Sort for particles bigger than FinalNeuronSize
         EdgeBinary = bwareaopen(EdgeBinary,FinalNeuronSize);
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
         % Perform brdige again to connect objects
         EdgeBinary = bwmorph(EdgeBinary,'bridge');
         % Repeat a second time to fill new holes
         BW = imfill(EdgeBinary,'holes');
         BW2 = BW - EdgeBinary;
         BW2 = bwmorph(BW2,'open');
         Holes = bwareaopen(BW2,HoleArea);
         EdgeBinary = BW - Holes;
         % Perform majority to fill small one pixel missing holes
         EdgeBinary = bwmorph(EdgeBinary,'majority');
         % Sort again for big particles
         EdgeBinary = bwareaopen(EdgeBinary,2000);
          % Converti image to 8-bit
         EdgeBinary = uint8(EdgeBinary);
         %Convert image to logical
         EdgeBinary = logical(EdgeBinary);
         %Now we bridge and close small holes again!
         EdgeBinary = bwmorph(EdgeBinary,'bridge');
         %Now we identify holes again
         Holes = imfill(EdgeBinary,'holes');
         Holes = Holes - EdgeBinary;
         HolesBig = bwareaopen(Holes,HoleArea);
         Holes = Holes - HolesBig;
         EdgeBinary = EdgeBinary + Holes;
         %Now we identify holes again
         EdgeBinary = uint8(EdgeBinary);
         %Convert image to logical
         EdgeBinary = logical(EdgeBinary);
         
          
        
         
        % This function is designed to split Neurons which are touching with a single process!
        % Idea we will label the binary particles from EdgeBinary and will
        % then check how many particles from binaryNeuriteFill are within
        % this area. The number of those particles will determine how often
        % we can cut!
        EdgeBinaryLabel = bwlabel(EdgeBinary);
        binaryNeuriteFill = logical(binaryNeuriteFill);
        binaryNeuriteFill = bwareaopen(binaryNeuriteFill,2000);
        NPartic = max(max(EdgeBinaryLabel));
        CutParticles = zeros(size(EdgeBinary));
        for(p=1:NPartic)
            EdgeBinaryOne = EdgeBinaryLabel;
            ii = EdgeBinaryOne == p;
            EdgeBinaryOne(ii) = 255;
            EdgeBinaryOne = EdgeBinaryOne -250;
            EdgeBinaryOne = uint8(EdgeBinaryOne);
            EdgeBinaryOne = logical(EdgeBinaryOne);
            ParticleFill = EdgeBinaryOne + binaryNeuriteFill;
            ParticleFill = ParticleFill - 1;
            ParticleFill = uint8(ParticleFill);
            ParticleFill = bwareaopen(ParticleFill,2000);
            ParticleFill = bwlabel(ParticleFill);
            ParticleFillN = max(max(ParticleFill));
            if(ParticleFillN > 1)
                ParticleFat = bwmorph(ParticleFill,'thicken',5);
                ParticleFat = bwmorph(ParticleFat,'bridge');
                ParticleFatTest = bwlabel(ParticleFat);
                ParticleFatTest = max(max(ParticleFatTest));
                      if(ParticleFatTest<ParticleFillN)
                          ParticleFat = ParticleFill;
                      end 
                CutTest = EdgeBinaryOne - ParticleFat;
                CutTest = uint8(CutTest);
                CutTest = logical(CutTest);
                CutTest = bwareaopen(CutTest,50);
                CutParticle = bwlabel(CutTest);
                NCutParticle = max(max(CutParticle));
                for(y=1:NCutParticle)
                    CutParticleOne = CutParticle; 
                    ii = CutParticleOne == y;
                    CutParticleOne(ii) = 255;
                    CutParticleOne = CutParticleOne - 250;
                    CutParticleOne = uint8(CutParticleOne);
                    CutParticleOne = logical(CutParticleOne);
                    TestSeparation = EdgeBinaryOne - CutParticleOne;
                    TestSeparation = uint8(TestSeparation);
                    TestSeparation = logical(TestSeparation);
                    TestSeparation = bwareaopen(TestSeparation,500);
                    TestSeparationLabel = bwlabel(TestSeparation);
                    NNeurons = max(max(TestSeparationLabel));
                    if(NNeurons >1)
                        CutParticleOneSkel = bwmorph(CutParticleOne,'thin',inf);
                        CutParticleOneSkel = bwmorph(CutParticleOneSkel,'spur',inf);
                        CutParticleOneSkel = bwmorph(CutParticleOneSkel,'thicken',10);
                        CutParticleOneSkel = bwmorph(CutParticleOneSkel,'bridge');
                        CutParticleOneSkel = bwmorph(CutParticleOneSkel,'majority');
                        if(sum(sum(CutParticleOneSkel))>750)
                            CutParticleOneSkel = CutParticleOne;
                        end  
                        CutParticles = CutParticles + CutParticleOneSkel;
                        
                    end
                end   
            end    
        end    
         
        EdgeBinary = EdgeBinary -  CutParticles;
        EdgeBinary = uint8(EdgeBinary);
        EdgeBinary = logical(EdgeBinary);
        EdgeBinary = bwareaopen(EdgeBinary,2000);
        %After this, we repeat the same procedure, to further cut
        %structures!
        EdgeBinaryLabel = bwlabel(EdgeBinary);
        binaryNeuriteFill = logical(binaryNeuriteFill);
        binaryNeuriteFill = bwareaopen(binaryNeuriteFill,2000);
        NPartic = max(max(EdgeBinaryLabel));
        CutParticles = zeros(size(EdgeBinary));
        for(p=1:NPartic)
            EdgeBinaryOne = EdgeBinaryLabel;
            ii = EdgeBinaryOne == p;
            EdgeBinaryOne(ii) = 255;
            EdgeBinaryOne = EdgeBinaryOne -250;
            EdgeBinaryOne = uint8(EdgeBinaryOne);
            EdgeBinaryOne = logical(EdgeBinaryOne);
            ParticleFill = EdgeBinaryOne + binaryNeuriteFill;
            ParticleFill = ParticleFill - 1;
            ParticleFill = uint8(ParticleFill);
            ParticleFill = bwareaopen(ParticleFill,2000);
            ParticleFill = bwlabel(ParticleFill);
            ParticleFillN = max(max(ParticleFill));
            if(ParticleFillN > 1)
                ParticleFat = bwmorph(ParticleFill,'thicken',5);
                ParticleFat = bwmorph(ParticleFat,'bridge');
                ParticleFatTest = bwlabel(ParticleFat);
                ParticleFatTest = max(max(ParticleFatTest));
                      if(ParticleFatTest<ParticleFillN)
                          ParticleFat = ParticleFill;
                      end    
                CutTest = EdgeBinaryOne - ParticleFat;
                CutTest = uint8(CutTest);
                CutTest = logical(CutTest);
                CutTest = bwareaopen(CutTest,50);
                CutParticle = bwlabel(CutTest);
                NCutParticle = max(max(CutParticle));
                for(y=1:NCutParticle)
                    CutParticleOne = CutParticle; 
                    ii = CutParticleOne == y;
                    CutParticleOne(ii) = 255;
                    CutParticleOne = CutParticleOne - 250;
                    CutParticleOne = uint8(CutParticleOne);
                    CutParticleOne = logical(CutParticleOne);
                    TestSeparation = EdgeBinaryOne - CutParticleOne;
                    TestSeparation = uint8(TestSeparation);
                    TestSeparation = logical(TestSeparation);
                    TestSeparation = bwareaopen(TestSeparation,500);
                    TestSeparationLabel = bwlabel(TestSeparation);
                    NNeurons = max(max(TestSeparationLabel));
                    if(NNeurons >1)
                        CutParticleOneSkel = bwmorph(CutParticleOne,'thin',inf);
                        CutParticleOneSkel = bwmorph(CutParticleOneSkel,'spur',inf);
                        CutParticleOneSkel = bwmorph(CutParticleOneSkel,'thicken',10);
                        CutParticleOneSkel = bwmorph(CutParticleOneSkel,'bridge');
                        CutParticleOneSkel = bwmorph(CutParticleOneSkel,'majority');
                        if(sum(sum(CutParticleOneSkel))>750)
                            CutParticleOneSkel = CutParticleOne;
                        end    
                        CutParticles = CutParticles + CutParticleOneSkel;
                    end
                end   
            end    
        end    
         
        EdgeBinary = EdgeBinary -  CutParticles;
        EdgeBinary = uint8(EdgeBinary);
        EdgeBinary = logical(EdgeBinary);
        % Here we can try to eliminate single points connecting two
        % neuronal areas
        PointsDel = zeros(size(EdgeBinary));
        NOrig = bwlabel(EdgeBinary);
        NOrig = max(max(NOrig));
        Singleton = bwmorph(EdgeBinary,'majority');
        Singleton = EdgeBinary - Singleton;
        Singleton = uint8(Singleton);
        SingletonLabel = bwlabel(Singleton);
        SingletonN = max(max(SingletonLabel));
        for(r=1:SingletonN)
            SingletonOne = SingletonLabel;
            ii = SingletonOne == r;
            SingletonOne(ii) = 50000;
            SingletonOne = SingletonOne -49995;
            SingletonOne = uint8(SingletonOne);
            SingletonOne = logical(SingletonOne);
            TestGap = EdgeBinary - SingletonOne;
            TestGap = bwareaopen(TestGap,2000);
            TestGapLabel = bwlabel(TestGap);
            TestGapLabel = max(max(TestGapLabel));
            if(TestGapLabel>NOrig)
                PointsDel = PointsDel + SingletonOne;
            end
        end    
           
            
        
        EdgeBinary = bwareaopen(EdgeBinary,2000);
        
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
        
        % Now we will also delete too big neuron clusters:
        DelCluster = zeros(size(EdgeBinary));
        EdgeBinaryLabel = bwlabel(EdgeBinary);
        NNParticle = max(max(EdgeBinaryLabel));
        for(y=1:NNParticle)
            EdgeBinaryOne = EdgeBinaryLabel;
            ii = EdgeBinaryOne == y;
            EdgeBinaryOne(ii) = 255;
            EdgeBinaryOne = EdgeBinaryOne - 250;
            EdgeBinaryOne = uint8(EdgeBinaryOne);
            EdgeBinaryOne = logical(EdgeBinaryOne);
            BoxCluster = regionprops(EdgeBinaryOne,'BoundingBox');
            WidthCluster = BoxCluster(1).BoundingBox(3);
            LengthCluster = BoxCluster(1).BoundingBox(4);
            AreaCluster = WidthCluster * LengthCluster;
            if(AreaCluster>750000)
                DelCluster = DelCluster + EdgeBinaryOne;
            end    
        end    
        DelCluster = uint8(DelCluster);
        EdgeBinary = EdgeBinary - DelCluster;
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
            % Extract only pixle with index l from labeld image
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
        imgbinaryKeep = imgbinary;
        
        if(Auto==1)
        %         %Check for double structures:
%         %Here we will automatically delete any double structure to enhance
%         %speed and accuracy for later evaluation:
        %We will use a minimum size criteria of 6000 pixel to exclude small
        %neurons and to speed up evluation:
        BinaryTestLabel = bwlabel(imgbinary);
        BinaryTestLabelN = max(max(BinaryTestLabel));
        imgPresort = zeros(size(imgbinary));
        for(h=1:BinaryTestLabelN)
            BinaryTest = BinaryTestLabel;
            ii = BinaryTest == h;
            BinaryTest(ii) = 255;
            BinaryTest = BinaryTest - 250;
            BinaryTest = uint8(BinaryTest);
            BinaryTest = logical(BinaryTest);
        if(sum(sum(BinaryTest))<4000)
                
        else
            %Here will now also delete dead neurons:
            BinaryTestI = imcomplement(BinaryTest);
            BinaryTestI = uint8(BinaryTestI*255);
            NeuriteSingle = imgNeurite - (BinaryTestI);
            NeuronImgOne = uint8(NeuriteSingle);
            level = isodata(NeuronImgOne(NeuronImgOne>0));
            if(level>0.4)
            NucTest = im2bw(NeuronImgOne,level*0.5);
            NucTest = bwmorph(NucTest,'thin',2);
            NucTestSoma = bwmorph(NucTest,'endpoints');
            NucTestSoma = bwmorph(NucTestSoma,'bridge');
            NucTestSoma = imfill(NucTestSoma,'holes');
            NucTestSoma=bwareaopen(NucTestSoma,250);
            NucTestArea = sum(sum(NucTestSoma));
            NucRatio = NucTestArea/sum(sum(BinaryTest));
            if(sum(sum(BinaryTest))<10000)
                Criteria = 300;
            else
                Criteria = 400;
            end    
            if(NucTestArea>Criteria && level<0.49)
                %Check for number of adjactant primary processis:
                ProcessSoma = bwmorph(NucTest,'spur',500);
                ProcessSoma = bwmorph(NucTest,'open');
                ProcessSoma = bwareaopen(ProcessSoma,250);
                Process = NucTest - ProcessSoma;
                Process = bwareaopen(Process,2);
                Process = bwmorph(Process,'endpoints');
                ProcessSomaFat = bwmorph(ProcessSoma,'thicken',4);
                Processis = Process+   ProcessSomaFat;
                Processis = Processis -1;
                Processis = uint8(Processis);
                Processis = logical(Processis);
                ProcessisS = bwareaopen(Processis,2);
                Processis = Processis - ProcessisS;
                ProcessisN = sum(sum(Processis));
                if(sum(sum(BinaryTest))>25000)
                    CriteriaProcessis = 20;
                else
                    CriteriaProcessis = 16;
                end    
                if(ProcessisN<CriteriaProcessis)
                    imgPresort = imgPresort+BinaryTest; 
                end
                
            end
            else
            imgPresort = imgPresort+BinaryTest;   
        end
            
        end 
        end
        imgbinary = logical(imgPresort); 
        
        %Check for double structures:
        %Here we will automatically delete any double structure to enhance
        %speed and accuracy for later evaluation:
        imgbinaryFinal = zeros(size(imgbinary));
        BinaryTestLabel = bwlabel(imgbinary);
        BinaryTestLabelN = max(max(BinaryTestLabel));
        %Here we will generate a binary image of the original image:
        level = isodata(imgNeurite);
        for(h=1:BinaryTestLabelN)
        t = t+1;
        ExcludeDouble = 0;
        ExcludeDoubleF = 1;
        BinaryTest = BinaryTestLabel;
        ii = BinaryTest == h;
        BinaryTest(ii) = 255;
        BinaryTest = BinaryTest - 250;
        BinaryTest = uint8(BinaryTest);
        BinaryTest = logical(BinaryTest);
        BinaryTestI = imcomplement(BinaryTest);
        BinaryTestI = uint8(BinaryTestI*255);
        NeuriteSingle = imgNeurite - (BinaryTestI);
        NeuriteSingle = uint8(NeuriteSingle);
        SkeletonThin=0;
        
        
             
        
        SE = strel('disk',3);
        %We will generate the right threshold here, based on the intensity of the Neuron:
         if(level<0.07)
           if(mean(mean(NeuriteSingle(NeuriteSingle>0)))>30) 
           DoubleTest = im2bw( NeuriteSingle,level*2);
           else
           DoubleTest = im2bw( NeuriteSingle,level*1.5);
           end
           DoubleTest0 = DoubleTest;
           %DoubleTest = bwareaopen(DoubleTest,500);
           DoubleTestL = imerode(DoubleTest,SE);
           DoubleTestL = bwmorph(DoubleTestL,'open');
           DoubleTestL = bwareaopen(DoubleTestL,75);
            SE = strel('disk',5);
           DoubleTestL = imdilate(DoubleTestL,SE);
           DoubleTestL = bwmorph(DoubleTestL,'bridge');
           DoubleTestN = bwlabel(DoubleTestL);
           DoubleTestN = max(max(DoubleTestN));
              
           
               
        end
        if(level<0.1 && level>=0.07)
           DoubleTest = im2bw( NeuriteSingle,level*1.8);
           DoubleTest0 = DoubleTest;
           %DoubleTest = bwareaopen(DoubleTest,500);
           SE = strel('disk',3);
           DoubleTestL = imerode(DoubleTest,SE);
           DoubleTestL = bwmorph(DoubleTestL,'open');
           DoubleTestL = bwareaopen(DoubleTestL,75);
           SE = strel('disk',5);
           DoubleTestL = imdilate(DoubleTestL,SE);
           DoubleTestL = bwmorph(DoubleTestL,'bridge');
           DoubleTestN = bwlabel(DoubleTestL);
           DoubleTestN = max(max(DoubleTestN));
             
           
        end
        if(level<0.2 && level>=0.1)
           DoubleTest = im2bw( NeuriteSingle,level*1.5);
           DoubleTest0 = DoubleTest;
           %DoubleTest = bwareaopen(DoubleTest,500);
           SE = strel('disk',3);
           DoubleTestL = imerode(DoubleTest,SE);
           DoubleTestL = bwmorph(DoubleTestL,'open');
           DoubleTestL = bwareaopen(DoubleTestL,100);
            SE = strel('disk',5);
           DoubleTestL = imdilate(DoubleTestL,SE);
           DoubleTestL = bwmorph(DoubleTestL,'bridge');
           DoubleTestN = bwlabel(DoubleTestL);
           DoubleTestN = max(max(DoubleTestN));
                
         
        end
        if(level>0.2 && level >=0.2)
           DoubleTest = im2bw( NeuriteSingle,level*0.5);
           %DoubleTest = bwareaopen(DoubleTest,500);
           DoubleTest0 = DoubleTest;
           DoubleTestL = imerode(DoubleTest,SE);
           DoubleTestL = bwmorph(DoubleTestL,'open');
           DoubleTestL = bwareaopen(DoubleTestL,100);
            SE = strel('disk',5);
           DoubleTestL = imdilate(DoubleTestL,SE);
           DoubleTestL = bwmorph(DoubleTestL,'bridge');
           DoubleTestN = bwlabel(DoubleTestL);
           DoubleTestN = max(max(DoubleTestN));
           SkeletonThin = 1;
                 
           
        end
        DoubleTest = bwareaopen(DoubleTest,250);
        if(DoubleTestN==1)
        if(SkeletonThin==1)
            DoubleTest = bwmorph(DoubleTest,'thin',4);
        else    
            DoubleTest = bwmorph(DoubleTest,'thin',7);
        end    
        SkeletonCon = bwmorph(DoubleTest,'open');
        CellBodiesOne = SkeletonCon;
        CellBodiesOne = bwareaopen(CellBodiesOne,100);
        SmallBodies = SkeletonCon - CellBodiesOne; 
        SmallBodies = uint8(SmallBodies );
        SmallBodies  = logical(SmallBodies );
        SE = strel('disk',3);
        CellBodiesOneFat = imdilate(CellBodiesOne,SE);
        CellBodiesOneN = bwlabel(CellBodiesOneFat);
        CellBodiesOneN = max(max(CellBodiesOneN));
        if(CellBodiesOneN==0)
            DoubleTest = bwareaopen(DoubleTest0,250);
            DoubleTest = bwmorph(DoubleTest,'thin',5);
            SkeletonCon = bwmorph(DoubleTest,'open');
            CellBodiesOne = SkeletonCon;
            CellBodiesOne = bwareaopen(CellBodiesOne,30);
            SmallBodies = SkeletonCon - CellBodiesOne; 
            SmallBodies = uint8(SmallBodies );
            SmallBodies  = logical(SmallBodies );
            CellBodiesOneFat = imdilate(CellBodiesOne,SE);
            CellBodiesOneN = bwlabel(CellBodiesOneFat);
            CellBodiesOneN = max(max(CellBodiesOneN));
        end    
        SkeletonCon = DoubleTest - CellBodiesOne + SmallBodies;
        SkeletonCon = uint8(SkeletonCon);
        SkeletonCon = logical(SkeletonCon);
        %Skeleton repair:
        SkeletonCon = bwmorph(SkeletonCon,'thin',inf);
        SkeletonCon = SkeletonCon + CellBodiesOne;
        SkeletonCon = bwmorph( SkeletonCon,'bridge');
        SkeletonCon = SkeletonCon - CellBodiesOne;
        SkeletonCon = bwmorph(SkeletonCon,'thin',inf);
        SkeletonCon = uint8(SkeletonCon);
        SkeletonCon = logical(SkeletonCon);
        %SkeletonCon = bwmorph( SkeletonCon,'thin',inf);
        SkeletonConLabel = bwlabel(SkeletonCon);
        SkeletonConN = max(max(SkeletonConLabel));
        DoubleSkel = bwmorph(DoubleTest,'thin',inf); 
        DoubleSkel = bwlabel(DoubleSkel);
        DoubleSkel =max(max(DoubleSkel));
        %Check if double structure is part of a linear neurite, by
        %analyzing endpoints:
        
            if(DoubleSkel>1)
                SkeletonConBig = bwareaopen(SkeletonCon,250);
                SkeletonConSmall = SkeletonCon - SkeletonConBig;
                EndpointCheck = bwmorph(SkeletonConSmall,'endpoints');
                if(sum(sum(EndpointCheck))>3)
                ExcludeDouble=1;
                end
            end 
        else
            DoubleTest = bwareaopen(DoubleTest0,500);
            DoubleTest = bwmorph(DoubleTest,'thin',inf);
            DoubleSkel = bwareaopen(DoubleTest,100);
            DoubleSkel = DoubleSkel+DoubleTestL;
            DoubleSkel = logical(DoubleSkel);
            DoubleSkel = bwlabel(DoubleTest);
            DoubleSkel =max(max(DoubleSkel));
            if(DoubleSkel>1)
                ExcludeDouble=1;
            end    
            CellBodiesOne=DoubleTestL;
            CellBodiesOneN = DoubleTestN;
            SkeletonCon = DoubleTest - CellBodiesOne;
            SkeletonCon = uint8(SkeletonCon);
            SkeletonCon = logical(SkeletonCon);
            SkeletonCon = bwmorph( SkeletonCon,'thin',inf);
            SkeletonConLabel = bwlabel(SkeletonCon);
            SkeletonConN = max(max(SkeletonConLabel));
        end
          %Here we add our criteria for excluding cells with big somata:
          BigBody=0;
    if(sum(sum(CellBodiesOne))>5500)
        ExcludeDouble=1;
        BigBody =1;
       
    end    
    if(sum(sum(CellBodiesOne))>3000&&sum(sum(BinaryTest))<12000)
        ExcludeDouble=1;
        BigBody =1;
        
    end 
        
        if(DoubleSkel==1&&BigBody==0)
        if(CellBodiesOneN>1)
        Found=0;
        minLength = 10000000000;
        SkeletonConAll = zeros(size(CellBodiesOne));
        TestConKeep = zeros(size(CellBodiesOne));
        SkeletonConOneKeep = zeros(size(CellBodiesOne));
        for(w=1:SkeletonConN)
              SkeletonConOne = SkeletonConLabel;
              ii = SkeletonConOne == w;
              SkeletonConOne(ii) = 255;
              SkeletonConOne = SkeletonConOne -250;
              SkeletonConOne = uint8(SkeletonConOne);
              SkeletonConOne = logical(SkeletonConOne);
              TestCon = SkeletonConOne + CellBodiesOne;
              TestCon = logical(TestCon);
              TestConN = bwlabel(TestCon);
              TestConN = max(max(TestConN));
              %Here we will simply look for the shortest connection to
              %exclude cases in which the longer one is found first
              if(TestConN<CellBodiesOneN)
                if(sum(sum(SkeletonConOne))<minLength)
                    minLength = sum(sum(SkeletonConOne));
                    TestConKeep = TestCon;
                    SkeletonConOneKeep = SkeletonConOne;
                    
                end
                SkeletonConAll = SkeletonConAll + SkeletonConOne;
              end  
        end
       
         
        
       
       TestCon = TestConKeep;
       SkeletonConOne = SkeletonConOneKeep;
       TestConN = bwlabel(TestCon);
       TestConN = max(max(TestConN));
       if(sum(sum(    TestConKeep))==0)
           TestConN=10;
       end  
              if(TestConN==1)
                  LengthCon = sum(sum(SkeletonConOne));
                  BranchCon = bwmorph(SkeletonConOne,'Branch');
                  BranchCon = sum(sum(BranchCon));
                  %Test for direct connection:
                  %Get rid of loops in skeleton:
                  SkeletonConOneC = imfill(SkeletonConOne,'holes');
                  SkeletonConOneC = bwmorph( SkeletonConOneC,'thin',inf);
                  %LoopStructure = bwmorph(SkeletonConOne,'spur',500);
                  TestDirect = SkeletonConOneC+CellBodiesOne;
                  %Check for open structures!
                  TestDirectClose = bwlabel(TestDirect);
                  TestDirectClose = max(max(TestDirectClose));
                  if(TestDirectClose>1)
                      TestDirect=bwmorph(TestDirect,'bridge');
                  end    
                  %TestDirect = TestDirect- LoopStructure;
                  TestDirect = uint8(TestDirect);
                  TestDirect= logical(TestDirect);
                  TestDirect = bwmorph(TestDirect,'spur',500);
                  %Now we have to account for multiple connections to one
                  %cell nuclei:
                  EndpointsTestDirekt =  TestDirect - CellBodiesOne;
                  EndpointsTestDirektE = bwmorph(EndpointsTestDirekt,'endpoints');
                  CellBodiesOneL = bwlabel(CellBodiesOne);
                  CellBodiesOneN = max(max(CellBodiesOneL));
                  CountNuc = 0;
                  Matrix_1x =0;
                  Matrix_1y=0;
                  Matrix_2x=0;
                  Matrix_2y=0;
                  if(sum(sum(EndpointsTestDirektE))>2 && CellBodiesOneN==2)
                      for(d=1:CellBodiesOneN)
                          CountNuc = CountNuc+1;
                          CellBodiesOneOne = CellBodiesOneL;
                          ii = CellBodiesOneOne== d;
                          CellBodiesOneOne(ii)=255;
                          CellBodiesOneOne = CellBodiesOneOne-250;
                          CellBodiesOneOne =uint8(CellBodiesOneOne);
                          CellBodiesOneOne=logical(CellBodiesOneOne);
                          CellBodiesOneOneF = bwmorph(CellBodiesOneOne,'thicken',2);
                          EndpointsSingle = EndpointsTestDirektE+CellBodiesOneOneF;
                          EndpointsSingle = EndpointsSingle-1;
                          EndpointsSingle=uint8(EndpointsSingle);
                          EndpointsSingle=logical(EndpointsSingle);
                          EndpointsSingle = bwmorph(EndpointsSingle,'spur',10);
                          EndpointsSingleN = bwlabel(EndpointsSingle);
                          EndpointsSingleN = max(max(EndpointsSingleN));
                          EndpointsSingleC = regionprops(EndpointsSingle,'Centroid');
                          for(p=1:EndpointsSingleN)
                              if(CountNuc==1)
                                Matrix_1x(p,1) = EndpointsSingleC(p).Centroid(1);
                                Matrix_1y(p,1) = EndpointsSingleC(p).Centroid(2);
                              else
                                Matrix_2x(p,1) = EndpointsSingleC(p).Centroid(1);
                                Matrix_2y(p,1) = EndpointsSingleC(p).Centroid(2);  
                              end
                          end
                                                    
                      end
                     
                  LengthMatrix_1=length(Matrix_1x);
                  LengthMatrix_2=length(Matrix_2x);
                  DistanceECompare=10000000;
                  ClosestEndpoints = zeros(size(CellBodiesOne));
                  DistanceE=0;
                  for(j=1:LengthMatrix_1)
                      x_1 = Matrix_1x(j,1);
                      y_1 = Matrix_1y(j,1);
                      for(z=1:LengthMatrix_2)
                          x_2= Matrix_2x(z,1);
                          y_2= Matrix_2y(z,1);
                          DistanceE(z,1) = ((x_2-x_1)^2+(y_2-y_1)^2)^0.5;
                      
                      DistanceEmin = min(DistanceE);
                      if(DistanceEmin<DistanceECompare)
                          x_1Keep = round(x_1,0);
                          y_1Keep = round(y_1,0);
                          x_2Keep = round(x_2,0);
                          y_2Keep = round(y_2,0);
                      end
                      DistanceECompare=DistanceEmin;
                      end
                  end
                  if(x_1Keep>0&& y_1Keep>0&&x_2Keep>0&&y_2Keep>0)
                  ClosestEndpoints(y_1Keep,x_1Keep)=1;
                  ClosestEndpoints(y_2Keep,x_2Keep)=1;
                  ClosestEndpoints = logical(ClosestEndpoints);
                  EndpointsTestDirektE = logical(EndpointsTestDirektE);
                  ClosestEndpoints = logical(ClosestEndpoints);
                  TestDirect = logical(TestDirect);
                  EndpointsDel =  EndpointsTestDirektE - ClosestEndpoints;
                  EndpointsDel = bwmorph(EndpointsDel,'thicken',2);
                  TestDirect = TestDirect -  CellBodiesOne;
                  TestDirect = TestDirect - EndpointsDel;
                  TestDirect = uint8(TestDirect);
                  TestDirect = logical(TestDirect);
                  TestDirect = TestDirect+CellBodiesOne;
                  end
                  end
                  TestDirect = bwmorph(TestDirect,'spur',500);
                  TestDirect = bwareaopen(TestDirect,100);
                  TestDirect = bwmorph( TestDirect,'bridge');
                  TestDirect = TestDirect -CellBodiesOne ;
                  TestDirect = logical(TestDirect);
                  %TestDirect = bwareaopen(TestDirect,10);
                  TestDirectL = bwlabel(TestDirect);
                  TestDirectL = max(max(TestDirectL));
                  if(TestDirectL>1)
                      TestDirect=bwareaopen(TestDirect,2);
                      TestDirectL = bwlabel(TestDirect);
                      TestDirectL = max(max(TestDirectL));
                  end    
                  if(TestDirectL==1)
                      LengthCon = sum(sum(TestDirect));
                      BranchCon = bwmorph(TestDirect,'Branch');
                      BranchCon = sum(sum(BranchCon));
                      if(LengthCon>50)
                      if(LengthCon>50 &&LengthCon<120)
                          Connection = SkeletonConAll;
                          %Now we can test if on of the cell bodies did not
                          %has any other long processis, like a dirt particle attached to a neurite:
                          TestParticle = SkeletonCon+CellBodiesOne-Connection;
                          TestParticle = logical(TestParticle);
                          TestParticleL = bwlabel(TestParticle);
                          TestParticleN = max(max(TestParticleL));
                          Number=0;
                          for(l=1:TestParticleN)
                              TestParticleOne = TestParticleL;
                              ii = TestParticleOne == l;
                              TestParticleOne(ii)=255;
                              TestParticleOne = TestParticleOne -250;
                              TestParticleOne = uint8(TestParticleOne);
                              TestParticleOne = logical(TestParticleOne);
                              SkeletonParticle = TestParticleOne - CellBodiesOne;
                              SkeletonParticle = uint8(SkeletonParticle);
                              SkeletonParticle = logical(SkeletonParticle);
                              SkeletonParticleL = bwlabel(SkeletonParticle);
                              SkeletonParticleN = max(max(SkeletonParticleL));
                              SkeletonParticleEndpoints = bwmorph(SkeletonParticle,'endpoints');
                              CellBodiesOneFat = bwmorph(CellBodiesOne,'thicken',2);
                              SkeletonParticleEndpoints = SkeletonParticleEndpoints+ CellBodiesOneFat;
                              SkeletonParticleEndpoints = SkeletonParticleEndpoints - 1;
                              SkeletonParticleEndpoints = uint8(SkeletonParticleEndpoints);
                              SkeletonParticleEndpoints = logical(SkeletonParticleEndpoints);
                              SkeletonParticleEndpointsN = bwlabel(SkeletonParticleEndpoints);
                              SkeletonParticleEndpointsN = max(max(SkeletonParticleEndpointsN));
                              if(SkeletonParticleEndpointsN<3 && SkeletonParticleN<2&&sum(sum(SkeletonParticle))<15)
                                 Number =1;
                              end
                          end
                              if(Number==0)                           
                                    ExcludeDouble = 1;
                                    
                              end  
                              Found = 1;
                      else
                          Connection = SkeletonConAll;
                          %Now we can test if on of the cell bodies did not
                          %has any other long processis, like a dirt particle attached to a neurite:
                          TestParticle = SkeletonCon+CellBodiesOne-Connection;
                          TestParticle = logical(TestParticle);
                          TestParticleL = bwlabel(TestParticle);
                          TestParticleN = max(max(TestParticleL));
                          if(TestParticleN<CellBodiesOneN)
                              Connection =zeros(size(CellBodiesOne));
                                 for(w=1:SkeletonConN)
                                      SkeletonConOne = SkeletonConLabel;
                                      ii = SkeletonConOne == w;
                                      SkeletonConOne(ii) = 255;
                                      SkeletonConOne = SkeletonConOne -250;
                                      SkeletonConOne = uint8(SkeletonConOne);
                                      SkeletonConOne = logical(SkeletonConOne);
                                      TestCon = SkeletonConOne + CellBodiesOne;
                                      TestCon = logical(TestCon);
                                      TestConN = bwlabel(TestCon);
                                      TestConN = max(max(TestConN));
                                      %Here we will simply look for the shortest connection to
                                      %exclude cases in which the longer one is found first
                                      if(TestConN<CellBodiesOneN)
                                            Connection = SkeletonConOne+Connection;

                                      end
                                        
                            end 
                               TestParticle = SkeletonCon+CellBodiesOne-Connection;
                              TestParticle = logical(TestParticle);
                              TestParticleL = bwlabel(TestParticle);
                              TestParticleN = max(max(TestParticleL));
                          end      
                          Number=0;
                          CountNumber =0;
                          for(l=1:TestParticleN)
                              TestParticleOne = TestParticleL;
                              ii = TestParticleOne == l;
                              TestParticleOne(ii)=255;
                              TestParticleOne = TestParticleOne -250;
                              TestParticleOne = uint8(TestParticleOne);
                              TestParticleOne = logical(TestParticleOne);
                              SkeletonParticle = TestParticleOne - CellBodiesOne;
                              SkeletonParticle = uint8(SkeletonParticle);
                              SkeletonParticle = logical(SkeletonParticle);
                              SkeletonParticleL = bwlabel(SkeletonParticle);
                              SkeletonParticleN = max(max(SkeletonParticleL));
                              SkeletonParticleEndpoints = bwmorph(SkeletonParticle,'endpoints');
                              CellBodiesOneFat = bwmorph(CellBodiesOne,'thicken',2);
                              SkeletonParticleEndpoints = SkeletonParticleEndpoints+ CellBodiesOneFat;
                              SkeletonParticleEndpoints = SkeletonParticleEndpoints - 1;
                              SkeletonParticleEndpoints = uint8(SkeletonParticleEndpoints);
                              SkeletonParticleEndpoints = logical(SkeletonParticleEndpoints);
                              SkeletonParticleEndpointsN = bwlabel(SkeletonParticleEndpoints);
                              SkeletonParticleEndpointsN = max(max(SkeletonParticleEndpointsN));
                              if(SkeletonParticleEndpointsN<3 && sum(sum(SkeletonParticle))<40)
                                 Number =1;
                              else
                                  CountNumber = CountNumber+1;
                              end
                          end
                          NumberConnections = bwlabel(Connection);
                          NumberConnections = max(max(NumberConnections));
                          if(NumberConnections==1)
                          if(Number==0)
                              
                                    ExcludeDouble = 1;
                                    
                          end  
                              Found = 1;
                          else
                              ExcludeDouble = 1;
                              Found = 1;
                          end    
                          
                      end
                      else
                          if(sum(sum(BinaryTest))>10000);
                          Found = 1;
                          ExcludeDouble = 0;
                          else
                              Found = 1;
                              ExcludeDouble = 1;
                          end    
                              
                      end    
                  else    
                  if(LengthCon>400 |BranchCon>4)
                          Connection = SkeletonConOne;
                          %Now we can test if on of the cell bodies did not
                          %has any other long processis, like a dirt particle attached to a neurite:
                          TestParticle = SkeletonCon+CellBodiesOne-Connection;
                          TestParticle = logical(TestParticle);
                          TestParticleL = bwlabel(TestParticle);
                          TestParticleN = max(max(TestParticleL));
                          Number=0;
                          for(l=1:TestParticleN)
                              TestParticleOne = TestParticleL;
                              ii = TestParticleOne == l;
                              TestParticleOne(ii)=255;
                              TestParticleOne = TestParticleOne -250;
                              TestParticleOne = uint8(TestParticleOne);
                              TestParticleOne = logical(TestParticleOne);
                              SkeletonParticle = TestParticleOne - CellBodiesOne;
                              SkeletonParticle = uint8(SkeletonParticle);
                              SkeletonParticle = logical(SkeletonParticle);
                              SkeletonParticleL = bwlabel(SkeletonParticle);
                              SkeletonParticleN = max(max(SkeletonParticleL));
                              if(SkeletonParticleN<2)
                                 Number =1;
                              end
                          end
                          if(Number==0)
                          ExcludeDouble = 1;  
                          end
                      Found = 1;
                  else
                      if(TestDirectL==0)
                          Found=1;
                          ExcludeDouble = 1;
                      else
                          Found = 1;
                          ExcludeDouble = 0;
                      end
                  end
                  
                  
                  end        
              end  
         if(Found==0)
            ExcludeDouble=1;
        end 
        end
        
        end
        if(CellBodiesOneN>2)
            ExcludeDouble=1;
        end  
               
            
          
        
        
        
        if(ExcludeDouble<1)
        imgbinaryFinal =imgbinaryFinal+BinaryTest;
        end
        end
        imgbinary = imgbinaryFinal;
        imageBinaryFirst = imgbinary;
        
        
        %Now we will do second seletion round:
        
        %         %Check for double structures:
%         %Here we will automatically delete any double structure to enhance
%         %speed and accuracy for later evaluation:
        Bouble = edge(imgNeurite,'canny');
        Bouble1 = edge(imgNeurite,'prewitt');
        Bouble2 = edge(imgNeurite,'sobel');
        Bouble3 = edge(imgNeurite,'log');
        Bouble= Bouble+Bouble1+Bouble2+Bouble3;
        imgbinayMaskFat = bwmorph(imgbinary,'thicken',3);
        Bouble = Bouble + logical(imgbinayMaskFat);
        Bouble = imfill(Bouble,'holes');
        Bouble = bwareaopen(Bouble,500);
        Bouble= uint8(Bouble);
        Bouble = logical(Bouble);
        SE = strel('disk',3);
        Bouble = imdilate(Bouble,SE);
        Bouble = imfill(Bouble,'holes');
        Bouble2nd = Bouble;
        Bouble = bwmorph(Bouble,'thin',100);
        BoublePart = bwmorph(Bouble2nd,'thin',10);
        BoublePart = BoublePart - logical(imgbinary);
        BoublePart = uint8(BoublePart);
        BoublePart = logical(BoublePart);
        BoublePart = bwareaopen(BoublePart,100);
        Lines = bwmorph(BoublePart,'open');
        BoublePart = BoublePart - Lines;
        BoublePart = bwareaopen(BoublePart,100);
        BoublePartLabel = bwlabel(BoublePart);
        BoublePartN = max(max(BoublePartLabel));
        BoublePartKeep = logical(zeros(size(Bouble)));
        if(sum(sum(BoublePart)))
        for(f=1:BoublePartN)
            BoublePartOne = BoublePartLabel;
            ii = BoublePartOne ==f;
            BoublePartOne(ii) = 65500;
            BoublePartOne = BoublePartOne -65000;
            BoublePartOne = uint8(BoublePartOne);
            BoublePartOne = logical(BoublePartOne);
            BoublePartOneP = bwmorph(BoublePartOne,'spur',20);
            BoubleEnd = bwmorph(BoublePartOneP,'Endpoint');
            BoubleEnd = sum(sum(BoubleEnd));
            AxisLength = regionprops(BoublePartOneP,'MajorAxis');
            AxisLength = AxisLength.MajorAxisLength;
            if(AxisLength>200&& BoubleEnd<3)
                BoublePartKeep = BoublePartKeep + BoublePartOne ;
            end    
        end    
        end    
        Bouble = bwmorph(Bouble,'open');
        Bouble = bwareaopen(Bouble,15000);
        Bouble = logical(Bouble);
        BoubleLabel = bwlabel(Bouble);
        BoubleN = max(max(BoubleLabel));
        BoubleKeep = zeros(size(Bouble));
        for(q=1:BoubleN);
            BoubleOne = BoubleLabel;
            ii = BoubleOne == q;
            BoubleOne(ii) = 255;
            BoubleOne = BoubleOne -250;
            BoubleOne = uint8(BoubleOne);
            BoubleOne = logical(BoubleOne);
            Ecc = regionprops(BoubleOne,'Perimeter');
            perimeter = Ecc(1).Perimeter;
            area = sum(sum(BoubleOne));
            metric = 4*pi*area/perimeter^2;
            if(metric>0.8)
                BoubleKeep = BoubleKeep+BoubleOne;
            end
        end
        Bouble = BoubleKeep;
          % Here we check for dead neurons! Idea, just look fora ares with in 255, since dead Neurons are very bright:
          imgNeuriteBrighth = imgNeurite - 254;
          imgNeuriteBright = uint8(imgNeuriteBrighth);
          imgNeuriteBrighth = logical(imgNeuriteBrighth);
          AreaBright = sum(sum(imgNeuriteBrighth));
        imgbinaryFinal = zeros(size(imgbinary));
        BinaryTestLabel = bwlabel(imgbinary);
        BinaryTestLabelN = max(max(BinaryTestLabel));
        for(h=1:BinaryTestLabelN)
        t = t+1;
        BinaryTest = BinaryTestLabel;
        ii = BinaryTest == h;
        BinaryTest(ii) = 255;
        BinaryTest = BinaryTest - 250;
        BinaryTest = uint8(BinaryTest);
        BinaryTest = logical(BinaryTest);
        FillTest = imfill(BinaryTest,'holes');
        FillTest = FillTest - BinaryTest;
        FillTest = uint8(FillTest);
        FillTest = logical(FillTest);
        FillTestLabel = bwlabel(FillTest);
        FillTestN = max(max(FillTestLabel));
        FillRatio = sum(sum(FillTest))/sum(sum(BinaryTest));
        if(FillRatio>1 && FillTestN>5)
            SE =strel('disk',5);
            FillTestArea = imdilate(FillTest,SE);
            FillTestArea = logical(FillTestArea);
            FillTestAreaN = regionprops(FillTestArea);
            FillTestAreaN = length(FillTestAreaN);
            if(FillTestAreaN==1)
                FillRatio = 0.5;
            end
        end  
        %Here we test for bubbles:
        if(sum(sum(Bouble>0)))
        BoubleTest = Bouble+BinaryTest;
        BoubleTest = BoubleTest-1;
        BoubleTest = uint8(BoubleTest);
        BoubleTest = logical(BoubleTest);
        if(sum(sum(BoubleTest))>0)
            FillTestN =100;
        end
        end
        if(sum(sum(BoublePartKeep))>0)
        %Now also for bubble parts:
        BoubleTestFat = bwmorph(BinaryTest,'thicken',15);
        BoublePartTest = BoubleTestFat + BoublePartKeep;
        BoublePartTest = BoublePartTest -1;
        BoublePartTest = uint8(BoublePartTest);
        BoublePartTest = logical(BoublePartTest);
        if(sum(sum(BoublePartTest))>0)
            FillTestN =100;
        end
        end
        BinaryTestFat = bwmorph(BinaryTest,'thicken',2);
        BoubleTest = Bouble + BinaryTestFat;
        BoubleTest = bwareaopen(BoubleTest,sum(sum(BinaryTest)));
        BoubleTest = BoubleTest-BinaryTestFat;
        BoubleTest = uint8(BoubleTest);
        BoubleTest = logical(BoubleTest);
        BoubleTest = bwareaopen(BoubleTest,50);
        BoubleTestLabel = bwlabel(BoubleTest);
        BoubleTestLabelN = max(max(BoubleTestLabel));
        if(sum(sum(BoubleTest))>0)
        for(l=1:BoubleTestLabelN);
            BoubleTestOne = BoubleTestLabel;
            ii = BoubleTestOne == l;
            BoubleTestOne(ii) = 255;
            BoubleTestOne = BoubleTestOne- 250;
            BoubleTestOne = uint8( BoubleTestOne);
            BoubleTestOne = logical( BoubleTestOne);
            MajorAxisB = regionprops(BoubleTestOne,'MajorAxis');
            lengthA = MajorAxisB(1).MajorAxisLength;
            EndpointsBubble = bwmorph(BoubleTestOne,'Endpoints');
            EndpointsBubble = sum(sum(EndpointsBubble));
            if(lengthA>250 && EndpointsBubble<3)
                FillTestN =100;
            end
        end
        end
%           %Dead Neuron test:
%           DeadNeuronOne = BinaryTest + imgNeuriteBrighth;
%           DeadNeuronOne = DeadNeuronOne -1;
%           DeadNeuronOne = uint8(DeadNeuronOne);
%           DeadNeuronOne = logical(DeadNeuronOne);
%           DeadNeuronOne = bwareaopen(DeadNeuronOne,3);
%           AreaNeuronOne = sum(sum(DeadNeuronOne));
%           RatioNeuronDead = AreaNeuronOne/sum(sum(BinaryTest));
%           DeadNeuronOneLabel = bwlabel(DeadNeuronOne);
%           DeadNeuronOneN = max(max(DeadNeuronOneLabel));
%           if(AreaNeuronOne>300&& DeadNeuronOneN>1)
%               FillTestN =100;
%           end    
        if(FillTestN<100)
         BinaryTestOrig = BinaryTest;
        %Now we generate the inverse of this image to use it
        %as a mask:
        BinaryTest = imcomplement(BinaryTest)*255;
        BinaryTest = uint8(BinaryTest);
        %Subtract from original image:
        NeuronImgOne = imgNeurite - BinaryTest;
        %Now we threshold and see on how many nuclei we will
        %obtain:
        level = isodata(imgNeurite);
        if(level*2>0.99)
            level = 0.49;
        end    
        NucTest = im2bw(NeuronImgOne,level*2);
        NucTestOrig = NucTest;
                    if(sum(sum(NucTestOrig))>2000)
                        if((level*4)<1)
                        NucTest = im2bw(NeuronImgOne,level*4);
                        if(sum(sum(NucTest))>2000)
                            NucTest = bwareaopen(NucTest,1000);
                        else    
                        if(sum(sum(NucTest))>300)
                            NucTest = bwareaopen(NucTest,250);
                        else    
                            NucTest = bwareaopen(NucTest,50);
                        end
                        end
                        Reduce = (sum(sum(NucTest)))/(sum(sum(NucTestOrig)));
                        if(Reduce <0.2)
                           NucTest = im2bw(NeuronImgOne,level*3);
                                if(sum(sum(NucTest))>1000)
                                    NucTest = bwareaopen(NucTest,250);
                                else    
                                    NucTest = bwareaopen(NucTest,50);
                                end 
                                if(sum(sum(NucTest)) > 2000)
                                    NucSep = bwmorph(NucTest,'thin',3);
                                    NucSep = bwmorph(NucSep,'open');
                                    NucTest = NucSep;
                                end
                        end        
                        else
                            if((level*3)<1)
                                NucTest = im2bw(NeuronImgOne,level*3);
                                NucTest = bwareaopen(NucTest,50);
                                if(sum(sum(NucTest)) > 2000)
                                    NucSep = bwmorph(NucTest,'thin',3);
                                    NucSep = bwmorph(NucSep,'open');
                                    NucTest = NucSep;
                                end    
                            end    
                        end
                    else    
                        NucTest = bwareaopen(NucTest,100);
                    end
                    NucTestLabel = bwlabel(NucTest);
                    NNucTest = max(max(NucTestLabel));
                    if(sum(sum(NucTest<1000)) &&NNucTest>1)
                        NucTest = im2bw(NeuronImgOne,level*1.5);
                        NucTest = bwareaopen(NucTest,100);
                    end    
                    NucTest = bwareaopen(NucTest,50);
                    if(sum(sum(NucTest)) == 0)
                        NucTest = im2bw(NeuronImgOne,level);
                        NucTest = bwareaopen(NucTest,50);
                    end    
                    SE = strel('disk',3);
                    NucTest = imdilate(NucTest,SE);
                    NucTest = bwmorph(NucTest,'bridge');
                    NucTestLabel = bwlabel(NucTest);
                    NNucTest = max(max(NucTestLabel));
                    if(NNucTest == 0)
                        NucTest = im2bw(NeuronImgOne,level*2);
                        NucTest = bwareaopen(NucTest,10);
                        NucTestLabel = bwlabel(NucTest);
                        NNucTest = max(max(NucTestLabel));
                    end
        %in some instances, we will have a dim neuron adjectant to a bright one. In this case we will simply use a lower threshold and check which particles were overlapping with the higher one:
        if(NNucTest == 1)
        NNucTestLow = im2bw(NeuronImgOne,level);
        NNucTestLow = bwareaopen(NNucTestLow,50);
        if(sum(sum(NNucTestLow))>1000)
            NNucTestLow = bwareaopen(NNucTestLow,250);
        end  
        if(sum(sum(NNucTestLow))>15000)
            NNucTestLow = bwareaopen(NNucTestLow,1000);
        end
        SE = strel('disk',3);
        NNucTestLow = imdilate(NNucTestLow,SE);;
        NNucTestLow = bwmorph(NNucTestLow,'bridge');
        NNucTestLow = bwareaopen(NNucTestLow,50);
        NNucTestLowLabel = bwlabel(NNucTestLow);
        NNucTestLowN = max(max(NNucTestLowLabel));
        for(k=1:NNucTestLowN)
            NNucTestLowOne = NNucTestLowLabel;
            ii = NNucTestLowOne == k;
            NNucTestLowOne(ii) = 255;
            NNucTestLowOne = NNucTestLowOne - 250;
            NNucTestLowOne = uint8(NNucTestLowOne);
            NNucTestLowOne = logical(NNucTestLowOne);
            OverlappNNuc  = NNucTestLowOne + NucTest;
            OverlappNNuc = OverlappNNuc - 1;
            OverlappNNuc = uint8(OverlappNNuc);
            OverlappNNuc = logical(OverlappNNuc);
            if(sum(sum(OverlappNNuc))==0)
                NNucTest = NNucTest + 1;
            end
               
        end
        if(NNucTest == 1)
            SeperateNeuronsValue = max(max(NeuronImgOne));
            SeperateNeuronsValue = SeperateNeuronsValue/4;
            SeperateNeurons = NeuronImgOne - SeperateNeuronsValue;
            SeperateNeurons  = uint8(SeperateNeurons);
            SeperateNeurons = logical(SeperateNeurons);
            SE = strel('disk',3);
            SeperateNeurons = bwareaopen(SeperateNeurons,500);
            SeperateNeurons = imdilate(SeperateNeurons,SE);
            SeperateNeuronsLabel = bwlabel(SeperateNeurons);
            SeperateNeuronsLabelN = max(max(SeperateNeuronsLabel));
            if(SeperateNeuronsLabelN>1)
                NNucTest = SeperateNeuronsLabelN;
            end    
        end    
        end
        if(NNucTest == 1)
            imgbinaryFinal = imgbinaryFinal + BinaryTestOrig;
        end
        end
        
        
        end
        imgbinary = imgbinaryKeep;
        end
        
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
            % We have to make sure holes are not too big like we did above
            Holes = NeuronNucOneProcess6 - NeuronNucOneProcess5;
            Holes = bwareaopen(Holes,HoleArea*2);
            NeuronNucOneProcess6 = NeuronNucOneProcess6 - Holes;
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
                if(BoxSomaArea <500)
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
       if(Auto==1)
       %Here we will check whether the somata a located on one of the
       %candidate Neurons:
       %Since some neurons have holes in the soma, we will add the soma to
       %the final image and then remove single somata without processis:
       imgbinaryFinal = imgbinaryFinal + NeuronNucOneProcess11;
       imgbinaryFinal = logical(imgbinaryFinal);
       imgbinaryFinal = bwareaopen(imgbinaryFinal,4000);
       NeuronNucOneProcess11Label = bwlabel(NeuronNucOneProcess11);
       NucLabel = max(max(NeuronNucOneProcess11Label));
       Candidates = zeros(size(imgbinary));
       for(l=1:NucLabel);
           NucCandidate = NeuronNucOneProcess11Label;
           ii = NucCandidate == l;
           NucCandidate(ii) =255;
           NucCandidate = NucCandidate -250;
           NucCandidate = uint8(NucCandidate);
           NucCandidate = logical(NucCandidate);
           NucCandidateTest = NucCandidate + imgbinaryFinal;
           NucCandidateTest = NucCandidateTest -1;
           NucCandidateTest = uint8(NucCandidateTest);
           NucCandidateTest = logical(NucCandidateTest);
           if(sum(sum(NucCandidateTest))>0)
               Candidates = Candidates + NucCandidateTest;
           end    
       end
       end
       % Now we need to reassess the centroid coordinates from
       % NeuronNucOneProcess11:
        NeuronNucOneProcess11 = logical(NeuronNucOneProcess11);
        if(Auto==1)
        Candidates = logical(Candidates);
        end
        % Create empty matrix for nucleus matrix
        NucleusM = sparse(zeros(size(imgNeurite)));
        NeuronsM = sparse(zeros(size(imgNeurite)));
        if(sum(sum(imgbinary)) >0)
        % If yes determine number of particles
        sR = regionprops(NeuronNucOneProcess11, 'Centroid');
        if(Auto==1)
        sRC = regionprops(Candidates, 'Centroid');
        end
        sRdimensions = size(sR);
        if(Auto==1)
        sRCdimensions = size(sRC);
        end
        else
        % If no, set number to 0
        sRdimensions = 0;
        if(Auto==1)
        sRCdimensions = 0;
        end
        end
        
        % For each particle assess the centroid (x,y coordinates and fill
        % the NucleusM
        for(k=1:sRdimensions)
            x_centroid = sR(k).Centroid(1);
            y_centroid = sR(k).Centroid(2);
            NucleusM(uint16(y_centroid),uint16(x_centroid)) = 1;
        end
        if(Auto==1)
        for(l=1:sRCdimensions)
            x_centroidC = sRC(l).Centroid(1);
            y_centroidC = sRC(l).Centroid(2);
            NeuronsM(uint16(y_centroidC),uint16(x_centroidC)) = 1;
        end
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
        

imgbinary=logical(imgbinary);
binaryimgNeuriteOrig = imgbinary;
%Algorithm to smooth edges of primary binary image from load cortical neurons:
%1) Find edges using canny method:
E = edge(imgbinary,'canny');
% Convert image to logical
E = logical(E);
%2) Add edge image to normal image
imgbinary = imgbinary + E;
%3) This image will have holes which have to be filled. Since crossing
%Neurite will create 'big' holes, the small 'holes' have to be extracted!
% Fill holes in image
Holes = imfill(imgbinary,'holes');
% Get only filled areas
Holes = Holes - imgbinary;
% Extract 'big' holes
HolesBig = bwareaopen(Holes, 20);
% Get small holes, by supstraction the big holes from all holes
Holes = Holes - HolesBig;
% Fill small holes:
imgbinary= imgbinary + Holes;
%Try to close small gaps with bridge function:
binaryimgNeuriteTry = bwmorph(imgbinary,'majority');
binaryimgNeuriteTry = bwareaopen(binaryimgNeuriteTry,4000);
% In some instances the thickening of the binary image will combine two
% binary images resulting in two nuclei but only one binary image.
% Therefore we need to check here, wheter we reduce the number of binary
% particles and if yes to seperate them again:
%1) Get number of particles:
binaryimgNeuriteOrigLabel = bwlabel(binaryimgNeuriteOrig);
binaryimgNeuriteTryLabel = bwlabel(binaryimgNeuriteTry);
nParticlebinaryOrig = max(max(binaryimgNeuriteOrigLabel));
nParticlebinary = max(max(binaryimgNeuriteTryLabel));
if(nParticlebinary>nParticlebinaryOrig)
    imgbinary= bwmorph(imgbinary,'majority');
    imgbinary = bwareaopen(imgbinary,4000);
    
else    
if(nParticlebinary==nParticlebinaryOrig)
    %Create image difference
    imgbinary= bwmorph(imgbinary,'majority');
    imgbinary = bwareaopen(imgbinary,4000);
       
else
    % In some imnstance already the addition of the Edge (E) will cause an
    % overlapp. We have not to seperate the combined areas. We could
    % perform the Edge additon one particle at once to figure out, when the
    % overlapp occurs:
    AreaSep = zeros(size(imgbinary));
    for(y=1:nParticlebinaryOrig)
        binaryimgNeuriteOrigOne = binaryimgNeuriteOrigLabel;
        ii = binaryimgNeuriteOrigOne == y;
        binaryimgNeuriteOrigOne(ii) = 255;
        binaryimgNeuriteOrigOne = binaryimgNeuriteOrigOne -250;
        binaryimgNeuriteOrigOne = uint8(binaryimgNeuriteOrigOne);
        binaryimgNeuriteOrigOne = logical(binaryimgNeuriteOrigOne);
        %Now perform Edge only on one Neuron:
        E = edge(binaryimgNeuriteOrigOne,'canny');
        E = logical(E);
        binaryimgNeuriteOneE = binaryimgNeuriteOrigOne + E;
        Holes = imfill(binaryimgNeuriteOneE,'holes');
        Holes = Holes - binaryimgNeuriteOneE;
        HolesBig = bwareaopen(Holes, 20);
        Holes = Holes - HolesBig;
        binaryimgNeuriteOneE= binaryimgNeuriteOneE + Holes;
        binaryimgNeuriteOneE = bwmorph(binaryimgNeuriteOneE,'majority');
        binaryimgNeuriteOneE = bwareaopen(binaryimgNeuriteOneE,4000);
        % Here we simply thicken the object artifically 
        binaryimgNeuriteOneEFat = bwmorph(binaryimgNeuriteOneE,'thicken',5);
        binaryimgNeuriteOneEFat = bwmorph(binaryimgNeuriteOneEFat,'bridge');
        binaryimgNeuriteOneEFat = bwmorph(binaryimgNeuriteOneEFat,'majority');
        % Now we check of pixels were overlapping with the original binary:
        TestOverlapp = binaryimgNeuriteOrig - binaryimgNeuriteOrigOne;
        TestOverlapp = binaryimgNeuriteOneEFat +  TestOverlapp;
        if(max(max(TestOverlapp))>1)
            TestOverlapp = TestOverlapp - 1;
            TestOverlapp = uint8(TestOverlapp);
            TestOverlapp = logical(TestOverlapp);
            AreaSep = TestOverlapp + AreaSep;
        end  
    end 
    % Now we will subtract the Area Sep from the original binary image and
    % will perform the normal Edge detection!
    binaryimgNeuriteOrig = binaryimgNeuriteOrig - AreaSep;
    binaryimgNeuriteOrig = uint8(binaryimgNeuriteOrig);
    binaryimgNeuriteOrig = logical(binaryimgNeuriteOrig);
    imgbinary = binaryimgNeuriteOrig;
    %Algorithm to smooth edges of primary binary image from load cortical neurons:
    %1) Find edges using canny method:
    E = edge(imgbinary,'canny');
    % Convert image to logical
    E = logical(E);
    %2) Add edge image to normal image
    imgbinary = imgbinary + E;
    %3) This image will have holes which have to be filled. Since crossing
    %Neurite will create 'big' holes, the small 'holes' have to be extracted!
    % Fill holes in image
    Holes = imfill(imgbinary,'holes');
    % Get only filled areas
    Holes = Holes - imgbinary;
    % Extract 'big' holes
    HolesBig = bwareaopen(Holes, 20);
    % Get small holes, by supstraction the big holes from all holes
    Holes = Holes - HolesBig;
    % Fill small holes:
    imgbinary= imgbinary + Holes;
    %Try to close small gaps with bridge function:
    imgbinary = bwmorph(imgbinary,'majority');
    imgbinary = bwareaopen(imgbinary,4000);    
end
end
%In some instances we will split off areas bigger than 4000, which do not
%contain a nucleus and will therefore lead to an error in NeuronTansRat!,
%Therefore we will delete these Areas. We will just run this, if the number of areas of neurons is bigger than number of nuclei:
nParticlePost = bwlabel(imgbinary);
nParticlePost = max(max(nParticlePost));
if(nParticlePost>nParticlebinaryOrig)
KeepNeurons = imgbinary;
imgbinary = zeros(size(imgbinary));
KeepNeuronsLabel = bwlabel(KeepNeurons);
nucleusPicBinary = imgNucleusWatershed;
nucleusPicBinary = logical(nucleusPicBinary);
for(y=1:nParticlePost)
    KeepNeuronOne = KeepNeuronsLabel;
    ii = KeepNeuronOne == y;
    KeepNeuronOne(ii) = 255;
    KeepNeuronOne = KeepNeuronOne -250;
    KeepNeuronOne = uint8(KeepNeuronOne);
    KeepNeuronOne = logical(KeepNeuronOne);
    TestOverFin = KeepNeuronOne + nucleusPicBinary;
    if(max(max(TestOverFin)) == 2)
      imgbinary = imgbinary +   KeepNeuronOne;
    end  
end
end
%Create skeleton by infinite thinning of the binary image
SkeletonImage = bwmorph(imgbinary,'thin',inf);
% Spur removes the last 4 pixel from endpoints. Thereby very small
% processes are eliminated
SkeletonImage = bwmorph(SkeletonImage,'spur',4);
%Convert image to 8-bit
SkeletonImage = uint8(SkeletonImage);
        
        
        
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
         newFile11 = regexprep(newfile, 'space','Skeleton.png','ignorecase');
          newFile12 = regexprep(newfile, 'space','SkeletonOrig.png','ignorecase');
        % Saving of Images
        imwrite(imgNucleus,[foldername1 '/' newFile9]);
        imwrite(imgNucleus,[foldername1 '/' newFile1]);
        imwrite(imgNeurite,[foldername1 '/' newFile2]);
        imwrite(imgNeurite,[foldername1 '/' newFile3]);
        imwrite(imgbinary,[foldername1 '/' newFile4]);
        imwrite(imgNucleusWatershed,[foldername1 '/' newFile5]);
        imwrite(SkeletonImage,[foldername1 '/' newFile11]);
        imwrite(SkeletonImage,[foldername1 '/' newFile12]);
        % Well name gets name from filename
        wellname = (newFile1(1:3));
        % Welllist is filled with well names
        wellList{n} = wellname;
        % Cell centroids are send to csvHandler for Omnisphero GUI
        csvHandler.CellPosMatrix(wellname) = NucleusM;
        csvHandler.ManualPositions1(wellname) = NeuronsM;
        % ScalingFactor is set to 1
        ScalingFactor = 1;
        clear newfile newfile1 newfile2 newfile3 newfile4 newfile5 newfile6 newfile7 newfile8 
        
    end
    
end

% Comments: It would be good to include a suitable filter for bubble rims