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

function [summary Summary2 Summary3 Summary4 Summary5, ExpName, csvHandler,neuronHandler, PercentInd] = MeasurementTransRat ( foldername, csvHandler, neuronHandler, Stats, ReAnalysis)
%Load fitting function for statistical analysis:
if(Stats==1)
sheet = 1;
FitFunctionTable = xlsread('True.xlsx',sheet);
[r c] = size(FitFunctionTable);
b=0;
t=0;
d=0;
k=0;
g=0;
for(l=1:c)
    if(isnan(FitFunctionTable(1,l))==0)
       if(l>1) 
       b= b+1; 
       MaxGraphVal = max(FitFunctionTable(1:31,l));
       FitFunctionTableOne = FitFunctionTable(1:31,l)/MaxGraphVal;
       FitFunctionTableCor(1:31,b) =  FitFunctionTableOne;
       else
           b= b+1;
        FitFunctionTableCor(1:31,b) =  FitFunctionTable(1:31,l);
       end
    end   
end  
% Now we isolate X- Values and y-Values:
XFit = FitFunctionTableCor(1:31,1);
YFit = FitFunctionTableCor(1:31,2:end);
YFit = mean(YFit,2);
%Now we can generate the curve fit:
f=fit(XFit,YFit,'poly9');
%Now we can also get R^2:
[curvefit,gof,output] = fit(XFit,YFit,f,'normalize','on');
sseFit = gof.sse;
RSquareFit = gof.rsquare;
AdRSquareFit = gof.adjrsquare;
xi= (0:10:300);
CurveFitY = f(xi);
Confidence = predint(f,xi);
ConfidenceL = Confidence(1:end,1);
ConfidenceU = Confidence(1:end,2);
%Here we will cross validate our training data set with our fitting model
sheet = 1;
ValidationData = xlsread('True.xlsx',sheet);
XVal = XFit;
YTableVal = ValidationData(1:31,2:end);



[r1 c1] = size(YTableVal);
for(l=1:c1)
    if(isnan(YTableVal(1,l))==0)
       k= k+1; 
       YTableCorVal(1:31,k) =  YTableVal(1:31,l);
    else
        d= d+1; 
       YTableCorVal(1:31,d) = zeros(size(1:31,1)); 
    end   
end  
[r c] = size(YTableCorVal);
e=0;
for(h=1:c)
    YVal = YTableCorVal(1:31,h);
    YOrigVal = YVal;
    %Normalization of Y-coordinates:
    YMaxVal = max(YVal);
    YVal = YVal/YMaxVal;
    sseVal = sum((CurveFitY-YVal).^2);
    RSquareVal = 1-(sseVal/(sum((YVal-mean(YVal)).^2)));
    Summary2Val(2,h) = sseVal;
    Summary2Val(3,h) = RSquareVal;
    if(isnan(sseVal))
    else
       
    g=g+1;
    sseVecVal(g,1) = sseVal;
    
    end
    
    

end
sseMedian = mean(sseVecVal);
sseSTDEV = std(sseVecVal);
sseCrit  = sseMedian + sseSTDEV;
sseCritStrict = sseMedian + (sseSTDEV*0.75);
sseCrit = round(sseCrit,1);
sseCritStrict = round(sseCritStrict,1);
%Now we will check how many of the training objects we will find back:
%Here we will calcualte the sse criteria based on the List:
for(h=1:c)    
     YVal = YTableCorVal(1:31,h);
    YOrigVal = YVal;
    %Normalization of Y-coordinates:
    YMaxVal = max(YVal);
    YVal = YVal/YMaxVal;
    sseVal = sum((CurveFitY-YVal).^2);
     if(sseVal<sseCrit)
        e=e+1;
        SumVal(1:31,e)=YOrigVal(1:31);
        SumNormVal(1:31,e)=YVal(1:31); 
        TabelReval(1:31,h) = YOrigVal(1:31);
        TabelReval(105,h) = sseVal;
        TabelReval(106,h) = RSquareVal;
        TabelValFinal(1:31,e)=YOrigVal(1:31);
     else    
        TabelReval(106,h) = zeros(size(1:31,1)); 
     end 
end    
%Percentage calculation of Identified objects from traing set:
PercentInd = length(TabelValFinal)/length(YTableCorVal)*100;
else
    PercentInd = 0;
end



% Create forder directory
mkdir(foldername);
% Get name of the original file by first generating a directory of the
% parent folder
allfilesOrig = dir(foldername);
% Go throigh all files in parent directory
for(j=1:numel(allfilesOrig))
    % Check of parent folder contains png imageso
    indO = strfind([foldername '/' allfilesOrig(j).name],'.png');
    % if yes get name from this file and keep as Experiment name
if(numel(indO) > 0)
ExpName = allfilesOrig(j).name; 
end
end
% The artifical well naming A01- is removed to only obtain experimental
% name
newSub = ExpName(1:4);
ExpName = regexprep(ExpName,newSub,'','ignorecase');
% foldername2 is defined
foldername2 = [foldername '/ConvertedCellomics'];
% Create forder directory
mkdir(foldername2);
% List of all files in foldername2
allfiles = dir(foldername2);
n=0;
NN = 1;
HH=1;
% This variable defines the thickness of circles used for sholl analysis
RingDistance = 10/0.65;
for(i=1:numel(allfiles))
    %Check if file ends with Skeleton
    ind = strfind([foldername2 '/ConvertedCellomics' allfiles(i).name],'SkeletonOrig');
    if(numel(ind) > 0)%
         n = n + 1
         mm = n + 1;
         % Read images
         SkeletonImage = imread([foldername2 '/' allfiles(i).name]);
         newFileBinary = regexprep(allfiles(i).name,'SkeletonOrig','Binary','ignorecase');  
         binaryimgNeurite = imread([foldername2 '/' newFileBinary]);
         %We have to keep the original image for the neurons:
         binaryimgNeuriteOrig =binaryimgNeurite;
         newFileNeuronNuc = regexprep(allfiles(i).name,'SkeletonOrig','NucleusBigWatershed','ignorecase');
         NeuronNuc = imread([foldername2 '/' newFileNeuronNuc]);
         newFileNeuronImg = regexprep(allfiles(i).name,'SkeletonOrig','NeuriteBig','ignorecase');
         NeuronImg = imread([foldername2 '/' newFileNeuronImg]);
         newFileNucImg = regexprep(allfiles(i).name,'SkeletonOrig','NucleusBig','ignorecase');
         NucleusImg = imread([foldername2 '/' newFileNucImg]);
         newFileNeuriteImg = regexprep(allfiles(i).name,'SkeletonOrig','NeuriteBig','ignorecase');
         NeuriteImg = imread([foldername2 '/' newFileNeuriteImg]);
         %Here we will also load the NeuronCompositeMatrix!
         %Now we have to generate the well name:
         if(n<10)
            a{1} = 'A0';
         else
            a{1} = 'A';
         end 
         a{2} = num2str(n);
         wellName = [a{1} a{2}];
         Dummy = zeros(size(SkeletonImage));
         try
             Recheck = csvHandler.ManualPositions2(wellName);
             Recheck = Recheck + Dummy;
             Recheck = logical(Recheck);
         catch
             Recheck=0;
         end    
         if(ReAnalysis == 1)
             NeuronCandidates = Recheck;
         else    
         NeuronCandidates = csvHandler.ManualPositions1(wellName);
         NeuronCandidates = Dummy+NeuronCandidates;
         NeuronCandidates = logical(NeuronCandidates);
         end
         %We will also read in the original Neuron image
         % Create image names
         if(n<10)
         a{1} = 'A0';
         else
         a{1} = 'A';
         end
         a{2} = int2str(n);
         a{3} = 'space';
         a{4} = '-';
         a{5} = '-++';
         newfile = [a{1} a{2} a{3}];
         newfile1 = [a{4} a{1} a{2} a{5} ];
         newFileC = [a{4} a{1} a{2} ];
         newFileImage = regexprep(ExpName, '.png',newfile1,'ignorecase');
         newFileImage1 = regexprep(ExpName, '.png',newFileC,'ignorecase');
         newFile2 = regexprep(newfile, 'space','NucleusSmall.png','ignorecase');
         wellname = (newFile2(1:3))
         
         % Here we determine new the number of centroids of the remaining
         % nuclei!
        NucleusM = csvHandler.CellPosMatrix(wellname);
        NucleusM = NucleusM + Dummy;
        NucleusM = logical(NucleusM);
        NeuronNuc = logical(NeuronNuc);
        %Here we will delete coordinates outside of cell somata:
        NucleusM = NucleusM+NeuronNuc;
        NucleusM = NucleusM-1;
        NucleusM = uint8(NucleusM);
        NucleusM = logical(NucleusM);
    
        
        %Sometimes we will have remaining structures without a cell nuclei,
        %we will delete them here:
        imgBinaryRemain = zeros(size(SkeletonImage));
        binaryimgNeuriteL = bwlabel(binaryimgNeurite);
        binaryimgNeuriteN = max(max(binaryimgNeuriteL));
        for(h=1:binaryimgNeuriteN);
            binaryimgNeuriteOne = binaryimgNeuriteL;
            ii= binaryimgNeuriteOne == h;
            binaryimgNeuriteOne(ii) = 255;
            binaryimgNeuriteOne = binaryimgNeuriteOne - 250;
            binaryimgNeuriteOne = uint8(binaryimgNeuriteOne);
            binaryimgNeuriteOne = logical(binaryimgNeuriteOne);     
            TestCor = binaryimgNeuriteOne + NucleusM;
            TestCor = TestCor -1;
            TestCor = uint8(TestCor);
            TestCor = logical(TestCor);
            if(sum(sum(TestCor))>0)
               HH = HH + 1;
               imgBinaryRemain = imgBinaryRemain +  binaryimgNeuriteOne;
               NeuriteMass = sum(sum(binaryimgNeuriteOne));
               Summary5{HH,5} = NeuriteMass;
            end
        end
        if(binaryimgNeuriteN==0)
            HH = HH+1;
            Summary5{HH,5}=[];
        end    
        binaryimgNeuriteOrig = logical(imgBinaryRemain);
        binaryimgNeurite = logical(imgBinaryRemain);
        
        %Now do the same for NeuronNuc
        binaryimgNeuriteNucKeep = zeros(size(binaryimgNeurite));
        NeuronNucLabel = bwlabel(NeuronNuc);
        NeuronNucN = max(max(NeuronNucLabel));
        for(r=1:NeuronNucN)
            imgbinaryOne = NeuronNucLabel ;
            ii = imgbinaryOne == r;
            imgbinaryOne(ii)=255;
            imgbinaryOne = imgbinaryOne -250;
            imgbinaryOne = uint8(imgbinaryOne);
            imgbinaryOne = logical(imgbinaryOne);
            TestNuc = imgbinaryOne + logical(NucleusM);
            TestNuc = TestNuc -1;
            TestNuc = uint8(TestNuc);
            TestNuc = logical(TestNuc);
            TestNuc = sum(sum(TestNuc));
            if(TestNuc>0)
                binaryimgNeuriteNucKeep= binaryimgNeuriteNucKeep + imgbinaryOne;
            end
        end
        NeuronNuc = logical(binaryimgNeuriteNucKeep);
        
        %Correct NucImage:
        [r c] = size(NeuronNuc);
        for(e=1:r)
            for(f=1:c)
                if(NeuronNuc(e,f) == 0)
                    NucleusImg(e,f) = NeuriteImg(e,f);
                else
                    NeuronNuc(e,f) = 1;
                    NucleusImg(e,f) = (NeuriteImg(e,f))*5;
                end
            end
        end
        
        
        
        
        
        %Now delete skeletons:
        
        SkeletonImage = logical(SkeletonImage) + binaryimgNeurite;
        SkeletonImage = SkeletonImage -1;
        SkeletonImage = uint8(SkeletonImage);
        SkeletonImage = logical(SkeletonImage);
         % Get number of Neurons
         NumberofNeurons=sum(sum(NucleusM));
        % Remove small particles from skeleton
        SkeletonImage = bwareaopen(SkeletonImage,200);
        % Convert image to 8-bit
        SkeletonImage = uint8(SkeletonImage);
        % Convert image to 8-bit
        NeuronNuc = uint8(NeuronNuc);
        % We will artifically increase the neuron soma to get a better cut
        % from the binary image!
        NeuronNucBig = bwmorph(NeuronNuc,'thicken',2);
        NeuronNucBig = bwmorph(NeuronNucBig,'bridge');
        % Convert image to 8-bit
        NeuronNucBig = uint8(NeuronNucBig);
        %Convert image to logical
        binaryimgNeurite = logical(binaryimgNeurite);
        % Convert image to 8-bit
        binaryimgNeurite = uint8(binaryimgNeurite);
        % The cell soma is removed from the binary image of the neuron for
        % later assessment of neurite mass
        binaryimgNeurite = binaryimgNeurite - NeuronNucBig;
        % Convert image to 8-bit
        binaryimgNeurite = uint8(binaryimgNeurite);
        %Convert image to logical
        binaryimgNeurite = logical(binaryimgNeurite);
        % Convert image to logical       
        SkeletonImage = logical(SkeletonImage);
        % Convert image to 8-bit
        SkeletonImage = uint8(SkeletonImage);
        % Create Raw skeleton without cell soma
        SkeletonWoNuc = SkeletonImage - NeuronNuc;
        % Convert image to 8-bit
        SkeletonWoNuc = uint8(SkeletonWoNuc);
        % Convert image to logical 
        SkeletonWoNuc = logical(SkeletonWoNuc);
        % Convert image to logical 
        SkeletonImage = logical(SkeletonImage);
        
       
        
        % Postprocessing the skeleton. In many cases we will derive artifical branching points as neighbors of endpoints. To circumvent this, this function eleminate this branches as well as endpoints for further counting! 
        % Gets endpoints of SkeletonWoNuc
        % We need to take the original image, since otherwise branching
        % points in the inner part will be removed through cuting of the
        % nucleus
        % Get image will all branching points
        eB = bwmorph(SkeletonImage,'branch');
        % Get image will all end points
        eE = bwmorph(SkeletonImage,'endpoints');
        % For a unknown reason the finding endpoints creates entpoints not
        % located on the skeleton, which we have to eliminate:
        CorrecteE = eE + SkeletonImage;
        CorrecteE = CorrecteE -1;
        CorrecteE = uint8(CorrecteE);
        CorrecteE = logical(CorrecteE);
        CorrecteE = eE - CorrecteE;
        CorrecteE = uint8(CorrecteE);
        CorrecteE = logical(CorrecteE);
        eE = eE - CorrecteE; 
        %In case of very short processess the endpoint will be very close
        %to an branching point
        % Therefore we add end and branching points
        Delete = eE + eB;
        %Perfrom 'bridge' to connect pixles with distance of pixel
        Delete = bwmorph(Delete,'bridge');
        % Now we filter for areas bigger than 2 (End Branch or Branch
        % Branch points close together)
        Delete = bwareaopen(Delete,2);
        % Now we subtract those areas from the branching points to get rid
        % of these artifical points
        eB = eB - Delete;
        % Convert image to 8-bit
        eB = uint8(eB);
        % Convert image to logical
        eB = logical(eB);
        % Now we subtract those areas from the endpoints to get rid
        % of these artifical points
        eE = eE - Delete;
        % Convert image to 8-bit
        eE = uint8(eE);
        % Convert image to logical
        eE = logical(eE);
        % To only obtain the branches and endpoints of the SkeletonWo Nuc
        % we add the skeleton and than substract 1. Only branches within
        % the skeleton will remain:
        eB = eB + SkeletonWoNuc;
        eB = eB -1;
        eB = uint8(eB);
        eB = logical(eB);
        eE = eE + SkeletonWoNuc;
        eE = eE -1;
        eE = uint8(eE);
        eE = logical(eE);
        
        %This function will check for small processis which will be removed
        % Now we generate a skeleton wo branches 
        SkeletonWoBranch = SkeletonImage - eB;
        % Convert image to 8-bit
        SkeletonWoBranch = uint8(SkeletonWoBranch);
        %Convert image to logical;
        SkeletonWoBranch = logical(SkeletonWoBranch);
        % Check if after removal of branches there are new formed
        AdditionalBranch = bwmorph(SkeletonWoBranch,'branch');
        % We keep the additional branches and add those to the skeleton, in
        % order to not produce gaps:
        eBRestore = eB + AdditionalBranch;
        % Also remove those for further fragmentation
        SkeletonWoBranch = SkeletonWoBranch - AdditionalBranch;
        % Convert image to 8-bit
        SkeletonWoBranch = uint8(SkeletonWoBranch);
        %Convert image to logical;
        SkeletonWoBranch = logical(SkeletonWoBranch);
        % Filter for every line bigger than 10 pixels
        ProcessSmall = bwareaopen(SkeletonWoBranch,10);
        % Obtain only small processes
        ProcessSmall = SkeletonWoBranch - ProcessSmall;
        %Index small processes
        ProcessSmallLabelAll = bwlabel(ProcessSmall);
        % Get maximal index
        nP = max(max(ProcessSmallLabelAll));
        % Generate empty image
        ProcessDelete = zeros(size(SkeletonImage));
        %Check for all processes, whether the contain an endpoint and only
        %delete those containing one, in order to not delete inner parts of
        %the skeleton
        for(y=1:nP)
            % Get only one process with pixel value y from the image
            ProcessSmallLabel = ProcessSmallLabelAll;
            ii=ProcessSmallLabel==y;
            ProcessSmallLabel(ii) = 255;
            ProcessSmallLabel = ProcessSmallLabel - 250;
            ProcessSmallLabel = uint8(ProcessSmallLabel);
            ProcessSmallLabel = logical(ProcessSmallLabel);
            %Add endpoints to process
            ProcessSmallLabelB = ProcessSmallLabel + eE;
            %Subtract one from image, so that only endpoints located on the
            %process will remain
            ProcessSmallLabelB = ProcessSmallLabelB -1;
            ProcessSmallLabelB = uint8(ProcessSmallLabelB);
            ProcessSmallLabelB = logical(ProcessSmallLabelB);
            % Check if number is bigger than 0
            if(sum(sum(ProcessSmallLabelB))>0)
               % If yes add process to deletion 
               % However, we have also to check, if we can close the gap
               % after removal of the process, to not generate big gaps!
               TestConnectivity = SkeletonWoBranch - ProcessSmallLabel;
               TestConnectivity = TestConnectivity + eBRestore;
               TestConnectivity = logical(TestConnectivity);
               % Here we need to delete branches which are now 
               % located outside the skeleton:
               TestConnectivity = bwareaopen(TestConnectivity,3);
               TestConnectivity = regionprops(TestConnectivity,'Area');
               TestConnectivity = size(TestConnectivity);
               TestConnectivity = TestConnectivity(1);
               if(TestConnectivity>NumberofNeurons)
               else    
               ProcessDelete = ProcessDelete +  ProcessSmallLabel;
               end
            end
        end
        % Remove small extensions
        SkeletonWoBranch = SkeletonImage - ProcessDelete;
        SkeletonWoBranch = uint8(SkeletonWoBranch);
        SkeletonWoBranch = logical(SkeletonWoBranch);
        % Filter remaining small paricles
        SkeletonWoBranch = bwareaopen(SkeletonWoBranch,10);
        
        % In some instances correction of the skeleton leads to gaps in
        % skeleton. This functions will reconect those.
        % Get number of particles
        testDisconnect = regionprops(SkeletonWoBranch,'area');
        testDisconnect = size(testDisconnect);
        testDisconnect = testDisconnect(1);
        % Check if number of particles exceed number of neurons, meaning
        % fragmentation of the skeleton
        if(testDisconnect>NumberofNeurons)
            % Perform 'bridge to reconnect lines with a gap of one pixel
            SkelRep = bwmorph(SkeletonWoBranch,'bridge');
            % Subtract SkeletonWoBranch (SkeletonImage) from the bridged image 
            SkelRep = SkelRep - SkeletonWoBranch;
            % Index all areas
            SkelRep = bwlabel(SkelRep);
            % Get maximal index
            sKR = max(max(SkelRep));
            %Test for all areas...
            for(y=1:sKR)
                %Only extract pixels of index y from the image
                SkelRepOne = SkelRep;
                ii = SkelRepOne == y;
                SkelRepOne(ii) = 255;
                SkelRepOne = SkelRepOne -250;
                SkelRepOne = uint8(SkelRepOne);
                SkelRepOne = logical(SkelRepOne);
                % Add area to skeleton image
                testRep = SkeletonWoBranch + SkelRepOne;
                % Get number of resulting areas
                testRep = regionprops(testRep, 'Area');
                testRep = size(testRep);
                testRep = testRep(1);
                % Test if introduction of bridge is decreasing the number o
                if(testRep <= testDisconnect)
                    % If it did not, the bridge point is used to fill the
                    % gap
                    SkeletonWoBranch = SkeletonWoBranch +  SkelRepOne;
                end
            end
        end    
        
        
        
        % After correction of the skeleton, we will now need to perform
        % assessment of branching and endpoints again, like above! And get
        % rid of artifial branch end endpoints
        % Get number of branching points from refined skeleton
        eB = bwmorph(SkeletonWoBranch,'branch');
        % Get number of endpoints from refined skeleton
        eE = bwmorph(SkeletonWoBranch,'endpoints');
        % For a unknown reason the finding endpoints creates entpoints not
        % located on the skeleton, which we have to eliminate:
        CorrecteE = eE + SkeletonWoBranch;
        CorrecteE = CorrecteE -1;
        CorrecteE = uint8(CorrecteE);
        CorrecteE = logical(CorrecteE);
        CorrecteE = eE - CorrecteE;
        CorrecteE = uint8(CorrecteE);
        CorrecteE = logical(CorrecteE);
        eE = eE - CorrecteE; 
        % Branches directly next to each other are eliminated
        DeleteB = bwareaopen(eB,2);
        eB = eB - DeleteB;
        Delete = eE + eB;
        Delete = bwmorph(Delete,'bridge');
        Delete = bwareaopen(Delete,2);
        % Get branches on the 2 or more pixel component
        DeleteB = eB + Delete;
        DeleteB = DeleteB -1;
        DeleteB = uint8(DeleteB);
        DeleteB = logical(DeleteB);
        % Index all branches on Delete
        DeleteLabelB = bwlabel(DeleteB);
        DeleteLabel = bwlabel(Delete);
        % Get max index
        sDel = max(max(DeleteLabel));
        % Check for each branch, ......
        %We will create a matrix containing all endpoints wich will be
        %deleted during the process!
        eEDeleteFinal = logical(zeros(size(eE)));
        eERescue = logical(zeros(size(eE)));
        eBRescue = logical(zeros(size(eE)));
        % We also keep the endpoints before post processing to figure out,
        % which we lose:
        eEOrig = eE;
        for(k=1:sDel)
        % Get only pixel from DeleteLabel and DeleteLabelB with index k    
        DeleteOne = DeleteLabel;
        DeleteBOne = DeleteLabelB;
        ii = DeleteOne ==k;
        DeleteOne(ii) = 255;
        ii = DeleteBOne ==k;
        DeleteBOne(ii) = 255;
        DeleteOne = DeleteOne -250;
        DeleteOne= uint8(DeleteOne);
        DeleteOne=logical(DeleteOne);
        DeleteBOne = DeleteBOne -250;
        DeleteBOne= uint8(DeleteBOne);
        DeleteBOne=logical(DeleteBOne);
        %After removing the branch we need to check, whethere there is a
        %new endpoints (just an artifact in the skel line) or if there is
        %no new endpoint (short branch!):
        % We remove the Branching point
        SkelTest = SkeletonWoBranch - DeleteBOne;
        % Filter for only particles bigger than 100
        SkelTest = bwareaopen(SkelTest,10);
        % Get endpoints from filtered Skeleton
        eETest = bwmorph(SkelTest,'endpoints');
        % Check if this removal of particles leads to a different number of
        % endpoints (Would e.g. happen if we delete somethin in the line
        % points, creating a gap presenting two new endpoints)
        TestNewEnd = eETest - eE;
        TestNewEnd = uint8(TestNewEnd);
        TestNewEnd = logical(TestNewEnd);
        % Get number of endpoints
        TestNewEnd = sum(sum(TestNewEnd));
        % Case difference
        %1) If there are 2 new endpoints, this means the introduction of a
        %gap
        if(TestNewEnd==2)
            % In this case the skeleton is repaired
            SkeletonWoBranch = SkelTest;
            SkeletonWoBranch = SkeletonWoBranch + DeleteBOne;
            % However, this still means an artifuical branch and therefore
            % it is deleted from the list of branching points
            eB = eB - DeleteOne;
            eB = uint8(eB);
            eB = logical(eB);
            % Same hold true for endpoints
            eE = eE - DeleteOne;
            eE = uint8(eE);
            eE = logical(eE);
            eEDelete = eEOrig - eE;
            eEDeleteFinal = eEDeleteFinal + eEDelete;
        else
            % If number of endpoints if unchaged, this means that a short
            % process is eliminated. Skeleton is not repaired
            if(TestNewEnd==0)
            SkeletonWoBranch = SkelTest; 
            % Again artifical endpoints and branches are eliminated from
            % the list
            eB = eB - DeleteOne;
            eB = uint8(eB);
            eB = logical(eB);
            eE = eE - DeleteOne;
            eE = uint8(eE);
            eE = logical(eE);
            eEDelete = eEOrig - eE;
            eEDeleteFinal = eEDeleteFinal + eEDelete;
            % In some instances we will create new endpoints. We have to
            % identify them and add them to eE!
            NeweE = bwmorph(SkeletonWoBranch,'endpoints');
            TestNewE = NeweE + DeleteOne;
            TestNewE = TestNewE -1;
            TestNewE = uint8(TestNewE);
            TestNewE = logical(TestNewE);
            if(sum(sum(TestNewE))>0)
                eERescue = eERescue + TestNewE;
            end    
            else
            % If number of endpoints is not altered, this means a branch is
            % too close to the endpoint and only the branch is eliminated
            eB = eB - DeleteBOne;
            eB = uint8(eB);
            eB = logical(eB);
            %However, we can also obtain a value of one for a gap, if there was a small branch of one pixel counted as an endpoint.
            %Therefore we will check if the removed endpoint is adjactend
            %to any endpoint, resulting in an area of 2 pixel. If this is
            %the case we will delete the endpoint from both the endpoint
            %image as well as from the skeleton!
            testEndBranch = eE + DeleteBOne;
            testEndBranchKeep = testEndBranch;
            testEndBranch = bwareaopen(testEndBranch,2);
            testEndBranch = testEndBranch - DeleteBOne;
            testEndBranch = uint8(testEndBranch);
            testEndBranch = logical(testEndBranch);
                if(sum(sum(testEndBranch))>0)
                    eE = eE - testEndBranch;
                    eE = uint8(eE);
                    eE = logical(eE);
                    eB = eB - testEndBranch;
                    eB = uint8(eB);
                    eB = logical(eB);
                    eEDelete = eEOrig - eE;
                    eEDeleteFinal = eEDeleteFinal + eEDelete;
                    SkeletonWoBranch = SkeletonWoBranch - testEndBranch;
                    % Here we need to chekc if removal of the endpoint
                    % creates an new branch or if there was a real branch,
                    % which got deleted:
                    eBTest = bwmorph(SkeletonWoBranch,'branch');
                    eBTest = eBTest + testEndBranchKeep;
                    eBTest = eBTest - 1;
                    eBTest = uint8(eBTest);
                    eBTest = logical(eBTest);
                    if(sum(sum(eBTest))>0)
                        eBRescue = eBRescue + eBTest;
                    end    
                end 
            end
            end
        end
        %Finnally we have to check if there are new endpoints created after
        %cutting processess and add those. We will however delete again all
        %endpoints deleted in the above postprocessing:
        eE = bwmorph(SkeletonWoBranch,'endpoints');
        eE = eE - eEDeleteFinal;
        eE = uint8(eE);
        eE = logical(eE);
        eE = eE +eERescue;
        eB = eB + eBRescue;
        % Now we finaly prune small processes. The idea of this algorithm
        % is to use the spur function to shorten the skeleton from each
        % endpoint, and to keep only those endpoints which can be shortened
        % 10 pixles by this procedure:
        %First do this for very small branches (This is important if we have two endpoints close to each other:
        SkelSpur = bwmorph(SkeletonWoBranch,'spur',3);
        ParticleSpur = SkeletonWoBranch - SkelSpur;
        ParticleSpurBig = bwareaopen(ParticleSpur,3);
        ParticleSpur = ParticleSpur - ParticleSpurBig;
        SkeletonWoBranch = SkeletonWoBranch - ParticleSpur;
        %Then repeat for bigger branch
        SkelSpur = bwmorph(SkeletonWoBranch,'spur',10);
        ParticleSpur = SkeletonWoBranch - SkelSpur;
        ParticleSpurBig = bwareaopen(ParticleSpur,10);
        ParticleSpur = ParticleSpur - ParticleSpurBig;
        SkeletonWoBranch = SkeletonWoBranch - ParticleSpur;
        NeuronNuc = logical(NeuronNuc);
        %In case single pixel remain we will finally filter them here
        SkeletonWoBranch = bwareaopen(SkeletonWoBranch,2);
        % Now we check a last time for new branches uiquly added this time
        % and add them eB:
        eBFinal = bwmorph(SkeletonWoBranch,'branch');
        eBFinalDif = eBFinal + eB;
        eBFinalDif = eBFinalDif -1;
        eBFinalDif = uint8(eBFinalDif);
        eBFinalDif = logical(eBFinalDif);
        eBFinal = eBFinal - eBFinalDif;
        eBFinalBig = bwareaopen(eBFinal,2);
        eBFinal = eBFinal - eBFinalBig;
        eBFinal = uint8(eBFinal);
        eBFinal = logical(eBFinal);
        eB = eBFinal + eB;
        SkeletonWoBranchWoNuc = SkeletonWoBranch - NeuronNuc;
        SkeletonWoBranchWoNuc = uint8(SkeletonWoBranchWoNuc);
        SkeletonWoBranchWoNuc = logical(SkeletonWoBranchWoNuc);
        % This procedure might lead to branches now located on the skelton
        % line (process spured away!). Therefore we have to identify those
        % and delete them from eB. Therefore we will first assess branches
        % in the new Skeleton and check if there are some in eB not present
        % in the new one and delte those!
        eBAfter = bwmorph(SkeletonWoBranch,'branch');
        %This return common branches!
        eBTest = eB + eBAfter;
        eBTest = eBTest -1;
        eBTest = uint8(eBTest);
        eBTest = logical(eBTest);
        % Now we can get all other than common as artifical branches:
        eBSub = eB - eBTest;
        eB= eB - eBSub;
        % To only obtain the branches and endpoints of the SkeletonWo Nuc
        % we add the skeleton and than substract 1. Only branches within
        % the skeleton will remain:
        eB = eB + SkeletonWoBranchWoNuc;
        eB = eB -1;
        eB = uint8(eB);
        eB = logical(eB);
        % In a final step we have to check if branches, located close to
        % each other were deleted as the result of deletion of end and
        % branches close together:
        
        
        
        
        % Branches were thickend to be better visible
        BranchesSave = bwmorph(eB,'thicken', 2);
        % The inner pixel of the circle is removed
        BranchesSave = bwmorph(BranchesSave,'remove', 1);
        BranchesSave = uint8(BranchesSave*255);
        eE = eE + SkeletonWoBranchWoNuc;
        eE = eE -1;
        eE = uint8(eE);
        eE = logical(eE);
        % Endpoints were thickend to be better visible
        EndpointsSave = bwmorph(eE,'thicken', 2); 
        EndpointsSave = uint8(EndpointsSave*255);
        % Sum endpoints and branchpoints
        EndpointBranchSave = BranchesSave + EndpointsSave;
        % Save branch and endpoint image
        newFileEndBranch = regexprep(allfiles(i).name,'SkeletonOrig','OligoSmall','ignorecase');
        newFileEndBranch2 = regexprep(allfiles(i).name,'SkeletonOrig','OligoBig','ignorecase');
        newFileSkel3 = regexprep(allfiles(i).name,'SkeletonOrig','Skeleton','ignorecase');
        imwrite(EndpointBranchSave,[foldername2 '/' newFileEndBranch]);
        imwrite(EndpointBranchSave,[foldername2 '/' newFileEndBranch2]);
        % Save final Skeleton Image
        imwrite(SkeletonWoBranchWoNuc,[foldername2 '/' newFileSkel3]);
        
        %This function calculates the number of primary tips going of from
        %the the soma:
        % Here we will simply add the SkeletonWoBranch to our binary
        % Nucleus image and subtract 1 to only obtain the skeleton within
        % this Area:
        PrimaryProcess = SkeletonWoBranch + NeuronNuc;
        PrimaryProcess = PrimaryProcess -1 ;
        PrimaryProcess = uint8(PrimaryProcess);
        PrimaryProcess = logical(PrimaryProcess);
        % Now we calculate the number of Endpoints on this skeletons:
        PrimaryProcess = bwmorph(PrimaryProcess,'endpoints');
        % Now we summ all endpoints of all Neuron somas
        PrimaryProcess = sum(sum(PrimaryProcess));
        
        
        
        
        %Morphological extracttion
        %1) Neurite Mass:
        NeuriteMass = (sum(sum(binaryimgNeurite)))/NumberofNeurons;
        %2) Total Neurite length:
        TotalNeuritelength = ((sum(sum(SkeletonWoBranchWoNuc)))/NumberofNeurons);
        %3) Number of Neurites:
        NumberofNeurites = PrimaryProcess/NumberofNeurons;
        %4) Average Neurte length:
        AverageNeuritelength = TotalNeuritelength/NumberofNeurites;
        %5) Number of branching points:
        NumberBranches = sum(sum(eB))/NumberofNeurons;
        %6) Number of endpoints:
        endpointSkel = sum(sum(eE))/NumberofNeurons;
        %7) Neurite Mass Axons:
        
        
        
        %Generate matrix for neurons identified with Stats apporach
        if(Stats == 1)
            NeuronStats = zeros(size(NeuronNuc));
        else
            NeuronStats = NeuronCandidates;
        end    
        %Performing sholl analysis
        NeuronNuc = logical(NeuronNuc);
        NeuronStatsRe = zeros(size(NeuronNuc));
        % Get number of neuron nuclei
        nNuc = regionprops(NeuronNuc,'Area');
        nNuc = size(nNuc);
        nNuc = nNuc(1);
        % Check if we have more particles then neuron nuclei
        if(nNuc>NumberofNeurons);
            NeuronNuc = logical(NeuronNuc);
            % Sometimes there will be more particles resulting from the
            % cutting of the soma in the preprocessing! We will filter very
            % small particles out 
            NeuronNuc = bwareaopen(NeuronNuc,5);
            NeuronNuc = bwlabel(NeuronNuc);
            % Get number of neuron nuclei
            nNuc = max(max(NeuronNuc));
        else
            NeuronNuc = bwlabel(NeuronNuc);
        end 
        %We also have to make sure, the SkeletonWoBranch does not contain
        %more particles than NumberofNeurons
        LabelSkel = bwlabel(SkeletonWoBranch);
        NSkel= max(max(LabelSkel));
        if(NSkel>NumberofNeurons);
           %Idea: We will subtract SkeletonWoBranch from SkeletonImage and will chekc for all remaining points, whether they can reconnect the skeleton!
           Reconnect = SkeletonImage - SkeletonWoBranch ;
           Reconnect = uint8(Reconnect);
           Reconnect = logical(Reconnect);
           %Now we need to label all areas in Reconnect and need to chekc
           %individually whether they can reconnect two skeleton parts!
           ReconnectL = bwlabel(Reconnect);
           ReconnectN = max(max(ReconnectL));
           for(q=1:ReconnectN)
               ReconnectOne = ReconnectL;
               ii = ReconnectOne == q;
               ReconnectOne(ii) = 60000;
               ReconnectOne = ReconnectOne -59995;
               ReconnectOne = uint8(ReconnectOne);
               ReconnectOne = logical(ReconnectOne);
               TestRecon = SkeletonWoBranch + ReconnectOne;
               TestRecon = logical(TestRecon);
               TestReconL = bwlabel(TestRecon);
               TestReconN = max(max(TestReconL));
               if(TestReconN<NSkel)
                   SkeletonWoBranch = SkeletonWoBranch + ReconnectOne;
               end    
           end    
        end  
        
        
        
        
        % We add here the skeleton, since the indexing is dependant on the
        % first x-coordinate which might differe between skeletons and
        % neurite  somas!
        LabelNuc = NeuronNuc + SkeletonWoBranch;
        % In some instances erosion of binary image to obtain nucleus
        % results in a nucleus which is not overlapping the skeleton. In
        % this case we will just increase nucleus size by thickening it:
        for(d=1:nNuc)
            LabelNucTest = NeuronNuc;
            ii = LabelNucTest == d;
            LabelNucTest(ii) = 255;
            LabelNucTest = LabelNucTest -250;
            LabelNucTest = uint8(LabelNucTest);
            LabelNucTest = logical(LabelNucTest);
            % We Overlapp each LabelNuc with the original skeleton
            TestOverlapp = LabelNucTest + SkeletonWoBranch;
            % One is subtracted. If the nuclei was loacted on the skeleton
            % there shoud be pixle left
            TestOverlapp = TestOverlapp -1;
            TestOverlapp = uint8(TestOverlapp);
            TestOverlapp = logical(TestOverlapp);
            if(sum(sum(TestOverlapp)) == 0)
                % If the area of the overlapp is 0, we thicken the
                % NeuronNuc
                LabelNucTest = bwmorph(LabelNucTest,'thicken',10);
                LabelNuc = LabelNuc + LabelNucTest;
                NeuronNuc = NeuronNuc + LabelNucTest;
            end
        end   
          
        LabelNuc = uint8(LabelNuc);
        LabelNuc = logical(LabelNuc);
        % Index all neuron somas
        LabelNuc = bwlabel(LabelNuc);
        % Label all Skeletons
        LabelNeuron = SkeletonWoBranch;
        LabelNeuron = uint8(LabelNeuron);
        LabelNeuron = bwlabel(LabelNeuron);
        %Generate label for Neurons:
        LabelbinaryimgNeuriteOrig = bwlabel(binaryimgNeuriteOrig);
        [r c] = size(NeuronNuc);
        % At the centroind number of intersections are 0
        TotalNumberIntersections = 0;
        AverageIntersections = 0;
        %Defines the maximal ring number in sholl analysis
        Ringnumber = 100;
        % Empty images are created
        RingLable = zeros(size(SkeletonImage));
        RingLable = uint8(RingLable);
        RingMask = zeros(size(SkeletonImage));
        RingMask = uint8(RingMask);
        %Cut our inner circle   
       for(h=1:NumberofNeurons)
           %Now we define the inner Ring, Which again need the the complete
           %Neurons with soma!
           InnerCircle = LabelNeuron;
           % Now we extract only one skeleton
           ii = InnerCircle == h;
           InnerCircle(ii) =255;
           InnerCircle = InnerCircle -250;
           InnerCircle = uint8(InnerCircle);
           InnerCircle = logical(InnerCircle);
           SkeletonOne = InnerCircle;
           LabelNucOne = LabelNuc;
           ii=LabelNucOne==h;
           LabelNucOne(ii)=255;
           % Here we add the Neuron Nuc, so that the neuron soma contained
           % in the indexed LabekNucOne has a pixel value of 256
           LabelNucOne = LabelNucOne + NeuronNuc;
           % Subtracting 255 results in only the area of the skeleton within the neuron nuclei
           LabelNucOne = LabelNucOne - 255;
           LabelNucOne = uint8(LabelNucOne);
           LabelNucOne = logical(LabelNucOne);
           %Get number centroid of LabelNucOne         
           s = regionprops(LabelNucOne,'Centroid');
           % Get closest suroounding box arround LabelNucOne
           sAS = regionprops(LabelNucOne,'BoundingBox');
           % Getting the axis of the InnerCircle used for the radius of the
           % inner circle
           Axis = sAS(1).BoundingBox(3);
           Axis2 = sAS(1).BoundingBox(4);
           if(Axis>Axis2)
               Axis = Axis;
           else
               Axis = Axis2;
           end 
           Axis = RingDistance*2;
           %Now we cut the inner circle (Cell soma)
           % Here we create a mask only containing the inner circle, so we do not add any further ring!
           [r c] = size(InnerCircle);
           %Get x and y coordinates from the centroid
            x_centroid = s(1).Centroid(1);
            y_centroid = s(1).Centroid(2);
           radIni = (Axis/2);
            for(x=1:r)
                for(y=1:c)
                if((x-x_centroid)^2  +(y-y_centroid)^2 > radIni^2)
                    InnerCircle(y,x) = 0;
                end 
                end
            end
            SkeletonOne = SkeletonOne - InnerCircle; 
           
           % Cretae header for all summaries
            NeuronN = int2str(h);
            NN = NN + 1;
            newFileSingleNeuron = regexprep(newFileImage, '++',NeuronN,'ignorecase');
            %Summary2
            Summary2{1,1} = 'Radius';
            Summary2{1,NN} = newFileSingleNeuron;
            %Summary3
            Summary3{1,1} = 'Radius';
            Summary3{1,NN} = newFileSingleNeuron;
            %Summary4
            Summary4{1,1} = 'Radius';
            Summary4{1,NN} = newFileSingleNeuron;
            %Here we can save endpoints like total length, average length,
            %number of branches , number of tips etc...
            TotalLengthSingle = sum(sum(SkeletonOne));
            %The number of Tips can be assessed by adding the eE to
            %SkeletonOne and subtract one:
            EndpointsSingle = SkeletonOne + eE;
            EndpointsSingle = EndpointsSingle -1;
            EndpointsSingle = uint8(EndpointsSingle);
            EndpointsSingle = logical(EndpointsSingle);
            EndpointsSingle = sum(sum(EndpointsSingle));
            BranchSingle = SkeletonOne + eB;
            BranchSingle = BranchSingle -1;
            BranchSingle = uint8(BranchSingle);
            BranchSingle = logical(BranchSingle);
            BranchSingle = sum(sum(BranchSingle));
            Summary5{1,1} = 'Image Name';
            Summary5{1,2} = 'Total Length';
            Summary5{1,3} = 'Number of Tips';
            Summary5{1,4} = 'Number of Branches';
            Summary5{1,5} = 'Neurite Mass';
            Summary5{NN,1} = newFileSingleNeuron;
            if(TotalLengthSingle>0)
            Summary5{NN,2} = TotalLengthSingle;
            Summary5{NN,3} = EndpointsSingle;
            Summary5{NN,4} = BranchSingle;
            else
            Summary5{NN,2} = [];
            Summary5{NN,3} = [];
            Summary5{NN,4} = [];
            Summary5{NN,5} = [];
            end
            % This variable chekcs, if there are still intersections with
            % the next ring.
            Fill = 1;
            m = 0; 
        % This is used to cut out a ring from the skeleton wo nuc         
        CutDummy = SkeletonOne;  
        InnerCircle =uint8(InnerCircle);
        % Create an empty image for final heat map image (Being the inner
        % part 1 and then for each ring n+1 higher in intensity)
        RingNeuron = zeros(size(SkeletonImage));
        RingNeuron = uint8(RingNeuron);
        % The inner circle is added to the final heat map image
        RingNeuron = RingNeuron + InnerCircle;
        InnerCircle = logical(InnerCircle);    
        p = 1;    
            % We now iterate over all ring numbers 
            for(b=1:Ringnumber)
                b = b
                % We check if Fill is >0
            if(Fill>0)
            p = p + 1; 
            % The actual Ring is set to the CutDummy
            Ring = CutDummy;
            % This gives the maximal radius of the n ring (in the first
            % instane the inner circle is removed and now everything
            % outside this ring is set to 0, resulting in a ring image
            rad = (Axis/2)+(RingDistance*b);
            for(x=1:r)
                for(y=1:c)
                if((x-x_centroid)^2  +(y-y_centroid)^2 > rad^2)
                    Ring(y,x) = 0;
                end 
                end
            end
            % Create an empty image with white Rings:
            % Angels of the ring
            th = 0:pi/5000:2*pi;
            % Since we need also the inner ring we start with ringnumber -1
            l = b -1;
            radMask = (Axis/2)+(RingDistance*l);
            sth = length(th);
            % Draw white circles          
            for(t=1:sth)
                thSingle = th(t);
                xW = radMask * cos(thSingle) + x_centroid;
                yW = radMask * sin(thSingle) + y_centroid;
                xW = round(xW);
                yW = round(yW);
                if(xW < 1 | xW > r)
                    xW = 0;
                end
                if(yW < 1 | yW > c)
                    yW = 0;
                end
                if(xW > 0 && yW > 0)
                RingMask(yW,xW) = 1;
                end
            end    
            % This is only done for the first circle    
            if(b==1)
                % The number of intersections are the number of endpoints
                % of the inner circle
                Intersection = bwmorph(InnerCircle,'endpoints');
                % Now we add the Ring image to the innercircle
                IntersectionInter = Intersection + Ring;
                % We eliminate all particles smaller 2 to delte non
                % adjactant endpoints (Sometime processes are so small that
                % they did not corss the first ring)
                IntersectionInter = bwareaopen(IntersectionInter,2);
                % We subtract the ring image, so that only endpoints on
                % lines on the ring image remain, which have to be intersection points       
                Intersection = IntersectionInter - Ring;
                % We mark the intercetion points by making circles out of
                % the points
                IntersectionAdd = bwmorph(Intersection,'thicken',2);
                IntersectionAdd = bwmorph(IntersectionAdd,'remove');
                IntersectionAdd = logical(IntersectionAdd);
                % Add the points to the circle
                IntersectionAddDummy = zeros(size(InnerCircle));
                IntersectionAddDummy = IntersectionAdd(1:c,1:r);
                IntersectionAdd = IntersectionAddDummy;                
                RingMask = logical(RingMask);
                RingMask = RingMask + IntersectionAdd;
                % Get number of intersections
                Intersection = sum(sum(Intersection));
            else
                % For all later circles check if there are neurites within
                % the ring
                if(sum(sum(Ring))>0)
                % Here we always overly endpoints from the previous circle with the new one to verifiy intersections    
                Intersection = bwmorph(RingNeuron,'endpoints');
                IntersectionInter = Intersection + Ring;
                IntersectionInter = bwareaopen(IntersectionInter,2);
                Intersection = IntersectionInter - Ring;
                Intersection = uint8(Intersection);
                Intersection = logical(Intersection);
                IntersectionAdd = bwmorph(Intersection,'thicken',2);
                IntersectionAdd = bwmorph(IntersectionAdd,'remove');
                IntersectionAdd = logical(IntersectionAdd);
                IntersectionAddDummy = zeros(size(InnerCircle));
                IntersectionAddDummy = IntersectionAdd(1:c,1:r);
                IntersectionAdd = IntersectionAddDummy;
                RingMask = logical(RingMask);
                RingMask = RingMask + IntersectionAdd;
                Intersection = bwlabel(Intersection);
                Intersection = max(max(Intersection)) ;  
                %Sometimes a circle is touching a ring, not presenting an
                %endpoint. However, if this happens , we will count the
                %number of intersections as one to circumvent stops!
                if(Intersection == 0)
                    Intersection =1;
                end    
                else
                % If rings are empty we will set number of Intersections to
                % 0 and the loop will stop
                Intersection = 0;   
                end
            end
            % The ring will be cutted from Cut Dummy, and therefore giving already the inner ring for the next ring           
            CutDummy = CutDummy - Ring;
            Ring = uint8(Ring);          
            Ring = logical(Ring);
            % Get number of Intersection
            NumberofIntersections = Intersection;
            % The length added per ring:
            LengthAdd = sum(sum(Ring));
            % The branches added per ring
            Branches = bwmorph(Ring,'Branchpoints');
            NumberofSections = Intersection;
            % Sum up all intersection throughout all circles
            TotalNumberIntersections = TotalNumberIntersections + NumberofIntersections;
            % Number of all branches
            NumberofBranches = sum(sum(Branches));
            % Sets Fill to NumberofIntersections
            Fill = NumberofIntersections;
           % This fills values in matrix (for number of intersections we
           % start with ring 0, which is inner circle
            RN = p -1;
            q = p +1;
            %Summary 2
            Summary2{2,1} = 0;
            Summary2{q,1} = RN*10;
            Summary2{p,NN} = NumberofIntersections(1);
            %Summary 3
            Summary3{2,1} = 0;
            Summary3{q,1} = RN*10;
            Summary3{2,NN} = 0;
            Summary3{q,NN} = LengthAdd;
            %Summary 4
            Summary4{2,1} = 0;
            Summary4{q,1} = RN*10;
            Summary4{2,NN} = 0;
            Summary4{q,NN} = NumberofBranches(1);
            % Ring is added to the final heat map image
            Ring = uint8(Ring);
            RingNeuron = uint8(RingNeuron);
            RingNeuron = RingNeuron + (Ring*(b+1));
            end
            end
            % All heat mapped skeletons are combined in one image
            RingLable = RingLable + RingNeuron;
            if(Stats==1)
            %Here we test if neuron is similar to the model function:
            %Therefore we extract the Neuron from Summary2:
            NeuronVal = Summary2(2:end,NN);
            NeuronVal = cell2mat(NeuronVal);
            INTNumber = sum(NeuronVal);
            NeuronValOrig = NeuronVal;
            try
            NeuronValRest = NeuronValOrig(32:end);
            MaxRest = max(NeuronValRest);
            catch
            MaxRest = 0;
            end
            %Normalization of Sholl plot;
            YMaxNeuronSingle = max(NeuronVal);
            NeuronVal = NeuronVal/YMaxNeuronSingle;
            LN = length(NeuronVal);
            if(LN>31)
                NeuronVal = NeuronVal(1:31,1);
                LN = 31;
            end    
            DummyV = zeros(31,1);
            DummyV(1:LN) = NeuronVal;
            NeuronVal = DummyV;
            sse = sum((CurveFitY-NeuronVal).^2);
            end
            %Now we have to check whether, the identified Neurons is
            %one of those chosen as candidates:
            TestCandidate = NeuronCandidates + LabelNucOne;
            TestCandidate = TestCandidate -1;
            TestCandidate = uint8(TestCandidate);
            TestCandidate = logical(TestCandidate);
            %For large structures we do not take the sse into account:
            if(Stats==1)
            if(INTNumber>165)
                sse = 0.1;
            end    
            if(INTNumber<100)
                sseCritCriteria = sseCritStrict;
            else
                sseCritCriteria = sseCrit;
            end
            end
            if(sum(sum(TestCandidate))>0)
            if(Stats==1)    
            if(sse<sseCritCriteria)
                if(YMaxNeuronSingle>5&&INTNumber>74)
                %Now we need to generate new matrix for Neurons! Manual 1!
                NeuronStatsRe = NeuronStatsRe + TestCandidate;
%                 SummarySSE(105,NN) = sse;
%                 SummarySSE(106,NN) = YMaxNeuronSingle;
%                 SummarySSE(107,NN) = INTNumber;
%                 SummarySSE(108,NN) = 1;
                else
                Summary2(:,NN) = [];
                Summary3(:,NN) = [];
                Summary4(:,NN) = [];
                Summary5{NN,2} = [];
                Summary5{NN,3} = [];
                Summary5{NN,4} = [];
                Summary5{NN,5} = [];
%                 SummarySSE(105,NN) = sse;
%                 SummarySSE(106,NN) = YMaxNeuronSingle;
%                 SummarySSE(107,NN) = INTNumber;
%                 SummarySSE(108,NN) = 1;
                %SummarySSE(35,NN) = [];
                end
                
            else
                
                %Here we will delete unsimilar Neurons
                Summary2(:,NN) = [];
                Summary3(:,NN) = [];
                Summary4(:,NN) = [];
                Summary5{NN,2} = [];
                Summary5{NN,3} = [];
                Summary5{NN,4} = []; 
                Summary5{NN,5} = [];
                %SummarySSE(35,NN) = [];
%                 SummarySSE(105,NN) = sse;
%                 SummarySSE(106,NN) = YMaxNeuronSingle;
%                 SummarySSE(107,NN) = INTNumber;
%                 SummarySSE(108,NN) = 1;
            end
            else
                NeuronStatsRe = NeuronStatsRe + TestCandidate;
            end    
            else
                Summary2(:,NN) = [];
                Summary3(:,NN) = [];
                Summary4(:,NN) = [];
                Summary5{NN,2} = [];
                Summary5{NN,3} = [];
                Summary5{NN,4} = []; 
                Summary5{NN,5} = [];
%                 SummarySSE(105,NN) = sse;
%                 SummarySSE(106,NN) = YMaxNeuronSingle;
%                 SummarySSE(107,NN) = INTNumber;
%                 SummarySSE(108,NN) = 0;
                %SummarySSE(35,NN) = [];
            end    
            end
            
       
       % Save ring mask     
       newfileRing = regexprep(allfiles(i).name,'SkeletonOrig','AstroBig','ignorecase');
       newfileRing2 = regexprep(allfiles(i).name,'SkeletonOrig','AstroSmall','ignorecase');
       imwrite(RingMask,[foldername2 '/' newfileRing]);
       imwrite(RingMask,[foldername2 '/' newfileRing2]);
       %We also need to rewrite
       %Binary,NucleusBig,NucleusSmall,NucleusBigWatershed, 
       imwrite(NucleusImg,[foldername2 '/' newFileNucImg]);
       newFileNucImgSmall = regexprep(allfiles(i).name,'SkeletonOrig','NucleusSmall','ignorecase');
       imwrite(NucleusImg,[foldername2 '/' newFileNucImgSmall]);
       imwrite(NucleusImg,[foldername2 '/' newFileNucImg]);
       imwrite(NeuronNuc,[foldername2 '/' newFileNeuronNuc]);
       imwrite(binaryimgNeuriteOrig,[foldername2 '/' newFileBinary]);
              
       % If there are no neurons, the Average number of intersections needs
       % to be excluded since we would divide by zero
       AverageIntersections = AverageIntersections + TotalNumberIntersections;
       if(NumberofNeurons>0)
           AverageIntersections = AverageIntersections/NumberofNeurons;
       else
           AverageIntersections =[];
       end   
       %Fill matrix if there are no neurons with empty cells
          if(NumberofNeurons==0)
            mkdir(foldername);
            NeuronN = int2str(0);
            NN = NN + 1;
            newFileSingleNeuron = regexprep(newFileImage, '++',NeuronN,'ignorecase');
            %Summary 2
            Summary2{1,1} = 'Radius';
            Summary2{1,NN} = newFileSingleNeuron;
            %Summary 3
            Summary3{1,1} = 'Radius';
            Summary3{1,NN} = newFileSingleNeuron;
            %Summary 4
            Summary4{1,1} = 'Radius';
            Summary4{1,NN} = newFileSingleNeuron;  
                p = 1; 
                for (b=1:Ringnumber)
                p = p + 1;
                q = p +1;
                RN = p -1;
                %Summary 2              
                Summary2{2,1} = 0;     
                Summary2{q,1} = RN*10;
                Summary2{p,NN} = [];
                %Summary 3
                Summary3{2,1} = 0;
                Summary3{q,1} = RN*10;
                Summary3{2,NN} = [];
                Summary3{p,NN} =[];
                %Summary 4
                Summary4{2,1} = 0;
                Summary4{q,1} = RN*10;
                Summary4{2,NN} = [];
                Summary4{p,NN} = [];
                end
        end        
        %Create summary  Matrix
        summary{1,1} =  'Image Name';
        summary{1,2} =  'Neurite Mass';
        summary{1,3} =  'Total Neurite Length';
        summary{1,4} =  'Number of Neurites';
        summary{1,5} =  'Average Neurite Length';
        summary{1,6} =  'Number of Branching Points';
        summary{1,7} =  'Number of tips'
        summary{1,8} =  'Total Number Intersections';
        summary{1,9} =  'Number of Neurons'
        % Add values
        summary{mm,1} =  newFileImage1;
        summary{mm,2} =  NeuriteMass;
        summary{mm,3} =  TotalNeuritelength;
        summary{mm,4} =  NumberofNeurites;
        summary{mm,5} =  AverageNeuritelength;
        summary{mm,6} =  NumberBranches;
        summary{mm,7} =  endpointSkel;
        summary{mm,8} =  AverageIntersections;
        summary{mm,9} = NumberofNeurons;
        
        NeuronStats = NeuronStatsRe; 
        wellname = (newfileRing(1:3));
        wellList{n} = wellname;
        csvHandler.CellPosMatrix(wellname) = NucleusM;
        csvHandler.ManualPositions2(wellname) = NeuronStats;
    end
    
end
         