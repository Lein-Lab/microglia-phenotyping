classdef Neurite<handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cutXPosStart;
        cutYPosStart;
        cutXPosEnd;
        cutYPosEnd;
        cutWidth;
        cutHeight;
        neuriteLength;
        image;
        xBranchingPoints;
        yBranchingPoints;
        xCrossingPoints;
        yCrossingPoints;
        
        %Bais Parameter
        skeletonImage;
        skeletonImageOrig;
        xVertices;
        yVertices;
        xEndpoints;
        yEndpoints;
        xShortedEndpoints;
        yShortedEndpoints;
        savedXEndpoints;
        savedYEndpoints;
        xContour;
        yContour;        
        
        ConnectionList;
        NeuriteAreaSize;
        StrongThresholdedPixels;
    end
    
    methods
        
        function CalculateBaiSkeleton(self,no_vertices)
            [bw,I0,x,y,x1,y1,aa,bb]=div_skeleton_new(4,1,self.image,no_vertices);
            self.skeletonImage = parsiSkel(bw);
            self.skeletonImageOrig = self.skeletonImage;
            self.yVertices = x;
            self.xVertices = y;
            self.yEndpoints = x1;
            self.xEndpoints = y1;
            self.savedXEndpoints = y1;
            self.savedYEndpoints = x1;
            self.yContour = aa;
            self.xContour = bb;   
            self.xShortedEndpoints = zeros(numel(self.savedXEndpoints),1);
            self.yShortedEndpoints = zeros(numel(self.savedXEndpoints),1);
            %Trace every endpoint 15 px back
            for i=1:numel(self.xEndpoints)
                [yTraced,xTraced] = self.TraceBackSkeleton(self.yEndpoints(i),self.xEndpoints(i));
                self.xShortedEndpoints(i) = xTraced;
                self.yShortedEndpoints(i) = yTraced;
            end
        end
        
        function PlotSkeleton(self,~)
            figure(1);
            hFig = figure(1);
            [height width] = size(self.image);
            set(hFig, 'Position', [300 300 width height])
            imshow(self.skeletonImage+self.image);
           % hold on;
           % plot(self.yContour, self.xContour, '.r');
           % plot(self.yEndpoints, self.xEndpoints, 'og');
           % plot(self.yVertices, self.xVertices, '.g');
        end
        
        function neighbours=GetNeighbourCount(self,y,x)
            [whitePointsY, whitePointsX] = find(self.skeletonImage);
            [ySizeSkel xSizeSkel] = size(self.skeletonImage);
            %Check every pixel in 8 neighbourhood.
            neighbourcount = 0;
            if(y-1 > 0 && x-1 > 0 && y+1<size(self.skeletonImage,1) && x+1 < size(self.skeletonImage,2))
                neighbours=self.skeletonImage(y-1:y+1,x-1:x+1);
                neighbours=nnz(neighbours)-1;
            else
                neighbours=0;
            end
        end

        
        function FindBranchingCrossingPoints(self)            
            self.xBranchingPoints = zeros(0);
            self.yBranchingPoints = zeros(0);
            self.xCrossingPoints = zeros(0);
            self.yCrossingPoints = zeros(0);


            %If yes, try to differentiate both Neurites by angle.
            %Detection of different Neurite types:
            %Available: EndPoints, get Branching and Crossing Points from Skeleton
            [whitePointsY, whitePointsX] = find(self.skeletonImage);
            [ySizeSkel xSizeSkel] = size(self.skeletonImage);
            for(i=1:numel(whitePointsY))
                %Check every pixel in 8 neighbourhood.
                neighbourcount = 0;
                for(j=-1:1)
                    for(k=-1:1)
                        if(whitePointsY(i)+j > 0 && whitePointsY(i)+j <= size(self.skeletonImage,1) && whitePointsX(i)+k > 0 && whitePointsX(i)+k <= size(self.skeletonImage,2) && self.skeletonImage(whitePointsY(i)+j,whitePointsX(i)+k) == 1)
                            neighbourcount = neighbourcount+1;
                        end
                    end
                end
                %Two pixels: Endpoint
                if(neighbourcount == 2)
                elseif(neighbourcount == 3)
                %Three pixels in neighbourhood: normal
                elseif(neighbourcount==4)
                %Four pixels: Branching point
                    self.xBranchingPoints = [self.xBranchingPoints whitePointsX(i)];
                    self.yBranchingPoints = [self.yBranchingPoints whitePointsY(i)];
                elseif(neighbourcount==5)
                %Five pixels: Crossing point
                    self.xCrossingPoints = [self.xCrossingPoints whitePointsX(i)];
                    self.yCrossingPoints = [self.yCrossingPoints whitePointsY(i)];
                end
            end
        end
        
        function [y,x] = TraceBackSkeleton(self,currentY,currentX)
        %Trace Route
        center = regionprops(self.image,'Centroid');
        nothingFound=0;
        pathLength=0;
        tracedX = [currentX];
        tracedY = [currentY];
        while(~nothingFound && pathLength<10)
                %ToDo: Check that initial direction is correct.
                %Allow only direction to balance point of Neurite/Skeleton
                if(currentX > center.Centroid(1))
                    xStart=-1;
                    xStop=0;
                else
                    xStart = 0;
                    xStop=1;
                end
                if(currentY > center.Centroid(2))
                    yStart=-1;
                    yStop=0;
                else
                    yStart = 0;
                    yStop=1;
                end                
                breakOuter=0;
                nothingFound=1;                
                for(j=xStart:xStop)
                    for(k=yStart:yStop)
                       tracingMatrix = reshape([tracedY tracedX],numel(tracedX),2);
                       if( ~ismember([currentY+k currentX+j],tracingMatrix,'rows') && currentY+k>0 && currentX+j>0 && (currentY+k<=size(self.skeletonImageOrig,1)) && (currentX+j<=size(self.skeletonImageOrig,2))  && (self.skeletonImageOrig(currentY+k,currentX+j) == 1) && ~(j==0 && k==0))
                        %ToDo: Check if branching point or crossing point
                        %hit
                        %If yes: Save that point and go along random path.
                        %If hit end of neurite before path length, go back to previous branching point and go other direction
                        currentX = currentX+j;
                        currentY = currentY+k;
                        tracedX = [tracedX currentX];
                        tracedY = [tracedY currentY];
                        %Exclude predecessor and same element
                        pathLength = pathLength + 1;
                        nothingFound=0;
                        breakOuter=1;                        
                        break;
                       end 
                    end
                    if(breakOuter==1)
                        break;
                    end
                end 
                
                %If nothing found try all different directions
                if(nothingFound==1)
                    for(j=-1:1)
                        for(k=-1:1)
                           tracingMatrix = reshape([tracedY tracedX],numel(tracedX),2);
                           if( ~ismember([currentY+k currentX+j],tracingMatrix,'rows') && currentY+k>0 && currentX+j>0 && (currentY+k<=size(self.skeletonImageOrig,1)) && (currentX+j<=size(self.skeletonImageOrig,2))  && (self.skeletonImageOrig(currentY+k,currentX+j) == 1) && ~(j==0 && k==0))
                            %ToDo: Check if branching point or crossing point
                            %hit
                            %If yes: Save that point and go along random path.
                            %If hit end of neurite before path length, go back to previous branching point and go other direction
                            currentX = currentX+j;
                            currentY = currentY+k;
                            tracedX = [tracedX currentX];
                            tracedY = [tracedY currentY];
                            %Exclude predecessor and same element
                            pathLength = pathLength + 1;
                            nothingFound=0;
                            breakOuter=1;                        
                            break;
                           end 
                        end
                        if(breakOuter==1)
                            break;
                        end
                    end
                end
                
                %If still nothing found try 2 steps in all directions
                if(nothingFound==1)
                    for(j=-2:2)
                        for(k=-2:2)
                           tracingMatrix = reshape([tracedY tracedX],numel(tracedX),2);
                           if( ~ismember([currentY+k currentX+j],tracingMatrix,'rows') && currentY+k>0 && currentX+j>0 && (currentY+k<=size(self.skeletonImageOrig,1)) && (currentX+j<=size(self.skeletonImageOrig,2))  && (self.skeletonImageOrig(currentY+k,currentX+j) == 1) && ~(j==0 && k==0))
                            %ToDo: Check if branching point or crossing point
                            %hit
                            %If yes: Save that point and go along random path.
                            %If hit end of neurite before path length, go back to previous branching point and go other direction
                            currentX = currentX+j;
                            currentY = currentY+k;
                            tracedX = [tracedX currentX];
                            tracedY = [tracedY currentY];
                            %Exclude predecessor and same element
                            pathLength = pathLength + 1;
                            nothingFound=0;
                            breakOuter=1;                        
                            break;
                           end 
                        end
                        if(breakOuter==1)
                            break;
                        end
                    end
                end
            end
            y=currentY;
            x=currentX;
        end
        
        function KillShortBranchesOLD(self)
            %Check shortest distances on path from every EndPoint to next Branching point
            %If shortest distance < thresh -> Kill EndPoint, Branching Point and Path to Branching point
            pathLengthThreshold = 12;
            deletionMarkedEndpointsX = zeros(0);
            deletionMarkedEndpointsY = zeros(0);
            for(i=1:numel(self.xEndpoints))
                currentX = 0;
                currentY = 0;
                tracedX = [self.xEndpoints(i)];
                tracedY = [self.yEndpoints(i)];
                pathLength=0;
                %check all neighbours. If one neighbour is branching point: Done
                %Go step in right direction
                %[neighboursY neighboursX] = zeros(ySkelSize, xSkelSize);
                for(j=-1:1)
                    for(k=-1:1)
                        if(self.skeletonImage(self.yEndpoints(i)+j,self.xEndpoints(i)+k) == 1 && ~(j==0 && k==0))
                            currentX = self.xEndpoints(i)+k;
                            currentY = self.yEndpoints(i)+j;
                            tracedX = [tracedX currentX];
                            tracedY = [tracedY currentY];
                        end
                    end
                end
                pathLength=pathLength+1;
                dx = currentX - self.xEndpoints(i);
                dy = currentY - self.yEndpoints(i);
                %While no branching, crossing or end point reached
                branchingMatrix = reshape([self.yBranchingPoints self.xBranchingPoints],numel(self.yBranchingPoints),2);
                crossingMatrix = reshape([self.yCrossingPoints self.xCrossingPoints],numel(self.yCrossingPoints),2);
                endpointMatrix = reshape([self.yEndpoints self.xEndpoints],numel(self.yEndpoints),2);
                nothingFound=0;
                while(~nothingFound && ~(pathLength>100 || ismember([currentY currentX],branchingMatrix,'rows') || ismember([currentY currentX],crossingMatrix,'rows') || ismember([currentY currentX],endpointMatrix,'rows')))
                    %Try to get next neighbour. Check direction.
                    %Direction can't be the opposite:
                    %dx or dy mustn't be -dx or -dy
                    %Get all neighbours of current Point
                    if(dx==-1)
                        xStart = -1;
                        xStop = 1;
                    elseif(dx==1)
                        xStart=-1;
                        xStop=1;
                    else
                        xStart = -1;
                        xStop=1;
                    end
                    if(dy==-1)
                        yStart = -1;
                        yStop = 1;
                    elseif(dy==1)
                        yStart=-1;
                        yStop=1;
                    else
                        yStart = -1;
                        yStop=1;
                    end
                    breakOuter=0;
                    nothingFound=1;
                    for(j=xStart:xStop)
                        for(k=yStart:yStop)
                           tracingMatrix = reshape([tracedY tracedX],numel(tracedX),2);
                           if( ~ismember([currentY+k currentX+j],tracingMatrix,'rows') && currentY+k>0 && currentX+j>0 && (currentY+k<=size(self.skeletonImage,1)) && (currentX+j<=size(self.skeletonImage,2))  && (self.skeletonImage(currentY+k,currentX+j) == 1) && ~(j==0 && k==0))
                            currentX = currentX+j;
                            currentY = currentY+k;
                            tracedX = [tracedX currentX];
                            tracedY = [tracedY currentY];
                            %Exclude predecessor and same element
                            pathLength = pathLength + 1;
                            breakOuter=1;
                            nothingFound=0;
                            break;
                           end
                        end
                        if(breakOuter==1)
                            break;
                        end
                    end
                end
                if(pathLength < pathLengthThreshold)
                    %Kill path
                    for(m=1:numel(tracedX)-1)
                        self.skeletonImage(tracedY(m),tracedX(m)) = 0;
                    end
                    %Kill branching point
                    [ismem branchpos] = ismember([currentY currentX],branchingMatrix,'rows');
                    if(ismem)
                        branchingMatrix(branchpos,:) = [];
                        self.xBranchingPoints(branchpos) = [];
                        self.yBranchingPoints(branchpos) = [];
                        deletionMarkedEndpointsX = [deletionMarkedEndpointsX i];
                    end
                    %Make crossing point to branching point
                    [ismem crosspos] = ismember([currentY currentX],crossingMatrix,'rows');
                    if(ismem)
                        crossingMatrix(crosspos,:) = [];
                        self.xCrossingPoints(crosspos)=[];
                        self.yCrossingPoints(crosspos)=[];
                        self.yBranchingPoints = [self.yBranchingPoints currentY];
                        self.xBranchingPoints = [self.xBranchingPoints currentX];
                        branchingMatrix = [branchingMatrix;currentY currentX];
                        deletionMarkedEndpointsX = [deletionMarkedEndpointsX i];
                    end
                end
            end
            deletionMarkedEndpointsX = sort(deletionMarkedEndpointsX,'descend');
            for(i=1:numel(deletionMarkedEndpointsX))
                self.xEndpoints(deletionMarkedEndpointsX(i)) = [];
                self.yEndpoints(deletionMarkedEndpointsX(i)) = [];
            end
        end
        
        function KillShortBranches(self)
            %Check shortest distances on path from every EndPoint to next Branching point
            %If shortest distance < thresh -> Kill EndPoint, Branching Point and Path to Branching point
            pathLengthThreshold = 12;
            deletionMarkedEndpointsX = zeros(0);
            deletionMarkedEndpointsY = zeros(0);
            i=1;
            while(i<=numel(self.xEndpoints))
            %for(i=1:numel(self.xEndpoints))
                currentX = 0;
                currentY = 0;
                tracedX = [self.xEndpoints(i)];
                tracedY = [self.yEndpoints(i)];
                pathLength=0;
                %check all neighbours. If one neighbour is branching point: Done
                %Go step in right direction
                %[neighboursY neighboursX] = zeros(ySkelSize, xSkelSize);
                for(j=-1:1)
                    for(k=-1:1)            
                        if(self.skeletonImage(self.yEndpoints(i)+j,self.xEndpoints(i)+k) == 1 && ~(j==0 && k==0))
                            currentX = self.xEndpoints(i)+k;
                            currentY = self.yEndpoints(i)+j;
                            tracedX = [tracedX currentX];
                            tracedY = [tracedY currentY];
                        end
                    end
                end
                pathLength=pathLength+1;
                dx = currentX - self.xEndpoints(i);
                dy = currentY - self.yEndpoints(i);
                %While no branching, crossing or end point reached
                branchingMatrix = reshape([self.yBranchingPoints self.xBranchingPoints],numel(self.yBranchingPoints),2);
                crossingMatrix = reshape([self.yCrossingPoints self.xCrossingPoints],numel(self.yCrossingPoints),2);
                endpointMatrix = reshape([self.yEndpoints self.xEndpoints],numel(self.yEndpoints),2);
                nothingFound=0;
                while(~nothingFound && ~(pathLength>100 || ismember([currentY currentX],branchingMatrix,'rows') || ismember([currentY currentX],crossingMatrix,'rows') || ismember([currentY currentX],endpointMatrix,'rows')))
                    %Try to get next neighbour. Check direction.
                    %Direction can't be the opposite:
                    %dx or dy mustn't be -dx or -dy
                    %Get all neighbours of current Point
                    if(dx==-1)
                        xStart = -1;
                        xStop = 1;
                    elseif(dx==1)
                        xStart=-1;
                        xStop=1;
                    else
                        xStart = -1;
                        xStop=1;
                    end
                    if(dy==-1)
                        yStart = -1;
                        yStop = 1;
                    elseif(dy==1)
                        yStart=-1;
                        yStop=1;
                    else
                        yStart = -1;
                        yStop=1;
                    end
                    breakOuter=0;
                    nothingFound=1;
                    for(j=xStart:xStop)
                        for(k=yStart:yStop)
                           tracingMatrix = reshape([tracedY tracedX],numel(tracedX),2);
                           if( ~ismember([currentY+k currentX+j],tracingMatrix,'rows') && currentY+k>0 && currentX+j>0 && (currentY+k<=size(self.skeletonImage,1)) && (currentX+j<=size(self.skeletonImage,2))  && (self.skeletonImage(currentY+k,currentX+j) == 1) && ~(j==0 && k==0))
                            currentX = currentX+j;
                            currentY = currentY+k;
                            tracedX = [tracedX currentX];
                            tracedY = [tracedY currentY];
                            %Exclude predecessor and same element
                            pathLength = pathLength + 1;
                            breakOuter=1;
                            nothingFound=0;
                            break;
                           end 
                        end
                        if(breakOuter==1)
                            break;
                        end
                    end
                end
                if(pathLength < pathLengthThreshold)
                    %Kill path
                    for(m=1:numel(tracedX)-1)
                        self.skeletonImage(tracedY(m),tracedX(m)) = 0;
                    end     
                    %Check if last one has only one neighbour: Then convert
                    %it to endpoint:
                    %if(self.GetNeighbourCount(tracedY(numel(tracedY)),tracedX(numel(tracedX))) <= 1)
                    %        self.xEndpoints = [self.xEndpoints;tracedX(numel(tracedX))];
                    %        self.yEndpoints = [self.yEndpoints;tracedY(numel(tracedY))];
                    %end
                    
                    %Kill branching point
                    [ismem branchpos] = ismember([currentY currentX],branchingMatrix,'rows');
                    if(ismem)
                        branchingMatrix(branchpos,:) = [];
                        if(self.GetNeighbourCount(self.yBranchingPoints(branchpos),self.xBranchingPoints(branchpos)) <= 1)
                            self.xEndpoints = [self.xEndpoints;self.xBranchingPoints(branchpos)];
                            self.yEndpoints = [self.yEndpoints;self.yBranchingPoints(branchpos)];
                            [yTraced,xTraced] = self.TraceBackSkeleton(self.yBranchingPoints(branchpos),self.xBranchingPoints(branchpos));
                            self.xShortedEndpoints = [self.xShortedEndpoints;xTraced];
                            self.yShortedEndpoints = [self.yShortedEndpoints;yTraced];
                        end
                        self.xBranchingPoints(branchpos) = [];
                        self.yBranchingPoints(branchpos) = [];                        
                        deletionMarkedEndpointsX = [deletionMarkedEndpointsX i];                       
                    end
                    %Make crossing point to branching point
                    [ismem crosspos] = ismember([currentY currentX],crossingMatrix,'rows');
                    if(ismem)
                        crossingMatrix(crosspos,:) = [];
                        self.xCrossingPoints(crosspos)=[];
                        self.yCrossingPoints(crosspos)=[];
                        self.yBranchingPoints = [self.yBranchingPoints currentY];
                        self.xBranchingPoints = [self.xBranchingPoints currentX];
                        branchingMatrix = [branchingMatrix;currentY currentX];
                        deletionMarkedEndpointsX = [deletionMarkedEndpointsX i];
                    end
                end
                i=i+1;
            end
            deletionMarkedEndpointsX = sort(deletionMarkedEndpointsX,'descend');
            for(i=1:numel(deletionMarkedEndpointsX))
                self.xEndpoints(deletionMarkedEndpointsX(i)) = [];
                self.yEndpoints(deletionMarkedEndpointsX(i)) = [];
                self.xShortedEndpoints(deletionMarkedEndpointsX(i)) = [];
                self.yShortedEndpoints(deletionMarkedEndpointsX(i)) = [];
            end
            %Final ToDo:
            %Check number of neighbours of every branching point and their neighbours recursively.
            %If only one neighbour: Convert them to EndPoints
            z=1;
            while(z<=numel(self.xBranchingPoints))
                if(self.GetNeighbourCount(self.yBranchingPoints(z),self.xBranchingPoints(z)) <= 1)
                            self.xEndpoints = [self.xEndpoints;self.xBranchingPoints(z)];
                            self.yEndpoints = [self.yEndpoints;self.yBranchingPoints(z)];
                            %Add shortened Endpoint
                            [yTraced,xTraced] = self.TraceBackSkeleton(self.yBranchingPoints(z),self.xBranchingPoints(z));
                            self.xShortedEndpoints = [self.xShortedEndpoints;xTraced];
                            self.yShortedEndpoints = [self.yShortedEndpoints;yTraced];                            
                            self.xBranchingPoints(z) = [];
                            self.yBranchingPoints(z) = [];
                else
                    finished=0;
                    for(j=-1:1)
                        for(k=-1:1)
                            if(self.GetNeighbourCount(self.yBranchingPoints(z)+j,self.xBranchingPoints(z)+k) <= 1)
                               self.xEndpoints = [self.xEndpoints;self.xBranchingPoints(z)+j];
                               self.yEndpoints = [self.yEndpoints;self.yBranchingPoints(z)+k];
                               
                               [yTraced,xTraced] = self.TraceBackSkeleton(self.yBranchingPoints(z)+k,self.xBranchingPoints(z)+j);
                               self.xShortedEndpoints = [self.xShortedEndpoints;xTraced];
                               self.yShortedEndpoints = [self.yShortedEndpoints;yTraced];    
                               
                               self.xBranchingPoints(z) = [];
                               self.yBranchingPoints(z) = [];
                               finished=1;
                               break;
                            end
                        end
                        if(finished==1)
                            break;
                        end
                    end
                end
                z=z+1;
            end
        end
        
        function [bw]=GetEndPoints(self,ro,T1,I0,no_vertice)
            %Iterate over image and try to find white points with just one
            %neighbour
        end
        
        function SplitUp(self)
        end
        
        function CalculateNeuriteLength(self)
            %Neurite length is sum of FIXED! Skeleton Image
            self.neuriteLength = numel(find(self.skeletonImageOrig));
        end
        
        function SplitNeurites(self)
            %If branching points available:
            %Look from branching point for angles to EndPoints
            %When two EndPoints have an angle around 180�, they are
            %probably belong together.
            %Endpoints with low angle difference belong to the regarding other Endpoint 
        end
        
    end
    
end

