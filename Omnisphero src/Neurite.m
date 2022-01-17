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
        
        nucleusImage;
        numberNuclei=0;
        rgbImage;
    end
    
    methods

        function [yBranchingPoints xBranchingPoints] = KillCloseBranches(self,yBranchingPoints,xBranchingPoints)
%Again: Kill branches which are close together.
xMid = zeros(0);
yMid = zeros(0);
killIndicesList = zeros(0);
newPointsListY = zeros(0);
newPointsListX = zeros(0);

    %Solange es zwei Branching Points mit pdist < 30 gibt, fï¿½hre folgenden
    %Algorithmus durch:
    %Fasse beide Branches zu einem zusammen.
    for(u=1:3)
        killIndicesList=zeros(0);
        newPointsListY=zeros(0);
        newPointsListX=zeros(0);
        for(a=1:numel(yBranchingPoints))
            for(b=a:numel(yBranchingPoints))
                if(a~=b && pdist([yBranchingPoints(a) xBranchingPoints(a);yBranchingPoints(b) xBranchingPoints(b)]) <= 6 && ~ismember(a,killIndicesList) && ~ismember(b,killIndicesList))
                    %Fasse beide Branches zu einem zusammen!
                    killIndicesList=[killIndicesList a];
                    killIndicesList=[killIndicesList b];
                    xMid = mean([xBranchingPoints(a) xBranchingPoints(b)]);
                    yMid = mean([yBranchingPoints(a) yBranchingPoints(b)]);
                    newPointsListY = [newPointsListY;yMid];
                    newPointsListX = [newPointsListX;xMid];
                end
            end
        end

        yBranchingPoints(killIndicesList) = [];
        xBranchingPoints(killIndicesList) = [];
        yBranchingPoints = [yBranchingPoints;newPointsListY];
        xBranchingPoints = [xBranchingPoints;newPointsListX];
    end
end
        
        function FindBranchingCrossingPoints2(self)
            %  Idee:
            % 1. Lösche Branching Points
            %for(z=1:2)
            branchImage = bwmorph(self.skeletonImage, 'branchpoints');
            [sizeY sizeX] = size(branchImage);
            subImage = self.skeletonImage - branchImage;
            subImage(subImage<0)=0;
            subImage=logical(subImage);
            
            % 2. Lösche alle connected Areas mit <= 11 Pixeln
            labeledNeurite = bwlabel(subImage);
            for(i=1:max(labeledNeurite(:)))
                if(numel(find(labeledNeurite==i)) <= 4)
                    labeledNeurite(labeledNeurite==i)=0;
                end
            end
            % 3. Verbinde verbliebene Regionen wieder miteinander, aber so, dass neuer Verbindungspunkt nur 1-2 Verbindungen in seiner 8-Nachbarschaft hat.
            
            endpointImage = bwmorph(labeledNeurite, 'endpoints');
            %For each endpoint:
            %Check shortest distance to next endpoint with another label.
            %For each label: Make one connection to another endpoint
            pathSave = logical(zeros(sizeY,sizeX));
            for(i=1:max(labeledNeurite(:)))
                %Iterate over every endpoint of current label
                %Get pixels which are 1 in endpointimage and which are
                %==i in labeledNeurite
                [endpointsCurrentLabelY endpointsCurrentLabelX] = find(endpointImage & labeledNeurite==i);
                subImageIndices = find(labeledNeurite~=i & labeledNeurite>0 & endpointImage);
                otherNeuriteImageIndices = find(labeledNeurite~=i & labeledNeurite>0);
                %Then get all other endpoints in endpointimage, where
                %i~= labeledNeurite
                endpointsDifferentLabel = find(endpointImage & labeledNeurite ~=i);
                
                subImage = logical(zeros(sizeY,sizeX));
                subImage(subImageIndices)=1;
                
                otherNeuriteImage = logical(zeros(sizeY,sizeX));
                otherNeuriteImage(otherNeuriteImageIndices)=1;
                
                [D IDX] = bwdist(subImage);
                maxDistToNextPic=9999;
                maxEndpointIndex=0;                
                for(j=1:numel(endpointsCurrentLabelY))
                    %Get shortest path between points in
                    %endpointsCurrentLabel and next pixel in
                    %subImage, which doesn't create new branching points
                    currentDist = D(endpointsCurrentLabelY(j),endpointsCurrentLabelX(j));
                    if(currentDist<maxDistToNextPic)
                        maxDistToNextPic=currentDist;
                        maxEndpointIndex=j;
                    end                    
                end                
            
                %Draw line between endpointCurrentLabelXY(maxEndpointIndex) and
                %D(endpointCurrentLabelXY(maxEndpointIndex))
                if(maxEndpointIndex>0)
                    y1=endpointsCurrentLabelY(maxEndpointIndex);
                    x1=endpointsCurrentLabelX(maxEndpointIndex);
                    [y2 x2] = ind2sub([sizeY sizeX],IDX(y1,x1));
                    %Draw line between xy1 and xy2
                    distY = y2-y1;
                    distX = x2-x1;
                    stepX = distX/10;
                    stepY = distY/10;
                    currentX=x1;
                    currentY=y1;
                    neighbourReached=0;
                    
                    %maze are indices where new branching points would
                    %occur
                    maze=zeros(sizeY,sizeX);
                    maze(labeledNeurite>0)=1;
                    %Remove regarding endpoints from maze
                    maze(y1,x1)=0;
                    maze(y2,x2)=0;
                    %Add neighboured points to maze
                    [D2 IDX2] = bwdist(maze);
                    maze(D2==1)=1;
                    maze(y1,x1)=0;
                    maze(y2,x2)=0;
                    %Connect with BFS
                    %path = solve_maze(maze, [x1 y1], [x2 y2]);
                    path = solve_maze2(self,maze,y1,x1,y2,x2,sizeY,sizeX,15);
                    pathSave(find(path))=1;
%                     while(neighbourReached==0 && (currentY ~= y2 || currentX ~= x2))
%                         oldX = currentX;
%                         oldY = currentY;
%                         currentX=currentX+stepX;
%                         currentY=currentY+stepY;
%                         %Check if only two neighbours in 8-neighbourhood of
%                         %labeledNeuriteImage.
%                         %If there are three neighbours go back one step and
%                         %check neighboured directions
%                         %e.g.
%                         self.skeletonImage(round(currentY),round(currentX)) = 1;
%                         %Check if one neighbour is already touched and
%                         %ensure that no new branching point is created
%                         if(nnz(otherNeuriteImage(currentY-1:currentY+1,currentX-1:currentX+1)) > 0)
%                             neighbourReached=1;
%                         end
%                     end
                   %labeledNeurite(labeledNeurite>1)=1;                   
                end
                self.skeletonImage=logical(labeledNeurite);
                self.skeletonImage(find(pathSave))=1;
                self.skeletonImage = parsiSkel(self.skeletonImage);
           end
            
            
            % Dazu: Finde Endpunkte in existierenden Linien.
            % 4. Sobald jede Region mit einer anderen verbunden wurde, es aber weiterhin mehrere Regionen gibt, werden auch diese auf dem kürzesten Weg miteinander verbunden.
            % 5. Suche erneut nach Branching Points
            
            

            %3 Endpoints: 1 branch
            %4 or more endpoints: Many branches
            
            % 6. Repeat
            %end           
            
            %Connect remaining branches independently of endpoints with
            %maze strategy!
            
            labeledNeurite = bwlabel(self.skeletonImage);
            %maze are indices where new branching points would
            %occur
            while(max(labeledNeurite(:)) > 1)
                labelOne = zeros(sizeY,sizeX);
                labelTwo = zeros(sizeY,sizeX);
                labelOne(labeledNeurite==1)=1;
                labelTwo(labeledNeurite==2)=1;
                [D IDX] = bwdist(labelOne);
                %Find min distance between One and Two
                D(labelOne==1) = 999;
                D(labelTwo~=1)=999;
                [minVal minInd]=min(D(:));
                target=IDX(minInd);
                [y1 x1] = ind2sub([sizeY sizeX], minInd);
                [y2 x2] = ind2sub([sizeY sizeX], target);
                maze=zeros(sizeY,sizeX);
                             
                %Connect with BFS
                %path = solve_maze(maze, [x1 y1], [x2 y2]);
                path = solve_maze2(self,maze,y1,x1,y2,x2,sizeY,sizeX,0);
                self.skeletonImage(find(path))=1;
                labeledNeurite = bwlabel(self.skeletonImage);
            end
            
            
            branchImage = logical(bwmorph(self.skeletonImage, 'branchpoints'));
            endpointImage = logical(bwmorph(self.skeletonImage, 'endpoints'));
            %2 Endpoints: No branches, maybe just circles
            self.xBranchingPoints = zeros(0);
            self.yBranchingPoints = zeros(0);
            [yBranchingPoints xBranchingPoints] = find(branchImage==1);
            self.yBranchingPoints = yBranchingPoints;
            self.xBranchingPoints = xBranchingPoints;
            
            [yEndPoints xEndPoints] = find(endpointImage==1);
            self.yEndpoints = yEndPoints;
            self.xEndpoints = xEndPoints;
            
           if(nnz(endpointImage) <= 2)
                self.xBranchingPoints = [];
                self.yBranchingPoints = [];
            elseif(nnz(endpointImage) == 3)
                %Calculate MidPoint of branchImage
                xMid = zeros(0);
                yMid = zeros(0);
                for(l=1:nnz(self.xBranchingPoints))
                    xMid = [xMid self.xBranchingPoints];
                    yMid = [yMid self.yBranchingPoints];
                end
                self.xBranchingPoints = mean(xMid);
                self.yBranchingPoints = mean(yMid);                
           end

        %If two branches are clother together than 6 pixels, merge them to one
        %branching point
        
        xMid = zeros(0);
        yMid = zeros(0);
        killIndicesList = zeros(0);
        newPointsListY = zeros(0);
        newPointsListX = zeros(0);
        
        %Solange es zwei Branching Points mit pdist < 30 gibt, führe folgenden
        %Algorithmus durch:
        %Fasse beide Branches zu einem zusammen.
        for(u=1:3)
            killIndicesList=zeros(0);
            newPointsListY=zeros(0);
            newPointsListX=zeros(0);
            for(a=1:numel(self.yBranchingPoints))
                for(b=a:numel(self.yBranchingPoints))
                    if(a~=b && pdist([self.yBranchingPoints(a) self.xBranchingPoints(a);self.yBranchingPoints(b) self.xBranchingPoints(b)]) <= 6 && ~ismember(a,killIndicesList) && ~ismember(b,killIndicesList))
                        %Fasse beide Branches zu einem zusammen!
                        killIndicesList=[killIndicesList a];
                        killIndicesList=[killIndicesList b];
                        xMid = mean([self.xBranchingPoints(a) self.xBranchingPoints(b)]);
                        yMid = mean([self.yBranchingPoints(a) self.yBranchingPoints(b)]);
                        newPointsListY = [newPointsListY;yMid];
                        newPointsListX = [newPointsListX;xMid];
                    end
                end
            end

            self.yBranchingPoints(killIndicesList) = [];
            self.xBranchingPoints(killIndicesList) = [];
            self.yBranchingPoints = [self.yBranchingPoints;newPointsListY];
            self.xBranchingPoints = [self.xBranchingPoints;newPointsListX];
        end
    end
        
        
        function CalculateNeuriteLength(self)
            %Neurite length is sum of FIXED! Skeleton Image
            self.neuriteLength = numel(find(self.skeletonImageOrig));
            %ToDo:
            %1. Cut off nuclei
            %2. Evaluate how much Neuron-Nuclei are within the neurite area
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
        
        function [y,x] = TraceBackSkeleton(self,currentY,currentX)
        %Trace Route
        center = regionprops(double(self.image),'Centroid');
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
        
        function CalculateBaiSkeleton(self,no_vertices)
            [bw,I0,x,y,x1,y1,aa,bb]=div_skeleton_new(4,1,self.image,no_vertices);
            self.skeletonImage = parsiSkel(bw);
            self.skeletonImageOrig = self.skeletonImage;
            self.nucleusImage = zeros(size(self.skeletonImage));
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
        
        
        
        
        function path = solve_maze(self, maze, start, finish)
%           %% Init data
%           img = imread(img_file);
%           img = rgb2gray(img);
%           maze = img > 0;
%           start = [985 398];
%           finish = [26 399];

          % Init BFS
          n = numel(maze);
          Q = zeros(n, 2);
          M = zeros([size(maze) 2]);
          front = 0;
          back = 1;

          function push(p, d)
            q = p + d;
            if maze(q(1), q(2)) && M(q(1), q(2), 1) == 0
              front = front + 1;
              Q(front, :) = q;
              M(q(1), q(2), :) = reshape(p, [1 1 2]);
            end
          end

          push(start, [0 0]);

          d = [0 1; 0 -1; 1 0; -1 0];

          % Run BFS
          while back <= front
            p = Q(back, :);
            back = back + 1;
            for i = 1:4
              push(p, d(i, :));
            end
          end

          % Extracting path
          path = finish;
          while true
            q = path(end, :);
            p = reshape(M(q(1), q(2), :), 1, 2);
            path(end + 1, :) = p;
            if isequal(p, start) 
              break;
            end
          end
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
        
    
        
        function [bw]=GetEndPoints(self,ro,T1,I0,no_vertice)
            %Iterate over image and try to find white points with just one
            %neighbour
        end
        
        function SplitUp(self)
        end       
        
        
        function SplitNeurites(self)
            %If branching points available:
            %Look from branching point for angles to EndPoints
            %When two EndPoints have an angle around 180ï¿½, they are
            %probably belong together.
            %Endpoints with low angle difference belong to the regarding other Endpoint 
        end
        
    
%The following parameters should be set on Neurite:
        %NucleusImage
        %skeleton Image as Raw Skeletonization from BAI
        function CompleteNeuritePostprocessing(self)
        [sizeYskel sizeXskel] = size(self.nucleusImage);
        %Same Postprocessing as in Skeletonization Algorithm
        oldSkel = self.skeletonImage;            
      
        bothIndices = find(logical(self.skeletonImage) & logical(self.nucleusImage));

        %Doppelt gemoppelt hï¿½lt besser :-D
        %if(ignoreNeurite==0)
        self.FindBranchingCrossingPoints2();           
        self.FindBranchingCrossingPoints2(); 
        %Finde ï¿½berdeckungslï¿½nge zwischen Nucleus und Neuritenbild

        %minNeuriteLength = optionHandler.SkeletonMinNeuriteLength - numel(bothIndices);
        %if(minNeuriteLength < 11)
            minNeuriteLength=11;
        %end
        
        subtractImage = self.skeletonImage; %- logical(self.nucleusImage);
        subtractImage(subtractImage<0)=0;
        subtractImage=logical(subtractImage);
        labelSubtract = bwlabel(subtractImage);

        maxIndex=0;
        maxValue=0;
        for(z=1:max(labelSubtract(:)))
            %Cut off all neurite parts with length < Min Neurite Length
            %If all are below Min Neurite Length: Take the longest one!
            if(numel(find(labelSubtract==z)) > maxValue)
                maxValue = numel(find(labelSubtract==z));
                maxIndex=z;
            end
            if(numel(find(labelSubtract==z)) <= minNeuriteLength)
                subtractImage(labelSubtract==z) = 0;
            end
        end
        self.skeletonImage(labelSubtract==maxIndex)=1;
        subtractImage(labelSubtract==maxIndex) = 1;
        %Ensure that at least one subtractImageLabel remains and 
        %ensure that bothIndices are counted before finding new
        %BranchingCrossingPoints.


        %Do this only if there is at least no Nucleus Point, touching a
        %Neurite Point

        [DSubtract,IDX] = bwdist(subtractImage);
        [DNuc,IDXNuc] = bwdist(self.nucleusImage);
        %neuriteCounter=neuriteCounter+1;

        %Get minimum of DSubtract on Indices where
        %newNeuriteBig.nucleusImage is 1
        
        minDist = min(DSubtract(self.nucleusImage==1));
        if(minDist >= 2)
            %1. Get mid of Nucleus and find shortest path to next Neurite
            [Ay Ax]= find(self.nucleusImage>=1);
            if(numel(Ay) > 0)
                midNucleusY = round(mean(Ay));
                midNucleusX = round(mean(Ax));

                rimPointNeurite=IDX(midNucleusY,midNucleusX);
                %Analog: find rimPointNucleus
                rimPointNucleus=IDXNuc(rimPointNeurite);
                [yRimNeurite xRimNeurite] = ind2sub([sizeYskel sizeXskel],rimPointNeurite);            
                [yRimNucleus xRimNucleus] = ind2sub([sizeYskel sizeXskel],rimPointNucleus);
            else
                yRimNeurite=0;
                xRimNeurite=0;
                yRimNucleus=0;
                xRimNucleus=0;
            end
            
            %If rimPointNeurite and rimPointNucleus are identical,
            %everything is ok, otherwise add line between rimPointNeurite
            %and rimPointNucleus to skeleton Image. Max length: 20
            %pixels
            stepSizeY=yRimNeurite-yRimNucleus;
            stepSizeX=xRimNeurite-xRimNucleus;
            
            while(stepSizeX > 1 || stepSizeY > 1 || stepSizeX < -1 || stepSizeY < -1)
                stepSizeX=stepSizeX/2;
                stepSizeY=stepSizeY/2;
            end
            yIndex = yRimNucleus;
            xIndex = xRimNucleus;
            if(pdist([yIndex xIndex;yRimNeurite xRimNeurite]) <= 20)
                while(round(yIndex) ~= yRimNeurite || round(xIndex) ~= xRimNeurite)
                    yIndex=yIndex+stepSizeY;
                    xIndex=xIndex+stepSizeX;
                    subtractImage(round(yIndex),round(xIndex))=1;
                end
            else
                %Ignore whole neurite
                ignoreNeurite=1;
            end
        end
        self.skeletonImage = subtractImage;
        %newNeuriteBig.FindBranchingCrossingPoints2();           
        %end


        %Create figure with overlayed Nucleus, Skeleton and original
        %Neurite in /SkeletonExport/selectedWell/ID.png
        
        currentObject = zeros(15,20);
        if(numel(find(self.skeletonImage==0))==0)
             self.skeletonImage=zeros(sizeYskel,sizeXskel);
             ignoreNeurite=1;
        end
        %if(ignoreNeurite==0)
            self.skeletonImage = parsiSkel(self.skeletonImage);
            branchImage = logical(bwmorph(self.skeletonImage, 'branchpoints'));
            endpointImage = logical(bwmorph(self.skeletonImage, 'endpoints'));
            %2 Endpoints: No branches, maybe just circles
            self.xBranchingPoints = zeros(0);
            self.yBranchingPoints = zeros(0);
            [yBranchingPoints xBranchingPoints] = find(branchImage==1);
            [yBranchingPoints xBranchingPoints] = self.KillCloseBranches(yBranchingPoints,xBranchingPoints);
            self.yBranchingPoints = yBranchingPoints;
            self.xBranchingPoints = xBranchingPoints;
            
            %Reset skeleton image if it is empty.
            
               
            
            self.skeletonImage = parsiSkel(self.skeletonImage);


            %self.nucleusImage = AreaResultSmall;
            %figure(3);
            %imshow(imfuse(newNeuriteBig.skeletonImage,newNeuriteBig.nucleusImage));
            exportRGB = zeros(sizeYskel,sizeXskel,3);
            exportRGB(:,:,3) = 1-self.image;
            exportRGB(:,:,1) = self.nucleusImage;
            exportRGB(:,:,2) = self.skeletonImage;
            %figure(1);
            %imshow(exportRGB);
            
            %positionID has to be mid of Neurite in big picture
            self.rgbImage = exportRGB;
            %positionID = sub2ind([sizeY sizeX],cutYPosStartRef, cutXPosStartRef);
            %Check if Neurite is already available in Dictionary            
            %connections=SplitUpNeuronsAndNuclei(SubNucM,newNeuriteBig.nucleusImage,newNeuriteBig.skeletonImage,newNeuriteBig,sizeY,sizeX);
            
        
        end
    
    end
end

