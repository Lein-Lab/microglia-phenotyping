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


%Check shortest distances on path from every EndPoint to next Branching point
%If shortest distance < thresh -> Kill EndPoint, Branching Point and Path to Branching point
function bw=EliminateShortBranches(bw,pathLengthThreshold)
    %bw: image
    %pathLengthThreshold: maximal length from Endpoint to Branching point
    origBw = bw;
    deletionMarkedEndpointsX = zeros(0);
    deletionMarkedEndpointsY = zeros(0);
    endpointImage = bwmorph(bw, 'endpoints');
    branchImage = bwmorph(bw, 'branchpoints');
    [yEndpoints xEndpoints] = find(endpointImage==1);
    [yBranchingPoints xBranchingPoints] = find(endpointImage==1);
    i=1;
    while(i<=numel(xEndpoints))
        currentX = 0;
        currentY = 0;
        tracedX = [xEndpoints(i)];
        tracedY = [yEndpoints(i)];
        pathLength=0;
        %check all neighbours. If one neighbour is branching point: Done
        %Go step in right direction
        %[neighboursY neighboursX] = zeros(ySkelSize, xSkelSize);
        for(j=-1:1)
            for(k=-1:1)            
                if(bw(yEndpoints(i)+j,xEndpoints(i)+k) == 1 && ~(j==0 && k==0))
                    currentX = xEndpoints(i)+k;
                    currentY = yEndpoints(i)+j;
                    tracedX = [tracedX currentX];
                    tracedY = [tracedY currentY];
                end
            end
        end
        pathLength=pathLength+1;
        dx = currentX - xEndpoints(i);
        dy = currentY - yEndpoints(i);
        %While no branching, crossing or end point reached
        branchingMatrix = reshape([yBranchingPoints xBranchingPoints],numel(yBranchingPoints),2);
        endpointMatrix = reshape([yEndpoints xEndpoints],numel(yEndpoints),2);
        nothingFound=0;
        while(~nothingFound && ~(pathLength>100 || ismember([currentY currentX],branchingMatrix,'rows') || ismember([currentY currentX],endpointMatrix,'rows')))
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
                   if( ~ismember([currentY+k currentX+j],tracingMatrix,'rows') && currentY+k>0 && currentX+j>0 && (currentY+k<=size(bw,1)) && (currentX+j<=size(bw,2))  && (bw(currentY+k,currentX+j) == 1) && ~(j==0 && k==0))
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
                bw(tracedY(m),tracedX(m)) = 0;
            end
            %Kill branching point
            [ismem branchpos] = ismember([currentY currentX],branchingMatrix,'rows');
            if(ismem)
                branchingMatrix(branchpos,:) = [];
                xBranchingPoints(branchpos) = [];
                yBranchingPoints(branchpos) = [];
                deletionMarkedEndpointsX = [deletionMarkedEndpointsX i];
            end
            %Make crossing point to branching point            
        end
        i=i+1;
    end
    deletionMarkedEndpointsX = sort(deletionMarkedEndpointsX,'descend');
    for(i=1:numel(deletionMarkedEndpointsX))
        xEndpoints(deletionMarkedEndpointsX(i)) = [];
        yEndpoints(deletionMarkedEndpointsX(i)) = [];
    end
    %Final ToDo:
    %Check number of neighbours of every branching point and their neighbours recursively.
    %If only one neighbour: Convert them to EndPoints
    z=1;
    while(z<=numel(xBranchingPoints))
        if(GetNeighbourCount(bw,yBranchingPoints(z),xBranchingPoints(z)) <= 1)
                    xEndpoints = [xEndpoints;xBranchingPoints(z)];
                    yEndpoints = [yEndpoints;yBranchingPoints(z)];
                    %Add shortened Endpoint
                    [yTraced,xTraced] = TraceBackSkeleton(origBw,yBranchingPoints(z),xBranchingPoints(z));                             
                    xBranchingPoints(z) = [];
                    yBranchingPoints(z) = [];
        else
            finished=0;
            for(j=-1:1)
                for(k=-1:1)
                    if(GetNeighbourCount(yBranchingPoints(z)+j,xBranchingPoints(z)+k) <= 1)
                       xEndpoints = [xEndpoints;xBranchingPoints(z)+j];
                       yEndpoints = [yEndpoints;yBranchingPoints(z)+k];

                       [yTraced,xTraced] = TraceBackSkeleton(origBw,yBranchingPoints(z)+k,xBranchingPoints(z)+j);                      
                       xBranchingPoints(z) = [];
                       yBranchingPoints(z) = [];
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

function neighbours=GetNeighbourCount(bw,y,x)
    [whitePointsY, whitePointsX] = find(bw);
    [ySizeSkel xSizeSkel] = size(bw);
    %Check every pixel in 8 neighbourhood.
    neighbourcount = 0;
    if(y-1 > 0 && x-1 > 0 && y+1<size(bw,1) && x+1 < size(bw,2))
        neighbours=bw(y-1:y+1,x-1:x+1);
        neighbours=nnz(neighbours)-1;
    else
        neighbours=0;
    end
end



function [y,x] = TraceBackSkeleton(skeletonImageOrig,currentY,currentX)
%Trace Route
center = regionprops(skeletonImageOrig,'Centroid');
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
               if( ~ismember([currentY+k currentX+j],tracingMatrix,'rows') && currentY+k>0 && currentX+j>0 && (currentY+k<=size(skeletonImageOrig,1)) && (currentX+j<=size(skeletonImageOrig,2))  && (skeletonImageOrig(currentY+k,currentX+j) == 1) && ~(j==0 && k==0))
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
                   if( ~ismember([currentY+k currentX+j],tracingMatrix,'rows') && currentY+k>0 && currentX+j>0 && (currentY+k<=size(skeletonImageOrig,1)) && (currentX+j<=size(skeletonImageOrig,2))  && (skeletonImageOrig(currentY+k,currentX+j) == 1) && ~(j==0 && k==0))
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
                   if( ~ismember([currentY+k currentX+j],tracingMatrix,'rows') && currentY+k>0 && currentX+j>0 && (currentY+k<=size(skeletonImageOrig,1)) && (currentX+j<=size(skeletonImageOrig,2))  && (skeletonImageOrig(currentY+k,currentX+j) == 1) && ~(j==0 && k==0))
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