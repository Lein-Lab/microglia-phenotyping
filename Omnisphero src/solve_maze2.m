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

function path = solve_maze2(self,maze,startY,startX,finishY,finishX,sizeY,sizeX,pathMax)
            Q = java.util.LinkedList();
            Q.add([startY startX]);
            maze(startY,startX)=1;
            parent=zeros(sizeY,sizeX,3);
            parent(1:sizeY, 1:sizeX, 1) = maze;
            found=0;
            while(Q.size > 0)
                node = Q.remove();
                if(node(1) == finishY && node(2) == finishX)
                    found=1;
                    break;
                end
                %Check for each neighbour if already visited.
                if(node(1)-1 > 0 && parent(node(1)-1,node(2),1) == 0)
                    Q.add([node(1)-1 node(2)]);
                    parent(node(1)-1, node(2), 1) = 1;
                    parent(node(1)-1, node(2), 2) = node(1);
                    parent(node(1)-1, node(2), 3) = node(2);
                end
                if(node(1)-1 > 0 && node(2)-1 > 0 && parent(node(1)-1,node(2)-1,1) == 0)
                    Q.add([node(1)-1 node(2)-1]);
                    parent(node(1)-1, node(2)-1, 1) = 1;
                    parent(node(1)-1, node(2)-1, 2) = node(1);
                    parent(node(1)-1, node(2)-1, 3) = node(2);
                end
                if(node(2)-1 > 0 && parent(node(1),node(2)-1,1) == 0)
                    Q.add([node(1) node(2)-1]);
                    parent(node(1), node(2)-1, 1) = 1;
                    parent(node(1), node(2)-1, 2) = node(1);
                    parent(node(1), node(2)-1, 3) = node(2);
                end
                if(node(1)+1 <= sizeY && node(2)-1 > 0 && parent(node(1)+1,node(2)-1,1) == 0)
                    Q.add([node(1)+1 node(2)-1]);
                    parent(node(1)+1, node(2)-1, 1) = 1;
                    parent(node(1)+1, node(2)-1, 2) = node(1);
                    parent(node(1)+1, node(2)-1, 3) = node(2);
                end
                if(node(1)-1 > 0 && node(2)+1 <= sizeX && parent(node(1)-1,node(2)+1,1) == 0)
                    Q.add([node(1)-1 node(2)+1]);
                    parent(node(1)-1, node(2)+1, 1) = 1;
                    parent(node(1)-1, node(2)+1, 2) = node(1);
                    parent(node(1)-1, node(2)+1, 3) = node(2);
                end
                if(node(1)+1 <= sizeY && parent(node(1)+1,node(2),1) == 0)
                    Q.add([node(1)+1 node(2)]);
                    parent(node(1)+1, node(2), 1) = 1;
                    parent(node(1)+1, node(2), 2) = node(1);
                    parent(node(1)+1, node(2), 3) = node(2);
                end
                if(node(2)+1 <= sizeX && parent(node(1),node(2)+1,1) == 0)
                    Q.add([node(1) node(2)+1]);
                    parent(node(1), node(2)+1, 1) = 1;
                    parent(node(1), node(2)+1, 2) = node(1);
                    parent(node(1), node(2)+1, 3) = node(2);
                end
                if(node(1)+1 <= sizeY && node(2)+1 <= sizeX && parent(node(1)+1,node(2)+1,1) == 0)
                    Q.add([node(1)+1, node(2)+1]);
                    parent(node(1)+1, node(2)+1, 1) = 1;
                    parent(node(1)+1, node(2)+1, 2) = node(1);
                    parent(node(1)+1, node(2)+1, 3) = node(2);
                end
            end       
            path=logical(zeros(sizeY,sizeX));
            currentX = node(2);
            currentY = node(1);
            while(currentX ~= startX | currentY ~= startY)
                if(found==1)
                    path(currentY,currentX)=1;
                end
                currentYNew = parent(currentY, currentX, 2);
                currentXNew = parent(currentY, currentX, 3);
                currentY=currentYNew;
                currentX=currentXNew;
            end         
            if(pathMax ~=0 && nnz(path)>pathMax)
                path=logical(zeros(sizeY,sizeX));
            end
        end