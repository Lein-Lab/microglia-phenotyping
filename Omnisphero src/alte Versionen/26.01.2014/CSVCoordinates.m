classdef CSVCoordinates < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        OffsetSize = 512;
    end
    
    properties
        Coordinates;
        CellPosMap;
        %Cellomics Positions
        CellPosMatrix;
        Well96Offsets;
        Well96Offsets484;
        Well96Offsets16Bit;
        Well48Offsets;
        Well48Offsets16Bit;
        WellWholeOTOffsets;
        WellWholeOTOffsets16Bit;
        WellSingleChamberOffsets;
        WellSingleChamberOffsets16Bit;
        %NeuronNucleusMapMatrix;
        %Manual Neuron Positions:
        ManualNeuronPositions;
        ManualNeuronPositionsSparse;
        %Automatic Neuron Positions by Edge Composite
        NeuronPositionsEdgeComposite;
        
        NeuronPositionsEdgeFill;
        NeuronPositionsEdgeFillSecond;
        
        %Data structure: Hash with xPos, yPos as Key
        %Values: Array with param.
        %1. Param: Skel marked (1=Skeleton found from Edge Fill, 2=Skeleton
        %found from EdgeFill Second, 3=Skeleton found from endpoint of Skeleton)
        %2. Param: Edge fill marked
        %3. Param: Edge fill marked second
        %4. Param: Manual marked
        %5. Param: Deleted by Skel (1=Deleted because too less Neurite
        %length or Distance between Skeleton endpoints, 2=Deleted because
        %another Neuron on same Skeleton and both Neurons are too near
        %together
        %6. Param: Edge fill overlap
        %7. Param: Neurite length
        %8. Param: Max dist between Skeleton Endpoints
        
        NeuronStatMatrix;
        
        NeuronPositionsSkelDeleted;
        NeuronPositionsSkeletonization;
        
        NeuronPositionsEdgeFillNeurite;
        NeuriteLengthMatrix;
        
        ManualPositions1;
        ManualPositions2;
        ManualPositions3;
        ManualPositions4;
        
        ManualCountingMode;
        ManualDeletingMode;
        StretchlimResult;
        IsodataResult;
        MinimumResult;
        ThresholdingLevel;
        ThresholdingLevelStrong;
        
        AreaDictionary;
    end

    %Singleton Pattern
    methods (Access = public)
        %Constructor
        function e = CSVCoordinates()
            %Fill Offsets for Whole OT
            structureMatrixOTWhole = zeros(20,43);            
            currentPos = [1,1];
            direction = Direction.Right;
            for i=860:-1:1
                e.WellWholeOTOffsets{i} = {e.OffsetSize * (currentPos(1)-1), e.OffsetSize * (currentPos(2)-1)};              
                structureMatrixOTWhole(currentPos(1),currentPos(2)) = 1;
                if(direction == Direction.Right)
                    if(currentPos(1) + 1 <= 20 && structureMatrixOTWhole(currentPos(1)+1, currentPos(2)) == 0)
                        currentPos(1) = currentPos(1) + 1;
                    else
                        currentPos(2) = currentPos(2) + 1;
                        direction = Direction.Down;
                    end %if
                elseif(direction == Direction.Down)
                    if(currentPos(2) + 1 <= 43 && structureMatrixOTWhole(currentPos(1), currentPos(2)+1) == 0)
                        currentPos(2) = currentPos(2) + 1;
                    else
                        currentPos(1) = currentPos(1) - 1;
                        direction = Direction.Left;
                    end %if
                elseif(direction == Direction.Left)
                    if(currentPos(1) - 1 >= 1 && structureMatrixOTWhole(currentPos(1)-1, currentPos(2)) == 0)
                        currentPos(1) = currentPos(1) - 1;
                    else
                        currentPos(2) = currentPos(2) - 1;
                        direction = Direction.Up;
                    end %if
                elseif(direction == Direction.Up)
                    if(currentPos(2) - 1 >= 1 && structureMatrixOTWhole(currentPos(1), currentPos(2)-1) == 0)
                        currentPos(2) = currentPos(2) - 1;
                    else
                        currentPos(1) = currentPos(1) + 1;
                        direction = Direction.Right;
                    end %if
                end %if
            end %for
            
            
            
            %Fill Offsets for OT Single Chamber
            structureMatrixOTSingle = zeros(21,27);
            currentPos = [1,1];
            direction = Direction.Right;
            for i=567:-1:1
                e.WellSingleChamberOffsets{i} = {e.OffsetSize * (currentPos(1)-1), e.OffsetSize * (currentPos(2)-1)};              
                structureMatrixOTSingle(currentPos(1),currentPos(2)) = 1;
                if(direction == Direction.Right)
                    if(currentPos(1) + 1 <= 21 && structureMatrixOTSingle(currentPos(1)+1, currentPos(2)) == 0)
                        currentPos(1) = currentPos(1) + 1;
                    else
                        currentPos(2) = currentPos(2) + 1;
                        direction = Direction.Down;
                    end %if
                elseif(direction == Direction.Down)
                    if(currentPos(2) + 1 <= 27 && structureMatrixOTSingle(currentPos(1), currentPos(2)+1) == 0)
                        currentPos(2) = currentPos(2) + 1;
                    else
                        currentPos(1) = currentPos(1) - 1;
                        direction = Direction.Left;
                    end %if
                elseif(direction == Direction.Left)
                    if(currentPos(1) - 1 >= 1 && structureMatrixOTSingle(currentPos(1)-1, currentPos(2)) == 0)
                        currentPos(1) = currentPos(1) - 1;
                    else
                        currentPos(2) = currentPos(2) - 1;
                        direction = Direction.Up;
                    end %if
                elseif(direction == Direction.Up)
                    if(currentPos(2) - 1 >= 1 && structureMatrixOTSingle(currentPos(1), currentPos(2)-1) == 0)
                        currentPos(2) = currentPos(2) - 1;
                    else
                        currentPos(1) = currentPos(1) + 1;
                        direction = Direction.Right;
                    end %if
                end %if
            end %for
            
            
            %Fill Offsets for OT Single Chamber in 16 Bit pattern
            structureMatrixOTSingle = zeros(21,27);
            currentPos = [1,1];
            direction = Direction.Right;
            for i=567:-1:1
                e.WellSingleChamberOffsets16Bit{i} = {e.OffsetSize * (currentPos(1)-1), e.OffsetSize * (currentPos(2)-1)};              
                structureMatrixOTSingle(currentPos(1),currentPos(2)) = 1;
                if(direction == Direction.Right)
                    if(currentPos(1) + 1 <= 21 && structureMatrixOTSingle(currentPos(1)+1, currentPos(2)) == 0)
                        currentPos(1) = currentPos(1) + 1;
                    else
                        currentPos(2) = currentPos(2) + 1;
                        direction = Direction.Down;
                    end %if
                elseif(direction == Direction.Down)
                    if(currentPos(2) + 1 <= 27 && structureMatrixOTSingle(currentPos(1), currentPos(2)+1) == 0)
                        currentPos(2) = currentPos(2) + 1;
                    else
                        currentPos(1) = currentPos(1) - 1;
                        direction = Direction.Left;
                    end %if
                elseif(direction == Direction.Left)
                    if(currentPos(1) - 1 >= 1 && structureMatrixOTSingle(currentPos(1)-1, currentPos(2)) == 0)
                        currentPos(1) = currentPos(1) - 1;
                    else
                        currentPos(2) = currentPos(2) - 1;
                        direction = Direction.Up;
                    end %if
                elseif(direction == Direction.Up)
                    if(currentPos(2) - 1 >= 1 && structureMatrixOTSingle(currentPos(1), currentPos(2)-1) == 0)
                        currentPos(2) = currentPos(2) - 1;
                    else
                        currentPos(1) = currentPos(1) + 1;
                        direction = Direction.Right;
                    end %if
                end %if
            end %for
            
            %Fill Offsets for 48er
            structureMatrix48 = zeros(36,36);
            currentPos = [1,1];
            direction = Direction.Right;
            for i=1296:-1:1
                e.Well48Offsets{i} = {e.OffsetSize * (currentPos(1)-1), e.OffsetSize * (currentPos(2)-1)};              
                structureMatrix48(currentPos(1),currentPos(2)) = 1;
                if(direction == Direction.Right)
                    if(currentPos(1) + 1 <= 36 && structureMatrix48(currentPos(1)+1, currentPos(2)) == 0)
                        currentPos(1) = currentPos(1) + 1;
                    else
                        currentPos(2) = currentPos(2) + 1;
                        direction = Direction.Down;
                    end %if
                elseif(direction == Direction.Down)
                    if(currentPos(2) + 1 <= 36 && structureMatrix48(currentPos(1), currentPos(2)+1) == 0)
                        currentPos(2) = currentPos(2) + 1;
                    else
                        currentPos(1) = currentPos(1) - 1;
                        direction = Direction.Left;
                    end %if
                elseif(direction == Direction.Left)
                    if(currentPos(1) - 1 >= 1 && structureMatrix48(currentPos(1)-1, currentPos(2)) == 0)
                        currentPos(1) = currentPos(1) - 1;
                    else
                        currentPos(2) = currentPos(2) - 1;
                        direction = Direction.Up;
                    end %if
                elseif(direction == Direction.Up)
                    if(currentPos(2) - 1 >= 1 && structureMatrix48(currentPos(1), currentPos(2)-1) == 0)
                        currentPos(2) = currentPos(2) - 1;
                    else
                        currentPos(1) = currentPos(1) + 1;
                        direction = Direction.Right;
                    end %if
                end %if
            end %for
                        
            
            %Fill Offsets for 96er
            structureMatrix96 = zeros(14,14);
            currentPos = [1,1];
            direction = Direction.Right;
            for i=196:-1:1
                e.Well96Offsets{i} = {e.OffsetSize * (currentPos(1)-1), e.OffsetSize * (currentPos(2)-1)};         
                structureMatrix96(currentPos(1),currentPos(2)) = 1;
                if(direction == Direction.Right)
                    if(currentPos(1) + 1 <= 14 && structureMatrix96(currentPos(1)+1, currentPos(2)) == 0)
                        currentPos(1) = currentPos(1) + 1;
                    else
                        currentPos(2) = currentPos(2) + 1;
                        direction = Direction.Down;
                    end %if
                elseif(direction == Direction.Down)
                    if(currentPos(2) + 1 <= 14 && structureMatrix96(currentPos(1), currentPos(2)+1) == 0)
                        currentPos(2) = currentPos(2) + 1;
                    else
                        currentPos(1) = currentPos(1) - 1;
                        direction = Direction.Left;
                    end %if
                elseif(direction == Direction.Left)
                    if(currentPos(1) - 1 >= 1 && structureMatrix96(currentPos(1)-1, currentPos(2)) == 0)
                        currentPos(1) = currentPos(1) - 1;
                    else
                        currentPos(2) = currentPos(2) - 1;
                        direction = Direction.Up;
                    end %if
                elseif(direction == Direction.Up)
                    if(currentPos(2) - 1 >= 1 && structureMatrix96(currentPos(1), currentPos(2)-1) == 0)
                        currentPos(2) = currentPos(2) - 1;
                    else
                        currentPos(1) = currentPos(1) + 1;
                        direction = Direction.Right;
                    end %if
                end %if
            end %for
            
            
            
            %Fill Offsets for 96er 484 Pics
            structureMatrix96484 = zeros(22,22);
            currentPos = [1,1];
            direction = Direction.Right;
            for i=484:-1:1
                e.Well96Offsets484{i} = {e.OffsetSize * (currentPos(1)-1), e.OffsetSize * (currentPos(2)-1)};         
                structureMatrix96484(currentPos(1),currentPos(2)) = 1;
                if(direction == Direction.Right)
                    if(currentPos(1) + 1 <= 22 && structureMatrix96484(currentPos(1)+1, currentPos(2)) == 0)
                        currentPos(1) = currentPos(1) + 1;
                    else
                        currentPos(2) = currentPos(2) + 1;
                        direction = Direction.Down;
                    end %if
                elseif(direction == Direction.Down)
                    if(currentPos(2) + 1 <= 22 && structureMatrix96484(currentPos(1), currentPos(2)+1) == 0)
                        currentPos(2) = currentPos(2) + 1;
                    else
                        currentPos(1) = currentPos(1) - 1;
                        direction = Direction.Left;
                    end %if
                elseif(direction == Direction.Left)
                    if(currentPos(1) - 1 >= 1 && structureMatrix96484(currentPos(1)-1, currentPos(2)) == 0)
                        currentPos(1) = currentPos(1) - 1;
                    else
                        currentPos(2) = currentPos(2) - 1;
                        direction = Direction.Up;
                    end %if
                elseif(direction == Direction.Up)
                    if(currentPos(2) - 1 >= 1 && structureMatrix96484(currentPos(1), currentPos(2)-1) == 0)
                        currentPos(2) = currentPos(2) - 1;
                    else
                        currentPos(1) = currentPos(1) + 1;
                        direction = Direction.Right;
                    end %if
                end %if
            end %for
            
            
            %Fill Offsets for 96er in 16 Bit pattern
            structureMatrix96 = zeros(14,14);
            currentPos = [1,14];
            direction = Direction.Right;
            for i=196:-1:1
                e.Well96Offsets16Bit{i} = {e.OffsetSize * (currentPos(1)-1), e.OffsetSize * (currentPos(2)-1)};         
                structureMatrix96(currentPos(1),currentPos(2)) = 1;
                if(direction == Direction.Right)
                    if(currentPos(1) + 1 <= 14 && structureMatrix96(currentPos(1)+1, currentPos(2)) == 0)
                        currentPos(1) = currentPos(1) + 1;
                    else
                        currentPos(2) = currentPos(2) - 1;
                        direction = Direction.Up;
                    end %if
                elseif(direction == Direction.Down)
                    if(currentPos(2) + 1 <= 14 && structureMatrix96(currentPos(1), currentPos(2)+1) == 0)
                        currentPos(2) = currentPos(2) + 1;
                    else
                        currentPos(1) = currentPos(1) + 1;
                        direction = Direction.Right;
                    end %if
                elseif(direction == Direction.Left)
                    if(currentPos(1) - 1 >= 1 && structureMatrix96(currentPos(1)-1, currentPos(2)) == 0)
                        currentPos(1) = currentPos(1) - 1;
                    else
                        currentPos(2) = currentPos(2) + 1;
                        direction = Direction.Down;
                    end %if
                elseif(direction == Direction.Up)
                    if(currentPos(2) - 1 >= 1 && structureMatrix96(currentPos(1), currentPos(2)-1) == 0)
                        currentPos(2) = currentPos(2) - 1;
                    else
                        currentPos(1) = currentPos(1) - 1;
                        direction = Direction.Left;
                    end %if
                end %if
            end %for
        end
    
        function i=GetPictureNumberFromPosition(self, xPos, yPos, wellType)
            if(wellType == WellType.Well96)
                dict = self.Well96Offsets;
            elseif(wellType == WellType.Well48)
                dict = self.Well48Offsets;
            elseif(wellType == WellType.Well96484)
                dict = self.Well96Offsets484;
            elseif(wellType == WellType.OTSingle)
                dict = self.WellSingleChamberOffsets;
            elseif(wellType == WellType.OTWhole)
                dict = self.WellWholeOTOffsets;
            elseif(wellType == WellType.Well4816)
                dict = self.Well48Offsets16Bit;
            elseif(wellType == WellType.Well9616)
                 dict = self.Well96Offsets16Bit;
            elseif(wellType == WellType.OTSingle16)
                dict = self.WellSingleChamberOffsets16Bit;
            elseif(wellType == WellType.OTWhole16)
                dict = self.WellWholeOTOffsets16Bit;
            end
            %Iterate over dictionary and check, where values are = Position
            for j=1:numel(dict)
                currentValue = dict{j};
                if(xPos == currentValue{1} && yPos == currentValue{2})
                    i=j-1;
                    break;
                end
            end
        end

       function [x y] = CalculateGlobalCoordinatesCSV(self, field, leftIndex, topIndex, wellType)
           %Formel für richtige Position:
           % Sqrt(ObjectArea) = Kantenlänge
           % Versatz = Sqrt(Kantenlänge^2+Kantenlänge^2) / 2
           % y = y-Versatz
           % x = x+Versatz
           
           if(wellType == WellType.Well96)
                x=self.Well96Offsets{field}{1} + leftIndex + 1;
                y=self.Well96Offsets{field}{2} + (512 - topIndex) + 1;
           elseif(wellType == WellType.Well9616)
                x=self.Well96Offsets{field}{1} + leftIndex + 1;
                y=self.Well96Offsets{field}{2} + (512 - topIndex) + 1;
           elseif(wellType == WellType.Well96484)
                x=self.Well96Offsets484{field}{1} + leftIndex + 1;
                y=self.Well96Offsets484{field}{2} + (512 - topIndex) + 1;
            elseif(wellType == WellType.Well48)
                x=self.Well48Offsets{field}{1} + leftIndex + 1;
                y=self.Well48Offsets{field}{2} + (512-topIndex) + 1;
            elseif(wellType == WellType.Well4816)
                x=self.Well48Offsets{field}{1} + leftIndex + 1;
                y=self.Well48Offsets{field}{2} + (512-topIndex) + 1;
            elseif(wellType == WellType.OTSingle)
                x=self.WellSingleChamberOffsets{field}{1} + leftIndex + 1;
                y=self.WellSingleChamberOffsets{field}{2} + (512-topIndex) + 1;
            elseif(wellType == WellType.OTSingle16)
                x=self.WellSingleChamberOffsets{field}{1} + leftIndex + 1;
                y=self.WellSingleChamberOffsets{field}{2} + (512-topIndex) + 1;
            elseif(wellType == WellType.OTWhole)
                x=self.WellWholeOTOffsets{field}{1} + leftIndex + 1;
                y=self.WellWholeOTOffsets{field}{2} + (512-topIndex) + 1;
            elseif(wellType == WellType.OTWhole16)
                x=self.WellWholeOTOffsets{field}{1} + leftIndex + 1;
                y=self.WellWholeOTOffsets{field}{2} + (512-topIndex) + 1;
            end %if
        end %function CalculateGlobalCoordinates
        
        function [x y] = CalculateGlobalCoordinatesImage(self, field, leftIndex, topIndex, wellType)
           %Formel für richtige Position:
           % Sqrt(ObjectArea) = Kantenlänge
           % Versatz = Sqrt(Kantenlänge^2+Kantenlänge^2) / 2
           % y = y-Versatz
           % x = x+Versatz
           
           if(wellType == WellType.Well96)
                x=self.Well96Offsets{field}{1} + leftIndex;
                y=self.Well96Offsets{field}{2} + (topIndex);
           elseif(wellType == WellType.Well9616)
                x=self.Well96Offsets16Bit{field}{1} + leftIndex;
                y=self.Well96Offsets16Bit{field}{2} + (topIndex);
                
           elseif(wellType == WellType.Well96484)
                x=self.Well96Offsets484{field}{1} + leftIndex;
                y=self.Well96Offsets484{field}{2} + (topIndex);
            elseif(wellType == WellType.Well48)
                x=self.Well48Offsets{field}{1} + leftIndex;
                y=self.Well48Offsets{field}{2} + (topIndex);
            elseif(wellType == WellType.Well4816)
                x=self.Well48Offsets16Bit{field}{1} + leftIndex;
                y=self.Well48Offsets16Bit{field}{2} + (topIndex);
            elseif(wellType == WellType.OTSingle)
                x=self.WellSingleChamberOffsets{field}{1} + leftIndex;
                y=self.WellSingleChamberOffsets{field}{2} + (topIndex);
            elseif(wellType == WellType.OTSingle16)
                x=self.WellSingleChamberOffsets16Bit{field}{1} + leftIndex;
                y=self.WellSingleChamberOffsets16Bit{field}{2} + (topIndex);
            elseif(wellType == WellType.OTWhole)
                x=self.WellWholeOTOffsets{field}{1} + leftIndex;
                y=self.WellWholeOTOffsets{field}{2} + (topIndex);
            elseif(wellType == WellType.OTWhole16)
                x=self.WellWholeOTOffsets16Bit{field}{1} + leftIndex;
                y=self.WellWholeOTOffsets16Bit{field}{2} + (topIndex);
            end %if
        end %function CalculateGlobalCoordinates
        
        
        function ReadCSVFile(self, datafile, wellType, varargin)
            waitbarHandle = waitbar(0,'Please wait. Your CSV data will soon be ready.');
            lineCounter = 0;
            fieldIndex = 0;
            topIndex = 0;
            leftIndex = 0;            
            nVarargs = length(varargin);
            if(nVarargs == 1)
                saveValues = varargin{1}; 
            elseif(nVarargs == 2)
                saveValues = varargin{1};
                filter = varargin{2};
            end
            self.CellPosMatrix = containers.Map(); 
            self.ManualNeuronPositionsSparse = containers.Map();
            headerMap = containers.Map();
            if exist(datafile, 'file') > 0
                filePointer = fopen(datafile,'r');
                fseek(filePointer, 0,'eof');
                filelength = ftell(filePointer);
                fseek(filePointer, 0,'bof');
                %Get indices of Columns "Field", "Top" and "Left"
                if(~feof(filePointer))
                    lineCounter = lineCounter+1;
                    line = fgetl(filePointer);
                    %Check if CSV is comma or semicolon separated
                    semiOccur = findstr(line,';');
                    commaOccur = findstr(line,',');
                    if(length(semiOccur) > length(commaOccur))
                        delimiter = ';';
                    else
                        delimiter = ',';
                    end
                    cellsplit = regexp(line,delimiter,'split');
                    for i=1:length(cellsplit)
                        headerMap(char(cellsplit(i))) = i;
                    end
                end
                counter = 0;
                while (~feof(filePointer))
                    counter = counter + 1;
                    waitbar((counter/filelength) * 213,waitbarHandle);
                    if(mod(counter,5000) == 0)
                        logMessage=['Read CSV Line number',' ', num2str(counter)];
                        disp(logMessage);
                    end                    
                    lineCounter = lineCounter+1;
                    line = fgetl(filePointer);
                    cellsplit = regexp(line,delimiter,'split');
                    
                    %if CellPosMatrix(Well) isn't initialized -> initialize
                    %it!
                    currentWellName = cellsplit{headerMap('Well')};
                    if(~isKey(self.CellPosMatrix, currentWellName))
                        if(wellType == WellType.Well96 || wellType == WellType.Well9616)
                            self.CellPosMatrix(currentWellName) = sparse(7168, 7168);
                            self.ManualNeuronPositionsSparse(currentWellName) = sparse(7168, 7168);
                        elseif(wellType == WellType.Well96484)
                            self.CellPosMatrix(currentWellName) = sparse(11264,11264);
                            self.ManualNeuronPositionsSparse(currentWellName) = sparse(11264,11264);
                        elseif(wellType == WellType.Well48 ||wellType == WellType.Well4816)
                            self.CellPosMatrix(currentWellName) = sparse(18432, 18432);
                            self.ManualNeuronPositionsSparse(currentWellName) = sparse(18432, 18432);
                        elseif(wellType == WellType.OTSingle || wellType == WellType.OTSingle16)
                            self.CellPosMatrix(currentWellName) = sparse(13824, 10752);
                            self.ManualNeuronPositionsSparse(currentWellName) = sparse(13824, 10752);
                        elseif(wellType == WellType.OTWhole || wellType == WellType.OTWhole16)
                            self.CellPosMatrix(currentWellName) = sparse(22016, 10240);
                            self.ManualNeuronPositionsSparse(currentWellName) = sparse(22016, 10240);
                        end %if
                    end
                    
                    
                    %Create new CellPosition Object.
                    %newCellPos = CellPosition();                    
                    %Filter Format: Array of 2D Array
                    %1st Dimension of 2D Array: Name of Dimension
                    %2st Dimension of 2D Array: Operator
                    %3nd Dimension of 2D Array: Value
                    allMatch = 1;
                    if(exist('filter', 'var'))
                        for i=1:numel(filter)
                            currentFilter = filter{i};
                            attributeName = currentFilter{1};
                            operator = currentFilter{2};
                            operand = currentFilter{3};
                            value = cellsplit{headerMap(attributeName)};
                            if(~isempty(value) && all(ismember(value,'0123456789')))
                                value = str2num(value);
                            elseif(numel(value) ~= numel(operand))
                                    allMatch = 0;
                            end 
                            if(operator == '=' && allMatch == 1)
                                if (value ~= operand)
                                    allMatch = 0;
                                end
                            elseif(operator == '>')
                                if (value <= operand)
                                    allMatch=0;
                                end
                            elseif(operator == '<')
                                if(value >= operand)
                                    allMatch=0;
                                end
                            end
                        end
                    end %if
                    
                    %Name, its operator and expected value
                    
                    %newCellPos.EventTypeProfile = cellsplit{headerMap('EventTypeProfile')};
                    %newCellPos.EventType1Status = cellsplit{headerMap('EventType1Status')};
                    %newCellPos.EventType2Status = cellsplit{headerMap('EventType2Status')};
                    %newCellPos.EventType3Status = cellsplit{headerMap('EventType3Status')};
                    %newCellPos.ObjectAreaCh1 = cellsplit{headerMap('ObjectAreaCh1')};
                    %newCellPos.ObjectAreaCh1Status = cellsplit{headerMap('ObjectAreaCh1Status')};
                    %newCellPos.ObjectShapeP2ACh1 = cellsplit{headerMap('ObjectShapeP2ACh1')};
                    %newCellPos.ObjectShapeP2ACh1Status = cellsplit{headerMap('ObjectShapeP2ACh1Status')};
                    %newCellPos.ObjectShapeLWRCh1 = cellsplit{headerMap('ObjectShapeLWRCh1')};
                    %newCellPos.ObjectShapeLWRCh1Status = cellsplit{headerMap('ObjectShapeLWRCh1Status')};
                    %newCellPos.ObjectTotalIntenCh1 = cellsplit{headerMap('ObjectTotalIntenCh1')};
                    %newCellPos.ObjectTotalIntenCh1Status = cellsplit{headerMap('ObjectTotalIntenCh1Status')};
                    %newCellPos.ObjectAvgIntenCh1 = cellsplit{headerMap('ObjectAvgIntenCh1')};
                    %newCellPos.ObjectAvgIntenCh1Status = cellsplit{headerMap('ObjectAvgIntenCh1Status')};
                    %newCellPos.ObjectVarIntenCh1 = cellsplit{headerMap('ObjectVarIntenCh1')};
                    %newCellPos.ObjectVarIntenCh1Status = cellsplit{headerMap('ObjectVarIntenCh1Status')};
                                      
                    if(allMatch == 1)
                        
                        fieldIndexCell = str2num((cellsplit{headerMap('Field')}));
                        if(isKey(headerMap, 'XCentroid'))
                            leftIndexCell = strrep(cellsplit{headerMap('XCentroid')},',','.');
                        else
                            leftIndexCell = strrep(cellsplit{headerMap('Left')},',','.');
                        end
                        if(isKey(headerMap, 'YCentroid'))
                            topIndexCell = strrep(cellsplit{headerMap('YCentroid')},',','.');  
                        else
                            topIndexCell = strrep(cellsplit{headerMap('Top')},',','.');  
                        end
                                          
                        leftIndexCell = round(str2num(leftIndexCell));
                        topIndexCell = round(str2num(topIndexCell));
                        %Cellomics CSV Coordinates are 0 based. Matlab is 1
                        %based! Increment X and Y by 1.
                        leftIndexCell = leftIndexCell + 1;
                        topIndexCell = topIndexCell + 1;
                        %newCellPos.Field = fieldIndexCell;
                        [x y] = self.CalculateGlobalCoordinatesCSV(fieldIndexCell, leftIndexCell, topIndexCell, wellType);
                   
                        %TODO: Filter
                        %Create Mask of Cell Positions 
                        currentMatrix = self.CellPosMatrix(currentWellName);
                        currentMatrix(y,x) = 1;
                        self.CellPosMatrix(currentWellName) = currentMatrix;
                        if(exist('saveValues', 'var'))
                            %newCellPos.AttributeList = cell(numel(saveValues));
                            cellPosAttributes = containers.Map();
                            for i=1:numel(saveValues)
                                attributeName = saveValues{i};
                                attributeValue = cellsplit{headerMap(attributeName)};
                                cellPosAttributes(attributeName) = attributeValue;
                            end
                            %newCellPos.AttributeList = cellPosAttributes;
                        end
                    end
                   %Comment out for complex Data Structure
                   %if(isKey(self.CellPosMap,num2str(newCellPos.Field)))
                   %    self.CellPosMap(num2str(newCellPos.Field)) = [self.CellPosMap(num2str(newCellPos.Field)) newCellPos];
                   %else
                   %    self.CellPosMap(num2str(newCellPos.Field)) = [newCellPos];
                   %end %if isKey                   
                  
                end
            else
                return;
            end %if exist
            waitbar(100,waitbarHandle);
            close(waitbarHandle);
        end %function ReadCSVFile
        
        %Get For NeuronPosition col, row nearest Nucleus in NucleusM
        function [nucCol nucRow euclidDist] = FindNucleusForNeuron(self, neuCol, neuRow, NucleusM, maxDist)
           nucCol = -1;
           nucRow = -1;
           found = 0;
           iterator=0;
           euclidDist=99999;
           x=neuCol;
           y=neuRow;
           [sizeY, sizeX] = size(NucleusM);
           while(found < 1 && iterator<maxDist)
               iterator=iterator+1;
               for i=y-iterator:y+iterator
               for j=x-iterator:x+iterator
                   %Check if point is in bounds.
                   if(i>0 && j<=sizeX && i<=sizeY && j>0)
                        %Check if match                                       
                        if(NucleusM(i,j) == 1)
                            currentDistVec = [i j;y x];         
                            dist = pdist(currentDistVec, 'euclidean');
                            if(dist < euclidDist)                                
                                nucCol = j;
                                nucRow = i;
                                euclidDist = dist;
                                found = 1;
                            end
                        end %if
                   end %if
               end %for
               end %for
          end %while
        end %function FindNucleusForNeuron   
        
        function GetSize(this)
          props = properties(this);
          totSize = 0;
          for ii=1:length(props)
            currentProperty = getfield(this, char(props(ii)));
            s = whos('currentProperty');
            totSize = totSize + s.bytes;
          end
          fprintf(1, '%d bytes\n', totSize);
        end
        
    end %methods    
end %classdef

