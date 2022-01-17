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

function varargout = NeuronDetectorGUI(varargin)
% NEURONDETECTORGUI MATLAB code for NeuronDetectorGUI.fig
%      NEURONDETECTORGUI, by itself, creates a new NEURONDETECTORGUI or raises the existing
%      singleton*.
%
%      H = NEURONDETECTORGUI returns the handle to a new NEURONDETECTORGUI or the handle to
%      the existing singleton*.
%
%      NEURONDETECTORGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEURONDETECTORGUI.M with the given input arguments.
%
%      NEURONDETECTORGUI('Property','Value',...) creates a new NEURONDETECTORGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NeuronDetectorGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NeuronDetectorGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NeuronDetectorGUI

% Last Modified by GUIDE v2.5 27-May-2013 17:27:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NeuronDetectorGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @NeuronDetectorGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before NeuronDetectorGUI is made visible.
function NeuronDetectorGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NeuronDetectorGUI (see VARARGIN)

% Choose default command line output for NeuronDetectorGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NeuronDetectorGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NeuronDetectorGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function LoadCSV_Callback(hObject, eventdata, handles)
% hObject    handle to LoadCSV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.csv', 'Select CSV File with Cellomics Data');
handles = guidata(handles.figure1);
if (~isempty(filename) && ~isempty(pathname))
    csvHandler = CSVCoordinates();
    filepath = strcat(pathname,filename);
    if(get(handles.rb48,'Value'))
        csvHandler.ReadCSVFile(filepath, WellType.Well48);
    elseif(get(handles.rb4816,'Value'))
        csvHandler.ReadCSVFile(filepath, WellType.Well4816);
    elseif(get(handles.rb96,'Value'))
        csvHandler.ReadCSVFile(filepath, WellType.Well96);
    elseif(get(handles.rb9616,'Value'))
        csvHandler.ReadCSVFile(filepath, WellType.Well9616);
    elseif(get(handles.rbWholeOt,'Value'))
        csvHandler.ReadCSVFile(filepath, WellType.OTWhole);
    elseif(get(handles.rbWholeOt16,'Value'))
        csvHandler.ReadCSVFile(filepath, WellType.OTWhole16);
    elseif(get(handles.rbOtSingleChamber,'Value'))
        csvHandler.ReadCSVFile(filepath, WellType.OTSingle);
    elseif(get(handles.rbOtSingleChamber16,'Value'))
        csvHandler.ReadCSVFile(filepath, WellType.OTSingle16);
    end    
    handles.CSVCoordinates = csvHandler;
    
    guidata(handles.figure1, handles)
end% --------------------------------------------------------------------


function LoadNeuronCSV_Callback(hObject, eventdata, handles)
% hObject    handle to LoadNeuronCSV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.csv', 'Select CSV File with Cellomics Data');
handles = guidata(handles.figure1);
if (~isempty(filename) && ~isempty(pathname))
    neuronHandler = CSVCoordinates();
    filepath = strcat(pathname,filename);
    filterETP = cell(3,1);
    filterETP{1} = 'EventTypeProfile';
    filterETP{2} = '=';
    filterETP{3} = 1;
    filterArray = cell(1,1);    
    filterArray{1} = filterETP;
    saveValuesArray = cell(1,1);
    saveValuesArray{1} = 'EventTypeProfile';
    if(get(handles.rb48,'Value'))
        neuronHandler.ReadCSVFile(filepath, WellType.Well48, saveValuesArray, filterArray);
    elseif(get(handles.rb4816,'Value'))
        neuronHandler.ReadCSVFile(filepath, WellType.Well4816, saveValuesArray, filterArray);
    elseif(get(handles.rb96,'Value'))
        neuronHandler.ReadCSVFile(filepath, WellType.Well96, saveValuesArray, filterArray);
    elseif(get(handles.rb9616,'Value'))
        neuronHandler.ReadCSVFile(filepath, WellType.Well9616, saveValuesArray, filterArray);
    elseif(get(handles.rbWholeOt,'Value'))
        neuronHandler.ReadCSVFile(filepath, WellType.OTWhole, saveValuesArray, filterArray);
    elseif(get(handles.rbWholeOt16,'Value'))
        neuronHandler.ReadCSVFile(filepath, WellType.OTWhole16, saveValuesArray, filterArray);
    elseif(get(handles.rbOtSingleChamber,'Value'))
        neuronHandler.ReadCSVFile(filepath, WellType.OTSingle, saveValuesArray, filterArray);
    elseif(get(handles.rbOtSingleChamber16,'Value'))
        neuronHandler.ReadCSVFile(filepath, WellType.OTSingle16, saveValuesArray, filterArray);
    end    
    neuronHandler.ManualNeuronPositions = containers.Map();
    handles.NeuronCoordinates = neuronHandler;
    guidata(handles.figure1, handles)
end

% --------------------------------------------------------------------
% Calculates for every point in NeuronCoordinates its regarding point in
% CSVCoordinates and delivers medium distance
function CalcNucleusNeuronDist_Callback(hObject, eventdata, handles)
% hObject    handle to CalcNucleusNeuronDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(handles.figure1);
if(isfield(handles,'NeuronCoordinates') && isfield(handles,'CSVCoordinates'))
           %Map Neuron Positions to Nucleus Positions
    selectedWellNumber = get(handles.lbWell,'Value');
    wellList = get(handles.lbWell, 'string');
    selectedWell = wellList{selectedWellNumber};
    csvHandler = handles.CSVCoordinates;
    neuronHandler = handles.NeuronCoordinates;
    WellNeuronDict = neuronHandler.CellPosMatrix;
    WellNucleusDict = csvHandler.CellPosMatrix;
    waitbarHandle = waitbar(0,'Please wait. I am currently mapping Neurons to Nuclei.');
    for i=1:numel(wellList)
      waitbar((i/numel(wellList)),waitbarHandle);
      selectedWell = wellList{i};  
      if(isKey(WellNeuronDict,selectedWell))
       NeuronM = WellNeuronDict(selectedWell);
       NucleusM = WellNucleusDict(selectedWell);
       distsum=0;
       [neuronRows, neuronCols] =  find(NeuronM);
       [nucleusRows, nucleusCols] = find(NucleusM);
    
        
       %For every Pixel in NeuronM: Get regarding Pixel in NucleusM with
       %minimum Euclidian Distance.
       %Save that in Mapping Object.
       %Data Structure:
       %Dictionary with NucleusM Matrix Positions as Key and NeuronM
       %Position as value. Save this Martix as object in neuronHandler.
       %neuronHandler.NeuronNucleusMapMatrix = containers.Map();
       
       for i=1:numel(neuronRows)
          row = neuronRows(i);
          col = neuronCols(i);

          [xMapped,yMapped,distance] = neuronHandler.FindNucleusForNeuron(col,row,NucleusM);
          distsum=distsum+distance;
          NeuronM(row,col) = 0;
          NeuronM(yMapped,xMapped) = 1;
       end
       WellNeuronDict(selectedWell) = NeuronM;
       neuronHandler.CellPosMatrix = WellNeuronDict;
       distsum=distsum/numel(nucleusRows);
       disp('Medium distance is ');
       disp(distsum);   
      end
    end
    close(waitbarHandle);
else
    disp('ERROR. Neuron and Nucles CSV need both to be loaded.');
end
guidata(handles.figure1, handles)


% --------------------------------------------------------------------
function LoadImageFolder_Callback(hObject, eventdata, handles)
% hObject    handle to LoadImageFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Restructure this method.
% Currently: Image Folder is saved. Only current Well is saved in memory.
% New: Save Image Folder. Create Pics for all Wells and save them in own
% SubFolder of Folder.

foldername = uigetdir;
handles = guidata(handles.figure1);
if (~isempty(foldername))
  if(isfield(handles,'ImageHandler'))        
        imageHandler = handles.ImageHandler;
  else
        imageHandler = ImageHandler();
        imageHandler.NegativeFilters = containers.Map();
        imageHandler.PositiveFilters = containers.Map();
  end
  if(isfield(handles,'CSVCoordinates'))
        csvHandler = handles.CSVCoordinates;
  else
        csvHandler = CSVCoordinates();
  end
  handles.CSVCoordinates = csvHandler;
  guidata(handles.figure1, handles)
  imageHandler.Foldername = foldername;
  maxValueWholeImage=0;
  maxFilterList = zeros(0);
  
    %Check if Subfolder with files already exists.
    subfoldername = [foldername '/ConvertedCellomics'];
    if(~exist(subfoldername, 'dir') || ~numel(dir(subfoldername)) > 10)
        waitbarHandle = waitbar(0,'Please wait. All Cellomics images will be preprocessed. This could take a while. Go and have a coffee!');
        selectedWellNumber = get(handles.lbWell,'Value');
        wellList = get(handles.lbWell, 'string');
        selectedWell = wellList{selectedWellNumber};
    
        % Load all images in folder.
        allFiles = dir(foldername); 
        %Check if Subfolder for own files exists. If not: Create
        mkdir(foldername,'ConvertedCellomics');
        subfoldername = [foldername '\ConvertedCellomics'];
    
        %Outer Loop: Iterate over all Wells.
        for currentWellIndex=1:numel(wellList)
         imageHandler.NucleusPicArray = cell(2500);
         imageHandler.NeuritePicArray = cell(2500);
         waitbar((currentWellIndex/numel(wellList)),waitbarHandle);
        % 1. Load all Pics for current well.
          currentWell = wellList{currentWellIndex};
          currentWellWithoutNull = currentWell;
          if(length(currentWell) == 2)
             currentWell = [currentWell(1) '0' currentWell(2)];
          end
          oneFileFound = 0;
          for i=1:numel(allFiles)
            nucleusRegexp = strcat('ARRAY.*\d*_(',currentWell,')f(\d+)d0');
            neuriteRegexp = strcat('ARRAY.*\d*_(',currentWell,')f(\d+)d1');
            tokensNucleus = regexpi(allFiles(i).name, nucleusRegexp, 'tokens');
            tokensNeurite = regexpi(allFiles(i).name, neuriteRegexp, 'tokens');
        
            if(length(tokensNucleus) > 0)
                tokensNucleus = tokensNucleus{1};
                fileWell = tokensNucleus{1};
                oneFileFound = 1;
            elseif(length(tokensNeurite) > 0)
                tokensNeurite = tokensNeurite{1};
                fileWell = tokensNeurite{1};
                oneFileFound = 1;
            end
            if((length(tokensNucleus) > 1) && strcmp(fileWell, currentWell))
                %Nucleus Pic. Get index number and add to nucleus Dictionary
                index = tokensNucleus{2};
                imageHandler.NucleusPicArray{str2num(index) + 1} = allFiles(i).name;
            elseif((length(tokensNeurite) > 1) && strcmp(fileWell, currentWell))
                %Neurite Pic. Get index number and add to neurite Dictionary
                index = tokensNeurite{2};
                imageHandler.NeuritePicArray{str2num(index) + 1} = allFiles(i).name;
            end %if
        end %for
    %Check if there is at least one file for well. Otherwise continue with
    %next well.
    if(oneFileFound==0)
        continue;
    end

    % 2. Load all Neurite Pics in the right order to Axes
    %Differentiate between different Well Types
    
    
    %Get Snail Pictures in correct order
    %Idea: Calculate for each Field Index global X and Y Coordinate
    %Save in 2 Matrices: Keys: X,Y-Position ; Values: Field Index
    %Sort this HashMap and Fill Axes with correct pictures
    
    %KeyArray is a Matrix. 1. Row: Field Index. 2. Row: X Values. 3. Row: Y
    % Values
       
        wellType = WellType.Undefined;
        if(get(handles.rb48,'Value'))
         maxIndex = numel(csvHandler.Well48Offsets);
         wellType = WellType.Well48;
         KeyArray = zeros(maxIndex,3);
        elseif(get(handles.rb4816,'Value'))
         maxIndex = numel(csvHandler.Well48Offsets16Bit);
         wellType = WellType.Well4816;
         KeyArray = zeros(maxIndex,3);
         
        elseif(get(handles.rb96,'Value'))
         maxIndex = numel(csvHandler.Well96Offsets);
         wellType = WellType.Well96;
         KeyArray = zeros(maxIndex,3);
        elseif(get(handles.rb9616,'Value'))
         maxIndex = numel(csvHandler.Well96Offsets16Bit);
         wellType = WellType.Well9616;
         KeyArray = zeros(maxIndex,3);
        elseif(get(handles.rbWholeOt,'Value'))
         maxIndex = numel(csvHandler.WellWholeOTOffsets);
         wellType = WellType.OTWhole;
         KeyArray = zeros(maxIndex,3);
        elseif(get(handles.rbWholeOt16,'Value'))
         maxIndex = numel(csvHandler.WellWholeOTOffsets16Bit);
         wellType = WellType.OTWhole16;
         KeyArray = zeros(maxIndex,3);
        elseif(get(handles.rbOtSingleChamber,'Value'))
         maxIndex = numel(csvHandler.WellSingleChamberOffsets);
         wellType = WellType.OTSingle;
         KeyArray = zeros(maxIndex,3);
        elseif(get(handles.rbOtSingleChamber16,'Value'))
         maxIndex = numel(csvHandler.WellSingleChamberOffsets16Bit);
         wellType = WellType.OTSingle16;
         KeyArray = zeros(maxIndex,3);
        end
    
        for i=1:maxIndex;
            [x y] = csvHandler.CalculateGlobalCoordinatesImage(i,0,0,wellType);
            %Save Index, because Matrix will be sorted.
            KeyArray(i,1) = i;
            KeyArray(i,2) = x;
            KeyArray(i,3) = y;
        end
    
        %Sort KeyArray
        SortedKeyArray = sortrows(KeyArray,[3,2]);
    
        y=SortedKeyArray(1,3);
        x=SortedKeyArray(1,2);
        previousY = y;
        i=1;
        neuriteImage = zeros(0,0);
        nucleusImage = zeros(0,0);
        sixteenBit = 0;
    
        %Fill Image Row by Row, with SortedKeyArray
        breaked = 0;
        highestNucleusValue=0;
        highestNeuriteValue=0;
        maxValueWhole=0;
        
        while(i<maxIndex)
            neuriteImageLine = zeros(0,0);
            nucleusImageLine = zeros(0,0);
            while(y==previousY && i<=maxIndex)
              x = SortedKeyArray(i,2);
              y = SortedKeyArray(i,3);
              imageNumber = SortedKeyArray(i,1);
              filenamenuc = imageHandler.NucleusPicArray(imageNumber);
              filenameneu = imageHandler.NeuritePicArray(imageNumber);
              filenamenuc = filenamenuc{1};
              filenameneu = filenameneu{1};
             % if((length(filenamenuc) < 2) || (length(filenameneu) < 2))
             %     %Only if ONE images is available for the well
             %     breaked = 1;
             %     break;
             % end
              imageString = strcat(foldername, '\', filenamenuc);
              if(exist(imageString,'file') && ~isempty(filenamenuc))
                currentNucleusPic = imread(imageString);
              else
                currentNucleusPic = zeros(512,512);


              end


              maxValueNucleus = max(currentNucleusPic(:));
              imageString = strcat(foldername, '\', filenameneu);
              if(exist(imageString,'file') && ~ isempty(filenameneu))
                  currentNeuritePic = imread(imageString);
              else
                currentNeuritePic = zeros(512,512);
              end
              maxValueNeurite = max(currentNeuritePic(:));
              maxValueWhole = max(maxValueNucleus,maxValueNeurite);
              maxValueWholeImage = max(maxValueWholeImage,maxValueWhole);
              %if(sixteenBit==0 && maxValueWhole>255)
              %  sixteenBit=1;
              %end
              %Check if in area x to x+512 and y to y+512 is an Nucleus. If
              %not, check for highest value for high pass filter.
              %if(sixteenBit == 1)
              %r=isNucleusInArea(eventdata,handles,currentWellWithoutNull,x+1,x+512,y+1,y+512);
          %Wofür?
              %    if(r==0)

          %         maxFilterList = [maxFilterList maxValueWhole];
          %    end

              
              neuriteImageLine = cat(2,neuriteImageLine,currentNeuritePic);
              nucleusImageLine = cat(2,nucleusImageLine,currentNucleusPic);
              i = i + 1;
              if(i<=maxIndex)
                  x = SortedKeyArray(i,2);
                  y = SortedKeyArray(i,3);
              end
            end
            if(breaked==1)
              break;
            end
            previousY = y;
            %neuriteImage = cat(1,neuriteImage,neuriteImageLine);
            %nucleusImage = cat(1,nucleusImage,nucleusImageLine);
            
            if(get(handles.rbOtSingleChamber16,'Value'))
                neuriteImageLine = flipdim(neuriteImageLine,1);
                nucleusImageLine = flipdim(nucleusImageLine,1);
            end
            neuriteImage = cat(1,neuriteImage,neuriteImageLine);
            nucleusImage = cat(1,nucleusImage,nucleusImageLine);
            end
            if(breaked==1)
                continue;
            end
            
            %If 16 Bit image, invert y Achses on image
    if(get(handles.rb4816,'Value') || get(handles.rb9616,'Value') || get(handles.rbWholeOt16,'Value')) %|| get(handles.rbOtSingleChamber16,'Value'))
        neuriteImage = flipdim(neuriteImage,1);
        nucleusImage = flipdim(nucleusImage,1);
    end
            imageHandler.ResizedNeuriteImage = imresize(neuriteImage,0.1);
            imageHandler.ResizedNucleusImage = imresize(nucleusImage,0.1);    
            imageHandler.NeuriteImage = neuriteImage;
            imageHandler.NucleusImage = nucleusImage;
    
            %For every well: Save Original Pics in Subfolder ConvertedCellomics
            filepath = [subfoldername '/' currentWell 'NeuriteBig.tif'];
            imwrite(neuriteImage,filepath);
            filepath = [subfoldername '/' currentWell 'NucleusBig.tif'];
            imwrite(nucleusImage,filepath);
            filepath = [subfoldername '/' currentWell 'NeuriteSmall.tif'];
            imwrite(imageHandler.ResizedNeuriteImage,filepath);
            filepath = [subfoldername '/' currentWell 'NucleusSmall.tif'];
            imwrite(imageHandler.ResizedNucleusImage,filepath);
        end  
        %Wofür?
        %standardDeviationMaxFilter = std(double(maxFilterList));
        %meanMaxFilter = mean(maxFilterList);
        %maxValue = meanMaxFilter + standardDeviationMaxFilter;
        close(waitbarHandle);
    end
    handles.ImageHandler = imageHandler;
    handles.CSVCoordinates = csvHandler;
    guidata(handles.figure1, handles) 
    
end

function r=isNucleusInArea(hObject, handles, selectedWell, xMin, xMax, yMin, yMax)
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
WellNucleusDict = csvHandler.CellPosMatrix;
NucleusM = WellNucleusDict(selectedWell);
SubNucleusM = NucleusM(yMin:yMax,xMin:xMax);
r=max(SubNucleusM(:));

% --- Executes on mouse press over axes background.
function executeOnImageClick(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
neuronHandler = handles.NeuronCoordinates;
csvHandler = handles.CSVCoordinates;
imageHandler = handles.ImageHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
hPixelInfoChildren = get(imageHandler.MousePosition,'Children');
stringPix = get(hPixelInfoChildren(1),'String');
%Extract X and Y Coordinate from stringPix
tokenMatrix = regexpi(stringPix,'\((\d+),\s(\d+)\)','tokens');
x = tokenMatrix{1}(1);
y = tokenMatrix{1}(2);
x=str2num(x{1});
y=str2num(y{1});
if(imageHandler.CheckSquareMode)
    %Map current Position to big picture and then to regarding Square
    %Mark Squares in 512er Steps
    xOffset = mod(x*10,512);
    yOffset = mod(y*10,512);
    xStart = x*10-xOffset;
    yStart = y*10-yOffset;
    h = imrect(gca,[xStart/10,yStart/10,512/10, 512/10]);
    setResizable(h, 0);
    filter = logical(createMask(h));
    imageHandler.CurrentSquareMask = logical(imageHandler.CurrentSquareMask | filter);    
elseif(neuronHandler.ManualCountingMode)
    manualPositionList = neuronHandler.ManualNeuronPositions(selectedWell);
    %Map current Position to big picture and then to regarding Nucleus
    xMin = imageHandler.ZoomState(1);
    yMin = imageHandler.ZoomState(2);
    xMapped = int32(x + xMin);
    yMapped = int32(y + yMin);
    WellNucleusDict = csvHandler.CellPosMatrix;
    NucleusM = WellNucleusDict(selectedWell);
    [xNuc,yNuc,distance] = csvHandler.FindNucleusForNeuron(xMapped,yMapped,NucleusM);
    disp(distance);
    %Save current Position
    manualPositionList{numel(manualPositionList)+1} = [xNuc yNuc];
    neuronHandler.ManualNeuronPositions(selectedWell) = manualPositionList;
    currentSparse = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
    currentSparse(yNuc,xNuc) = 1;
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = currentSparse;
    %Calculate position back for marker
    markerX = xNuc-xMin;
    markerY = yNuc-yMin;    
    hold on;
    plot(markerX,markerY,'Linestyle','none','Marker','.','Markersize',20);
    hold off;
    handles.NeuronCoordinates = neuronHandler;
    guidata(handles.figure1, handles) 
elseif(neuronHandler.ManualDeletingMode)
    manualPositionList = neuronHandler.ManualNeuronPositions(selectedWell);
    %Map current Position to big picture and then to regarding Nucleus
    xMin = imageHandler.ZoomState(1);
    yMin = imageHandler.ZoomState(2);
    xMapped = int32(x + xMin);
    yMapped = int32(y + yMin);
    WellNucleusDict = csvHandler.CellPosMatrix;
    NucleusM = WellNucleusDict(selectedWell);
    [xNuc,yNuc,distance] = neuronHandler.FindNucleusForNeuron(xMapped,yMapped,NucleusM);
    %Delete current Position, if available.
    %Find index of current Position
    deleteList = cell(0);
    filterEntries = 0;
    manListCopy = cell(0);
    for i=1:numel(manualPositionList)
        manList = manualPositionList(i);
        manList = manList{1};
        if(manList(1) == xNuc && manList(2) == yNuc)
            deleteList{filterEntries+1} = i;
            filterEntries = filterEntries + 1;
        else
            manListCopy(i-filterEntries) = manualPositionList(i);
        end
    end    
    
    neuronHandler.ManualNeuronPositions(selectedWell) = manListCopy;
    currentSparse = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
    currentSparse(yNuc,xNuc) = 0;
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = currentSparse;
    %Calculate position back for marker
    markerX = xNuc-xMin;
    markerY = yNuc-yMin;    
    hold on;
    plot(markerX,markerY,'Linestyle','none','Marker','.','Markersize',20,'Color','red');
    %set(p,'Color','red');
    hold off;
    handles.NeuronCoordinates = neuronHandler;
    guidata(handles.figure1, handles) 
elseif(imageHandler.MeasureMigDistMode)
    imageHandler.MigDistPoints = [imageHandler.MigDistPoints; [x*10 y*10]];
    if(mod(numel(imageHandler.MigDistPoints)/2, 2) == 0)
        set(handles.lbFirstMigDistPoint,'Visible','off');
    else
        set(handles.lbFirstMigDistPoint,'Visible','on');
    end
end


% --- Executes on button press in pbExtractSphere.
function pbExtractSphere_Callback(hObject, eventdata, handles)
% hObject    handle to pbExtractSphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(isfield(handles,'ImageHandler'))
    handles = guidata(handles.figure1);
    imageHandler = handles.ImageHandler;
    [x, y, BW, xi, yi] = roipoly;
    handles.ImageHandler = imageHandler;
    guidata(handles.figure1, handles) 
end

function LoadManualNeuronPositionsToGUI(handles)
handles = guidata(handles.figure1);
neuronHandler = handles.NeuronCoordinates;
csvHandler = handles.CSVCoordinates;
imageHandler = handles.ImageHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
zoomState = imageHandler.ZoomState;
if(isKey(neuronHandler.ManualNeuronPositions,selectedWell) && numel(imageHandler.ZoomState) > 1)
xMin = imageHandler.ZoomState(1);
yMin = imageHandler.ZoomState(2);
%Iterate over all selected Points in current Well

manualPositionList = neuronHandler.ManualNeuronPositions(selectedWell);
for i=1:numel(manualPositionList)
    %Extract X and Y Coordinate
    pos = manualPositionList(i);
    pos = pos{1};
    xNuc = pos(1);
    yNuc = pos(2);
    %Calculate position back for marker
    markerX = xNuc-xMin;
    markerY = yNuc-yMin;    
    hold on;
    plot(markerX,markerY,'Linestyle','none','Marker','.','Markersize',20);
    hold off;
end
end


% --- Executes on button press in pbZoomIn.
function pbZoomIn_Callback(hObject, eventdata, handles)
% hObject    handle to pbZoomIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
if(imageHandler.ZoomState)
    filterList = imageHandler.PositiveFilters;
    filterNumber = get(handles.lbFilter, 'string');
    if(class(filterNumber) == 'cell')
        filterNumber = filterNumber{1};
    end
    filterNumber = str2num(filterNumber);
    % Check, if ZoomState is already set.
    % If yes, make additional zooming available

    zoomOnZoomedImage = getrect;
    xmin = imageHandler.ZoomState(1);
    ymin = imageHandler.ZoomState(2); 
    width = imageHandler.ZoomState(3);
    height = imageHandler.ZoomState(4);
    
    subNucImage = imageHandler.ShowNucleusImage(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
    subNeuImage = imageHandler.ShowNeuriteImage(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
    if(~numel(imageHandler.NucleusImage) == 0)              
            mat = imfuse(subNucImage,subNeuImage);
            hImage = imshow(mat);        
            set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
            hPixelInfo = impixelinfo;        
            imageHandler.MousePosition = hPixelInfo;
    end
    %Set new ZoomState
    xmin = xmin + zoomOnZoomedImage(1);
    ymin = ymin + zoomOnZoomedImage(2);
    width = zoomOnZoomedImage(3);
    height = zoomOnZoomedImage(4);
    imageHandler.ZoomState = [xmin ymin width height];
    handles.ImageHandler = imageHandler;
    guidata(handles.figure1, handles)
    LoadManualNeuronPositionsToGUI(handles);
else
    csvHandler = handles.CSVCoordinates;
    neuronHandler = handles.NeuronCoordinates;
    %optionHandler = handles.OptionHandler;
    zoomOnZoomedImage = getrect;
    imageHandler.ZoomState = [zoomOnZoomedImage(1)*10 zoomOnZoomedImage(2)*10 zoomOnZoomedImage(3)*10 zoomOnZoomedImage(4)*10];
    currentPicNuc = imageHandler.NucleusImage(zoomOnZoomedImage(2)*10:zoomOnZoomedImage(2)*10+zoomOnZoomedImage(4)*10, zoomOnZoomedImage(1)*10:zoomOnZoomedImage(1)*10+zoomOnZoomedImage(3)*10);
    currentPicNeu = imageHandler.NeuriteImage(zoomOnZoomedImage(2)*10:zoomOnZoomedImage(2)*10+zoomOnZoomedImage(4)*10, zoomOnZoomedImage(1)*10:zoomOnZoomedImage(1)*10+zoomOnZoomedImage(3)*10);
    %minBrightness = optionHandler.HistogramMinNeurite;
    %maxBrightness = optionHandler.HistogramMaxNeurite;
    %currentPicNeu = double(currentPicNeu)./3996;
    %currentPicNeu = imadjust(currentPicNeu, [double(minBrightness)/3996;double(maxBrightness)/3996], [0;1]);
    %currentPicNeu = uint8(currentPicNeu.*255);   
    %minBrightness = optionHandler.HistogramMinNucleus;
    %maxBrightness = optionHandler.HistogramMaxNucleus;
    %currentPicNuc = double(currentPicNuc)./3996;
    %currentPicNuc = imadjust(currentPicNuc, [double(minBrightness)/3996;double(maxBrightness)/3996], [0;1]);
    %currentPicNuc = uint8(currentPicNuc.*255); 
    mat = imfuse(currentPicNuc,currentPicNeu);
    hImage = imshow(mat);
    cbManualNeuronMatrixChecked = get(handles.cbManualNeuronMatrix ,'Value');
    cbCellomicsNeuronMatrixChecked = get(handles.cbCellomicsNeuronMatrix ,'Value');
    cbNucleusMatrixChecked = get(handles.cbNucleusMatrix ,'Value');
    %Load selected Matrices to GUI
    if(cbNucleusMatrixChecked)
        WellNucleusDict = csvHandler.CellPosMatrix;
        NucleusM = WellNucleusDict(selectedWell);
        NucleusM = NucleusM(zoomOnZoomedImage(2)*10:zoomOnZoomedImage(2)*10+zoomOnZoomedImage(4)*10, zoomOnZoomedImage(1)*10:zoomOnZoomedImage(1)*10+zoomOnZoomedImage(3)*10);
        [nucleusRows, nucleusCols] = find(NucleusM);   
        hold on;
        plot((nucleusCols),(nucleusRows),'Linestyle','none','Marker','.','Markersize',15,'Color',[1 .5 0]);       
        hold off;
    end
    if(cbCellomicsNeuronMatrixChecked)
        WellNeuronDict = neuronHandler.CellPosMatrix;
        NeuronM = WellNeuronDict(selectedWell); 
        NeuronM = NeuronM(zoomOnZoomedImage(2)*10:zoomOnZoomedImage(2)*10+zoomOnZoomedImage(4)*10, zoomOnZoomedImage(1)*10:zoomOnZoomedImage(1)*10+zoomOnZoomedImage(3)*10);
        [neuronRows, neuronCols] = find(NeuronM);   
        hold on;
        plot((neuronCols),(neuronRows),'Linestyle','none','Marker','.','Markersize',20,'Color','red');       
        hold off;
    end
    if(cbManualNeuronMatrixChecked)
        WellNeuronDict = neuronHandler.ManualNeuronPositionsSparse;
        selectedWellNumber = get(handles.lbWell,'Value');
        NeuronM = WellNeuronDict(selectedWell);
        NeuronM = NeuronM(zoomOnZoomedImage(2)*10:zoomOnZoomedImage(2)*10+zoomOnZoomedImage(4)*10, zoomOnZoomedImage(1)*10:zoomOnZoomedImage(1)*10+zoomOnZoomedImage(3)*10);
        [neuronRows, neuronCols] =  find(NeuronM);    
        hold on;
        plot((neuronCols),(neuronRows),'Linestyle','none','Marker','.','Markersize',20,'Color','green');       
        hold off;
    end
end



% --- Executes when selected object is changed in uipanel6.
function uipanel6_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel6 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
% 
handles = guidata(handles.figure1);
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'rbNeuritePic'
        if(isfield(handles,'ImageHandler'))
            imageHandler = handles.ImageHandler;
            imageHandler.ZoomState = zeros(1);
            if(~numel(imageHandler.NucleusImage) == 0 && ~numel(imageHandler.ResizedNeuriteImage) == 0)       
                hImage = imshow(imageHandler.ResizedNeuriteImage);
                set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
                hPixelInfo = impixelinfo;        
                imageHandler.MousePosition = hPixelInfo;
            end            
        end        
    case 'rbNucleusPic'
        if(isfield(handles,'ImageHandler'))
            imageHandler = handles.ImageHandler;
            imageHandler.ZoomState = zeros(1);
            if(~numel(imageHandler.NucleusImage) == 0 && ~numel(imageHandler.ResizedNeuriteImage) == 0)                     
                hImage = imshow(imageHandler.ResizedNucleusImage);  
                set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
                hPixelInfo = impixelinfo;        
                imageHandler.MousePosition = hPixelInfo;
            end
        end
    case 'rbNucleusNeuritePic'
        if(isfield(handles,'ImageHandler'))
            imageHandler = handles.ImageHandler;
            imageHandler.ZoomState = zeros(1);
            if(~numel(imageHandler.NucleusImage) == 0 && ~numel(imageHandler.ResizedNeuriteImage) == 0)              
                mat = imfuse(imageHandler.ResizedNucleusImage,imageHandler.ResizedNeuriteImage);
                hImage = imshow(mat);        
                set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
                hPixelInfo = impixelinfo;        
                imageHandler.MousePosition = hPixelInfo;
            end
        end
    case 'rbNucleusMatrix'
        if(isfield(handles,'CSVCoordinates'))
            csvHandler = handles.CSVCoordinates;
            WellNucleusDict = csvHandler.CellPosMatrix;
            selectedWellNumber = get(handles.lbWell,'Value');
            wellList = get(handles.lbWell, 'string');
            selectedWell = wellList{selectedWellNumber};
            NucleusM = WellNucleusDict(selectedWell);
            NucleusMGUI = WellNucleusDict(selectedWell);                   
            [nucleusRows, nucleusCols] = find(NucleusM);
            [sizeY, sizeX] = size(NucleusM);          
            hold on;
            plot((nucleusCols./10),(nucleusRows./10),'Linestyle','none','Marker','.','Markersize',2,'Color',[1 .5 0]);       
            hold off;
            %spy(NucleusMGUI);
            hImage = gca;
            set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
        end
    case 'rbNeuronMatrix'
        if(isfield(handles,'NeuronCoordinates'))
            neuronHandler = handles.NeuronCoordinates;
            WellNeuronDict = neuronHandler.CellPosMatrix;
            selectedWellNumber = get(handles.lbWell,'Value');
            wellList = get(handles.lbWell, 'string');
            selectedWell = wellList{selectedWellNumber};
            NeuronM = WellNeuronDict(selectedWell);
            NeuronMGUI = WellNeuronDict(selectedWell);
            [neuronRows, neuronCols] =  find(NeuronM);
            [sizeY, sizeX] = size(NeuronM);
            hold on;
            plot((neuronCols./10),(neuronRows./10),'Linestyle','none','Marker','.','Markersize',10,'Color','red');       
            hold off;
            hImage = gca;
            set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
        end
    case 'rbManualNeuron'
        if(isfield(handles,'NeuronCoordinates'))
            neuronHandler = handles.NeuronCoordinates;
            WellNeuronDict = neuronHandler.ManualNeuronPositionsSparse;
            selectedWellNumber = get(handles.lbWell,'Value');
            wellList = get(handles.lbWell, 'string');
            selectedWell = wellList{selectedWellNumber};
            NeuronM = WellNeuronDict(selectedWell);
            NeuronMGUI = WellNeuronDict(selectedWell);
            [neuronRows, neuronCols] =  find(NeuronM);
            [sizeY, sizeX] = size(NeuronM);
            for i=1:numel(neuronRows)
                row = neuronRows(i);
                col = neuronCols(i);
                if(row > 100 && row+100 < sizeY && col > 100 && col+100 < sizeX)
                    for j=-100;100;
                        for k=-100;100;
                            NeuronMGUI(row+k, col+k) = 1;
                        end
                    end
                end %if
            end %for
            spy(NeuronMGUI);
            hImage = gca;
            set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
        end
    otherwise
        disp('Error. No valid GUI Radiobutton selected');
end


% --- Executes on selection change in lbWell.
function lbWell_Callback(hObject, eventdata, handles)
% hObject    handle to lbWell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
set(handles.lbFilter, 'Value',1);
set(handles.lbFilter, 'String', '');
if(isfield(handles,'ImageHandler'))   
    %Empty filter list
    imageHandler = handles.ImageHandler;
    imageHandler.ZoomState = zeros(1);
    imageHandler.PositiveFilters = 0;
    imageHandler.NegativeFilters = 0;
    imageHandler.CheckSquareMode=0;
    selectedWellNumber = get(handles.lbWell,'Value');
    wellList = get(handles.lbWell, 'string');
    selectedWell = wellList{selectedWellNumber};
    
    if(imageHandler.MeasureMigDistMode==1)
        imageHandler.MeasureMigDistMode=0;
        %Calculate and save measured migration distance for current well
        mediumMigDist=0;
        count = 0;
        for i=1:numel(imageHandler.MigDistPoints)/2
            if(mod(i,2) ~= 0)
                pointOne = [imageHandler.MigDistPoints(i,1) imageHandler.MigDistPoints(i,2)];
                pointTwo = [imageHandler.MigDistPoints(i+1,1) imageHandler.MigDistPoints(i+1,2)];
                mediumMigDist = mediumMigDist + pdist([pointOne;pointTwo]);
                count = count+1;
            end
        end
        mediumMigDist = mediumMigDist / count; 
        imageHandler.MigDistWellMapping(imageHandler.PreviousWell) = mediumMigDist;
        disp(strcat('Manual Migration distance is ',num2str(mediumMigDist),' pixels'));
    end
    
    imageHandler.PreviousWell = selectedWell;
    %Load filters for current well
    path = [imageHandler.Foldername '/ConvertedCellomics/' selectedWell '.mat'];
    if(exist(path,'file'))
        load(path);
        imageHandler.PositiveFilters = FMask.PositiveFilters;
        imageHandler.NegativeFilters = FMask.NegativeFilters;
        filterEntries = numel(imageHandler.PositiveFilters);
        currentFilter = imageHandler.PositiveFilters;

        %For every positive entry in Filter List: Create entry in Listbox
        fulltext = '';
        for i=1:numel(currentFilter)
            fulltext = [fulltext;{i}];
        end
        set(handles.lbFilter, 'String', fulltext);
    end
    if(length(selectedWell) == 2)
      selectedWell = [selectedWell(1) '0' selectedWell(2)];
    end
    
    %Get images from image folder
    foldername = [imageHandler.Foldername '/ConvertedCellomics'];
    imagePathNucleusBig = [foldername '/' selectedWell 'NucleusBig.tif'];
    imagePathNucleusSmall = [foldername '/' selectedWell 'NucleusSmall.tif'];
    imagePathNeuriteBig = [foldername '/' selectedWell 'NeuriteBig.tif'];
    imagePathNeuriteSmall = [foldername '/' selectedWell 'NeuriteSmall.tif'];    
    imageHandler.ResizedNeuriteImage = imread(imagePathNeuriteSmall);
    imageHandler.ResizedNucleusImage = imread(imagePathNucleusSmall);    
    imageHandler.NeuriteImage = imread(imagePathNeuriteBig);
    imageHandler.NucleusImage = imread(imagePathNucleusBig);
    handles.ImageHandler = imageHandler;
    %Get selected GUI value
    guidata(handles.figure1, handles)   
    RefreshUnzoomedImage(handles);    
end
guidata(handles.figure1, handles)         


% --- Executes on button press in pbSetFilter.
function pbSetFilter_Callback(hObject, eventdata, handles)
% hObject    handle to pbSetFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
path = [imageHandler.Foldername '/ConvertedCellomics/' selectedWell '.mat'];
if(exist(path,'file'))
    load(path);
    imageHandler.PositiveFilters = FMask.PositiveFilters;
    imageHandler.NegativeFilters = FMask.NegativeFilters;
    FMask=0;
end
if(iscell(imageHandler.PositiveFilters) > 0)
  filterEntries = numel(imageHandler.PositiveFilters);
  currentFilter = imageHandler.PositiveFilters;
else
  filterEntries = 0;
  currentFilter = cell(0);
end
filterEntriesStart = filterEntries;
if(get(handles.rbRectangleFilter,'Value'))
    [rowNumber, colNumber]=size(imageHandler.NeuriteImage);
    A = get(handles.tbA,'String');
    B = get(handles.tbB,'String');
    if(imageHandler.ZoomState)
        h = imrect(gca,[1,1,(str2num(A)), str2num(B)]);
    else
        h = imrect(gca,[1,1,(str2num(A)/10),(str2num(B)/10)]);
    end
    setResizable(h, 0);
    A2 = get(handles.tbA,'String');
    B2 = get(handles.tbB,'String');
    if(imageHandler.ZoomState)
        h2 = imrect(gca,[1,1,(str2num(B)), str2num(A)])
    else
        h2 = imrect(gca,[1,1,(str2num(B)/10),(str2num(A)/10)]);    
    end
    setResizable(h2, 0);
    posH = wait(h);
    posH2 = wait(h2);
    if(imageHandler.ZoomState)
        xmin = imageHandler.ZoomState(1);
        ymin = imageHandler.ZoomState(2); 
        polyPos = getPosition(h);
        polyPos2 = getPosition(h2);
        %Wandle alle Positionen in Position für das große Bild um        
        currentFilter{filterEntries+1} = poly2mask(polyPos(:,1)+xmin,polyPos(:,2)+ymin,colNumber, rowNumber);
        currentFilter{filterEntries+2} = poly2mask(polyPos2(:,1)+xmin,polyPos2(:,2)+ymin,colNumber, rowNumber);
    else
        currentFilter{filterEntries+1} = createMask(h);
        currentFilter{filterEntries+2} = createMask(h2);
        currentFilter{filterEntries+1} = imresize(currentFilter{filterEntries+1},[rowNumber, colNumber]);
        currentFilter{filterEntries+2} = imresize(currentFilter{filterEntries+2},[rowNumber, colNumber]);
    end
    lbFilterString = get(handles.lbFilter,'String');
    %For every positive entry in Filter List: Create entry in Listbox
    filterEntries = filterEntries+1;
    lbFilterString = [lbFilterString;num2str(filterEntries)];
    imageHandler.PositiveFilters = currentFilter;
    set(handles.lbFilter, 'String', lbFilterString);
    filterEntries = filterEntries+1;
elseif(get(handles.rbPolygon,'Value'))
    h = impoly(gca);
    [rowNumber, colNumber]=size(imageHandler.NeuriteImage);
    if(imageHandler.ZoomState)
        xmin = imageHandler.ZoomState(1);
        ymin = imageHandler.ZoomState(2);  
        polyPos = getPosition(h);
        %Wandle alle Positionen in Position für das große Bild um        
        mask = poly2mask(polyPos(:,1)+xmin,polyPos(:,2)+ymin,colNumber, rowNumber);
    else        
        mask = createMask(h);
        %Map small mask to big picture
        mask = imresize(mask,[rowNumber, colNumber]);        
    end
    currentFilter{filterEntries+1} = mask;
    imageHandler.PositiveFilters = currentFilter;
    filterEntries = filterEntries+1;
elseif(get(handles.rbSquares,'Value'))
    %Set flag for CheckSquareMode
    if(imageHandler.CheckSquareMode~=1)
        [ySize xSize] = size(imageHandler.ResizedNucleusImage);
        imageHandler.CheckSquareMode=1;
        imageHandler.CurrentSquareMask = zeros(ySize,xSize);
    else
        [rowNumber, colNumber]=size(imageHandler.NeuriteImage);
        imageHandler.CheckSquareMode=0;
        imageHandler.CurrentSquareMask = imresize(imageHandler.CurrentSquareMask,[rowNumber, colNumber]);
        currentFilter{filterEntries+1} = imageHandler.CurrentSquareMask;
        imageHandler.CurrentSquareMask=0;
        imageHandler.PositiveFilters = currentFilter;
        filterEntries = filterEntries+1;
    end
elseif(get(handles.rbAutoFilter,'Value'))
        waitbarHandle = waitbar(0,'Set Filter for Sphere automatically.');
        selectedWellNumber = get(handles.lbWell,'Value');
        wellList = get(handles.lbWell, 'string');
        selectedWell = wellList{selectedWellNumber};
        waitbar((20/100),waitbarHandle);
        sizingFactorY = 64;
        sizingFactorX=64;
        
        [filterDistance nonFilterDistance SphereArea NucleusArea] = calculateMigrationDistance(selectedWell,sizingFactorX,sizingFactorY, handles);
        waitbar((60/100),waitbarHandle);
        [rowNumber, colNumber]=size(imageHandler.NeuriteImage);
        sizingFactorX2 = colNumber/sizingFactorX;
        sizingFactorY2 = rowNumber/sizingFactorY;     
        %Calculate Factor for preview Pic
        sizingFactorX2 = sizingFactorX2 / 10;
        sizingFactorY2 = sizingFactorY2 / 10;
        allowedDistanceNucleus = 3*((sizingFactorX2 + sizingFactorY2)/2);
        allowedDistanceSphere = 10*((sizingFactorX2 + sizingFactorY2)/2);
        rowColMatSorted = sortMatrixByNearestNeighbour(SphereArea,allowedDistanceSphere,sizingFactorX,sizingFactorY,handles);
        waitbar((70/100),waitbarHandle);
        %negativeRowColMatSorted = sortMatrixByNearestNeighbour(NucleusArea,allowedDistanceNucleus,sizingFactorX,sizingFactorY,handles);
        newNucleusMask = imresize(NucleusArea,[rowNumber,colNumber]);
        waitbar((80/100),waitbarHandle);
        h = impoly(handles.axes2,rowColMatSorted);
        cbNegativeFilterChecked = get(handles.cbNegativeFilter ,'Value');
        if(cbNegativeFilterChecked)
            %h2 = impoly(handles.axes2,negativeRowColMatSorted);            
            %setColor(h2,'red');
            %maskNegative = createMask(h2);
            %maskNegative = imresize(maskNegative,[rowNumber, colNumber]);
            if(iscell(imageHandler.NegativeFilters) > 0)
                negFilterEntries = numel(imageHandler.NegativeFilters);
                currentNegativeFilter = imageHandler.NegativeFilters;
            else
                negFilterEntries = 0;
                currentNegativeFilter = cell(0);
            end
            currentNegativeFilter{negFilterEntries+1} = newNucleusMask;
            imageHandler.NegativeFilters = currentNegativeFilter;
            negFilterEntries = negFilterEntries + 1;
        end
        mask = createMask(h);
        %Map small mask to big picture
        mask = imresize(mask,[rowNumber, colNumber]);
        currentFilter{filterEntries+1} = mask;
        imageHandler.PositiveFilters = currentFilter;
        filterEntries = filterEntries+1;
        waitbar((90/100),waitbarHandle);
        close(waitbarHandle);
end
if(filterEntries>0 && filterEntries ~= filterEntriesStart)
    %Write all Filters to ListBox
    lbFilterString = get(handles.lbFilter,'String');
    %For every positive entry in Filter List: Create entry in Listbox

    lbFilterString = [lbFilterString;num2str(filterEntries)];
    set(handles.lbFilter, 'String', lbFilterString);

    %Save to file
    foldername = imageHandler.Foldername;
    subfoldername = [foldername '/ConvertedCellomics'];
    FMask = FilterMask();
    FMask.PositiveFilters = imageHandler.PositiveFilters;
    FMask.NegativeFilters = imageHandler.NegativeFilters;
    save('-v7.3',strcat(subfoldername,'\',selectedWell),'FMask');
    FMask=0;
    imageHandler.PositiveFilters = 0;
    imageHandler.NegativeFilters = 0;
end
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles)


function rowColMatSorted = sortMatrixByNearestNeighbour(SphereArea,allowedMaxDistance,AreaSizeX, AreaSizeY,handles)
imageHandler=handles.ImageHandler;
[NonZeroRows NonZeroCols] = find(SphereArea);
[rowNumber, colNumber]=size(imageHandler.NeuriteImage);
sizingFactorX = colNumber/AreaSizeX;
sizingFactorY = rowNumber/AreaSizeY;
%Divide by 10 for Preview Pic
sizingFactorX = sizingFactorX / 10;
sizingFactorY = sizingFactorY / 10;
NonZeroCols = (NonZeroCols .* sizingFactorX);
NonZeroRows = (NonZeroRows .* sizingFactorY);
CalculatedVector = zeros(size(NonZeroRows));
%Sort Cols and rows by nearest neighbour
        
rowColMat = ([NonZeroCols NonZeroRows CalculatedVector]);
rowColMatSorted = zeros(0,0);
i=1;
while(i~=0)
    minDist=9999999;
    nearestJ=0;
    for j=1:size(NonZeroRows)
       %Wenn i und j unterschiedlich sind, und j noch nicht
       %bearbeitet wurde, berechne distanz.
       if(i~=j && rowColMat(j,3) == 0)
           dist = pdist([rowColMat(j,1) rowColMat(j,2); rowColMat(i,1) rowColMat(i,2)]);
           if (dist<minDist && dist<=allowedMaxDistance)
              minDist=dist;
              nearestJ=j;
           end
        end
    end
    if(nearestJ~=0)
        rowColMatSorted = [rowColMatSorted;NonZeroCols(nearestJ) NonZeroRows(nearestJ)];
        rowColMat(nearestJ,3) = 1;
        i=nearestJ;
    else
        i=0;
    end
end

% --- Executes on button press in pbExcludePoly.
function pbExcludePoly_Callback(hObject, eventdata, handles)
% hObject    handle to pbExcludePoly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};

if(iscell(imageHandler.NegativeFilters) > 0)
  filterEntries = numel(imageHandler.NegativeFilters);
  currentFilter = imageHandler.NegativeFilters;
else
  filterEntries = 0;
  currentFilter = cell(0);
end
if(get(handles.rbPolygon,'Value'))
    h = impoly(gca);
    setColor(h,'red');
    [rowNumber, colNumber]=size(imageHandler.NeuriteImage);
    if(imageHandler.ZoomState)        
        xmin = imageHandler.ZoomState(1);
        ymin = imageHandler.ZoomState(2);  
        polyPos = getPosition(h);
        %Wandle alle Positionen in Position für das große Bild um        
        mask = poly2mask(polyPos(:,1)+xmin,polyPos(:,2)+ymin,colNumber, rowNumber);
    else
        mask = createMask(h);
        %Map small mask to big picture
        mask = imresize(mask,[rowNumber, colNumber]);        
    end
    currentFilter{filterEntries+1} = mask;
    imageHandler.NegativeFilters = currentFilter;
    filterEntries = filterEntries + 1;
end

%Save to file
foldername = imageHandler.Foldername;
subfoldername = [foldername '/ConvertedCellomics'];
FMask = FilterMask();
FMask.PositiveFilters = imageHandler.PositiveFilters;
FMask.NegativeFilters = imageHandler.NegativeFilters;
save('-v7.3',strcat(subfoldername,'\',selectedWell),'FMask');
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pbZoomOut.
function pbZoomOut_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pbZoomOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function tbA_Callback(hObject, eventdata, handles)
% hObject    handle to tbA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tbA as text
%        str2double(get(hObject,'String')) returns contents of tbA as a double


% --- Executes during object creation, after setting all properties.
function tbA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tbA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tbB_Callback(hObject, eventdata, handles)
% hObject    handle to tbB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tbB as text
%        str2double(get(hObject,'String')) returns contents of tbB as a double


% --- Executes during object creation, after setting all properties.
function tbB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tbB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function NeuronCount_Callback(hObject, eventdata, handles)
% hObject    handle to NeuronCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Extract Rectangles from ImageHandler. 
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
filterList = get(handles.lbFilter, 'string');
filterNumber = get(handles.lbFilter,'Value');
%Iterate over all filters
filterList = imageHandler.PositiveFilters;
filterMask = filterList(filterNumber);
filterList = 0;
filterMask = filterMask{1};
%Remove excluded parts of picture
if(iscell(imageHandler.NegativeFilters) > 0)
    negativeMask = imageHandler.NegativeFilters;
    for i=1:numel(negativeMask)
      workMask = negativeMask{i};
      workMask = int8(workMask -1);
      workMask = logical(workMask .* (-1));
      workMask = logical(workMask);
      filterMask = logical((filterMask .* workMask));
    end
end
workMask = 0;
negativeMask = 0;
%Multiply filterMask with Neuron and Nucleus Matrices
NucleusM = csvHandler.CellPosMatrix(selectedWell);
NeuronM = neuronHandler.CellPosMatrix(selectedWell);
NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
[filterSizeY filterSizeX] = size(filterMask);
[NucleusMSizeY NucleusMSizeX] = size(NucleusM);
[NeuronMSizeY NeuronMSizeX] = size(NeuronM);
[NeuronManualMSizeY NeuronManualMSizeX] = size(NeuronManualM);
if(NucleusMSizeX > filterSizeX)
     NucleusM(:,filterSizeX+1:NucleusMSizeX) = [];
     csvHandler.CellPosMatrix(selectedWell) = NucleusM;
end
if(NeuronMSizeX > filterSizeX)
     NeuronM(:,filterSizeX+1:NeuronMSizeX) = [];
     neuronHandler.CellPosMatrix(selectedWell) = NeuronM;
end
if(NeuronManualMSizeX > filterSizeX)
     NeuronManualM(:,filterSizeX+1:NeuronManualMSizeX) = [];
     neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
end
if(NucleusMSizeY > filterSizeY)
     NucleusM(filterSizeY+1:NucleusMSizeY,:) = [];
     csvHandler.CellPosMatrix(selectedWell) = NucleusM;
end
if(NeuronMSizeY > filterSizeY)
     NeuronM(filterSizeY+1:NeuronMSizeY,:) = [];
     neuronHandler.CellPosMatrix(selectedWell) = NeuronM;
end
if(NeuronManualMSizeY > filterSizeY)
    NeuronManualM(filterSizeY+1:NeuronManualMSizeY,:) = [];
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
end

NucleusM = NucleusM .* filterMask;
numberNuclei = nnz(NucleusM);
NucleusM = 0;
NeuronM = NeuronM .* filterMask;   
numberNeurons = nnz(NeuronM);
NeuronM=0;
NeuronManualM = NeuronManualM .* filterMask;
numberManual = nnz(NeuronManualM);
%Add entry in list
table = get(handles.tblData,'Data');
k = size(table,1)+1;
if k == 1 % Data is converted the first time, but only needs one
    B = num2cell(table);
else
    B = table;
end
B{k,1} = selectedWell;
B{k,2} = filterNumber;
B{k,3} =numberNuclei;
B{k,4} =numberNeurons;
B{k,5} =numberManual;
set(handles.tblData,'Data',B)
handles.NeuronCoordinates = neuronHandler;
handles.CSVCoordinates = csvHandler;
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function DeleteNeuronsManual_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteNeuronsManual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
%If not available, create List with manual positions for current Well
if(isKey(neuronHandler.ManualNeuronPositions, selectedWell))
    manualPositionList = neuronHandler.ManualNeuronPositions(selectedWell);
else
    manualPositionList = cell(0);
end
if(neuronHandler.ManualDeletingMode)
    neuronHandler.ManualDeletingMode = 0;
    neuronHandler.ManualCountingMode = 1;
else
    neuronHandler.ManualDeletingMode = 1;
    neuronHandler.ManualCountingMode = 0;
end
neuronHandler.ManualNeuronPositions(selectedWell) = manualPositionList;
handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles)

% --------------------------------------------------------------------
function CountNeuronsManual_Callback(hObject, eventdata, handles)
% hObject    handle to CountNeuronsManual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
%If not available, create List with manual positions for current Well
if(isKey(neuronHandler.ManualNeuronPositions, selectedWell))
    manualPositionList = neuronHandler.ManualNeuronPositions(selectedWell);
else
    manualPositionList = cell(0);
end
if(neuronHandler.ManualCountingMode)
    neuronHandler.ManualCountingMode = 0;
else
    neuronHandler.ManualCountingMode = 1;
end
neuronHandler.ManualNeuronPositions(selectedWell) = manualPositionList;
handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles)

% --- Executes on selection change in lbFilter.
function lbFilter_Callback(hObject, eventdata, handles)
% hObject    handle to lbFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% If available: Apply Filter mask on combined image.
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
filterNumber = get(hObject,'Value');
path = [imageHandler.Foldername '/ConvertedCellomics/' selectedWell '.mat'];
if(exist(path,'file'))
        load(path);
        imageHandler.PositiveFilters = FMask.PositiveFilters;
        imageHandler.NegativeFilters = FMask.NegativeFilters;
end
FMask=0; 
    %Get Filter
    filterList = imageHandler.PositiveFilters;
    filterHandle = filterList(filterNumber);
    filterHandle = filterHandle{1};
    filterMask = filterHandle;
    %Remove excluded parts of picture
    if(iscell(imageHandler.NegativeFilters) > 0)
        negativeMask = imageHandler.NegativeFilters;
        for i=1:numel(negativeMask)
          workMask = negativeMask{i};
          workMask = int8(workMask -1);
          workMask = logical(workMask .* (-1));
          workMask = logical(workMask);
          filterMask = logical((filterMask .* workMask));
        end
    end
    workMask=0;
    %resize filter mask for pic
    filterMaskSmall = imresize(filterMask,0.1);
    %Multiply mask with image
    nucImage = uint8(imageHandler.NucleusImage).*uint8(filterMask);
    neuImage = uint8(imageHandler.NeuriteImage).*uint8(filterMask);
    filterMask=0;    
    %Transpose only resized image, because Transposing is very cost
    %intensive
    nucImageResized = double(imageHandler.ResizedNucleusImage).*double(filterMaskSmall);
    nucResizedTransposed = transpose(nucImageResized);
    %nucImageTransposed = transpose(nucImage);
    %Cut Matrix and set new Zoom State
    [maxPossibleRow, maxPossibleCol] = size(nucImage);
    [minRow, minCol] = find(nucImage,1,'first');
    [minColTrans, minRowTrans] = find(nucResizedTransposed,1,'first');
    minRowTrans = minRowTrans * 10;
    minColTrans = minColTrans * 10;
    minRow = min(minRow, minRowTrans);
    minCol = min(minCol, minColTrans);
    [maxRow, maxCol] = find(nucImage,1,'last');
    [maxColTrans, maxRowTrans] = find(nucResizedTransposed,1,'last');
    
    maxColTrans = maxColTrans * 10;
    
    maxRowTrans = maxRowTrans * 10;
    
    maxRow = max(maxRow,maxRowTrans);
    maxCol = max(maxCol,maxColTrans);
    maxRow = min(maxRow, maxPossibleRow);
    maxCol = min(maxCol, maxPossibleCol);
    %Set correct zoom State
    imageHandler.ShowNucleusImage = nucImage(minRow:maxRow,minCol:maxCol);
    imageHandler.ShowNeuriteImage = neuImage(minRow:maxRow,minCol:maxCol);
    nucImage=0;
    neuImage=0;
    if(~numel(imageHandler.NucleusImage) == 0 && ~numel(imageHandler.ResizedNeuriteImage) == 0)              
            mat = imfuse(imageHandler.ShowNucleusImage,imageHandler.ShowNeuriteImage);
            hImage = imshow(mat);        
            set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
            hPixelInfo = impixelinfo;        
            imageHandler.MousePosition = hPixelInfo;
    end
    imageHandler.ZoomState = [minCol minRow maxCol-minCol maxRow-minRow];
    imageHandler.PositiveFilters = 0;
    imageHandler.NegativeFilters = 0;
    FMask = 0;
    handles.ImageHandler = imageHandler;
    guidata(handles.figure1, handles);
    LoadManualNeuronPositionsToGUI(handles);
    guidata(handles.figure1, handles);


% --- Executes during object creation, after setting all properties.
function lbFilter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function SaveFixedNeuronPositions_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFixedNeuronPositions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Iterate over Neuron Coordinates and Save mapped data to CSV
handles = guidata(handles.figure1);
imageHandler = 0;
neuronHandler = 0;
csvHandler = 0;
optionHandler = 0;
if(isfield(handles,'ImageHandler'))
    imageHandler = handles.ImageHandler;
    %Set not needed Fields of imageHandler to 0, to reduce amount of needed
    %memory and prevent Out of Memory exception
    neuriteImage = imageHandler.NeuriteImage;
    nucleusImage = imageHandler.NucleusImage;
    resizedNeuriteImage = imageHandler.ResizedNeuriteImage;
    resizedNucleusImage = imageHandler.ResizedNucleusImage;
    if(isprop(imageHandler,'ShowNucleusImage'))
        showNucleusImage = imageHandler.ShowNucleusImage;
        showNeuriteImage = imageHandler.ShowNeuriteImage;
    else
        showNucleusImage=0;
        showNeuriteImage=0;
    end
    imageHandler.NeuriteImage=0;
    imageHandler.NucleusImage=0;
    imageHandler.ResizedNeuriteImage=0;
    imageHandler.ResizedNucleusImage=0;
    imageHandler.ShowNucleusImage=0;
    imageHandler.ShowNeuriteImage=0;
    imageHandler.NeuritePicArray = 0;
    imageHandler.NucleusPicArray = 0;
end
if(isfield(handles,'OptionHandler'))
    optionHandler = handles.OptionHandler;
end
if(isfield(handles,'NeuronCoordinates'))
    neuronHandler = handles.NeuronCoordinates;
end
if(isfield(handles,'CSVCoordinates'))
    csvHandler = handles.CSVCoordinates;
end
[FileName,PathName] = uiputfile;
save('-v7.3',strcat(PathName,FileName),'imageHandler','neuronHandler','csvHandler','optionHandler');
if(isfield(handles,'ImageHandler'))
    imageHandler.NeuriteImage = neuriteImage;
    imageHandler.NucleusImage = nucleusImage;
    imageHandler.ResizedNeuriteImage = resizedNeuriteImage;
    imageHandler.ResizedNucleusImage = resizedNucleusImage;
    imageHandler.ShowNucleusImage = showNucleusImage;
    imageHandler.ShowNeuriteImage = showNeuriteImage;
end
disp('Saved');


% --------------------------------------------------------------------
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
[filename, pathname] = uigetfile('*.mat', 'Select MAT File');
load(strcat(pathname,filename));
%New: One file per well for all Filters
%Load and save on every lbWellClick

%Make split at this place
if(exist('imageHandler') && imageHandler ~= 0)
  s = size(imageHandler.PositiveFilters);
  s = s(1);
  if(s > 1)
    keyList = imageHandler.PositiveFilters.keys;
    foldername = imageHandler.Foldername;
    subfoldername = [foldername '/ConvertedCellomics'];
    for i=1:numel(keyList)
        currentWell = keyList{i};
        %Save all filters to files      
        FMask = FilterMask();
        FMask.PositiveFilters = imageHandler.PositiveFilters(currentWell);
        FMask.NegativeFilters = imageHandler.NegativeFilters(currentWell);
        %Save File for filter
        save('-v7.3',strcat(subfoldername,'\',currentWell),'FMask');
        %Reset old Filter
        imageHandler.PositiveFilters(currentWell)=0;
        imageHandler.NegativeFilters(currentWell)=0;
    end
  end
  imageHandler.PositiveFilters=0;
  imageHandler.NegativeFilters=0;
  imageHandler.MeasureMigDistMode=0;
  imageHandler.MigDistWellMapping = containers.Map();
  handles.ImageHandler = imageHandler;
end
if(exist('neuronHandler'))
    neuronHandler.ManualCountingMode=0;
    handles.NeuronCoordinates = neuronHandler;
end
if(exist('csvHandler'))
    handles.CSVCoordinates = csvHandler;
    guidata(handles.figure1, handles)
    ReloadWellList(handles);
end
if(exist('optionHandler'))
    handles.OptionHandler = optionHandler;
end
guidata(handles.figure1, handles)



function edMinBrightness_Callback(hObject, eventdata, handles)
% hObject    handle to edMinBrightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edMinBrightness as text
%        str2double(get(hObject,'String')) returns contents of edMinBrightness as a double


% --- Executes during object creation, after setting all properties.
function edMinBrightness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edMinBrightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edMaxBrightness_Callback(hObject, eventdata, handles)
% hObject    handle to edMaxBrightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edMaxBrightness as text
%        str2double(get(hObject,'String')) returns contents of edMaxBrightness as a double


% --- Executes during object creation, after setting all properties.
function edMaxBrightness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edMaxBrightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function ConvertPicturesTo8Bit_Callback(hObject, eventdata, handles)
% hObject    handle to ConvertPicturesTo8Bit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Iterate over all pictures in folder
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
foldername = [imageHandler.Foldername '/ConvertedCellomics'];
allFiles = dir(foldername); 
waitbarHandle = waitbar(0,'Please wait. All images are converted to 8 Bit.');
for i=3:numel(allFiles)
    waitbar((i/numel(allFiles)),waitbarHandle);
    currentFile=allFiles(i).name;
    if(strfind(currentFile, 'Neurite'))
        minBrightness = optionHandler.HistogramMinNeurite;
        maxBrightness = optionHandler.HistogramMaxNeurite;
    else
        minBrightness = optionHandler.HistogramMinNucleus;
        maxBrightness = optionHandler.HistogramMaxNucleus;
    end
    imageString = strcat(foldername, '\', currentFile);
    if(length(findstr(imageString, '.mat')) == 0)
      currentPic = imread(imageString);
      deleteIndices = uint16(find(abs(currentPic) < minBrightness));
      currentPic(deleteIndices) = 0;
    
      currentPic = double(currentPic)./3996;
      currentPic = imadjust(currentPic, [double(minBrightness)/3996;double(maxBrightness)/3996], [0;1]);
      currentPic = uint8(currentPic.*255);    
      imwrite(currentPic,imageString);
    end
end
close(waitbarHandle);



function edMinBrightnessNuc_Callback(hObject, eventdata, handles)
% hObject    handle to edMinBrightnessNuc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edMinBrightnessNuc as text
%        str2double(get(hObject,'String')) returns contents of edMinBrightnessNuc as a double


% --- Executes during object creation, after setting all properties.
function edMinBrightnessNuc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edMinBrightnessNuc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edMaxBrightnessNuc_Callback(hObject, eventdata, handles)
% hObject    handle to edMaxBrightnessNuc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edMaxBrightnessNuc as text
%        str2double(get(hObject,'String')) returns contents of edMaxBrightnessNuc as a double


% --- Executes during object creation, after setting all properties.
function edMaxBrightnessNuc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edMaxBrightnessNuc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb8bitPreview.
function pb8bitPreview_Callback(hObject, eventdata, handles)
% hObject    handle to pb8bitPreview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
minBrightness = optionHandler.HistogramMinNeurite;
maxBrightness = optionHandler.HistogramMaxNeurite;
currentPicNeu = imageHandler.ResizedNeuriteImage;
currentPicNeu = double(currentPicNeu)./3996;
currentPicNeu = imadjust(currentPicNeu, [double(minBrightness)/3996;double(maxBrightness)/3996], [0;1]);
currentPicNeu = uint8(currentPicNeu.*255);   
minBrightness = optionHandler.HistogramMinNucleus;
maxBrightness = optionHandler.HistogramMaxNucleus;
currentPicNuc = imageHandler.ResizedNucleusImage;
currentPicNuc = double(currentPicNuc)./3996;
currentPicNuc = imadjust(currentPicNuc, [double(minBrightness)/3996;double(maxBrightness)/3996], [0;1]);
currentPicNuc = uint8(currentPicNuc.*255); 
imageHandler.ShowNucleusImage = currentPicNuc;
imageHandler.ShowNeuriteImage = currentPicNeu;
mat = imfuse(currentPicNuc,currentPicNeu);
hImage = imshow(mat);
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles)


% --------------------------------------------------------------------
function NucleiArea_Callback(hObject, eventdata, handles)
% hObject    handle to NucleiArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
imageHandler = handles.ImageHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
[result withoutFilterResult SphereArea NucleusArea]=calculateMigrationDistance(selectedWell,64,64,handles);
hPixelInfo = impixelinfo;        
imageHandler.MousePosition = hPixelInfo;
msgbox(['Migration distance is ' result ' pixel. Without filter it would be ' withoutFilterResult ' pixel.']); 


function result = FloodFillIterative(DensityM,AreaResult,threshold,above,maxX,maxY,stack,image)
[sizeY, sizeX] = size(image);
while(stack.isEmpty() == 0)
    vector = stack.pop();
    x=vector(1);
    y=vector(2);
    if(above==0)
        [checkX, checkY]=MapPoint(x,y,sizeX,sizeY,64,64,1);
        if(DensityM(y,x) <= threshold && AreaResult(y,x) == 0 && x-1 > 0 && y-1 > 0 && y+1<= maxY && x+1 <= maxX && image(checkY,checkX) > 250)
            AreaResult(y,x) = 1;
            stack.push([x y+1]);
            stack.push([x-1 y]);
            stack.push([x y-1]);
            stack.push([x+1 y]);
            %8 Neighbourhood
            stack.push([x+1 y+1]);
            stack.push([x-1 y-1]);
            stack.push([x-1 y+1]);
            stack.push([x+1 y-1]);
        elseif(DensityM(y,x) <= threshold && AreaResult(y,x) == 0)
            AreaResult(y,x)=1;
        end
    else
        if(DensityM(y,x) >= threshold && AreaResult(y,x) == 0 && x-1 > 0 && y-1 > 0 && y+1<= maxY && x+1 <= maxX)
            AreaResult(y,x) = 1;
            stack.push([x y+1]);
            stack.push([x-1 y]);
            stack.push([x y-1]);
            stack.push([x+1 y]);
            %8 Neighbourhood
            stack.push([x+1 y+1]);
            stack.push([x-1 y-1]);
            stack.push([x-1 y+1]);
            stack.push([x+1 y-1]);
        elseif(DensityM(y,x) >= threshold && AreaResult(y,x) == 0)
            AreaResult(y,x)=1;
        end
    end
end
result=AreaResult;

function [x y] = MapPoint(xOriginal,yOriginal,sizeXOld,sizeYOld,sizeXNew,sizeYNew, calculateOffset)
xFactor = sizeXOld/sizeXNew;
yFactor=sizeYOld/sizeYNew;
if(calculateOffset==1)
    x=uint16((double(xOriginal)*xFactor)-xFactor/2);
    y=uint16((double(yOriginal)*yFactor)-yFactor/2);
else
    x=uint16((double(xOriginal)*xFactor));
    y=uint16((double(yOriginal)*yFactor));
end

function RefreshUnzoomedImage(handles)
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
imageHandler = handles.ImageHandler;
neuronHandler = handles.NeuronCoordinates;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
%Check which checkboxes are checked.
cbNucleusPictureChecked = get(handles.cbNucleusPicture ,'Value');
cbNeuritePictureChecked = get(handles.cbNeuritePicture ,'Value');
cbNucleusMatrixChecked = get(handles.cbNucleusMatrix ,'Value');
cbCellomicsNeuronMatrixChecked = get(handles.cbCellomicsNeuronMatrix ,'Value');
cbManualNeuronMatrixChecked = get(handles.cbManualNeuronMatrix ,'Value');
imageHandler.ZoomState = zeros(1);
[sizeY, sizeX] = size(imageHandler.ResizedNeuriteImage);
m = zeros(sizeY,sizeX);
imshow(m);
if(cbNucleusPictureChecked && cbNeuritePictureChecked)    
    if(~numel(imageHandler.NucleusImage) == 0 && ~numel(imageHandler.ResizedNeuriteImage) == 0)              
         mat = imfuse(imageHandler.ResizedNucleusImage,imageHandler.ResizedNeuriteImage);
         hImage = imshow(mat);  
         set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});         
         hPixelInfo = impixelinfo;        
         imageHandler.MousePosition = hPixelInfo;
    end
elseif(cbNucleusPictureChecked)
    if(~numel(imageHandler.NucleusImage) == 0 && ~numel(imageHandler.ResizedNeuriteImage) == 0)              
         hImage = imshow(imageHandler.ResizedNucleusImage);
         set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
         hPixelInfo = impixelinfo;        
         imageHandler.MousePosition = hPixelInfo;
    end
elseif(cbNeuritePictureChecked)
    if(~numel(imageHandler.NucleusImage) == 0 && ~numel(imageHandler.ResizedNeuriteImage) == 0)              
         hImage = imshow(imageHandler.ResizedNeuriteImage);      
         set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
         hPixelInfo = impixelinfo;        
         imageHandler.MousePosition = hPixelInfo;
    end
else
    %Check if one of the Pictures is selected. If not prepare plot!
end
if(cbNucleusMatrixChecked)
    WellNucleusDict = csvHandler.CellPosMatrix;
    NucleusM = WellNucleusDict(selectedWell);               
    [nucleusRows, nucleusCols] = find(NucleusM);   
    hold on;
    plot((nucleusCols./10),(nucleusRows./10),'Linestyle','none','Marker','.','Markersize',4.5000001,'Color', [1 .5 0]);       
    hold off;
end
if(cbCellomicsNeuronMatrixChecked)
    WellNeuronDict = neuronHandler.CellPosMatrix;
    NeuronM = WellNeuronDict(selectedWell);               
    [neuronRows, neuronCols] = find(NeuronM);   
    hold on;
    plot((neuronCols./10),(neuronRows./10),'Linestyle','none','Marker','.','Markersize',10,'Color','red');       
    hold off;
end
if(cbManualNeuronMatrixChecked)
    WellNeuronDict = neuronHandler.ManualNeuronPositionsSparse;
    selectedWellNumber = get(handles.lbWell,'Value');
    NeuronM = WellNeuronDict(selectedWell);
    [neuronRows, neuronCols] =  find(NeuronM);    
    hold on;
    plot((neuronCols./10),(neuronRows./10),'Linestyle','none','Marker','.','Markersize',10,'Color','green');       
    hold off;
end
cbRings = get(handles.cbPolys ,'Value');
if(cbRings)
    path = [imageHandler.Foldername '/ConvertedCellomics/MarkerPointCoordinates-' selectedWell '.mat'];
    if(exist(path,'file'))
        optionHandler = handles.OptionHandler;
        ringNumber = optionHandler.DensityDistributionRingNumber;
        load(path);
        hInner = impoly(handles.axes2,double(markerPointCoordinates('0')./10));
        previousMask = createMask(hInner);
        for i=1:ringNumber  
          hCurrent = impoly(handles.axes2,double(markerPointCoordinates(num2str(i*10))./10));
        end
    end
end

function [DensityM] = getDensityDistributionsFromCSV(NucleusM,sizeXTarget,sizeYTarget, handles)
%csvHandler = handles.CSVCoordinates;
%neuronHandler = handles.NeuronCoordinates;
%NucleusM = csvHandler.CellPosMatrix;
%NeuronM = neuronHandler.CellPosMatrix;
%selectedWellLong=selectedWell;
%if(length(selectedWell) == 2)
%      selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
%end
%if(isKey(NucleusM,selectedWell))
%    NucleusM = NucleusM(selectedWell);
%    NucleusM = logical(full(NucleusM));
%    NeuronM = NeuronM(selectedWell);
%    NeuronM = logical(full(NeuronM));
    [sizeY, sizeX] = size(NucleusM);
    [row col] = find(NucleusM);
    %[neuronRow neuronCol] = find(NeuronM);
    %neuronRow = [neuronRow;sizeY];
    %neuronCol = [neuronCol;sizeX];
    row = [row;sizeY];
    col = [col;sizeX];
    %neuronRow = [neuronRow;0];
    %neuronCol = [neuronCol;0];
    row = [row;0];
    col = [col;0];
    DensityM = uint8(hist3([row col],[sizeYTarget sizeXTarget]));
    %DensityNeuron = uint8(hist3([neuronRow neuronCol],[sizeYTarget sizeXTarget]));
%Map density matrix to values between o and 1
    minBrightness=0;
    maxBrightness=max(DensityM(:));
    DensityM = double(DensityM)./255;   
    DensityM = imadjust(DensityM, [double(minBrightness)/255;double(maxBrightness)/255], [0;1]);
    DensityM = uint8(DensityM.*255); 
%    DensityNeuron = double(DensityNeuron)./255;   
%    DensityNeuron = imadjust(DensityNeuron, [double(minBrightness)/255;double(maxBrightness)/255], [0;1]);
%    DensityNeuron = uint8(DensityNeuron.*255);
%end

function [filterDistance nonFilterDistance SphereArea NucleusArea64 markerPointCoordinates result] = calculateDensityDistribution(selectedWell, SphereAreaSizeX, SphereAreaSizeY, handles)
markerPointCoordinates=-1;
filterDistance = -1;
[filterDistance nonFilterDistance SphereArea NucleusArea64 markerPointCoordinates] = calculateMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, handles);
%Migration Distance calculation gives distance from Sphere and Nucleus
%We need: LinePic
if(filterDistance == -1)
   result = -1;
elseif (markerPointCoordinates==0)
   result = -1;
elseif(numel(markerPointCoordinates('0')) == 0)
    result = -1;
else
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
%NucleusM = csvHandler.CellPosMatrix;
%NeuronM = neuronHandler.CellPosMatrix;
NucleusM = csvHandler.CellPosMatrix(selectedWell);
NeuronM = neuronHandler.CellPosMatrix(selectedWell);
if(~isKey(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = sparse(7168, 7168);
end
NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
ringNumber = optionHandler.DensityDistributionRingNumber;
[sizeY, sizeX] = size(imageHandler.NucleusImage);
mat = sparse(sizeY, sizeX);
%Fix size of CSVCoordinate Matrices
[NucleusMSizeY NucleusMSizeX] = size(NucleusM);
[NeuronMSizeY NeuronMSizeX] = size(NeuronM);
[NeuronManualMSizeY NeuronManualMSizeX] = size(NeuronManualM);
[filterSizeY filterSizeX] = size(imageHandler.NucleusImage);
if(NucleusMSizeX > filterSizeX)
    NucleusM(:,filterSizeX+1:NucleusMSizeX) = [];
    csvHandler.CellPosMatrix(selectedWell) = NucleusM;
end
if(NeuronMSizeX > filterSizeX)
   NeuronM(:,filterSizeX+1:NeuronMSizeX) = [];
   neuronHandler.CellPosMatrix(selectedWell) = NeuronM;
end
if(NeuronManualMSizeX > filterSizeX)
                NeuronManualM(:,filterSizeX+1:NeuronManualMSizeX) = [];
                neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
end
if(NucleusMSizeY > filterSizeY)
                NucleusM(filterSizeY+1:NucleusMSizeY,:) = [];
                csvHandler.CellPosMatrix(selectedWell) = NucleusM;
end
if(NeuronMSizeY > filterSizeY)
                NeuronM(filterSizeY+1:NeuronMSizeY,:) = [];
                neuronHandler.CellPosMatrix(selectedWell) = NeuronM;
end
if(NeuronManualMSizeY > filterSizeY)
                NeuronManualM(filterSizeY+1:NeuronManualMSizeY,:) = [];
                neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
end
result = zeros(6,7);
hInner = impoly(handles.axes2,double(markerPointCoordinates('0')./10));
previousMask = createMask(hInner);
cbRings = get(handles.cbPolys ,'Value');
if(~cbRings)
    delete(hInner);
end
%Save hInner to file
foldername = imageHandler.Foldername;
subfoldername = [foldername '/ConvertedCellomics'];
save('-v7.3',strcat(subfoldername,'\MarkerPointCoordinates-',selectedWell),'markerPointCoordinates');

for i=1:ringNumber
%Create masks from LinePic and cut them from CSVs
%hTest = impoly(handles.axes2, [1 1; 500 500]);   
    hCurrent = impoly(handles.axes2,double(markerPointCoordinates(num2str(i*10))./10));    
    currentMask = logical(createMask(hCurrent));
    currentRing = logical(currentMask-previousMask); 
    currentRing = imresize(currentRing, [sizeY, sizeX]);
    previousMask = currentMask;
    %Save hCurrent to file
    if(~cbRings)
        delete(hCurrent);
    end
%Filter each ring on CSVs to calculate Density            
    Nucleus2 = logical(logical(NucleusM) .* currentRing);
    NeuronCell2 = logical(logical(NeuronM) .* currentRing);
    NeuronManual2 = logical(logical(NeuronManualM) .* currentRing);
    result(i,1) = nnz(currentRing);
    result(i,2) = nnz(Nucleus2);
    result(i,3) = nnz(NeuronCell2);
    result(i,4) = nnz(NeuronManual2);
    result(i,5) = result(i,3) / result(i,2);
    result(i,6) = result(i,4) / result(i,2);
    result(i,7) = result(i,2) / result(i,1);
end
end


function [filterDistance nonFilterDistance SphereArea NucleusArea64 markerPointCoordinates] = calculateMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, handles)
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
NucleusM = csvHandler.CellPosMatrix;
NeuronM = neuronHandler.CellPosMatrix;
ringNumber = optionHandler.DensityDistributionRingNumber;
SphereArea=0;
NucleusArea=0;
filterDistance=-1;
nonFilterDistance=-1;
markerPointCoordinates = 0;
NucleusArea64=0;
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
      selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
elseif(length(selectedWell) == 4)
        selectedWellLong = [selectedWell(1) '0' selectedWell(2:length(selectedWell))];
end
    
    %Get images from image folder
foldername = [imageHandler.Foldername '/ConvertedCellomics'];
imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig.tif'];
imagePathNucleusSmall = [foldername '/' selectedWellLong 'NucleusSmall.tif'];
imagePathNeuriteBig = [foldername '/' selectedWellLong 'NeuriteBig.tif'];
imagePathNeuriteSmall = [foldername '/' selectedWellLong 'NeuriteSmall.tif'];
if(isKey(NucleusM,selectedWell) && exist(imagePathNucleusSmall,'file'))
    imageHandler.ResizedNeuriteImage = imread(imagePathNeuriteSmall);
    imageHandler.ResizedNucleusImage = imread(imagePathNucleusSmall);    
    imageHandler.NeuriteImage = imread(imagePathNeuriteBig);
    imageHandler.NucleusImage = imread(imagePathNucleusBig);
    [sizeY, sizeX] = size(imageHandler.NucleusImage);
    
    NucleusM = csvHandler.CellPosMatrix;
    NeuronM = neuronHandler.CellPosMatrix;
    selectedWellLong=selectedWell;
    if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    elseif(length(selectedWell) == 4)
        selectedWellLong = [selectedWell(1) '0' selectedWell(2:length(selectedWell))];
    end
    if(isKey(NucleusM,selectedWell))
        NucleusM = NucleusM(selectedWell);
        NucleusM = logical(full(NucleusM));
        NeuronM = NeuronM(selectedWell);
        NeuronM = logical(full(NeuronM));
    end

    [DensityM] = getDensityDistributionsFromCSV(NucleusM,32,32, handles);
    [DensityNeuron] = getDensityDistributionsFromCSV(NeuronM,32,32, handles);
    
    %Find circle in Density Edge
    %Get point of Maximum Density
    [C,rowmaxarray]=max(DensityM);
    [maxvalue,xMaxSphere]=max(C);
    yMaxSphere=rowmaxarray(xMaxSphere);
    densityNeurons = DensityNeuron(uint8(yMaxSphere/2), uint8(xMaxSphere/2));
    i=0;    
     thresholdNucImage = optionHandler.MigDistLowerNucleusThreshold;
     %New Calculation method for Nucleus Area
     [thresholdedRows thresholdedCols] = find(imageHandler.NucleusImage > thresholdNucImage);
     thresholdedRows = [thresholdedRows;sizeY];
     thresholdedCols = [thresholdedCols;sizeX];
     thresholdedRows = [thresholdedRows;0];
     thresholdedCols = [thresholdedCols;0];
     NucleusDensity = uint8(hist3([thresholdedRows,thresholdedCols],[256,256]));
     figure(1);
     imshow(NucleusDensity);
     thresholdDensityImage = optionHandler.MigDistLowerDensityThreshold;
     [thresholdedRows thresholdedCols] = find(NucleusDensity > thresholdDensityImage);
     meanY = uint8(mean(thresholdedRows));
     meanX = uint8(mean(thresholdedCols));
     %Get point of density 255,255 nearby meanX and meanY
     found=0;
     radius=0;
     stack = java.util.Stack();
     xMax=0;
     yMax=0;
     while(found==0 && meanX-radius > 0 && meanY-radius > 0 && meanX+radius < 300 && meanY+radius < 300);
         for i=meanX-radius:meanX+radius
             for j=meanY-radius:meanY+radius
                 if(NucleusDensity(j,i)==255)
                    found=1;
                    yMax=int16(j);
                    xMax=int16(i);
                    stack.push([xMax yMax]);
                 end
             end
         end
         radius = radius+1;
     end
     minimizingFactorX = 256/SphereAreaSizeX;
     minimizingFactorY = 256/SphereAreaSizeY;
    % NucleusArea=NucleusDensity;
     NucleusArea = logical(zeros(256,256));    
     floodFillThresh=optionHandler.MigDistLowerFloodFillThreshold;
     NucleusArea = FloodFillIterative(NucleusDensity,NucleusArea,floodFillThresh,1,256,256,stack,0);
     %Also exclude pixels with less than 2 neighbours
     filter = [1 1 1;1 0 1;1 1 1];
     %filteredPix = filter2(filter, NucleusArea);
     %[rows cols] = find(filteredPix < 2);
     %NucleusArea(rows,cols) = 0;
     %NucleusArea = edge(NucleusArea);
     [posRows, posCols] = find(NucleusArea);
     posRows = posRows ./ minimizingFactorY;
     posCols = posCols ./ minimizingFactorX;
     
     
          
     
     NucleusArea64 = zeros(SphereAreaSizeY,SphereAreaSizeX);
     for i=1:numel(posRows)
        NucleusArea64(uint8(posRows(i)),uint8(posCols(i))) = 1;
     end
     [row,col] = find(NucleusArea64);
     if(length(row) > 0 && length(col) > 0)
       sY = sum(row)/length(row);
       sX = sum(col)/length(col);
       %Map sY and sX to Big Picture
       [sXMapped sYMapped] = MapPoint(sX,sY,sizeX,sizeY,SphereAreaSizeX,SphereAreaSizeY,1);
       
       %Take out Nuclei from NucleusM, for Nuclei which are too far from sX
       %and sY
       %Dictionary Point -> Distance
       [yIndices,xIndices,values] = find(NucleusM);       
       %for every point in NucleusM
         %calculate pdist from sX, sY
       distanceList = zeros(numel(yIndices),1);
       for (i=1:numel(yIndices))
         distanceList(i) = pdist([yIndices(i) xIndices(i);sYMapped sXMapped]);
       end
       %Sort by Distance
       [sortedDistances, sortedDistancesIndices] = sort(distanceList);
       %for Top 90% in Sorted Dictionary
         %Copy point to NucleusMSorted
       NucleusMSorted = sparse(sizeY,sizeX);
       for (i=1:numel(yIndices) * 0.97)
         NucleusMSorted(yIndices(sortedDistancesIndices(i)), xIndices(sortedDistancesIndices(i))) = 1;
       end
      
       %Calculate Density with NucleusMSorted
       %Calculate SphereArea with that
       [DensityMSized] = getDensityDistributionsFromCSV(NucleusMSorted,SphereAreaSizeX, SphereAreaSizeY, handles);
    % NucleusArea64=imfill(NucleusArea64,'holes');
%Flood fill everything from max Density above Threshold
%Outer circle is area for Sphere
       stack = java.util.Stack();
    %Validate point of max density: It's only valid, if Cellomics Neurons
    %are also in region
      i=0;
      while( i<=5)
        DensityM(yMaxSphere,xMaxSphere) =0;
        [C,rowmaxarray]=max(DensityMSized);
        [maxvalue,xMaxSphere]=max(C);
        yMaxSphere=rowmaxarray(xMaxSphere);
        densityNeurons = DensityNeuron(uint8(yMaxSphere/2), uint8(xMaxSphere/2));
        i=i+1;
        stack.push([xMaxSphere yMaxSphere]);
      end
      xMax=xMax/3;
      yMax=yMax/3;
      SphereArea = logical(zeros(SphereAreaSizeY,SphereAreaSizeX));
      SphereArea = FloodFillIterative(DensityMSized,SphereArea,1,1,SphereAreaSizeX,SphereAreaSizeY,stack,0);
      SphereArea=imfill(SphereArea,'holes');
      SphereArea = edge(SphereArea);
      MigrationPic = SphereArea + NucleusArea64;
      figure(1);
      imshow(MigrationPic);
    %hPixelInfo = impixelinfo;        
    %imageHandler.MousePosition = hPixelInfo;
%Für weitere Verarbeitung:
% Berechne den Schwerpunkt aller weißen Pixel in NucleusArea
     
      minimizingFactorX = 32/SphereAreaSizeX;
      minimizingFactorY = 32/SphereAreaSizeY;
      mediumDistanceList = zeros(0,0);
      mediumDistanceListInInterval = zeros(0);
     
% Gehe jeweils im 10 Grad Winkel von dort aus los.
       successCount=0;
       %Declaration for density distribution calculation
       %Mapping: Distance in 5 steps
       markerPointCoordinates = containers.Map();
       for i=0:ringNumber
        markerPointCoordinates(num2str(i*10)) = zeros(0,0); 
       end
       counter=0;
       for i=pi/32:pi/32:2*pi
        
        deltaX = cos(i);
        deltaY = sin(i);
        xNucleus = 0;
        yNucleus = 0;
        xOuterSphere = 0;
        yOuterSphere = 0;    
        x=sX;
        y=sY;
        SphereArea(uint8(y),uint8(x)) =0;
        while(xOuterSphere ==0 && uint8(x) <= SphereAreaSizeX && uint8(y) <= SphereAreaSizeY && uint8(x) > 0 && uint8(y) > 0)
          [xBig yBig] = MapPoint(x,y,256,256,SphereAreaSizeX,SphereAreaSizeY,0);
            if(NucleusArea(uint8(yBig),uint8(xBig)) == 1)
            [xNucleus yNucleus]=MapPoint(xBig,yBig,sizeX,sizeY,256,256,0);
          end
          if(SphereArea(uint8(y),uint8(x)) == 1)
            [xOuterSphere yOuterSphere] = MapPoint(x,y,sizeX,sizeY,SphereAreaSizeX,SphereAreaSizeY,1);
          end
          x=x+deltaX;
          y=y+deltaY;                   
        end
        if(xNucleus == 0)
            [xNucleus yNucleus]= MapPoint(sX,sY,sizeX,sizeY,SphereAreaSizeX,SphereAreaSizeY,1);
        end
        
        if(xOuterSphere~=0 && yOuterSphere~=0 && xNucleus ~= 0 && yNucleus ~= 0)
          successCount = successCount+1;
          mediumDistanceList = [mediumDistanceList; pdist([xOuterSphere yOuterSphere;xNucleus yNucleus]) successCount];        
          %Create 5 points between xyOuterSphere and xyNucleus. Save these
          %points to markerPointCoordinates
          stepX = int16((int16(xOuterSphere) - int16(xNucleus)) / (ringNumber));
          stepY = int16((int16(yOuterSphere) - int16(yNucleus)) / (ringNumber));
          for k=0:ringNumber
            markerPointCoordinates(num2str(k*10)) = [markerPointCoordinates(num2str(k*10)); [int16(xNucleus)+k*stepX int16(yNucleus)+k*stepY]];
          end
        end
       end
     
     distancestddev = std(mediumDistanceList(:,1));
     mediumDistance = mean(mediumDistanceList(:,1));
     maxDistance = max(mediumDistanceList(:,1));
     for i=1:numel(mediumDistanceList(:,1))
      currentDistance = mediumDistanceList(i,1);
      if(i-1>0)
        leftNeighbourDistance = mediumDistanceList(i-1,1);
      else
          leftNeighbourDistance = mediumDistanceList(numel(mediumDistanceList(:,1)),1);
      end
      if(i+1<numel(mediumDistanceList(:,1)))
        rightNeighbourDistance = mediumDistanceList(i+1,1);
      else
          rightNeighbourDistance = 1;
      end
      threshold = (maxDistance-(0.65*maxDistance));
    %Check if current distance is in Interval. Check also if one of both
    %neighbourpoints is also in interval
      %if(currentDistance >= 150 && currentDistance >= (maxDistance-(2.0*distancestddev)))
      if(currentDistance >= threshold && (leftNeighbourDistance > threshold || rightNeighbourDistance > threshold))
        mediumDistanceListInInterval = [mediumDistanceListInInterval currentDistance];
      else
          %If not: Set points of MarkerPointCoordinates to inner point.
          currentCounter = mediumDistanceList(i,2);
          markerReferencePoint = markerPointCoordinates('0');
          markerReferencePoint = markerReferencePoint(currentCounter,:);
          for i=1:ringNumber
            markerPointList = markerPointCoordinates(num2str(i*10));
            markerPointList(currentCounter,:) = markerReferencePoint;
            markerPointCoordinates(num2str(i*10)) = markerPointList;
          end          
      end
      %Extra check: Spike prevention: At least two markerpoints on one ring
      %not set back.      
     end
     
     
     
%mediumDistance = mediumDistance / successCount;
     mediumDistance = num2str(mediumDistance,'%f');
     mediumDistanceSelected = num2str(mean(mediumDistanceListInInterval),'%f');
     disp('Medium distance total (in Pixel): ');
     disp(mediumDistance);
     disp('Medium distance selected (in Pixel): ');
     disp(mediumDistanceSelected);
     filterDistance = mediumDistanceSelected;
     nonFilterDistance = mediumDistance;
     end
    else
        filterDistance=-1;
        nonFilterDistance=-1;
    end
% --- Executes on button press in cbNucleusPicture.
function cbNucleusPicture_Callback(hObject, eventdata, handles)
% hObject    handle to cbNucleusPicture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles)
% Hint: get(hObject,'Value') returns toggle state of cbNucleusPicture


% --- Executes on button press in cbNeuritePicture.
function cbNeuritePicture_Callback(hObject, eventdata, handles)
% hObject    handle to cbNeuritePicture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles)
% Hint: get(hObject,'Value') returns toggle state of cbNeuritePicture


% --- Executes on button press in cbNucleusMatrix.
function cbNucleusMatrix_Callback(hObject, eventdata, handles)
% hObject    handle to cbNucleusMatrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles)
% Hint: get(hObject,'Value') returns toggle state of cbNucleusMatrix


% --- Executes on button press in cbCellomicsNeuronMatrix.
function cbCellomicsNeuronMatrix_Callback(hObject, eventdata, handles)
% hObject    handle to cbCellomicsNeuronMatrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles)
% Hint: get(hObject,'Value') returns toggle state of cbCellomicsNeuronMatrix


% --- Executes on button press in cbManualNeuronMatrix.
function cbManualNeuronMatrix_Callback(hObject, eventdata, handles)
% hObject    handle to cbManualNeuronMatrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles)
% Hint: get(hObject,'Value') returns toggle state of cbManualNeuronMatrix


% --------------------------------------------------------------------
function MigrationDistanceAll_Callback(hObject, eventdata, handles)
% hObject    handle to MigrationDistanceAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Outer Loop: Iterate over all Wells.
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
[FileName,PathName] = uiputfile('migdist.csv');
fileID = fopen(strcat(PathName,'\',FileName),'w');
handles = guidata(handles.figure1);
ringCount = optionHandler.DensityDistributionRingNumber;
densityWidth = optionHandler.MigrationDistanceDensityImageXSize;
densityHeight = optionHandler.MigrationDistanceDensityImageYSize
wellList = get(handles.lbWell, 'string');
waitbarHandle = waitbar(0,'Calculating and exporting Migration Distances for all wells');
selectedWellNumber = get(handles.lbWell,'Value');
exportText = sprintf('Well;Ring Number;Migration Distance(px);Manual measured Migration Distance (px);Migration Distance without filter (px);Number of Pixels;Number Nuclei;Cellomics Neurons;Manual Neurons;Cellomics Neurons per Nuclei;Manual Neurons per Nuclei;Nuclei per Pixels;\r\n');
for currentWellIndex=1:numel(wellList)
    currentWell = wellList(currentWellIndex);
    currentWell = currentWell{1};
    disp(currentWell);
    [distFilter, distWithoutFilter, SphereArea, NucleusArea, markerPointCoordinates densDist]=calculateDensityDistribution(currentWell,densityWidth,densityHeight,handles);
    if(numel(densDist) == 1)
        %Fill up densDist with dummy values
        densDist = zeros(ringCount,7);        
    end
    if(distFilter >-1)
        distFilter = num2str(distFilter,'%f');
        distFilter=strrep(distFilter,'.',',');
        distWithoutFilter = num2str(distWithoutFilter,'%f');
        distWithoutFilter=strrep(distWithoutFilter,'.',',');
        for i=1:optionHandler.DensityDistributionRingNumber        
          if(isKey(imageHandler.MigDistWellMapping,currentWell))
            manualDist = num2str(imageHandler.MigDistWellMapping(currentWell), '%f')
            manualDist = strrep(manualDist,'.',',');
            exportText = [exportText currentWell ';' num2str(i) ';' distFilter ';' manualDist ';' distWithoutFilter ';' strrep(num2str(densDist(i,1)),'.',',') ';' strrep(num2str(densDist(i,2)),'.',',') ';' strrep(num2str(densDist(i,3)),'.',',') ';' strrep(num2str(densDist(i,4)),'.',',') ';' strrep(num2str(densDist(i,5)),'.',',') ';' strrep(num2str(densDist(i,6)),'.',',') ';' strrep(num2str(densDist(i,7)),'.',',') sprintf('\r\n')];
          else
            exportText = [exportText currentWell ';' num2str(i) ';' distFilter ';;' distWithoutFilter ';' strrep(num2str(densDist(i,1)),'.',',') ';' strrep(num2str(densDist(i,2)),'.',',') ';' strrep(num2str(densDist(i,3)),'.',',') ';' strrep(num2str(densDist(i,4)),'.',',') ';' strrep(num2str(densDist(i,5)),'.',',') ';' strrep(num2str(densDist(i,6)),'.',',') ';' strrep(num2str(densDist(i,7)),'.',',') sprintf('\r\n')];
          end
        end
    end
    waitbar((currentWellIndex/numel(wellList)),waitbarHandle);
end
%Write exportText to File
fprintf(fileID,'%s',exportText);
close(waitbarHandle);


% --------------------------------------------------------------------
function Algorithms_Callback(hObject, eventdata, handles)
% hObject    handle to Algorithms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function EdgeCompositNeuronCount_Callback(hObject, eventdata, handles)
% hObject    handle to EdgeCompositNeuronCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Make overlay of both channels
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
%ToDo: Calculate density distribution to exclude too dense areas
neuriteEdges = edge(imageHandler.NeuriteImage,'nothinning');
nucleusEdges = edge(imageHandler.NucleusImage,'nothinning');
overlaypic = imfuse(nucleusEdges, neuriteEdges);
%Get white points (Neuron Points)
neuronPointIndices = overlaypic(:,:,1) > 250 & overlaypic(:,:,2) > 250 & overlaypic(:,:,3) > 250;
figure(2);
imshow(neuronPointIndices);
figure(3)
imshow(overlaypic);
%For every white point:
%Look for next nucleus in csvHandler
%Mark this nucleus as algorithm neuron

%overlaypic = edge(overlaypic);
%imshow(overlaypicsmall);


% --------------------------------------------------------------------
function StatisticsCSVExport_Callback(hObject, eventdata, handles)
% hObject    handle to StatisticsCSVExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Extract Rectangles from ImageHandler. 
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
[FileName,PathName] = uiputfile('WellStatistics.csv');
fileID = fopen(strcat(PathName,'\',FileName),'w');
waitbarHandle = waitbar(0,'Calculating and exporting Statistics for all wells');
%Iterate over all wells
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
exportText = sprintf('Well;Filter;#Nuclei by Cellomics;#Neurons Cellomics;#Neurons by Clicks;%%Neurons by Cellomics;%%Neurons by Clicks;Migration Distance calculated (px);Migration Distance manual (px)\r\n');
for i=1:numel(wellList)    
    selectedWell = wellList{i};
    disp(strcat('Current well: ',selectedWell)); 
    manualDist = '';
    if(isKey(imageHandler.MigDistWellMapping,selectedWell))
            manualDist = num2str(imageHandler.MigDistWellMapping(currentWell), '%f')
            manualDist = strrep(manualDist,'.',',');
            exportText = [exportText currentWell ';' distFilter ';' manualDist ';' distWithoutFilter sprintf('\r\n')];
    end
    [distFilter, distWithoutFilter, SphereArea, NucleusArea]=calculateMigrationDistance(selectedWell,64,64,handles);
    if(distFilter >-1)
        distFilter = num2str(distFilter,'%f');
        distFilter=strrep(distFilter,'.',',');
    else
        distFilter = '';
    end
    filterList = get(handles.lbFilter, 'string');
    
    %Load filters for current well
    filterpath = [imageHandler.Foldername '/ConvertedCellomics/' selectedWell '.mat'];
    if(exist(filterpath,'file'))
        load(filterpath);
        imageHandler.PositiveFilters = FMask.PositiveFilters;
        imageHandler.NegativeFilters = FMask.NegativeFilters;
    else
        imageHandler.PositiveFilters=0;
        imageHandler.NegativeFilters=0;
    end
    %Iterate over all filters
    filterList = imageHandler.PositiveFilters;
    for j=1:numel(filterList)
        if(iscell(filterList))
            filterMask = filterList(j);
            filterMask = filterMask{1};
            %Remove excluded parts of picture
            if(iscell(imageHandler.NegativeFilters) > 0)
                negativeMask = imageHandler.NegativeFilters;
                for i=1:numel(negativeMask)
                    workMask = negativeMask{i};
                    workMask = int8(workMask -1);
                    workMask = logical(workMask .* (-1));
                    workMask = logical(workMask);
                    filterMask = logical((filterMask .* workMask));
                end
            end
            workMask = 0;
            negativeMask = 0;
            %Multiply filterMask with Neuron and Nucleus Matrices
            NucleusM = csvHandler.CellPosMatrix(selectedWell);
            NeuronM = neuronHandler.CellPosMatrix(selectedWell);
            NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
           
            %Fix size of CSVCoordinate Matrices
            [NucleusMSizeY NucleusMSizeX] = size(NucleusM);
            [NeuronMSizeY NeuronMSizeX] = size(NeuronM);
            [NeuronManualMSizeY NeuronManualMSizeX] = size(NeuronManualM);
            [filterSizeY filterSizeX] = size(filterMask);
            if(NucleusMSizeX > filterSizeX)
                NucleusM(:,filterSizeX+1:NucleusMSizeX) = [];
                csvHandler.CellPosMatrix(selectedWell) = NucleusM;
            end
            if(NeuronMSizeX > filterSizeX)
                NeuronM(:,filterSizeX+1:NeuronMSizeX) = [];
                neuronHandler.CellPosMatrix(selectedWell) = NeuronM;
            end
            if(NeuronManualMSizeX > filterSizeX)
                NeuronManualM(:,filterSizeX+1:NeuronManualMSizeX) = [];
                neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
            end
            if(NucleusMSizeY > filterSizeY)
                NucleusM(filterSizeY+1:NucleusMSizeY,:) = [];
                csvHandler.CellPosMatrix(selectedWell) = NucleusM;
            end
            if(NeuronMSizeY > filterSizeY)
                NeuronM(filterSizeY+1:NeuronMSizeY,:) = [];
                neuronHandler.CellPosMatrix(selectedWell) = NeuronM;
            end
            if(NeuronManualMSizeY > filterSizeY)
                NeuronManualM(filterSizeY+1:NeuronManualMSizeY,:) = [];
                neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
            end
            
            NucleusM = NucleusM .* filterMask;
            numberNuclei = nnz(NucleusM);
            NucleusM = 0;
            NeuronM = NeuronM .* filterMask;   
            numberNeurons = nnz(NeuronM);
            NeuronM=0;
            NeuronManualM = NeuronManualM .* filterMask;
            numberManual = nnz(NeuronManualM);
            percentageCellomics = (numberNeurons/numberNuclei) * 100;
            percentageManual = (numberManual/numberNuclei) * 100;
            numberNuclei = num2str(numberNuclei,'%f');
            numberNuclei = strrep(numberNuclei, '.',',');
            numberNeurons = num2str(numberNeurons,'%f');
            numberNeurons = strrep(numberNeurons, '.',',');
            numberManual = num2str(numberManual,'%f');
            numberManual = strrep(numberManual, '.',',');
            percentageCellomics = num2str(percentageCellomics,'%f');
            percentageCellomics = strrep(percentageCellomics, '.',',');
            percentageManual = num2str(percentageManual,'%f');
            percentageManual = strrep(percentageManual, '.',',');
            %Add entry to CSV Exort string
            exportText = [exportText selectedWell ';' num2str(j) ';' numberNuclei ';' numberNeurons ';' numberManual ';' percentageCellomics ';' percentageManual ';' distFilter ';' manualDist sprintf('\r\n')];
        end
    end
    waitbar((i/numel(wellList)),waitbarHandle);
end
%Write exportText to File
fprintf(fileID,'%s',exportText);
close(waitbarHandle);
handles.NeuronCoordinates = neuronHandler;
handles.CSVCoordinates = csvHandler;
guidata(handles.figure1, handles)

% --- Executes on button press in cbNegativeFilter.
function cbNegativeFilter_Callback(hObject, eventdata, handles)
% hObject    handle to cbNegativeFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbNegativeFilter


% --------------------------------------------------------------------
function MeasureMigDist_Callback(hObject, eventdata, handles)
% hObject    handle to MeasureMigDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
if(imageHandler.MeasureMigDistMode==0)
    imageHandler.MeasureMigDistMode=1;
    neuronHandler.ManualCountingMode = 0;
    imageHandler.MigDistPoints = zeros(2,0);
else
    imageHandler.MeasureMigDistMode=0;
    %Calculate and save measured migration distance for current well
    mediumMigDist=0;
    count = 0;
    for i=1:numel(imageHandler.MigDistPoints)/2
        if(mod(i,2) ~= 0)
            pointOne = [imageHandler.MigDistPoints(i,1) imageHandler.MigDistPoints(i,2)];
            pointTwo = [imageHandler.MigDistPoints(i+1,1) imageHandler.MigDistPoints(i+1,2)];
            mediumMigDist = mediumMigDist + pdist([pointOne;pointTwo]);
            count = count+1;
        end
    end
    mediumMigDist = mediumMigDist / count; 
    imageHandler.MigDistWellMapping(selectedWell) = mediumMigDist;
    disp(strcat('Manual Migration distance is ',num2str(mediumMigDist),' pixels'));    
end
handles.ImageHandler = imageHandler;
handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles)

% --------------------------------------------------------------------
function Options_Callback(hObject, eventdata, handles)
% hObject    handle to Options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Create Input Dialogue
handles = guidata(handles.figure1);
if(isfield(handles,'OptionHandler'))
    optionHandler = handles.OptionHandler;
    if(optionHandler == 0)
    optionHandler = Option();
    optionHandler.HistogramMaxNucleus = 2995;
    optionHandler.HistogramMaxNeurite = 500;
    optionHandler.HistogramMinNucleus = 405;
    optionHandler.HistogramMinNeurite = 70;
    optionHandler.MigDistLowerNucleusThreshold = 250;
    optionHandler.MigDistLowerDensityThreshold = 250;
    optionHandler.MigDistLowerFloodFillThreshold = 175;
    optionHandler.DensityDistributionRingNumber = 10;
    optionHandler.MigrationDistanceDensityImageXSize = 32;
    optionHandler.MigrationDistanceDensityImageYSize = 32;
    end
else
    optionHandler = Option();
    optionHandler.HistogramMaxNucleus = 2995;
    optionHandler.HistogramMaxNeurite = 500;
    optionHandler.HistogramMinNucleus = 405;
    optionHandler.HistogramMinNeurite = 70;
    optionHandler.MigDistLowerNucleusThreshold = 250;
    optionHandler.MigDistLowerDensityThreshold = 250;
    optionHandler.MigDistLowerFloodFillThreshold = 175;
    optionHandler.DensityDistributionRingNumber = 10;
    optionHandler.MigrationDistanceDensityImageXSize = 32;
    optionHandler.MigrationDistanceDensityImageYSize = 32;
end
prompt = {'Histogramm Correction: Max Intensity Nucleus:','Histogramm Correction: Max Intensity Neurite:','Histogramm Correction: Min Intensity Nucleus:','Histogramm Correction: Min Intensity Neurite:', 'Migration Distance: Lower Threshold for Nucleus Pixels:', 'Migration Distance: Lower Threshold for Density Pixels:','Migration Distance: Lower Threshold for Floodfilling Nucleus Area:','Density Distribution: Number of Rings:','Density Distribution: Width of Density Picture','Density Distribution: Height of Density Picture'};
dlg_title = 'Options';
num_lines = 1;
def = {num2str(optionHandler.HistogramMaxNucleus),num2str(optionHandler.HistogramMaxNeurite),num2str(optionHandler.HistogramMinNucleus),num2str(optionHandler.HistogramMinNeurite),num2str(optionHandler.MigDistLowerNucleusThreshold),num2str(optionHandler.MigDistLowerDensityThreshold), num2str(optionHandler.MigDistLowerFloodFillThreshold), num2str(optionHandler.DensityDistributionRingNumber),num2str(optionHandler.MigrationDistanceDensityImageXSize), num2str(optionHandler.MigrationDistanceDensityImageYSize)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
optionHandler.HistogramMaxNucleus = str2num(answer{1});
optionHandler.HistogramMaxNeurite = str2num(answer{2});
optionHandler.HistogramMinNucleus = str2num(answer{3});
optionHandler.HistogramMinNeurite = str2num(answer{4});
optionHandler.MigDistLowerNucleusThreshold = str2num(answer{5});
optionHandler.MigDistLowerDensityThreshold = str2num(answer{6});
optionHandler.MigDistLowerFloodFillThreshold = str2num(answer{7});
optionHandler.DensityDistributionRingNumber = str2num(answer{8});
optionHandler.MigrationDistanceDensityImageXSize = str2num(answer{9});
optionHandler.MigrationDistanceDensityImageYSize = str2num(answer{10});
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles)


% --------------------------------------------------------------------
function DensDist_Callback(hObject, eventdata, handles)
% hObject    handle to DensDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
densityWidth = optionHandler.MigrationDistanceDensityImageXSize;
densityHeight = optionHandler.MigrationDistanceDensityImageYSize;
[filterDistance nonFilterDistance SphereArea NucleusArea64 markerPointCoordinates result] = calculateDensityDistribution(selectedWell, densityWidth, densityHeight, handles);
if(numel(result) > 1)
    disp('Ring Number;Pixels in Ring;Nuclei in Ring;Neurons by Cellomics in Ring;Manual counted Neurons in Ring; Cellomics Neurons per Nuclei;Manual Neurons per Nuclei; Nuclei per Pixels');
    for i=1:optionHandler.DensityDistributionRingNumber
        disp([num2str(i) ';' strrep(num2str(result(i,1)),'.',',') ';' strrep(num2str(result(i,2)),'.',',') ';' strrep(num2str(result(i,3)),'.',',') ';' strrep(num2str(result(i,4)),'.',',') ';' strrep(num2str(result(i,5)),'.',',') ';' strrep(num2str(result(i,6)),'.',',') ';' strrep(num2str(result(i,7)),'.',',')]);
    end
end


% --- Executes on button press in cbPolys.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to cbPolys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbPolys


% --- Executes on button press in cbPolys.
function cbPolys_Callback(hObject, eventdata, handles)
% hObject    handle to cbPolys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
ringNumber = optionHandler.DensityDistributionRingNumber;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
if(get(hObject,'Value') >0)
    path = [imageHandler.Foldername '/ConvertedCellomics/MarkerPointCoordinates-' selectedWell '.mat'];
    if(exist(path,'file'))
        load(path);
        hInner = impoly(handles.axes2,double(markerPointCoordinates('0')./10));
        previousMask = createMask(hInner);
        for i=1:ringNumber  
          hCurrent = impoly(handles.axes2,double(markerPointCoordinates(num2str(i*10))./10));
        end
    end
else
    RefreshUnzoomedImage(handles); 
end


% --------------------------------------------------------------------
function ConvertSingleChambTo96er_Callback(hObject, eventdata, handles)
% hObject    handle to ConvertSingleChambTo96er (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
      selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
elseif(length(selectedWell) == 4)
        selectedWellLong = [selectedWell(1) '0' selectedWell(2:length(selectedWell))];
end
%Get images from image folder
foldername = [imageHandler.Foldername '/ConvertedCellomics'];
imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig.tif'];
imagePathNucleusSmall = [foldername '/' selectedWellLong 'NucleusSmall.tif'];
imagePathNeuriteBig = [foldername '/' selectedWellLong 'NeuriteBig.tif'];
imagePathNeuriteSmall = [foldername '/' selectedWellLong 'NeuriteSmall.tif'];   

%Iteriere �ber Filter f�r das aktuelle Well
%For every positive entry in Filter List: Cut pictures

path = [imageHandler.Foldername '/ConvertedCellomics/' selectedWell '.mat'];

load(path);
imageHandler.PositiveFilters = FMask.PositiveFilters;
imageHandler.NegativeFilters = FMask.NegativeFilters;

filterCount = numel(imageHandler.PositiveFilters);
FMask=0;
imageHandler.PositiveFilters = 0;
imageHandler.NegativeFilters = 0;
for i=1:filterCount
    load(path);
    filterHandle = FMask.PositiveFilters(i);
    FMask=0;
    filterHandle = filterHandle{1};
    filterMask = logical(filterHandle);
    filterMaskSmall = imresize(filterMask,0.1);    
    nucImageResized = double(imageHandler.ResizedNucleusImage).*double(filterMaskSmall);
    %nucImage = zeros(size(imageHandler.NucleusImage,1),size(imageHandler.NucleusImage,2));
    %for i=1:size(nucImage,1)/4
    %    disp([num2str(i) 'of' num2str(size(nucImage,1))]);
    %    nucImage(i:i+3,:) = uint16(imageHandler.NucleusImage(i:i+3,:)) .* uint16(filterMask(i:i+3,:));
    %    i=i+3;
    %    for j=1:size(nucImage,2)
    %        nucImage(i,j) = uint16(imageHandler.NucleusImage(i,j)) * uint16(filterMask(i,j));
    %    end
    %end
    
    %nucImage = logical(imageHandler.NucleusImage).*logical(filterMask);
    %neuImage = uint16(imageHandler.NeuriteImage).*uint16(filterMask);
    
     nucResizedTransposed = transpose(nucImageResized);
     %nucImageTransposed = transpose(nucImage);
     %Cut Matrix and set new Zoom State
     [maxPossibleRow, maxPossibleCol] = (size(nucImageResized));
     maxPossibleRow = maxPossibleRow * 10;
     maxPossibleCol = maxPossibleCol * 10;
     [minRow, minCol] = (find(nucImageResized,1,'first'));
     minRow = minRow * 10;
     minCol = minCol * 10;
     [minColTrans, minRowTrans] = find(nucResizedTransposed,1,'first');
     minRowTrans = minRowTrans * 10;
     minColTrans = minColTrans * 10;
     minRow = min(minRow, minRowTrans);
     minCol = min(minCol, minColTrans);
     [maxRow, maxCol] = (find(nucImageResized,1,'last'));
     maxRow = maxRow * 10;
     maxCol = maxCol * 10;
     [maxColTrans, maxRowTrans] = find(nucResizedTransposed,1,'last');    
     maxColTrans = maxColTrans * 10;    
     maxRowTrans = maxRowTrans * 10;    
     maxRow = max(maxRow,maxRowTrans);
     maxCol = max(maxCol,maxColTrans);
     maxRow = min(maxRow, maxPossibleRow);
     maxCol = min(maxCol, maxPossibleCol);
     %resize filter mask for pic
     %filterMaskSmall = imresize(filterMask,0.1);
     %Multiply mask with image
     nucResizedTransposed=0;
     nucImageResized=0;
     if(maxRow>13824)
         maxRow=13824;
     end
     if(maxCol>10752)
         maxCol=10752;
     end
     
     nucImage = uint16(imageHandler.NucleusImage(minRow:maxRow,minCol:maxCol));
     neuImage = uint16(imageHandler.NeuriteImage(minRow:maxRow,minCol:maxCol));
     %if(size(nucImage,1) > 7168)
     %    nucImage(7169:end,:) = [];
     %elseif(size(nucImage,1) < 7168)
     
     NucleusM = logical(csvHandler.CellPosMatrix(selectedWell));
     NeuronM = logical(neuronHandler.CellPosMatrix(selectedWell));
    
    %Better: 
     currentWellNucMatrix = NucleusM(minRow:maxRow,minCol:maxCol);
     currentWellNeuMatrix = NeuronM(minRow:maxRow,minCol:maxCol);
     
     filterMaskResized = filterMask(minRow:maxRow,minCol:maxCol);
     
     currentWellNucMatrix = logical(logical(currentWellNucMatrix) .* logical(filterMaskResized));
     currentWellNeuMatrix = logical(logical(currentWellNeuMatrix) .* logical(filterMaskResized));
     nucImage = uint16(uint16(nucImage) .* uint16(filterMaskResized));
     neuImage = uint16(uint16(neuImage) .* uint16(filterMaskResized));
     
     nucImage = vertcat(nucImage,zeros(7168-size(nucImage,1),size(nucImage,2)));
     nucImage = horzcat(nucImage,zeros(size(nucImage,1),7168-size(nucImage,2)));
     neuImage = vertcat(neuImage,zeros(7168-size(neuImage,1),size(neuImage,2)));
     neuImage = horzcat(neuImage,zeros(size(neuImage,1),7168-size(neuImage,2)));
     
    %end
    %if(size(nucImage,2) > 7168)
    %    nucImage(:,7169:end) = [];
    %elseif(size(nucImage,2) < 7168)
    %end
    %nucImage(7169:end,7169:end) = [];
    %neuImage(7169:end,7169:end) = [];
    nucImageSmall = imresize(nucImage,0.1);
    neuImageSmall = imresize(neuImage,0.1);
    %Save new pics to ConvertedCellomics Folder
    imagePathNucleusBigNew = [foldername '/' selectedWellLong '_' num2str(i) 'NucleusBig.tif'];
    imagePathNucleusSmallNew = [foldername '/' selectedWellLong '_' num2str(i) 'NucleusSmall.tif'];
    imagePathNeuriteBigNew = [foldername '/' selectedWellLong '_' num2str(i) 'NeuriteBig.tif'];
    imagePathNeuriteSmallNew = [foldername '/' selectedWellLong '_' num2str(i) 'NeuriteSmall.tif']; 
    imwrite(nucImage,imagePathNucleusBigNew);
    imwrite(neuImage,imagePathNeuriteBigNew);
    imwrite(nucImageSmall,imagePathNucleusSmallNew);
    imwrite(neuImageSmall,imagePathNeuriteSmallNew);
    neuImage=0;
    nucImage=0;
    neuImageSmall=0;
    nucImageSmall=0;
    NucleusM = 0;
    NeuronM = 0;
    %Schneide CSVCoordinates zurecht
    






    currentWellNucMatrix = vertcat(currentWellNucMatrix,zeros(7168-size(currentWellNucMatrix,1),size(currentWellNucMatrix,2)));
    currentWellNucMatrix = horzcat(currentWellNucMatrix,zeros(size(currentWellNucMatrix,1),7168-size(currentWellNucMatrix,2)));
    currentWellNeuMatrix = vertcat(currentWellNeuMatrix,zeros(7168-size(currentWellNeuMatrix,1),size(currentWellNeuMatrix,2)));
    currentWellNeuMatrix = horzcat(currentWellNeuMatrix,zeros(size(currentWellNeuMatrix,1),7168-size(currentWellNeuMatrix,2)));
    csvHandler.CellPosMatrix([selectedWell '_' num2str(i)]) = currentWellNucMatrix;
    neuronHandler.CellPosMatrix([selectedWell '_' num2str(i)]) = currentWellNeuMatrix;    
    
    %Create sparse matrices for new Wells:
    neuronHandler.ManualNeuronPositionsSparse([selectedWell '_' num2str(i)]) = sparse(7168, 7168);    
end
remove(csvHandler.CellPosMatrix,selectedWell);
remove(neuronHandler.CellPosMatrix,selectedWell);
%ToDo:
%Delete old Files
delete(imagePathNucleusBig, imagePathNucleusSmall, imagePathNeuriteBig, imagePathNeuriteSmall);
handles.CSVCoordinates=csvHandler;
handles.NeuronCoordinates=neuronHandler;
guidata(handles.figure1, handles)
ReloadWellList(handles);


% --- Executes when selected object is changed in uipanel2.
function ReloadWellList(handles)
% hObject    handle to the selected object in uipanel2 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

%If csvHandler is available: Get all available Wells from CSVHandler in
%CellPosMatrix

handles = guidata(handles.figure1);
if(isfield(handles,'NeuronCoordinates') && isfield(handles,'CSVCoordinates'))
           %Map Neuron Positions to Nucleus Positions
    wellList = get(handles.lbWell, 'string');
    csvHandler = handles.CSVCoordinates;
    neuronHandler = handles.NeuronCoordinates;
    WellNeuronDict = neuronHandler.CellPosMatrix;
    WellNucleusDict = csvHandler.CellPosMatrix;
    wellList = keys(WellNeuronDict);
    %Remove all entries from lbWell
    set(handles.lbWell, 'String', '');    
    fulltext = '';
    for i=1:numel(wellList)
         fulltext = [fulltext;wellList(i)];
    end
    set(handles.lbWell, 'String', fulltext);
end

