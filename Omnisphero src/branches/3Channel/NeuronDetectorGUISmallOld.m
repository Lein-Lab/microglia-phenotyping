function varargout = NeuronDetectorGUISmall(varargin)
% NEURONDETECTORGUISMALL MATLAB code for NeuronDetectorGUISmall.fig
%      NEURONDETECTORGUISMALL, by itself, creates a new NEURONDETECTORGUISMALL or raises the existing
%      singleton*.
%
%      H = NEURONDETECTORGUISMALL returns the handle to a new NEURONDETECTORGUISMALL or the handle to
%      the existing singleton*.
%
%      NEURONDETECTORGUISMALL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEURONDETECTORGUISMALL.M with the given input arguments.
%
%      NEURONDETECTORGUISMALL('Property','Value',...) creates a new NEURONDETECTORGUISMALL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NeuronDetectorGUISmall_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NeuronDetectorGUISmall_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NeuronDetectorGUISmall

% Last Modified by GUIDE v2.5 05-Apr-2013 14:15:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NeuronDetectorGUISmall_OpeningFcn, ...
                   'gui_OutputFcn',  @NeuronDetectorGUISmall_OutputFcn, ...
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


% --- Executes just before NeuronDetectorGUISmall is made visible.
function NeuronDetectorGUISmall_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NeuronDetectorGUISmall (see VARARGIN)

% Choose default command line output for NeuronDetectorGUISmall
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NeuronDetectorGUISmall wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NeuronDetectorGUISmall_OutputFcn(hObject, eventdata, handles) 
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
    
        % Load all images in  folder.
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
          for i=1:numel(allFiles)
            nucleusRegexp = strcat('ARRAY.*\d*_(',currentWell,')f(\d+)d0');
            neuriteRegexp = strcat('ARRAY.*\d*_(',currentWell,')f(\d+)d1');
            tokensNucleus = regexpi(allFiles(i).name, nucleusRegexp, 'tokens');
            tokensNeurite = regexpi(allFiles(i).name, neuriteRegexp, 'tokens');
        
            if(length(tokensNucleus) > 0)
                tokensNucleus = tokensNucleus{1};
                fileWell = tokensNucleus{1};                   
            elseif(length(tokensNeurite) > 0)
                tokensNeurite = tokensNeurite{1};
                fileWell = tokensNeurite{1};           
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
              if((length(filenamenuc) < 2) || (length(filenameneu) < 2))
                  %Only if all images are available for the well
                  breaked = 1;
                  break;
              end
              imageString = strcat(foldername, '\', filenamenuc);
              currentNucleusPic = imread(imageString);
              maxValueNucleus = max(currentNucleusPic(:));
              imageString = strcat(foldername, '\', filenameneu);
              currentNeuritePic = imread(imageString);
              maxValueNeurite = max(currentNeuritePic(:));
              maxValueWhole = max(maxValueNucleus,maxValueNeurite);
              maxValueWholeImage = max(maxValueWholeImage,maxValueWhole);
              %if(sixteenBit==0 && maxValueWhole>255)
              %  sixteenBit=1;
              %end
              %Check if in area x to x+512 and y to y+512 is an Nucleus. If
              %not, check for highest value for high pass filter.
              %if(sixteenBit == 1)
              r=isNucleusInArea(eventdata,handles,currentWellWithoutNull,x+1,x+512,y+1,y+512);
              if(r==0)
                   maxFilterList = [maxFilterList maxValueWhole];
              end
              
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
            neuriteImage = cat(1,neuriteImage,neuriteImageLine);
            nucleusImage = cat(1,nucleusImage,nucleusImageLine);
            end
            if(breaked==1)
                continue;
            end
            
            %If 16 Bit image, invert y Achses on image
    if(get(handles.rb4816,'Value') || get(handles.rb9616,'Value') || get(handles.rbWholeOt16,'Value') || get(handles.rbOtSingleChamber16,'Value'))
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
        standardDeviationMaxFilter = std(double(maxFilterList));
        meanMaxFilter = mean(maxFilterList);
        maxValue = meanMaxFilter + standardDeviationMaxFilter;
        set(handles.edMinBrightness,'String',maxValue);
        set(handles.edMaxBrightness,'String',maxValueWholeImage);
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
if(imageHandler.CheckSquareMode)
    hPixelInfoChildren = get(imageHandler.MousePosition,'Children');
    stringPix = get(hPixelInfoChildren(1),'String')
    %Extract X and Y Coordinate from stringPix
    tokenMatrix = regexpi(stringPix,'\((\d+),\s(\d+)\)','tokens');
    x = tokenMatrix{1}(1);
    y = tokenMatrix{1}(2);
    x=str2num(x{1});
    y=str2num(y{1});
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
    hPixelInfoChildren = get(imageHandler.MousePosition,'Children');
    stringPix = get(hPixelInfoChildren(1),'String')
    %Extract X and Y Coordinate from stringPix
    tokenMatrix = regexpi(stringPix,'\((\d+),\s(\d+)\)','tokens');
    x = tokenMatrix{1}(1);
    y = tokenMatrix{1}(2);
    x=str2num(x{1});
    y=str2num(y{1});
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
    hPixelInfoChildren = get(imageHandler.MousePosition,'Children');
    stringPix = get(hPixelInfoChildren(1),'String')
    %Extract X and Y Coordinate from stringPix
    tokenMatrix = regexpi(stringPix,'\((\d+),\s(\d+)\)','tokens');
    x = tokenMatrix{1}(1);
    y = tokenMatrix{1}(2);
    x=str2num(x{1});
    y=str2num(y{1});
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
    zoomOnZoomedImage = getrect;
    currentPicNuc = imageHandler.NucleusImage(zoomOnZoomedImage(2)*10:zoomOnZoomedImage(2)*10+zoomOnZoomedImage(4)*10, zoomOnZoomedImage(1)*10:zoomOnZoomedImage(1)*10+zoomOnZoomedImage(3)*10);
    currentPicNeu = imageHandler.NeuriteImage(zoomOnZoomedImage(2)*10:zoomOnZoomedImage(2)*10+zoomOnZoomedImage(4)*10, zoomOnZoomedImage(1)*10:zoomOnZoomedImage(1)*10+zoomOnZoomedImage(3)*10);
    minBrightness = get(handles.edMinBrightness,'String');
    maxBrightness = get(handles.edMaxBrightness,'String');
    minBrightness = str2num(minBrightness);
    maxBrightness = str2num(maxBrightness);
    currentPicNeu = double(currentPicNeu)./3996;
    currentPicNeu = imadjust(currentPicNeu, [double(minBrightness)/3996;double(maxBrightness)/3996], [0;1]);
    currentPicNeu = uint8(currentPicNeu.*255);   
    minBrightness = get(handles.edMinBrightnessNuc,'String');
    maxBrightness = get(handles.edMaxBrightnessNuc,'String');
    minBrightness = str2num(minBrightness);
    maxBrightness = str2num(maxBrightness);
    currentPicNuc = double(currentPicNuc)./3996;
    currentPicNuc = imadjust(currentPicNuc, [double(minBrightness)/3996;double(maxBrightness)/3996], [0;1]);
    currentPicNuc = uint8(currentPicNuc.*255); 
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
        plot((nucleusCols),(nucleusRows),'Linestyle','none','Marker','.','Markersize',15,'Color','blue');       
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
        plot((neuronCols),(neuronRows),'Linestyle','none','Marker','.','Markersize',20,'Color',[0.8,0.3,0]);       
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
            plot((nucleusCols./10),(nucleusRows./10),'Linestyle','none','Marker','.','Markersize',2,'Color','blue');       
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
    currentFilter{filterEntries+1} = createMask(h);
    currentFilter{filterEntries+2} = createMask(h2);
    currentFilter{filterEntries+1} = imresize(currentFilter{filterEntries+1},[rowNumber, colNumber]);
    currentFilter{filterEntries+2} = imresize(currentFilter{filterEntries+2},[rowNumber, colNumber]);
    lbFilterString = get(handles.lbFilter,'String');
    %For every positive entry in Filter List: Create entry in Listbox
    filterEntries = filterEntries+1;
    lbFilterString = [lbFilterString;num2str(filterEntries)];
    imageHandler.PositiveFilters = currentFilter;
    set(handles.lbFilter, 'String', lbFilterString);
    filterEntries = filterEntries+1;
elseif(get(handles.rbPolygon,'Value'))
    if(imageHandler.ZoomState)
        %ToDo
    else
        h = impoly(gca);
        [rowNumber, colNumber]=size(imageHandler.NeuriteImage);
        mask = createMask(h);
        %Map small mask to big picture
        mask = imresize(mask,[rowNumber, colNumber]);
        currentFilter{filterEntries+1} = mask;
        imageHandler.PositiveFilters = currentFilter;
        filterEntries = filterEntries+1;
    end
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
        [filterDistance nonFilterDistance SphereArea NucleusArea] = calculateMigrationDistance(selectedWell, handles);
        waitbar((60/100),waitbarHandle);
        [rowNumber, colNumber]=size(imageHandler.NeuriteImage);
        sizingFactorX = colNumber/1000;
        sizingFactorY = rowNumber/1000;        
        allowedDistance = 15*((sizingFactorX + sizingFactorY)/2);
        rowColMatSorted = sortMatrixByNearestNeighbour(SphereArea,allowedDistance,handles);
        waitbar((70/100),waitbarHandle);
        negativeRowColMatSorted = sortMatrixByNearestNeighbour(NucleusArea,allowedDistance,handles);
        waitbar((80/100),waitbarHandle);
        h = impoly(handles.axes2,rowColMatSorted);
        cbNegativeFilterChecked = get(handles.cbNegativeFilter ,'Value');
        if(cbNegativeFilterChecked)
            h2 = impoly(handles.axes2,negativeRowColMatSorted);            
            setColor(h2,'red');
            maskNegative = createMask(h2);
            maskNegative = imresize(maskNegative,[rowNumber, colNumber]);
            if(iscell(imageHandler.NegativeFilters) > 0)
                negFilterEntries = numel(imageHandler.NegativeFilters);
                currentNegativeFilter = imageHandler.NegativeFilters;
            else
                negFilterEntries = 0;
                currentNegativeFilter = cell(0);
            end
            currentNegativeFilter{negFilterEntries+1} = maskNegative;
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
end
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles)


function rowColMatSorted = sortMatrixByNearestNeighbour(SphereArea,allowedMaxDistance,handles)
imageHandler=handles.ImageHandler;
[NonZeroRows NonZeroCols] = find(SphereArea);
[rowNumber, colNumber]=size(imageHandler.NeuriteImage);
sizingFactorX = colNumber/1000;
sizingFactorY = rowNumber/1000;
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
    if(imageHandler.ZoomState)
        %ToDo
        h = impoly(gca);
        %Create invisible copy of h
        setColor(h,'red');
        hPositionMatrix = getPosition(h);        
        %TODO
        %Add Offset Position to hcopy
        xmin = imageHandler.ZoomState(1);
        ymin = imageHandler.ZoomState(2);
        width = imageHandler.ZoomState(3);
        height = imageHandler.ZoomState(4);     
        hcopy = hcopy + [xmin ymin];
        %Map zoomed mask to big picture
        currentFilter{filterEntries+1} = mask;
        imageHandler.NegativeFilters(selectedWell) = currentFilter;
        filterEntries = filterEntries + 1;
    else
        h = impoly(gca);
        setColor(h,'red');
        [rowNumber, colNumber]=size(imageHandler.NeuriteImage);
        mask = createMask(h);
        %Map small mask to big picture
        mask = imresize(mask,[rowNumber, colNumber]);
        currentFilter{filterEntries+1} = mask;
        imageHandler.NegativeFilters = currentFilter;
        filterEntries = filterEntries + 1;
    end
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
if(isfield(handles,'NeuronCoordinates'))
    neuronHandler = handles.NeuronCoordinates;
end
if(isfield(handles,'CSVCoordinates'))
    csvHandler = handles.CSVCoordinates;
end
[FileName,PathName] = uiputfile;
save('-v7.3',strcat(PathName,FileName),'imageHandler','neuronHandler','csvHandler');
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
handles.ImageHandler = imageHandler;
handles.NeuronCoordinates = neuronHandler;
handles.CSVCoordinates = csvHandler;
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
foldername = [imageHandler.Foldername '/ConvertedCellomics'];
allFiles = dir(foldername); 
waitbarHandle = waitbar(0,'Please wait. All images are converted to 8 Bit.');
for i=3:numel(allFiles)
    waitbar((i/numel(allFiles)),waitbarHandle);
    currentFile=allFiles(i).name;
    if(strfind(currentFile, 'Neurite'))
        minBrightness = get(handles.edMinBrightness,'String');
        maxBrightness = get(handles.edMaxBrightness,'String');
        minBrightness = str2num(minBrightness);
        maxBrightness = str2num(maxBrightness);
    else
        minBrightness = get(handles.edMinBrightnessNuc,'String');
        maxBrightness = get(handles.edMaxBrightnessNuc,'String');
        minBrightness = str2num(minBrightness);
        maxBrightness = str2num(maxBrightness);
    end
    imageString = strcat(foldername, '\', currentFile);
    currentPic = imread(imageString);
    deleteIndices = uint16(find(abs(currentPic) < minBrightness));
    currentPic(deleteIndices) = 0;
    
    currentPic = double(currentPic)./3996;
    currentPic = imadjust(currentPic, [double(minBrightness)/3996;double(maxBrightness)/3996], [0;1]);
    currentPic = uint8(currentPic.*255);    
    imwrite(currentPic,imageString);
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
minBrightness = get(handles.edMinBrightness,'String');
maxBrightness = get(handles.edMaxBrightness,'String');
minBrightness = str2num(minBrightness);
maxBrightness = str2num(maxBrightness);
currentPicNeu = imageHandler.ResizedNeuriteImage;
currentPicNeu = double(currentPicNeu)./3996;
currentPicNeu = imadjust(currentPicNeu, [double(minBrightness)/3996;double(maxBrightness)/3996], [0;1]);
currentPicNeu = uint8(currentPicNeu.*255);   
minBrightness = get(handles.edMinBrightnessNuc,'String');
maxBrightness = get(handles.edMaxBrightnessNuc,'String');
minBrightness = str2num(minBrightness);
maxBrightness = str2num(maxBrightness);
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
[result withoutFilterResult SphereArea NucleusArea]=calculateMigrationDistance(selectedWell,handles);
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
        [checkX, checkY]=MapPoint(x,y,sizeX,sizeY,100,100);
        if(DensityM(y,x) <= threshold && AreaResult(y,x) == 0 && x-1 > 0 && y-1 > 0 && y+1<= maxY && x+1 <= maxX && image(checkY,checkX) > 250)
            AreaResult(y,x) = 1;
            stack.push([x y+1]);
            stack.push([x-1 y]);
            stack.push([x y-1]);
            stack.push([x+1 y]);
        end
    else
        if(DensityM(y,x) >= threshold && AreaResult(y,x) == 0 && x-1 > 0 && y-1 > 0 && y+1<= maxY && x+1 <= maxX)
            AreaResult(y,x) = 1;
            stack.push([x y+1]);
            stack.push([x-1 y]);
            stack.push([x y-1]);
            stack.push([x+1 y]);
        end
    end
end
result=AreaResult;

function [x y] = MapPoint(xOriginal,yOriginal,sizeXOld,sizeYOld,sizeXNew,sizeYNew)
xFactor = sizeXOld/sizeXNew;
yFactor=sizeYOld/sizeYNew;
x=uint16(double(xOriginal)*xFactor);
y=uint16(double(yOriginal)*yFactor);

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
    plot((nucleusCols./10),(nucleusRows./10),'Linestyle','none','Marker','.','Markersize',2,'Color','blue');       
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

function [DensityM, DensityNeuron] = getDensityDistributionsFromCSV(selectedWell, handles)
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
NucleusM = csvHandler.CellPosMatrix;
NeuronM = neuronHandler.CellPosMatrix;
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
      selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end
if(isKey(NucleusM,selectedWell))
    NucleusM = NucleusM(selectedWell);
    NucleusM = logical(full(NucleusM));
    NeuronM = NeuronM(selectedWell);
    NeuronM = logical(full(NeuronM));
    [row col] = find(NucleusM);
    [neuronRow neuronCol] = find(NeuronM);
    DensityM = uint8(hist3([row col],[100 100]));
    DensityNeuron = uint8(hist3([neuronRow neuronCol],[50 50]));
%Map density matrix to values between o and 1
    minBrightness=0;
    maxBrightness=max(DensityM(:));
    DensityM = double(DensityM)./255;   
    DensityM = imadjust(DensityM, [double(minBrightness)/255;double(maxBrightness)/255], [0;1]);
    DensityM = uint8(DensityM.*255); 
    DensityNeuron = double(DensityNeuron)./255;   
    DensityNeuron = imadjust(DensityNeuron, [double(minBrightness)/255;double(maxBrightness)/255], [0;1]);
    DensityNeuron = uint8(DensityNeuron.*255);
end

function [filterDistance nonFilterDistance SphereArea NucleusArea] = calculateMigrationDistance(selectedWell, handles)
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
NucleusM = csvHandler.CellPosMatrix;
NeuronM = neuronHandler.CellPosMatrix;
SphereArea=0;
NucleusArea=0;
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
      selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
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

    [DensityM DensityNeuron] = getDensityDistributionsFromCSV(selectedWell, handles);
    %Find circle in Density Edge
    %Get point of Maximum Density
    [C,rowmaxarray]=max(DensityM);
    [maxvalue,xMax]=max(C);
    yMax=rowmaxarray(xMax);
    densityNeurons = DensityNeuron(uint8(yMax/2), uint8(xMax/2));
    i=0;
    %Validate point of max density: It's only valid, if Cellomics Neurons
    %are also in region
    while(densityNeurons < 5 && i<=5)
        DensityM(yMax,xMax) =0;
        [C,rowmaxarray]=max(DensityM);
        [maxvalue,xMax]=max(C);
        yMax=rowmaxarray(xMax);
        densityNeurons = DensityNeuron(uint8(yMax/2), uint8(xMax/2));
        i=i+1;
    end
    
    %Confirm point of maximum density with Neuron Matrix. If there are no
    %neurons in area:
    %Look for another point of max Density
    
    
%Look for points of zero Density near Maximum Density in a radius of x.
%Just go to center of the pic
%pixels
    radius=0;
    euclidDist=9999;    
    centerdist = pdist([xMax yMax; 50 50]);
    zeroPoints = uint8(zeros(0));
    while(radius < 40)
      radius=radius+1;
      for i=yMax-radius:yMax+radius
        for j=xMax-radius:xMax+radius
        %Check if point is in bounds.
           if(i>0 && j<=100 && i<=100 && j>0)
              %Check if match and if distance to center has decreased.
              %Verify if regarding point in Neurite picture is white 
               currentDist = pdist([j i;50 50]);
               [picX picY] = MapPoint(j,i,sizeX,sizeY,100,100);                                        
               if(DensityM(i,j) == 0 && currentDist<(centerdist+40) && ~ismember([j i],zeroPoints, 'rows') && imageHandler.NucleusImage(picY,picX) > 250)
                  zeroPoints = [zeroPoints;j i];
               end %if
           end %if                                                                                                                                                
        end %for
      end %for
    end %while
%Flood fill low density points from zero point and count them as area for
%cell nucleus

    NucleusArea = logical(zeros(100,100));
    if(numel(zeroPoints) > 0)
     for i=1:numel(zeroPoints(:,1))
      stack = java.util.Stack();
      neuX = zeroPoints(i,1);
      neuY = zeroPoints(i,2);
      stack.push([neuX neuY]);
      NucleusArea = logical(logical(NucleusArea) + logical(FloodFillIterative(DensityM,NucleusArea,0,0,100,100,stack,imageHandler.NucleusImage)));
     end
     
     NucleusArea = NucleusArea-1;
     NucleusArea = NucleusArea.*(-1);
     NucleusArea = edge(NucleusArea);
     
     boundarieObjects = bwboundaries(NucleusArea);
     NucleusArea = logical(zeros(100,100));
%Take only boundarieobjects with more than points into account
%Draw white pixels on black image for relevant boundary objects
     for i=1:numel(boundarieObjects)
      currentBoundarieObject = boundarieObjects{i};
      if(numel(currentBoundarieObject) >= 20)
        for j=1:size(currentBoundarieObject,1)
            yBound = currentBoundarieObject(j,1);
            xBound = currentBoundarieObject(j,2);
            NucleusArea(yBound,xBound) = 1;
        end
      end
     end
%Flood fill everything from max Density above Threshold
%Outer circle is area for Sphere
     stack = java.util.Stack();
     stack.push([xMax yMax]);
     SphereArea = logical(zeros(100,100));
     SphereArea = FloodFillIterative(DensityM,SphereArea,1,1,100,100,stack,0);
     SphereArea=imfill(SphereArea,'holes');
     SphereArea = edge(SphereArea);
     MigrationPic = SphereArea + NucleusArea;
     figure(1);
     imshow(MigrationPic);
    %hPixelInfo = impixelinfo;        
    %imageHandler.MousePosition = hPixelInfo;
%Für weitere Verarbeitung:
% Berechne den Schwerpunkt aller weißen Pixel in NucleusArea
     [row,col] = find(NucleusArea);
     sY = sum(row)/length(row);
     sX = sum(col)/length(col);
% Gehe jeweils im 10 Grad Winkel von dort aus los.
     mediumDistanceList = zeros(0);
     mediumDistanceListInInterval = zeros(0);
     successCount=0;
     for i=pi/32:pi/32:2*pi
      deltaX = cos(i);
      deltaY = sin(i);
      xNucleus = 0;
      yNucleus = 0;
      xOuterSphere = 0;
      yOuterSphere = 0;    
      x=sX;
      y=sY;
      while(xOuterSphere ==0 && uint8(x) < 98 && uint8(y) < 98 && uint8(x) > 2 && uint8(y) > 2)
        if(NucleusArea(uint8(y),uint8(x)) == 1)
            [xNucleus yNucleus]=MapPoint(x,y,sizeX,sizeY,100,100);
        end
        if(SphereArea(uint8(y),uint8(x)) == 1)
            [xOuterSphere yOuterSphere] = MapPoint(x,y,sizeX,sizeY,100,100);
        end
        x=x+deltaX;
        y=y+deltaY;
      end
      if(xOuterSphere~=0 && yOuterSphere~=0 && xNucleus ~= 0 && yNucleus ~= 0)
        successCount = successCount+1;
        mediumDistanceList = [mediumDistanceList pdist([xOuterSphere yOuterSphere;xNucleus yNucleus])];        
        %mediumDistance = mediumDistance + pdist([xOuterSphere yOuterSphere;xNucleus yNucleus]);
      end
     end
     distancestddev = std(mediumDistanceList);
     mediumDistance = mean(mediumDistanceList);
     maxDistance = max(mediumDistanceList);
     for i=1:numel(mediumDistanceList)
      currentDistance = mediumDistanceList(i);
    %Check if current distance is in Interval
      if(currentDistance >= 150 && currentDistance >= (maxDistance-(2*distancestddev))) % (mediumDistance - (2 * distancestddev)) && currentDistance <= (mediumDistance +(2*distancestddev)))
        mediumDistanceListInInterval = [mediumDistanceListInInterval currentDistance];
      end
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
    else
        filterDistance=-1;
        nonFilterDistance=-1;
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
[FileName,PathName] = uiputfile('migdist.csv');
fileID = fopen(strcat(PathName,'\',FileName),'w');
handles = guidata(handles.figure1);
wellList = get(handles.lbWell, 'string');
waitbarHandle = waitbar(0,'Calculating and exporting Migration Distances for all wells');
selectedWellNumber = get(handles.lbWell,'Value');
exportText = sprintf('Well;Migration Distance(px);Migration Distance without filter (px)\r\n');
for currentWellIndex=1:numel(wellList)
    currentWell = wellList(currentWellIndex);
    currentWell = currentWell{1};
    disp(currentWell);
    [distFilter, distWithoutFilter, SphereArea, NucleusArea]=calculateMigrationDistance(currentWell,handles);
    if(distFilter >-1)
        distFilter = num2str(distFilter,'%f');
        distFilter=strrep(distFilter,'.',',');
        distWithoutFilter = num2str(distWithoutFilter,'%f');
        distWithoutFilter=strrep(distWithoutFilter,'.',',');
        exportText = [exportText currentWell ';' distFilter ';' distWithoutFilter sprintf('\r\n')];
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
exportText = sprintf('Well;Filter;#Nuclei by Cellomics;#Neurons Cellomics;#Neurons by Clicks;%%Neurons by Cellomics;%%Neurons by Clicks;Migration Distance (px)\r\n');
for i=1:numel(wellList)
    selectedWell = wellList{i};
    [distFilter, distWithoutFilter, SphereArea, NucleusArea]=calculateMigrationDistance(selectedWell,handles);
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
            exportText = [exportText selectedWell ';' num2str(j) ';' numberNuclei ';' numberNeurons ';' numberManual ';' percentageCellomics ';' percentageManual ';' distFilter  sprintf('\r\n')];
        end
    end
    waitbar((i/numel(wellList)),waitbarHandle);
end
%Write exportText to File
fprintf(fileID,'%s',exportText);
close(waitbarHandle);

% --- Executes on button press in cbNegativeFilter.
function cbNegativeFilter_Callback(hObject, eventdata, handles)
% hObject    handle to cbNegativeFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbNegativeFilter