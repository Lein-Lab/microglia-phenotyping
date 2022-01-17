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
%      existing singleton*.  Starting from the  left, property value pairs are
%      applied to the GUI before NeuronDetectorGUISmall_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NeuronDetectorGUISmall_OpeningFcn via varargin.
%
%      *See GUI Opt on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NeuronDetectorGUISmall

% Last Modified by GUIDE v2.5 21-Dec-2013 12:10:15

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
    elseif(get(handles.rb96484,'Value'))
        csvHandler.ReadCSVFile(filepath, WellType.Well96484);
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

% --------------------------------------------------------------------
function LoadNeuronCSV2_Callback(hObject, eventdata, handles)
% hObject    handle to LoadNeuronCSV2 (see GCBO)
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
    WellNeuronDict = CalcNucleusNeuronDist(handles, neuronHandler.CellPosMatrix);
    neuronHandler.CellPosMatrix = WellNeuronDict;
    handles.NeuronCoordinates2 = neuronHandler;
    guidata(handles.figure1, handles)
    
end

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
    WellNeuronDict = CalcNucleusNeuronDist(handles, neuronHandler.CellPosMatrix);
    neuronHandler.CellPosMatrix = WellNeuronDict;
    handles.NeuronCoordinates = neuronHandler;
    guidata(handles.figure1, handles)
end


% --------------------------------------------------------------------
function FixNeuronPositions_Callback(hObject, eventdata, handles)
% hObject    handle to FixNeuronPositions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
FixNeuronPositions(hObject, eventdata, handles);
guidata(handles.figure1, handles)

% --------------------------------------------------------------------
% Calculates for every point in NeuronCoordinates its regarding point in
% CSVCoordinates and delivers medium distance
function FixNeuronPositions(hObject, eventdata, handles)
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

          [xMapped,yMapped,distance] = neuronHandler.FindNucleusForNeuron(col,row,NucleusM, 100);
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
% Calculates for every point in NeuronCoordinates its regarding point in
% CSVCoordinates and delivers medium distance
function WellNeuronDict = CalcNucleusNeuronDist(handles, WellNeuronDict)
% hObject    handle to CalcNucleusNeuronDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(handles.figure1);
if(isfield(handles,'NeuronCoordinates') && isfield(handles,'CSVCoordinates'))
           %Map Neuron Positions to Nucleus Positions
    selectedWellNumber = get(handles.lbWell,'Value');
    wellList = get(handles.lbWell, 'string');
    csvHandler = handles.CSVCoordinates;
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

          [xMapped,yMapped,distance] = csvHandler.FindNucleusForNeuron(col,row,NucleusM,100);
          distsum=distsum+distance;
          NeuronM(row,col) = 0;
          NeuronM(yMapped,xMapped) = 1;
       end
       WellNeuronDict(selectedWell) = NeuronM;
       distsum=distsum/numel(nucleusRows);
       disp('Medium distance is ');
       disp(distsum);   
      end
    end
    neuronHandler.CellPosMatrix = WellNeuronDict;
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
         elseif(get(handles.rb96484,'Value'))
         maxIndex = numel(csvHandler.Well96Offsets484);
         wellType = WellType.Well96484;
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


              maxValueNucleus = uint16(max(currentNucleusPic(:)));
              imageString = strcat(foldername, '\', filenameneu);
              if(exist(imageString,'file') && ~ isempty(filenameneu))
                  currentNeuritePic = uint16(imread(imageString));
              else
                currentNeuritePic = uint16(zeros(512,512));
              end
              maxValueNeurite = uint16(max(currentNeuritePic(:)));
              maxValueWhole = uint16(max(maxValueNucleus,maxValueNeurite));
              maxValueWholeImage = uint16(max(uint16(maxValueWholeImage),maxValueWhole));
              %if(sixteenBit==0 && maxValueWhole>255)
              %  sixteenBit=1;
              %end
              %Check if in area x to x+512 and y to y+512 is an Nucleus. If
              %not, check for highest value for high pass filter.
              %if(sixteenBit == 1)
              %r=isNucleusInArea(eventdata,handles,currentWellWithoutNull,x+1,x+512,y+1,y+512);
          %Wof√ºr?
              %    if(r==0)

          %         maxFilterList = [maxFilterList maxValueWhole];
          %    end

              
              neuriteImageLine = uint16(cat(2,neuriteImageLine,currentNeuritePic));
              nucleusImageLine = uint16(cat(2,nucleusImageLine,currentNucleusPic));
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
                neuriteImageLine = uint16(flipdim(neuriteImageLine,1));
                nucleusImageLine = uint16(flipdim(nucleusImageLine,1));
            end
            neuriteImage = uint16(cat(1,neuriteImage,neuriteImageLine));
            nucleusImage = uint16(cat(1,nucleusImage,nucleusImageLine));
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
        %Wof√ºr?
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
    [xNuc,yNuc,distance] = csvHandler.FindNucleusForNeuron(xMapped,yMapped,NucleusM,30);
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
elseif(imageHandler.ZoomState)
    %Left click OK
    %Check if csvHandler.ManualPositions 1-4 exist
    clickType = get(gcf,'SelectionType');
    %Get state of 'ManualButtonGroup'
    if(numel(neuronHandler.ManualNeuronPositionsSparse) ==0) 
        neuronHandler.ManualNeuronPositionsSparse = containers.Map()
    end
    if(numel(csvHandler.ManualPositions1) ==0)        
        csvHandler.ManualPositions1 = containers.Map();
        csvHandler.ManualPositions2 = containers.Map();
        csvHandler.ManualPositions3 = containers.Map();
        csvHandler.ManualPositions4 = containers.Map();        
    end    
    if(isKey(csvHandler.ManualPositions1,selectedWell))
            currentNuclei1 = csvHandler.ManualPositions1(selectedWell);
            currentNuclei2 = csvHandler.ManualPositions2(selectedWell);
            currentNuclei3 = csvHandler.ManualPositions3(selectedWell);
            currentNuclei4 = csvHandler.ManualPositions4(selectedWell);
    else            
            [sizeY sizeX]=size(imageHandler.NeuriteImage);
            currentNuclei1 = sparse(sizeY,sizeX);
            currentNuclei2 = sparse(sizeY,sizeX);
            currentNuclei3 = sparse(sizeY,sizeX);
            currentNuclei4 = sparse(sizeY,sizeX);
            csvHandler.ManualPositions1(selectedWell) = currentNuclei1;
            csvHandler.ManualPositions2(selectedWell) = currentNuclei2;
            csvHandler.ManualPositions3(selectedWell) = currentNuclei3;
            csvHandler.ManualPositions4(selectedWell) = currentNuclei4;
    end  
    if(isKey(csvHandler.ManualNeuronPositionsSparse,selectedWell))
        currentNucleiMain = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
    else
        [sizeY sizeX]=size(imageHandler.NeuriteImage);
        currentNucleiMain = sparse(sizeY,sizeX);
        neuronHandler.ManualNeuronPositionsSparse(selectedWell) = currentNucleiMain;
    end
    if(get(handles.rbManual1,'Value')) 
        manualList=currentNuclei1;        
        %Add clicked position to manualList
    elseif(get(handles.rbManual2,'Value'))
        manualList=currentNuclei2;
    elseif(get(handles.rbManual3,'Value'))
        manualList=currentNuclei3;
    elseif(get(handles.rbManual4,'Value'))
        manualList=currentNuclei4;
    elseif(get(handles.rbManualM, 'Value'))
        manualList = currentNucleiMain;
    end
    xMin = imageHandler.ZoomState(1);
    yMin = imageHandler.ZoomState(2);
    xMapped = int32(x + xMin);
    yMapped = int32(y + yMin);
    WellNucleusDict = csvHandler.CellPosMatrix;
    NucleusM = WellNucleusDict(selectedWell);
    [xNuc,yNuc,distance] = csvHandler.FindNucleusForNeuron(xMapped,yMapped,NucleusM,30);
    
    if(numel(clickType)==6)
        %Save current Position
        manualList(yNuc,xNuc) = 1;
    else
        %Remove current Position
        manualList(yNuc,xNuc) = 0;
    end
    
    if(get(handles.rbManual1,'Value'))              
        csvHandler.ManualPositions1(selectedWell) = manualList;
        color = 'green';
        %Add clicked position to manualList
    elseif(get(handles.rbManual2,'Value'))
        csvHandler.ManualPositions2(selectedWell) = manualList;
        color='white';
    elseif(get(handles.rbManual3,'Value'))
        csvHandler.ManualPositions3(selectedWell) = manualList;
        color = 'cyan';
    elseif(get(handles.rbManual4,'Value'))
        csvHandler.ManualPositions4(selectedWell) = manualList;
        color = 'yellow';
    elseif(get(handles.rbManualM, 'Value'))
        neuronHandler.ManualNeuronPositionsSparse(selectedWell) = manualList;
        color = 'magenta';
    end
    markerX = xNuc-xMin;
    markerY = yNuc-yMin;    
    hold on;
    if(numel(clickType)==6)
        plot(markerX,markerY,'Linestyle','none','Marker','.','Markersize',20, 'Color', color);
    else
        plot(markerX,markerY,'Linestyle','none','Marker','.','Markersize',20,'Color','blue');
        %set(p,'Color','red');
    end
    hold off;
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
    newRect = getrect;
    xmin = imageHandler.ZoomState(1)+newRect(1);
    ymin = imageHandler.ZoomState(2)+newRect(2); 
    width = imageHandler.ZoomState(3)-(imageHandler.ZoomState(3)-newRect(3));
    height = imageHandler.ZoomState(4)-(imageHandler.ZoomState(4)-newRect(4));
    zoomOnZoomedImage = [xmin ymin width height];
    subNucImage = imageHandler.NucleusImage(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
    subNeuImage = imageHandler.NeuriteImage(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
    if(~numel(imageHandler.NucleusImage) == 0)              
            mat = imfuse(subNucImage,subNeuImage);
            hImage = imshow(mat);        
            set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
            hPixelInfo = impixelinfo;        
            imageHandler.MousePosition = hPixelInfo;
    end
    %Set new ZoomState
    imageHandler.ZoomState = zoomOnZoomedImage;
    handles.ImageHandler = imageHandler;
    guidata(handles.figure1, handles)
    RefreshZoomedImage(handles);   
else
    csvHandler = handles.CSVCoordinates;
    if(isfield(handles,'NeuronCoordinates'))
        neuronHandler = handles.NeuronCoordinates;
    end
    if(isfield(handles,'NeuronCoordinates2'))
        neuronHandler2 = handles.NeuronCoordinates2;
    end
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
    set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
    hPixelInfo = impixelinfo;        
    imageHandler.MousePosition = hPixelInfo;
    %Load selected Matrices to GUI
    RefreshZoomedImage(handles);
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
        imageHandler.PositiveFilters = 0;
        imageHandler.NegativeFilters = 0;
    end
    if(length(selectedWell) == 2)
      selectedWell = [selectedWell(1) '0' selectedWell(2)];
    elseif(length(selectedWell) == 4)
        selectedWell = [selectedWell(1) '0' selectedWell(2:length(selectedWell))];
    end
    
    %Get images from image folder
    foldername = [imageHandler.Foldername '/ConvertedCellomics'];
    imagePathNucleusBig = [foldername '/' selectedWell 'NucleusBig.tif'];
    imagePathNucleusSmall = [foldername '/' selectedWell 'NucleusSmall.tif'];
    imagePathNeuriteBig = [foldername '/' selectedWell 'NeuriteBig.tif'];
    imagePathNeuriteSmall = [foldername '/' selectedWell 'NeuriteSmall.tif'];    
    imagePathSkeleton = [foldername '/' selectedWell 'Skeleton.tif'];
    imagePathBinary = [foldername '/' selectedWell 'Binary.tif'];
    imageHandler.ResizedNeuriteImage = imread(imagePathNeuriteSmall);
    imageHandler.ResizedNucleusImage = imread(imagePathNucleusSmall);    
    imageHandler.NeuriteImage = imread(imagePathNeuriteBig);
    imageHandler.NucleusImage = imread(imagePathNucleusBig);
    
    if(exist(imagePathBinary,'file'))
        imageHandler.BinaryImage = imread(imagePathBinary);
        imageHandler.ResizedBinaryImage = imresize(imageHandler.BinaryImage,0.1);
    end
    if(exist(imagePathSkeleton,'file'))
        imageHandler.SkeletonImage = imread(imagePathSkeleton);
        imageHandler.ResizedSkeletonImage = imresize(imageHandler.SkeletonImage,0.1);
    end
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
    %A2 = get(handles.tbA,'String');
    %B2 = get(handles.tbB,'String');
    %if(imageHandler.ZoomState)
    %    h2 = imrect(gca,[1,1,(str2num(B)), str2num(A)])
    %else
    %    h2 = imrect(gca,[1,1,(str2num(B)/10),(str2num(A)/10)]);    
    %end
    %setResizable(h2, 0);
    posH = wait(h);
    %posH2 = wait(h2);
    if(imageHandler.ZoomState)
        xmin = imageHandler.ZoomState(1);
        ymin = imageHandler.ZoomState(2); 
        polyPos = getPosition(h);
        %polyPos2 = getPosition(h2);
        %Wandle alle Positionen in Position f¸r das groﬂe Bild um        
        currentFilter{filterEntries+1} = poly2mask(polyPos(:,1)+xmin,polyPos(:,2)+ymin,colNumber, rowNumber);
        %currentFilter{filterEntries+2} = poly2mask(polyPos2(:,1)+xmin,polyPos2(:,2)+ymin,colNumber, rowNumber);
    else
        currentFilter{filterEntries+1} = createMask(h);
        %currentFilter{filterEntries+2} = createMask(h2);
        currentFilter{filterEntries+1} = imresize(currentFilter{filterEntries+1},[rowNumber, colNumber]);
        %currentFilter{filterEntries+2} = imresize(currentFilter{filterEntries+2},[rowNumber, colNumber]);
    end
    %lbFilterString = get(handles.lbFilter,'String');
    %For every positive entry in Filter List: Create entry in Listbox
    filterEntries = filterEntries+1;
    %lbFilterString = [lbFilterString;num2str(filterEntries)];
    imageHandler.PositiveFilters = currentFilter;
    %set(handles.lbFilter, 'String', lbFilterString);
    %filterEntries = filterEntries+1;
elseif(get(handles.rbPolygon,'Value'))
    h = impoly(gca);
    [rowNumber, colNumber]=size(imageHandler.NeuriteImage);
    if(imageHandler.ZoomState)
        xmin = imageHandler.ZoomState(1);
        ymin = imageHandler.ZoomState(2);  
        polyPos = getPosition(h);
        %Wandle alle Positionen in Position f¸r das groﬂe Bild um        
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
        %Wandle alle Positionen in Position f¸r das groﬂe Bild um        
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
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
path = [imageHandler.Foldername '/ConvertedCellomics/' selectedWell '.mat'];
if(exist(path,'file'))
        load(path);
        imageHandler.PositiveFilters = FMask.PositiveFilters;
        imageHandler.NegativeFilters = FMask.NegativeFilters;
end
FMask=0;
if(isfield(handles,'NeuronCoordinates'))
    neuronHandler = handles.NeuronCoordinates;
end
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
filterList = get(handles.lbFilter, 'string');
filterNumber = get(handles.lbFilter,'Value');
%Iterate over all filters
filterList = imageHandler.PositiveFilters;
filterMask = filterList(filterNumber);
imageHandler.PositiveFilters = 0;
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
imageHandler.NegativeFilters=0;
negativeMask = 0;
[filterSizeY filterSizeX] = size(filterMask);
%Multiply filterMask with Neuron and Nucleus Matrices
NucleusM = csvHandler.CellPosMatrix(selectedWell);
if(isfield(handles,'NeuronCoordinates'))
    NeuronM = neuronHandler.CellPosMatrix(selectedWell);
    NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
    [NeuronMSizeY NeuronMSizeX] = size(NeuronM);
    [NeuronManualMSizeY NeuronManualMSizeX] = size(NeuronManualM);
    
    if(NeuronMSizeX > filterSizeX && exist('neuronHandler','var'))
     NeuronM(:,filterSizeX+1:NeuronMSizeX) = [];
     neuronHandler.CellPosMatrix(selectedWell) = NeuronM;
    end
    if(NeuronManualMSizeX > filterSizeX)
     NeuronManualM(:,filterSizeX+1:NeuronManualMSizeX) = [];
     neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
    end
    if(NeuronMSizeY > filterSizeY && exist('neuronHandler','var'))
     NeuronM(filterSizeY+1:NeuronMSizeY,:) = [];
     neuronHandler.CellPosMatrix(selectedWell) = NeuronM;
    end
    if(NeuronManualMSizeY > filterSizeY && exist('neuronHandler','var'))
     NeuronManualM(filterSizeY+1:NeuronManualMSizeY,:) = [];
     neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
    end
end



[NucleusMSizeY NucleusMSizeX] = size(NucleusM);
if(NucleusMSizeX > filterSizeX)
     NucleusM(:,filterSizeX+1:NucleusMSizeX) = [];
     csvHandler.CellPosMatrix(selectedWell) = NucleusM;
end

if(NucleusMSizeY > filterSizeY)
     NucleusM(filterSizeY+1:NucleusMSizeY,:) = [];
     csvHandler.CellPosMatrix(selectedWell) = NucleusM;
end


%%%PREVENT OUT OF MEMORY
filterMaskSmall = imresize(filterMask,0.1);    
     nucImageResized = double(imageHandler.ResizedNucleusImage).*double(filterMaskSmall);
    
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
     
    
    %Better: 
     currentWellNucMatrix = NucleusM(minRow:maxRow,minCol:maxCol);
     if (exist('neuronHandler','var'))  
        neuImage = uint16(imageHandler.NeuriteImage(minRow:maxRow,minCol:maxCol));
        currentWellNeuMatrix = NeuronM(minRow:maxRow,minCol:maxCol);
        currentWellNeuManualMatrix = NeuronManualM(minRow:maxRow,minCol:maxCol);
     end
     
     filterMaskResized = filterMask(minRow:maxRow,minCol:maxCol);
     
     currentWellNucMatrix = logical(logical(currentWellNucMatrix) .* logical(filterMaskResized));
     
     nucImage = uint16(uint16(nucImage) .* uint16(filterMaskResized));
     if (exist('neuronHandler','var'))  
        currentWellNeuMatrix = logical(logical(currentWellNeuMatrix) .* logical(filterMaskResized));
        neuImage = uint16(uint16(neuImage) .* uint16(filterMaskResized));
        neuImage = vertcat(neuImage,zeros(7168-size(neuImage,1),size(neuImage,2)));
        neuImage = horzcat(neuImage,zeros(size(neuImage,1),7168-size(neuImage,2)));
     end
     
     nucImage = vertcat(nucImage,zeros(7168-size(nucImage,1),size(nucImage,2)));
     nucImage = horzcat(nucImage,zeros(size(nucImage,1),7168-size(nucImage,2)));




numberNuclei = nnz(currentWellNucMatrix);
NucleusM = 0;
numberManual = 0;
numberNeurons=0;
if (exist('neuronHandler','var'))    
    numberNeurons = nnz(currentWellNeuMatrix);
    NeuronM=0;
    numberManual = nnz(currentWellNeuManualMatrix);
end
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
if (exist('neuronHandler','var'))  
    handles.NeuronCoordinates = neuronHandler;
end
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
function lbFilter_Callback(hObject, ~, handles)
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
    
      %currentPic = uint8(currentPic);
      currentPic = double(currentPic)./4095;
      currentPic = imadjust(currentPic, [double(minBrightness)/4095;double(maxBrightness)/4095], [0;1]);
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
currentPicNeu = double(currentPicNeu)./4095;
currentPicNeu = imadjust(currentPicNeu, [double(minBrightness)/4095;double(maxBrightness)/4095], [0;1]);
currentPicNeu = uint8(currentPicNeu.*255);   
minBrightness = optionHandler.HistogramMinNucleus;
maxBrightness = optionHandler.HistogramMaxNucleus;
currentPicNuc = imageHandler.ResizedNucleusImage;
currentPicNuc = double(currentPicNuc)./4095;
currentPicNuc = imadjust(currentPicNuc, [double(minBrightness)/4095;double(maxBrightness)/4095], [0;1]);
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
[result withoutFilterResult SphereArea NucleusArea]=calculateMigrationDistance(selectedWell,32,32,handles);
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

function RefreshZoomedImage(handles)
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
csvHandler = handles.CSVCoordinates;
imageHandler = handles.ImageHandler;
if(isfield(handles,'NeuronCoordinates'))
    neuronHandler = handles.NeuronCoordinates;
end
if(isfield(handles,'NeuronCoordinates2'))
    neuronHandler2 = handles.NeuronCoordinates2;
end
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
%Check which checkboxes are checked.
cbNucleusPictureChecked = get(handles.cbNucleusPicture ,'Value');
cbNeuritePictureChecked = get(handles.cbNeuritePicture ,'Value');
cbSkeletonPicChecked = get(handles.cbSkeletonPic ,'Value');
cbBinaryPicChecked = get(handles.cbBinaryPic , 'Value');
cbRemoveCoreChecked = get(handles.cbRemoveCore, 'Value');


%Analog to RefreshUnzoomed image but care about cutting Big Pictures before
zoomOnZoomedImage = imageHandler.ZoomState;
currentPicNuc = imageHandler.NucleusImage;
currentPicNeu = imageHandler.NeuriteImage;
binaryImage = imageHandler.BinaryImage;
skelImage = imageHandler.SkeletonImage;
if(cbRemoveCoreChecked)
    imMap = containers.Map();
    imMap('1') = currentPicNuc;
    imMap('2') = currentPicNeu;
    imMap('3') = binaryImage;
    %imMap('4') = skelImage;
    imMap = CutOutCircles(imMap,selectedWell,1, 0,handles);
    currentPicNuc = imMap('1');
    currentPicNeu = imMap('2');
    binaryImage = imMap('3');
    %skelImage = imMap('4');
end
currentPicNuc = currentPicNuc(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
currentPicNeu = currentPicNeu(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
if(numel(binaryImage)>0)
    binaryImage = binaryImage(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
end
if(numel(skelImage)>0)
    skelImage = skelImage(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
end
WellNucleusDict = csvHandler.CellPosMatrix;
NucleusM = WellNucleusDict(selectedWell);
NucleusM = NucleusM(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
WellNeuronDict = neuronHandler.CellPosMatrix;
NeuronM = WellNeuronDict(selectedWell);
NeuronM = NeuronM(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));


WellSkeletonDict = neuronHandler.NeuronPositionsSkeletonization;

if(numel(WellSkeletonDict) > 0 && isKey(WellSkeletonDict,selectedWell))
    SkeletonM = WellSkeletonDict(selectedWell);
    SkeletonM = SkeletonM(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
end
if(isfield(handles,'NeuronCoordinates2'))
 WellNeuronDict2 = neuronHandler2.CellPosMatrix;
 NeuronM2 = WellNeuronDict2(selectedWell);
 NeuronM2 = NeuronM2(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
end

WellManualDict = neuronHandler.ManualNeuronPositionsSparse;
ManualM = WellManualDict(selectedWell);
ManualM = ManualM(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
%WellEdgeCompDict = neuronHandler.NeuronPositionsEdgeComposite;
WellEdgeFilledDict = neuronHandler.NeuronPositionsEdgeFill;
%EdgeCompM = WellEdgeCompDict(selectedWell);
%EdgeCompM = EdgeCompM(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
EdgeFilledM = WellEdgeFilledDict(selectedWell);
EdgeFilledM = EdgeFilledM(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
WellDeletedDict = neuronHandler.NeuronPositionsSkelDeleted;

if(numel(WellDeletedDict) > 0 && isKey(WellDeletedDict,selectedWell))
    DeletedM = WellDeletedDict(selectedWell);
    DeletedM = DeletedM(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
end
if(cbSkeletonPicChecked)
    imshow(skelImage);
elseif(cbNeuritePictureChecked && cbBinaryPicChecked)
    if(~numel(imageHandler.NucleusImage) == 0 && ~numel(imageHandler.ResizedNeuriteImage) == 0)              
         mat = imfuse(imadjust(currentPicNeu,neuronHandler.StretchlimResult,[0 1]),binaryImage);
         hImage = imshow(mat);  
         set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});         
         hPixelInfo = impixelinfo;        
         imageHandler.MousePosition = hPixelInfo;
    end
elseif(cbNucleusPictureChecked && cbBinaryPicChecked)
    if(~numel(imageHandler.NucleusImage) == 0)              
         mat = imfuse(imadjust(currentPicNuc,neuronHandler.StretchlimResult,[0 1]),binaryImage);
         hImage = imshow(mat);  
         set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});         
         hPixelInfo = impixelinfo;        
         imageHandler.MousePosition = hPixelInfo;
    end
elseif(cbBinaryPicChecked)
    imshow(binaryImage);
elseif(cbNucleusPictureChecked && cbNeuritePictureChecked)    
    if(~numel(currentPicNuc) == 0 && ~numel(imageHandler.ResizedNeuriteImage) == 0)              
         mat = imfuse(currentPicNuc,currentPicNeu);
         hImage = imshow(mat);  
         set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});         
         hPixelInfo = impixelinfo;        
         imageHandler.MousePosition = hPixelInfo;
    end
elseif(cbNucleusPictureChecked)
    if(~numel(imageHandler.NucleusImage) == 0 && ~numel(imageHandler.ResizedNeuriteImage) == 0)              
         if(isfield(handles,'NeuronCoordinates') && ~isempty(neuronHandler.StretchlimResult))
             hImage = imshow(imadjust(currentPicNuc,neuronHandler.StretchlimResult,[0 1]));
         else
             hImage = imshow(currentPicNuc);
         end
         set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
         hPixelInfo = impixelinfo;        
         imageHandler.MousePosition = hPixelInfo;
    end
elseif(cbNeuritePictureChecked)
    if(~numel(imageHandler.NeuriteImage) == 0 && ~numel(imageHandler.ResizedNeuriteImage) == 0)              
         if(~isempty(neuronHandler.StretchlimResult))
            hImage = imshow(imadjust(currentPicNeu,neuronHandler.StretchlimResult,[0 1]));  
         else
             hImage = imshow(currentPicNeu); 
         end
         set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
         hPixelInfo = impixelinfo;        
         imageHandler.MousePosition = hPixelInfo;
    end
else
    %Check if one of the Pictures is selected. If not prepare plot!
end

%New:
listboxstrings = get(handles.popupPlus1);
listboxindex = get(handles.popupPlus1, 'value');
listboxstring = listboxstrings.String(listboxindex); 
MatrixPlus1 = GetMatrixFromListboxString(listboxstring,handles);
listboxindex = get(handles.popupPlus2, 'value');
listboxstring = listboxstrings.String(listboxindex); 
MatrixPlus2 = GetMatrixFromListboxString(listboxstring,handles);
listboxindex = get(handles.popupMinus1, 'value');
listboxstring = listboxstrings.String(listboxindex); 
MatrixMinus1 = GetMatrixFromListboxString(listboxstring,handles);
listboxindex = get(handles.popupMinus2, 'value');
listboxstring = listboxstrings.String(listboxindex); 
MatrixMinus2 = GetMatrixFromListboxString(listboxstring,handles);
listboxindex = get(handles.popupMinus3, 'value');
listboxstring = listboxstrings.String(listboxindex); 
MatrixMinus3 = GetMatrixFromListboxString(listboxstring,handles);
if(MatrixPlus1 ~= -1)
    MatrixPlus1 = MatrixPlus1(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
    if(MatrixMinus1 ~= -1)
        MatrixMinus1 = MatrixMinus1(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
        MatrixPlus1 = MatrixPlus1 - MatrixMinus1;
        ind = MatrixPlus1<0;
        MatrixPlus1(ind) = 0;
    end
    if(MatrixMinus2 ~= -1)
        MatrixMinus2 = MatrixMinus2(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
        MatrixPlus1 = MatrixPlus1 - MatrixMinus2;
        ind = MatrixPlus1<0;
        MatrixPlus1(ind) = 0;
    end
    if(MatrixMinus3 ~= -1)
        MatrixMinus3 = MatrixMinus3(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
        MatrixPlus1 = MatrixPlus1 - MatrixMinus3;
        ind = MatrixPlus1<0;
        MatrixPlus1(ind) = 0;
    end
    [nucleusRows, nucleusCols] = find(MatrixPlus1);   
     hold on;
     plot((nucleusCols),(nucleusRows),'Linestyle','none','Marker','.','Markersize',20,'Color','red');       
     hold off;
end
if(MatrixPlus2 ~= -1)
    MatrixPlus2 = MatrixPlus2(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));    
    [nucleusRows, nucleusCols] = find(MatrixPlus2);   
     hold on;
     plot((nucleusCols),(nucleusRows),'Linestyle','none','Marker','.','Markersize',20,'Color','blue');       
     hold off;
end


function mat=GetMatrixFromListboxString(listboxstring,handles)
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
listboxstring=listboxstring{1};
if(isfield(handles,'NeuronCoordinates'))
    neuronHandler = handles.NeuronCoordinates;
end
if(isfield(handles,'NeuronCoordinates2'))
    neuronHandler2 = handles.NeuronCoordinates2;
end
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
switch listboxstring
    case 'Empty'
        mat=-1;
    case 'Nucleus Matrix'
        mat = csvHandler.CellPosMatrix;         
        if(mat ~= -1)
            mat = mat(selectedWell);
        end
    case 'Cellomics Neurons'
        mat = neuronHandler.CellPosMatrix;
        if(mat ~= -1)
            mat = mat(selectedWell);
        end
    case 'Cellomics Neurons 2'
        mat = neuronHandler2.CellPosMatrix;
        if(mat ~= -1)
            mat = mat(selectedWell);
        end
    case 'Manual Neurons Main'
        mat = neuronHandler.ManualNeuronPositionsSparse;
        if(mat ~= -1)
            mat = mat(selectedWell);
        end
    case 'Manual Neurons 1'
        mat = csvHandler.ManualPositions1;
        if(mat ~= -1)
            mat = mat(selectedWell);
        end
    case 'Manual Neurons 2'
        mat = csvHandler.ManualPositions2;
        if(mat ~= -1)
            mat = mat(selectedWell);
        end
    case 'Manual Neurons 3'
        mat = csvHandler.ManualPositions3;
        if(mat ~= -1)
            mat = mat(selectedWell);
        end
    case 'Manual Neurons 4'
        mat = csvHandler.ManualPositions4;
        if(mat ~= -1)
            mat = mat(selectedWell);
        end
    case 'Edge Composit Neurons'
        mat = neuronHandler.NeuronPositionsEdgeComposite;
        if(mat ~= -1)
            mat = mat(selectedWell);
        end
    case 'Edge Fill Neurons'
        mat = neuronHandler.NeuronPositionsEdgeFill;
        if(mat ~= -1)
            mat = mat(selectedWell);
        end
    case 'Skeleton Neurons'
        del = neuronHandler.NeuronPositionsSkelDeleted(selectedWell);
        mat = neuronHandler.NeuronPositionsSkeletonization(selectedWell);
        mat=mat-del;
        ind = mat<0;
        mat(ind) = 0;
    otherwise
        mat=-1;        
end


    
function RefreshUnzoomedImage(handles)
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
if(numel(imageHandler.ZoomState) > 1)
    RefreshZoomedImage(handles)
else
if(isfield(handles,'NeuronCoordinates'))
    neuronHandler = handles.NeuronCoordinates;
end
if(isfield(handles,'NeuronCoordinates2'))
    neuronHandler2 = handles.NeuronCoordinates2;
end
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
%Check which checkboxes are checked.
cbNucleusPictureChecked = get(handles.cbNucleusPicture ,'Value');
cbNeuritePictureChecked = get(handles.cbNeuritePicture ,'Value');
cbBinaryPicChecked = get(handles.cbBinaryPic , 'Value');
cbSkeletonPicChecked = get(handles.cbSkeletonPic ,'Value');
cbRemoveCoreChecked = get(handles.cbRemoveCore, 'Value');
imageHandler.ZoomState = zeros(1);
[sizeY, sizeX] = size(imageHandler.ResizedNeuriteImage);
m = zeros(sizeY,sizeX);
imshow(m);


binaryImage = imageHandler.ResizedBinaryImage;    

%Analog to RefreshUnzoomed image but care about cutting Big Pictures before

if(cbRemoveCoreChecked)
    currentPicNuc = imageHandler.NucleusImage;
    currentPicNeu = imageHandler.NeuriteImage;
    binaryImage = imageHandler.BinaryImage;
    %skelImage = imageHandler.SkeletonImage;
    imMap = containers.Map();
    imMap('1') = currentPicNuc;
    imMap('2') = currentPicNeu;
    binAvailable=0;
    if(size(binaryImage,1) > 0)
        imMap('3') = binaryImage;
        binAvailable=1;
    end
    %imMap('4') = skelImage;
    imMap = CutOutCircles(imMap,selectedWell,1, 0,handles);
    currentPicNuc = imMap('1');
    currentPicNeu = imMap('2');
    if(binAvailable==1)
        binaryImage = imMap('3');
    end
    %skelImage = imMap('4');
    currentPicNuc = imresize(currentPicNuc, 0.1);
    currentPicNeu = imresize(currentPicNeu, 0.1);
    binaryImage = imresize(binaryImage, 0.1);
    %skelImage = imresize(skelImage, 0.1);
else
    currentPicNuc = imageHandler.ResizedNucleusImage;
    currentPicNeu = imageHandler.ResizedNeuriteImage;
    binaryImage = imageHandler.ResizedBinaryImage;    
end
skelImage = imageHandler.ResizedSkeletonImage;
skelImage = skelImage.*255;
if(cbSkeletonPicChecked)          
         imshow(skelImage);
         %set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
         %hPixelInfo = impixelinfo;        
         %imageHandler.MousePosition = hPixelInfo;
elseif(cbNeuritePictureChecked && cbBinaryPicChecked)
    if(~numel(imageHandler.NucleusImage) == 0 && ~numel(imageHandler.ResizedNeuriteImage) == 0)              
         mat = imfuse(imadjust(currentPicNeu,neuronHandler.StretchlimResult,[0 1]),binaryImage);
         hImage = imshow(mat);  
         set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});         
         hPixelInfo = impixelinfo;        
         imageHandler.MousePosition = hPixelInfo;
    end
elseif(cbNucleusPictureChecked && cbBinaryPicChecked)
    if(~numel(imageHandler.NucleusImage) == 0)              
         mat = imfuse(imadjust(currentPicNuc,neuronHandler.StretchlimResult,[0 1]),binaryImage);
         hImage = imshow(mat);  
         set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});         
         hPixelInfo = impixelinfo;        
         imageHandler.MousePosition = hPixelInfo;
    end
elseif(cbBinaryPicChecked)
    imshow(imageHandler.ResizedBinaryImage);
elseif(cbNucleusPictureChecked && cbNeuritePictureChecked)    
    if(~numel(imageHandler.NucleusImage) == 0 && ~numel(imageHandler.ResizedNeuriteImage) == 0)              
         mat = imfuse(currentPicNuc,currentPicNeu);
         hImage = imshow(mat);  
         set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});         
         hPixelInfo = impixelinfo;        
         imageHandler.MousePosition = hPixelInfo;
    end
elseif(cbNucleusPictureChecked)
    if(~numel(imageHandler.NucleusImage) == 0 && ~numel(imageHandler.ResizedNeuriteImage) == 0)              
         if(isfield(handles,'NeuronCoordinates') && ~isempty(neuronHandler.StretchlimResult))
             hImage = imshow(imadjust(currentPicNuc,neuronHandler.StretchlimResult,[0 1]));
         else
             hImage = imshow(currentPicNuc);
         end
         set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
         hPixelInfo = impixelinfo;        
         imageHandler.MousePosition = hPixelInfo;
    end
elseif(cbNeuritePictureChecked)
    if(~numel(imageHandler.NeuriteImage) == 0 && ~numel(imageHandler.ResizedNeuriteImage) == 0)              
         if(~isempty(neuronHandler.StretchlimResult))
            hImage = imshow(imadjust(currentPicNeu,neuronHandler.StretchlimResult,[0 1]));  
         else
             hImage = imshow(currentPicNeu); 
         end
         set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
         hPixelInfo = impixelinfo;        
         imageHandler.MousePosition = hPixelInfo;
    end
else
    %Check if one of the Pictures is selected. If not prepare plot!
end
%New: Add Matrices selected in Listboxes from listbox
listboxstrings = get(handles.popupPlus1);
listboxindex = get(handles.popupPlus1, 'value');
listboxstring = listboxstrings.String(listboxindex); 
MatrixPlus1 = GetMatrixFromListboxString(listboxstring,handles);
listboxindex = get(handles.popupPlus2, 'value');
listboxstring = listboxstrings.String(listboxindex); 
MatrixPlus2 = GetMatrixFromListboxString(listboxstring,handles);
listboxindex = get(handles.popupMinus1, 'value');
listboxstring = listboxstrings.String(listboxindex); 
MatrixMinus1 = GetMatrixFromListboxString(listboxstring,handles);
listboxindex = get(handles.popupMinus2, 'value');
listboxstring = listboxstrings.String(listboxindex); 
MatrixMinus2 = GetMatrixFromListboxString(listboxstring,handles);
listboxindex = get(handles.popupMinus3, 'value');
listboxstring = listboxstrings.String(listboxindex); 
MatrixMinus3 = GetMatrixFromListboxString(listboxstring,handles);
if(MatrixPlus1 ~= -1)
    if(MatrixMinus1 ~= -1)
        MatrixPlus1 = MatrixPlus1 - MatrixMinus1;
        ind = MatrixPlus1<0;
        MatrixPlus1(ind) = 0;
    end
    if(MatrixMinus2 ~= -1)
        MatrixPlus1 = MatrixPlus1 - MatrixMinus2;
        ind = MatrixPlus1<0;
        MatrixPlus1(ind) = 0;
    end
    if(MatrixMinus3 ~= -1)
        MatrixPlus1 = MatrixPlus1 - MatrixMinus3;
        ind = MatrixPlus1<0;
        MatrixPlus1(ind) = 0;
    end
    [nucleusRows, nucleusCols] = find(MatrixPlus1);   
     hold on;
     plot((nucleusCols./10),(nucleusRows./10),'Linestyle','none','Marker','.','Markersize',10,'Color','red');       
     hold off;
end
if(MatrixPlus2 ~= -1)
    %if(MatrixMinus2 ~= -1)
    %    MatrixPlus2 = MatrixPlus2 - MatrixMinus2;
    %   ind = MatrixPlus2<0;
    %    MatrixPlus2(ind) = 0;
    %end
    [nucleusRows, nucleusCols] = find(MatrixPlus2);   
     hold on;
     plot((nucleusCols./10),(nucleusRows./10),'Linestyle','none','Marker','.','Markersize',10,'Color','blue');       
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

function [filterDistance nonFilterDistance SphereArea markerPointCoordinates result] = calculateDensityDistribution(selectedWell, SphereAreaSizeX, SphereAreaSizeY, handles)
markerPointCoordinates=-1;
filterDistance = -1;
nonFilterDistance=-1;
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
SphereArea=-1;
optionHandler = handles.OptionHandler;
path = [imageHandler.Foldername '/ConvertedCellomics/MarkerPointCoordinates-' selectedWell '.mat'];
if(exist(path,'file'))    
    ringNumber = optionHandler.DensityDistributionRingNumber;
    load(path);
    filterDistance = str2double(filterDistance);
    nonFilterDistance = str2double(nonFilterDistance);
    %hInner = impoly(handles.axes2,double(markerPointCoordinates('0')./10));
    %previousMask = createMask(hInner);
    %for i=1:ringNumber  
    %  hCurrent = impoly(handles.axes2,double(markerPointCoordinates(num2str(i*10))./10));
    %end
  else
      
     [filterDistance nonFilterDistance SphereArea markerPointCoordinates] = calculateMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, handles); 
     %Save hInner to file
     foldername = imageHandler.Foldername;
     subfoldername = [foldername '/ConvertedCellomics'];
     %if(markerPointCoordinates~=0)
        save('-v7.3',strcat(subfoldername,'\MarkerPointCoordinates-',selectedWell),'markerPointCoordinates', 'filterDistance', 'nonFilterDistance');
     %end
  end
%[filterDistance nonFilterDistance SphereArea markerPointCoordinates] = calculateMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, handles);
%Migration Distance calculation gives distance from Sphere and Nucleus
%We need: LinePic
if(filterDistance == -1)
   result = zeros(ringNumber,9);
elseif (markerPointCoordinates==0)
   result = zeros(ringNumber,9);
elseif(numel(markerPointCoordinates('0')) == 0)
    result = zeros(ringNumber,9);
else

%NucleusM = csvHandler.CellPosMatrix;
%NeuronM = neuronHandler.CellPosMatrix;
NucleusM = csvHandler.CellPosMatrix(selectedWell);
NeuronM = neuronHandler.CellPosMatrix(selectedWell);
if(~isKey(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = sparse(7168, 7168);
end
NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
if(nnz(neuronHandler.NeuronPositionsSkeletonization) > 0 && isKey(neuronHandler.NeuronPositionsSkeletonization,selectedWell))
    NeuronSkelM = neuronHandler.NeuronPositionsSkeletonization(selectedWell);
    NeuronsDeleted = neuronHandler.NeuronPositionsSkelDeleted(selectedWell);

    NeuronSkelM = NeuronSkelM - (NeuronsDeleted);
    ind = NeuronSkelM<0;
    NeuronSkelM(ind) = 0;
end

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
result = zeros(ringNumber,9);
hInner = impoly(handles.axes2,double(markerPointCoordinates('0')./10));
previousMask = createMask(hInner);
cbRings = get(handles.cbPolys ,'Value');
if(~cbRings)
    delete(hInner);
end
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
    
    if(nnz(neuronHandler.NeuronPositionsSkeletonization) > 0 && isKey(neuronHandler.NeuronPositionsSkeletonization,selectedWell))
        NeuronSkel2 = logical(logical(NeuronSkelM) .* currentRing);
    else
        NeuronSkel2=0;
    end
    result(i,1) = nnz(currentRing);
    result(i,2) = nnz(Nucleus2);
    result(i,3) = nnz(NeuronCell2);
    result(i,4) = nnz(NeuronManual2);
    result(i,5) = nnz(NeuronSkel2);
    result(i,6) = result(i,3) / result(i,2);
    result(i,7) = result(i,4) / result(i,2);
    result(i,8) = result(i,5) / result(i,2);
    result(i,9) = result(i,2) / result(i,1);
end
end


function [filterDistance nonFilterDistance SphereArea markerPointCoordinates] = calculateMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, handles)
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
        NucleusM = logical(full(logical(csvHandler.CellPosMatrix(selectedWell))));
        %NucleusM = logical(full(NucleusM));
        NeuronM = logical(full(logical(neuronHandler.CellPosMatrix(selectedWell))));
        %NeuronM = logical(full(NeuronM));
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
     %figure(1);
     %imshow(NucleusDensity);
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
                    [xMax yMax] = MapPoint(xMax,yMax,sizeX,sizeY,256,256,1);
                    stack.push([xMax yMax]);
                 end
             end
         end
         radius = radius+1;
     end
     minimizingFactorX = 256/SphereAreaSizeX;
     minimizingFactorY = 256/SphereAreaSizeY;
     %New: Nucleus Area is calculated within full pic 
     NucleusArea = logical(zeros(sizeY,sizeX));
     
     floodFillThresh=optionHandler.MigDistLowerFloodFillThreshold;
     
     
     nucleusPic = imageHandler.NucleusImage;
    nucThreshold = optionHandler.NucleusThreshold;
    %[nucleusPic nucleusPicStrong] = ThresholdPic(nucleusPic,handles);
    nucleusPic(nucleusPic<nucThreshold) = 0;%double(double(imageHandler.NucleusImage) ./ 255);
    nucleusPic(nucleusPic>=nucThreshold) = 1;
    nucleusPic = logical(nucleusPic);
   
    LB = 3;
    UB = 3000;
    CuttedIndicesNuc = logical(nucleusPic);
    nucleusPic = xor(bwareaopen(nucleusPic,LB),  bwareaopen(nucleusPic,UB));
    CuttedIndicesNuc = logical(CuttedIndicesNuc - nucleusPic);
    nucleusPic=0;
    
    %Take from CuttedIndicesNuc only biggest area
    %Perform FloodFillOperation for all white pixels
    %Input: All white pixels
    %Return separate list of areas
    %Alle Randpixel weiﬂ!
    CuttedIndicesNuc(1,:)=0;
    CuttedIndicesNuc(:,1)=0;
    CuttedIndicesNuc(sizeY-1:sizeY,:)=0;
    CuttedIndicesNuc(:,sizeX-1:sizeX)=0;
    %Map points in stack to big coordinates
  %  while(numel(find(CuttedIndicesNuc)) > 0)
        [whiteRows whiteCols] = find(CuttedIndicesNuc);
        %floodStack = java.util.Stack();        
        %x=whiteCols(1);
        %y=whiteRows(1);
        %floodStack.push([x y]);
       % NucleusAreaTemp = logical(sparse(zeros(sizeY,sizeX)));
       i=0; 
       while(stack.isEmpty() == 0)
            vector = stack.pop();
            x=vector(1);
            y=vector(2);
            i=i+1;
            if(mod(i,10000) == 0)
                disp([num2str(nnz(NucleusArea)) ' of ' num2str(nnz(CuttedIndicesNuc))]);
            end
            if(CuttedIndicesNuc(y,x) == 1 && NucleusArea(y,x) == 0 && x-1 > 0 && y-1 > 0 && y+1<= sizeY && x+1 <= sizeX)
                NucleusArea(y,x) = 1;
                %CuttedIndicesNuc(y,x) = 0;
                stack.push([x y+1]);
                stack.push([x-1 y]);
                stack.push([x y-1]);
                stack.push([x+1 y]);
                %8 Neighbourhood
                %stack.push([x+1 y+1]);
                %stack.push([x-1 y-1]);
                %stack.push([x-1 y+1]);
                %stack.push([x+1 y-1]);
            elseif(NucleusArea(y,x) == 0)                
                NucleusArea(y,x) = 1;                
            end
        end
        
    %end
    CuttedIndicesNuc=0;
    NucleusAreaTemp=0;
    
    %NucleusArea is now filled with correct Nucleus and as big as picture.
    %NucleusAreaOld = logical(zeros(256,256));
    % NucleusAreaOld = FloodFillIterative(NucleusDensity,NucleusAreaOld,floodFillThresh,1,256,256,stack,0);
     %Also exclude pixels with less than 2 neighbours
    % filter = [1 1 1;1 0 1;1 1 1];
     %filteredPix = filter2(filter, NucleusArea);
     %[rows cols] = find(filteredPix < 2);
     %NucleusArea(rows,cols) = 0;
     %NucleusArea = edge(NucleusArea);
     %[posRows, posCols] = find(NucleusArea);
     %posRows = posRows ./ minimizingFactorY;
     %posCols = posCols ./ minimizingFactorX;
     
   %  NucleusArea64 = zeros(SphereAreaSizeY,SphereAreaSizeX);
   %  for i=1:numel(posRows)
   %     NucleusArea64(uint8(posRows(i)),uint8(posCols(i))) = 1;
   %  end
     [row,col] = find(NucleusArea);
     if(length(row) > 0 && length(col) > 0)
       
       sYMapped = sum(row)/length(row);
       sXMapped = sum(col)/length(col);
       %Map sY and sX to Big Picture
       [sX sY] = MapPoint(sXMapped,sYMapped,SphereAreaSizeX,SphereAreaSizeY,sizeX,sizeY,1);
       
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
      SE = [1 1 0;1 1 0; 0 0 0];
      % SE = strel('disk', 2);
      SphereArea = imdilate(SphereArea,SE);
      SphereArea = edge(SphereArea);
      %MigrationPic = SphereArea + NucleusArea64;
      %figure(1);
      %imshow(MigrationPic);
    %hPixelInfo = impixelinfo;        
    %imageHandler.MousePosition = hPixelInfo;
%F¸r weitere Verarbeitung:
% Berechne den Schwerpunkt aller weiﬂen Pixel in NucleusArea
     
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
        xBig=sXMapped;
        yBig=sYMapped;
        if(uint8(x)<1 | uint8(y)<1)
            continue;
        end
        SphereArea(uint8(y),uint8(x)) =0;
        while(xOuterSphere ==0 && uint8(x) <= SphereAreaSizeX && uint8(y) <= SphereAreaSizeY && uint8(x) > 0 && uint8(y) > 0 && uint8(xBig) > 0 && uint8(yBig) >0 && uint16(xBig) < sizeX && uint16(yBig) < sizeY)
          [x y] = MapPoint(xBig,yBig,SphereAreaSizeX,SphereAreaSizeY,sizeX,sizeY,1);
          if(uint8(x)<1 | uint8(y)<1)
              break;
          end
          %[xBig yBig] = MapPoint(x,y,sizeX,sizeY,SphereAreaSizeX,SphereAreaSizeY,0);
          if(NucleusArea(uint16(yBig),uint16(xBig)) == 1)
            xNucleus=xBig;
            yNucleus=yBig;
              %[xNucleus yNucleus]=[xBig yBig];%MapPoint(xBig,yBig,sizeX,sizeY,256,256,0);
          end
          
          if(SphereArea(uint8(y),uint8(x)) == 1)
            xOuterSphere=xBig;
            yOuterSphere=yBig;
              %[xOuterSphere yOuterSphere] = %MapPoint(x,y,sizeX,sizeY,SphereAreaSizeX,SphereAreaSizeY,1);
          end
          %x=x+deltaX;
          %y=y+deltaY;
          xBig=xBig+deltaX;
          yBig=yBig+deltaY;
        end
        if(xNucleus == 0)
            xNucleus=sXMapped;
            yNucleus=sYMapped;
            %[xNucleus yNucleus]= MapPoint(sX,sY,sizeX,sizeY,SphereAreaSizeX,SphereAreaSizeY,1);
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
        elseif(xNucleus ~=0 && yNucleus ~=0)
            %Set all points to Nucleus, if Outer Sphere is out of borders
            xOuterSphere = xNucleus;
            yOuterSphere = yNucleus;
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
     if(numel(mediumDistanceList)>0)
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
      threshold = (maxDistance-(0.85*maxDistance));
      %threshold = (maxDistance-(2*distancestddev));
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
     mediumDistance = num2str(mediumDistance,'%f');
     mediumDistanceSelected = num2str(mean(mediumDistanceListInInterval),'%f');
     disp('Medium distance total (in Pixel): ');
     disp(mediumDistance);
     disp('Medium distance selected (in Pixel): ');
     disp(mediumDistanceSelected);
     filterDistance = mediumDistanceSelected;
     nonFilterDistance = mediumDistance;

     else
        distancestddev = 0;
        mediumDistance = 0;
        mediumDistanceListInInterval=0;
        maxDistance = 0;
        filterDistance=-1;
        nonFilterDistance=-1;
     end
     
     
%mediumDistance = mediumDistance / successCount;
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
exportText = sprintf('Well;Ring Number;Migration Distance(px);Manual measured Migration Distance (px);Migration Distance without filter (px);Number of Pixels;Number Nuclei;Cellomics Neurons;Manual Neurons;Skeleton Neurons;Cellomics Neurons per Nuclei;Manual Neurons per Nuclei;Skeleton Neurons per Nuclei;Nuclei per Pixels;\r\n');
for currentWellIndex=1:numel(wellList)
    currentWell = wellList(currentWellIndex);
    currentWell = currentWell{1};
    disp(currentWell);
    [distFilter, distWithoutFilter, SphereArea, markerPointCoordinates densDist]=calculateDensityDistribution(currentWell,densityWidth,densityHeight,handles);
    if(numel(densDist) == 1)
        %Fill up densDist with dummy values
        densDist = zeros(ringCount,9);        
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
            exportText = [exportText currentWell ';' num2str(i) ';' distFilter ';' manualDist ';' distWithoutFilter ';' strrep(num2str(densDist(i,1)),'.',',') ';' strrep(num2str(densDist(i,2)),'.',',') ';' strrep(num2str(densDist(i,3)),'.',',') ';' strrep(num2str(densDist(i,4)),'.',',') ';' strrep(num2str(densDist(i,5)),'.',',') ';' strrep(num2str(densDist(i,6)),'.',',') ';' strrep(num2str(densDist(i,7)),'.',',') ';' strrep(num2str(densDist(i,8)),'.',',') ';' strrep(num2str(densDist(i,9)),'.',',') sprintf('\r\n')];
          else
            exportText = [exportText currentWell ';' num2str(i) ';' distFilter ';;' distWithoutFilter ';' strrep(num2str(densDist(i,1)),'.',',') ';' strrep(num2str(densDist(i,2)),'.',',') ';' strrep(num2str(densDist(i,3)),'.',',') ';' strrep(num2str(densDist(i,4)),'.',',') ';' strrep(num2str(densDist(i,5)),'.',',') ';' strrep(num2str(densDist(i,6)),'.',',') ';' strrep(num2str(densDist(i,7)),'.',',') ';' strrep(num2str(densDist(i,8)),'.',',') ';' strrep(num2str(densDist(i,9)),'.',',') sprintf('\r\n')];
          end
        end
    else
        for i=1:optionHandler.DensityDistributionRingNumber
            exportText = [exportText currentWell ';' num2str(i) ';' num2str(0) ';;' num2str(0) ';' strrep(num2str(0),'.',',') ';' strrep(num2str(0),'.',',') ';' strrep(num2str(0),'.',',') ';' strrep(num2str(0),'.',',') ';' strrep(num2str(0),'.',',') ';' strrep(num2str(0),'.',',') ';' strrep(num2str(0),'.',',') ';' strrep(num2str(0),'.',',') ';' strrep(num2str(0),'.',',') sprintf('\r\n')];
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


function EdgeCompositAlgorithm(selectedWell,handles)
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
NucleusM = csvHandler.CellPosMatrix(selectedWell);
densityWidth = optionHandler.MigrationDistanceDensityImageXSize;
densityHeight = optionHandler.MigrationDistanceDensityImageYSize;
[sizeY sizeX]= size(imageHandler.NeuriteImage);
ringNumber = optionHandler.DensityDistributionRingNumber;
%Get density distribution to exclude too dense areas
foldername = [imageHandler.Foldername '/ConvertedCellomics'];
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end
imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig.tif'];
imagePathNucleusSmall = [foldername '/' selectedWellLong 'NucleusSmall.tif'];
imagePathNeuriteBig = [foldername '/' selectedWellLong 'NeuriteBig.tif'];
imagePathNeuriteSmall = [foldername '/' selectedWellLong 'NeuriteSmall.tif'];

%Check if NeuronStatMatrix exists
if(~isa(neuronHandler.NeuronStatMatrix,'containers.Map'))
    neuronHandler.NeuronStatMatrix = containers.Map();
end
%if(isKey(neuronHandler.NeuronStatMatrix,selectedWell))
%    currentNeuronStat = neuronHandler.NeuronStatMatrix(selectedWell);
%else
%    currentNeuronStat = containers.Map();
%end

imageHandler.ResizedNeuriteImage = imread(imagePathNeuriteSmall);
imageHandler.ResizedNucleusImage = imread(imagePathNucleusSmall);    
imageHandler.NeuriteImage = imread(imagePathNeuriteBig);
imageHandler.NucleusImage = imread(imagePathNucleusBig);


path = [imageHandler.Foldername '/ConvertedCellomics/MarkerPointCoordinates-' selectedWell '.mat'];
if(exist(path,'file'))    
    ringNumber = optionHandler.DensityDistributionRingNumber;
    load(path);
    %hInner = impoly(handles.axes2,double(markerPointCoordinates('0')./10));
    %previousMask = createMask(hInner);
    %for i=1:ringNumber  
    %  hCurrent = impoly(handles.axes2,double(markerPointCoordinates(num2str(i*10))./10));
    %end
else
    [filterDistance nonFilterDistance SphereArea markerPointCoordinates] = calculateMigrationDistance(selectedWell, densityWidth, densityHeight, handles); 
    %Save hInner to file
     foldername = imageHandler.Foldername;
     subfoldername = [foldername '/ConvertedCellomics'];
     %if(markerPointCoordinates~=0)
        save('-v7.3',strcat(subfoldername,'\MarkerPointCoordinates-',selectedWell),'markerPointCoordinates', 'filterDistance', 'nonFilterDistance');
     %end
end
%MarkerPointCoordinates are available. Get second ring:
if(markerPointCoordinates ~= 0)
    hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./10));
    innerMask = logical(createMask(hInner));
    hOuter = impoly(handles.axes2,double(markerPointCoordinates(num2str(ringNumber * 10))./10));
    outerMask = logical(createMask(hOuter));
    currentRing = logical(outerMask-innerMask); 
    currentRing = imresize(currentRing, [sizeY, sizeX]);
    delete(hInner);
    delete(hOuter);
%if(strcmp(class(imageHandler.NeuriteImage),'uint16'))
%imageHandler.NeuriteImage = uint8(imageHandler.NeuriteImage);
%end
    brightNeurites = imageHandler.NeuriteImage .* uint8(currentRing);
else
    brightNeurites = imageHandler.NeuriteImage;
end
%Cut out dark pixels from NeuriteImage:     
%brightNeurites(brightNeurites > 150) = 0;
brightNeurites(brightNeurites < optionHandler.EdgeCompositNeuriteLowerThreshold) = 0;

sobelOperator = fspecial('sobel');
%overlaypicManuell = imfilter(imfuse(imageHandler.NucleusImage, brightNeurites),sobelOperator);
overlaypic = imfuse(edge(imageHandler.NucleusImage), edge(brightNeurites));
%Get white points (Neuron Points)
%neuronPointIndices = oimageHandler.NucleusImageverlaypicManuell(:,:,1) > 5 & overlaypicManuell(:,:,3) > 5 & overlaypicManuell(:,:,2) > 3;
neuronPointIndices = overlaypic(:,:,1) > 250 & overlaypic(:,:,2) > 250 & overlaypic(:,:,3) > 250;
figure(2);
imshow(neuronPointIndices);
figure(3);
imshow(overlaypic);
[sizeY sizeX] = size(NucleusM);
EdgeCompositPositionList = sparse(sizeY,sizeX);%neuronHandler.NeuronPositionsEdgeComposite(selectedWell);
%Look for circle structure in Neuron points:
%Execute FloodFill on neuronPointIndices and count number of Area Points
%If Area is big enough and not too big:
%Calculate focus of area and mark next Nucleus as Neuron

%Iterate over every nucleus and check if white point in overlayPic
[NonZeroY, NonZeroX] = find(NucleusM);
for i=1:numel(NonZeroY)    
    if(mod(i,100) == 0)
       disp([num2str(i) ' of ' num2str(numel(NonZeroY))]);
    end
    startX=0;
    startY=0;
    maxDistWhiteNuc = optionHandler.EdgeCompositeDistanceNucleusWhiteArea;
    for j=-maxDistWhiteNuc:maxDistWhiteNuc
        for k=-maxDistWhiteNuc:maxDistWhiteNuc             
            
            currentX = NonZeroX(i) + k;
            currentY = NonZeroY(i) + j;
            %If point or one neighbour is white in Overlaypic.
            if(currentY < sizeY && currentX < sizeX && currentY>0 && currentX>0 && neuronPointIndices(currentY,currentX) > 0)
                startX= currentX;
                startY = currentY;
            end
        end
    end
    if(startX > 0)
        %FloodFill white Area from current point.
        stack = java.util.Stack();
        
        
        %Idea to be more efficient:
        %Cut out NucleusArea from startX - 50 until startX + 50
        %Do FloodFill Operation only on this subset
        if(startY-50>0 && startX-50>0 && startX+50 < sizeX && startY + 50 < sizeY)
            SubArea = neuronPointIndices(startY-50:startY+50,startX-50:startX+50);
            WhiteArea2 = imfill(SubArea,[50 50]);         
            [yIndices xIndices] = find(WhiteArea2);
            yIndices = yIndices + (startY-50);
            xIndices = xIndices + (startX-50);
        else
            stack.push([startX startY]);
            WhiteArea = zeros(sizeY,sizeX);
            %WhiteArea2 = imfill(neuronPointIndices,[startY startX]);
            WhiteArea = FloodFillIterative(neuronPointIndices,WhiteArea,1,1,sizeY,sizeX,stack,0);
             %Schwerpunkt aller Punkte berechnen.
            [yIndices xIndices] = find(WhiteArea);
        end     
        %Check if Area has not too less and not too much pixels
        if(numel(xIndices) >= optionHandler.EdgeCompositMinArea)% && numel(xIndices) <=1200)
            xFocus = int16(mean(xIndices));
            yFocus = int16(mean(yIndices));
            %Mark nearest point as Neuron
            [xNuc,yNuc,distance] = csvHandler.FindNucleusForNeuron(xFocus,yFocus,NucleusM,75);
            EdgeCompositPositionList(yNuc,xNuc) = 1;
        end
        
    end
end
%Check if Neuronhandler already has Dictionary of NeuronPositionsEdgeComposite
if(numel(neuronHandler.NeuronPositionsEdgeComposite) ==0)
    neuronHandler.NeuronPositionsEdgeComposite = containers.Map();
end
neuronHandler.NeuronPositionsEdgeComposite(selectedWell) = EdgeCompositPositionList;
handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function EdgeCompositNeuronCount_Callback(hObject, eventdata, handles)
% hObject    handle to EdgeCompositNeuronCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Make overlay of both channels
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
EdgeCompositAlgorithm(selectedWell,handles);
% NucleusM = csvHandler.CellPosMatrix(selectedWell);
% densityWidth = optionHandler.MigrationDistanceDensityImageXSize;
% densityHeight = optionHandler.MigrationDistanceDensityImageYSize;
% [sizeY sizeX]= size(imageHandler.NeuriteImage);
% ringNumber = optionHandler.DensityDistributionRingNumber;
% %Get density distribution to exclude too dense areas
% path = [imageHandler.Foldername '/ConvertedCellomics/MarkerPointCoordinates-' selectedWell '.mat'];
% if(exist(path,'file'))    
%     ringNumber = optionHandler.DensityDistributionRingNumber;
%     load(path);
%     %hInner = impoly(handles.axes2,double(markerPointCoordinates('0')./10));
%     %previousMask = createMask(hInner);
%     %for i=1:ringNumber  
%     %  hCurrent = impoly(handles.axes2,double(markerPointCoordinates(num2str(i*10))./10));
%     %end
% else
%     [filterDistance nonFilterDistance SphereArea NucleusArea64 markerPointCoordinates] = calculateMigrationDistance(selectedWell, densityWidth, densityHeight, handles);    
% end
% %MarkerPointCoordinates are available. Get second ring:
% hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./10));
% innerMask = logical(createMask(hInner));
% hOuter = impoly(handles.axes2,double(markerPointCoordinates(num2str(ringNumber * 10))./10));
% outerMask = logical(createMask(hOuter));
% currentRing = logical(outerMask-innerMask); 
% currentRing = imresize(currentRing, [sizeY, sizeX]);
% delete(hInner);
% delete(hOuter);
% %if(strcmp(class(imageHandler.NeuriteImage),'uint16'))
% %imageHandler.NeuriteImage = uint8(imageHandler.NeuriteImage);
% %end
% brightNeurites = imageHandler.NeuriteImage .* uint8(currentRing);
% 
% %Cut out dark pixels from NeuriteImage:     
% %brightNeurites(brightNeurites > 150) = 0;
% brightNeurites(brightNeurites < optionHandler.EdgeCompositNeuriteLowerThreshold) = 0;
% 
% sobelOperator = fspecial('sobel');
% %overlaypicManuell = imfilter(imfuse(imageHandler.NucleusImage, brightNeurites),sobelOperator);
% overlaypic = imfuse(edge(imageHandler.NucleusImage), edge(brightNeurites));
% %Get white points (Neuron Points)
% %neuronPointIndices = overlaypicManuell(:,:,1) > 5 & overlaypicManuell(:,:,3) > 5 & overlaypicManuell(:,:,2) > 3;
% neuronPointIndices = overlaypic(:,:,1) > 250 & overlaypic(:,:,2) > 250 & overlaypic(:,:,3) > 250;
% figure(2);
% imshow(neuronPointIndices);
% figure(3);
% imshow(overlaypic);
% [sizeY sizeX] = size(NucleusM);
% EdgeCompositPositionList = sparse(sizeY,sizeX);%neuronHandler.NeuronPositionsEdgeComposite(selectedWell);
% %Look for circle structure in Neuron points:
% %Execute FloodFill on neuronPointIndices and count number of Area Points
% %If Area is big enough and not too big:
% %Calculate focus of area and mark next Nucleus as Neuron
% 
% %Iterate over every nucleus and check if white point in overlayPic
% [NonZeroY, NonZeroX] = find(NucleusM);
% for i=1:numel(NonZeroY)    
%     if(mod(i,100) == 0)
%        disp([num2str(i) ' of ' num2str(numel(NonZeroY))]);
%     end
%     startX=0;
%     startY=0;
%     maxDistWhiteNuc = optionHandler.EdgeCompositeDistanceNucleusWhiteArea;
%     for j=-maxDistWhiteNuc:maxDistWhiteNuc
%         for k=-maxDistWhiteNuc:maxDistWhiteNuc             
%             
%             currentX = NonZeroX(i) + k;
%             currentY = NonZeroY(i) + j;
%             %If point or one neighbour is white in Overlaypic.
%             if(currentY < sizeY && currentX < sizeX && currentY>0 && currentX>0 && neuronPointIndices(currentY,currentX) > 0)
%                 startX= currentX;
%                 startY = currentY;
%             end
%         end
%     end
%     if(startX > 0)
%         %FloodFill white Area from current point.
%         stack = java.util.Stack();
%         
%         
%         %Idea to be more efficient:
%         %Cut out NucleusArea from startX - 50 until startX + 50
%         %Do FloodFill Operation only on this subset
%         if(startY-50>0 && startX-50>0 && startX+50 < sizeX && startY + 50 < sizeY)
%             SubArea = neuronPointIndices(startY-50:startY+50,startX-50:startX+50);
%             WhiteArea2 = imfill(SubArea,[50 50]);         
%             [yIndices xIndices] = find(WhiteArea2);
%             yIndices = yIndices + (startY-50);
%             xIndices = xIndices + (startX-50);
%         else
%             stack.push([startX startY]);
%             WhiteArea = zeros(sizeY,sizeX);
%             %WhiteArea2 = imfill(neuronPointIndices,[startY startX]);
%             WhiteArea = FloodFillIterative(neuronPointIndices,WhiteArea,1,1,sizeY,sizeX,stack,0);
%              %Schwerpunkt aller Punkte berechnen.
%             [yIndices xIndices] = find(WhiteArea);
%         end     
%         %Check if Area has not too less and not too much pixels
%         if(numel(xIndices) >= optionHandler.EdgeCompositMinArea)% && numel(xIndices) <=1200)
%             xFocus = int16(mean(xIndices));
%             yFocus = int16(mean(yIndices));
%             %Mark nearest point as Neuron
%             [xNuc,yNuc,distance] = csvHandler.FindNucleusForNeuron(xFocus,yFocus,NucleusM,75);
%             EdgeCompositPositionList(yNuc,xNuc) = 1;
%         end
%         
%     end
% end
% %Check if Neuronhandler already has Dictionary of NeuronPositionsEdgeComposite
% if(numel(neuronHandler.NeuronPositionsEdgeComposite) ==0)
%     neuronHandler.NeuronPositionsEdgeComposite = containers.Map();
% end
% neuronHandler.NeuronPositionsEdgeComposite(selectedWell) = EdgeCompositPositionList;
% handles.NeuronCoordinates = neuronHandler;
% guidata(handles.figure1, handles);


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
function Opt_Callback(hObject, eventdata, handles)
% hObject    handle to Opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Create Input Dialogue



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
[filterDistance nonFilterDistance SphereArea markerPointCoordinates result] = calculateDensityDistribution(selectedWell, densityWidth, densityHeight, handles);
if(numel(result) > 1)
    disp('Ring Number;Pixels in Ring;Nuclei in Ring;Neurons by Cellomics in Ring;Manual counted Neurons in Ring;Skeletonized Neurons in Ring; Cellomics Neurons per Nuclei;Manual Neurons per Nuclei;Skeleton Neurons per Nuclei; Nuclei per Pixels');
    for i=1:optionHandler.DensityDistributionRingNumber
        disp([num2str(i) ';' strrep(num2str(result(i,1)),'.',',') ';' strrep(num2str(result(i,2)),'.',',') ';' strrep(num2str(result(i,3)),'.',',') ';' strrep(num2str(result(i,4)),'.',',') ';' strrep(num2str(result(i,5)),'.',',') ';' strrep(num2str(result(i,6)),'.',',') ';' strrep(num2str(result(i,7)),'.',',') ';' strrep(num2str(result(i,8)),'.',',') ';' strrep(num2str(result(i,9)),'.',',')]);
    end
end


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9


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
        if(numel(markerPointCoordinates('10') > 0))
          hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./10));
          previousMask = createMask(hInner);
          for i=2:ringNumber  
            hCurrent = impoly(handles.axes2,double(markerPointCoordinates(num2str(i*10))./10));
          end
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

%Iteriere ¸ber Filter f¸r das aktuelle Well
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
    if(isprop(neuronHandler, 'CellPosMatrix'))
        WellNeuronDict = neuronHandler.CellPosMatrix;
        wellList = keys(WellNeuronDict);
        %Remove all entries from lbWell
        set(handles.lbWell, 'String', '');    
        fulltext = '';
        for i=1:numel(wellList)
          fulltext = [fulltext;wellList(i)];
        end
        set(handles.lbWell, 'String', fulltext);
    end
    WellNucleusDict = csvHandler.CellPosMatrix;
end


% --- Executes on button press in cbEdgeCompositeNeurons.
function cbEdgeCompositeNeurons_Callback(hObject, eventdata, handles)
% hObject    handle to cbEdgeCompositeNeurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);


% --- Executes on button press in cbSkeletonNeurons.
function cbSkeletonNeurons_Callback(hObject, eventdata, handles)
% hObject    handle to cbSkeletonNeurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles)



% --------------------------------------------------------------------
function QualityCheckOld_Callback(hObject, eventdata, handles)
% hObject    handle to QualityCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
WellNeuronDict = neuronHandler.NeuronPositionsEdgeComposite;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};




densityWidth = 50;
densityHeight = 50;
[sizeY sizeX]= size(imageHandler.NeuriteImage);
ringNumber = optionHandler.DensityDistributionRingNumber;
%Get density distribution to exclude too dense areas
path = [imageHandler.Foldername '/ConvertedCellomics/MarkerPointCoordinates-' selectedWell '.mat'];
%if(exist(path,'file'))    
%    ringNumber = optionHandler.DensityDistributionRingNumber;
%    load(path);
%    hInner = impoly(handles.axes2,double(markerPointCoordinates('0')./10));
%    previousMask = createMask(hInner);
%    for i=1:ringNumber  
%      hCurrent = impoly(handles.axes2,double(markerPointCoordinates(num2str(i*10))./10));
%    end
%else
    [filterDistance nonFilterDistance SphereArea NucleusArea64 markerPointCoordinates] = calculateMigrationDistance(selectedWell, densityWidth, densityHeight, handles);    
%end
%MarkerPointCoordinates are available. Get second ring:
hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./10));
innerMask = logical(createMask(hInner));
%hOuter = impoly(handles.axes2,double(markerPointCoordinates(num2str(ringNumber * 10))./10));
%outerMask = logical(createMask(hOuter));
%if(CutInnerCircle == 1 && CutOuterCircle == 1)
%    currentRing = logical(outerMask-innerMask); 
%elseif(CutInnerCircle==1 && CutOuterCircle == 0)
    currentRing =(logical(ones(int32(sizeY/10),int32(sizeX/10)) - innerMask));
%else
%    currentRing = outerMask;
%end
currentRing = imresize(currentRing, [sizeY, sizeX]);
currentRing=sparse(currentRing);
delete(hInner);
%delete(hOuter);
innerMask=0;

EdgeCompNeurons = logical(WellNeuronDict(selectedWell));
area = optionHandler.EdgeFillNucleusAreaWithinNeurite;
[EdgeCompSizeY EdgeCompSizeX] = size(EdgeCompNeurons);
[sizeY sizeX] = size(imageHandler.NeuriteImage);
if(EdgeCompSizeX > sizeX)
      EdgeCompNeurons(:,sizeX+1:EdgeCompSizeX) = [];
      neuronHandler.NeuronPositionsEdgeComposite(selectedWell) = EdgeCompNeurons;
end
if(EdgeCompSizeY > sizeY)
      EdgeCompNeurons(sizeY+1:EdgeCompSizeY,:) = [];
      neuronHandler.NeuronPositionsEdgeComposite(selectedWell) = EdgeCompNeurons;
end


ManualNeurons = logical(neuronHandler.ManualNeuronPositionsSparse(selectedWell) .* currentRing);
filterDistance=0;
nonFilterDistance=0;
SphereArea=0;
NucleusArea64=0;
markerPointCoordinates=0;
ManualNeuronsInvert = ~ManualNeurons;
EdgeFilledNeurons = logical(neuronHandler.NeuronPositionsEdgeFill(selectedWell).* currentRing);

%EdgeFilledNeurons = EdgeFilledNeurons(:,:,2) > area;
CellomicsNeurons = neuronHandler.CellPosMatrix(selectedWell);
SkeletonNeurons = logical(neuronHandler.NeuronPositionsSkeletonization(selectedWell) .* currentRing);
currentRing=0;
NeuronsDeleted = neuronHandler.NeuronPositionsSkelDeleted(selectedWell);

SkeletonNeurons = SkeletonNeurons; %+ EdgeFilledNeurons;
SkeletonNeurons = SkeletonNeurons - (NeuronsDeleted);
ind = SkeletonNeurons<0;
SkeletonNeurons(ind) = 0;

if(isfield(handles,'NeuronCoordinates2'))
  neuronHandler2 = handles.NeuronCoordinates2;
  CellomicsNeuronsEnhanced = neuronHandler2.CellPosMatrix(selectedWell);
  bothManualCellomicsEnhanced = nnz(ManualNeurons .* CellomicsNeuronsEnhanced);
  
  CellomicsEnhancedInvert = ~CellomicsNeuronsEnhanced;
  AdditionalManualCellomicsEnhanced =  nnz(ManualNeurons .* CellomicsEnhancedInvert);
  AdditionalCellomicsEnhancedManual = nnz(CellomicsNeuronsEnhanced .* ManualNeuronsInvert);
else
    CellomicsNeuronsEnhanced = 0;
    bothManualCellomicsEnhanced = 0;
    CellomicsEnhancedInvert = 0;
    AdditionalManualCellomicsEnhanced = 0;
    AdditionalCellomicsEnhancedManual = 0;
end


disp('Neurons Manual;Neurons Cellomics;Neurons Edge Composite;Neurons Cellomics Enhanced;Neurons Skeleton;Neurons Edge Fill;Positions Manual & Edge;Additional Manual; Additional Edge Composite;;Positions Edge & Cellomics;Additional Cellomics; Additional Edge Composite;; Positions Manual & Cellomics; Additional Manual; Additional Cellomics;; Positions Manual & Enhanced Cellomics; Additional Manual; Additional Enhanced Cellomics;; Positions Manual & Edge Filled; Additioal Manual; Additional Edge Filled;; Positions Manual & Skeletonized; Additional Manual; Additional Skeletonization')
%disp('Neurons Manual;Neurons Cellomics;Neurons Edge Composite;Positions Manual & Edge;Additional Manual; Additional Edge Composite;;Positions Edge & Cellomics;Additional Cellomics; Additional Edge Composite;; Positions Manual & Cellomics; Additional Manual; Additional Cellomics;; Positions Manual & Edge Filled; Additioal Manual; Additional Edge Filled')
bothCount = nnz(EdgeCompNeurons .* ManualNeurons);
bothEdgeCellomics = nnz(EdgeCompNeurons .* CellomicsNeurons);
bothManualCellomics = nnz(ManualNeurons .* CellomicsNeurons);
bothManualEdgeFilled = nnz(EdgeFilledNeurons .* ManualNeurons);

bothManualSkeletonized = nnz(SkeletonNeurons .* ManualNeurons);






AdditionalEdge = nnz(EdgeCompNeurons .* ManualNeuronsInvert);
CellomicsNeuronsInvert = ~CellomicsNeurons;
AdditionalEdgeCellomics = nnz(EdgeCompNeurons .* CellomicsNeuronsInvert);
AdditionalManualCellomics = nnz(ManualNeurons .* CellomicsNeuronsInvert);
CellomicsNeuronsInvert=0;
AdditionalCellomicsManual = nnz(CellomicsNeurons .* ManualNeuronsInvert);
AdditionalSkeletonManual = nnz(SkeletonNeurons .* ManualNeuronsInvert);
SkeletonNeuronsInvert = ~SkeletonNeurons;
AdditionalManualSkeleton = nnz(ManualNeurons .* SkeletonNeuronsInvert);
SkeletonNeuronsInvert=0;

AdditionalEdgeFilledManual = nnz(EdgeFilledNeurons .* ManualNeuronsInvert);
ManualNeuronsInvert=0;
EdgeFilledNeuronsInvert = ~EdgeFilledNeurons;
AdditionalManualEdgeFilled = nnz(ManualNeurons .* EdgeFilledNeuronsInvert);
EdgeCompInvert = ~EdgeCompNeurons;
AdditionalCellomics = nnz(EdgeCompInvert .* CellomicsNeurons);
AdditionalManual = nnz(EdgeCompInvert .* ManualNeurons);

disp([num2str(nnz(ManualNeurons)) ';' num2str(nnz(CellomicsNeurons)) ';' num2str(nnz(EdgeCompNeurons)) ';' num2str(nnz(CellomicsNeuronsEnhanced)) ';' num2str(nnz(SkeletonNeurons)) ';' num2str(nnz(EdgeFilledNeurons)) ';' num2str(bothCount) ';' num2str(AdditionalManual) ';' num2str(AdditionalEdge) ';;' num2str(bothEdgeCellomics) ';' num2str(AdditionalCellomics) ';' num2str(AdditionalEdgeCellomics) ';;' num2str(bothManualCellomics) ';' num2str(AdditionalManualCellomics) ';' num2str(AdditionalCellomicsManual) ';;' num2str(bothManualCellomicsEnhanced) ';' num2str(AdditionalManualCellomicsEnhanced) ';' num2str(AdditionalCellomicsEnhancedManual) ';;' num2str(bothManualEdgeFilled) ';' num2str(AdditionalManualEdgeFilled) ';' num2str(AdditionalEdgeFilledManual) ';;' num2str(bothManualSkeletonized) ';' num2str(AdditionalManualSkeleton) ';' num2str(AdditionalSkeletonManual) ';;']);
%disp([num2str(nnz(ManualNeurons)) ';' num2str(nnz(CellomicsNeurons)) ';' num2str(nnz(EdgeCompNeurons)) ';' num2str(bothCount) ';' num2str(AdditionalManual) ';' num2str(AdditionalEdge) ';;' num2str(bothEdgeCellomics) ';' num2str(AdditionalCellomics) ';' num2str(AdditionalEdgeCellomics) ';;' num2str(bothManualCellomics) ';' num2str(AdditionalManualCellomics) ';' num2str(AdditionalCellomicsManual) ';;' num2str(bothManualEdgeFilled) ';' num2str(AdditionalManualEdgeFilled) ';' num2str(AdditionalEdgeFilledManual) ';;']);
handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function QualityCheck_Callback(hObject, eventdata, handles)
% hObject    handle to QualityCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
WellNeuronDict = neuronHandler.NeuronPositionsEdgeComposite;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};




densityWidth = 50;
densityHeight = 50;
[sizeY sizeX]= size(imageHandler.NeuriteImage);
ringNumber = optionHandler.DensityDistributionRingNumber;
%Get density distribution to exclude too dense areas
path = [imageHandler.Foldername '/ConvertedCellomics/MarkerPointCoordinates-' selectedWell '.mat'];
path = [imageHandler.Foldername '/ConvertedCellomics/MarkerPointCoordinates-' selectedWell '.mat'];
%Get density distribution to exclude too dense areas
%if(saveCircle)    
  if(exist(path,'file'))    
    ringNumber = optionHandler.DensityDistributionRingNumber;
    load(path);
    %hInner = impoly(handles.axes2,double(markerPointCoordinates('0')./10));
    %previousMask = createMask(hInner);
    %for i=1:ringNumber  
    %  hCurrent = impoly(handles.axes2,double(markerPointCoordinates(num2str(i*10))./10));
    %end
  else
     [filterDistance nonFilterDistance SphereArea markerPointCoordinates] = calculateMigrationDistance(selectedWell, densityWidth, densityHeight, handles); 
     %Save hInner to file
     foldername = imageHandler.Foldername;
     subfoldername = [foldername '/ConvertedCellomics'];
     save('-v7.3',strcat(subfoldername,'\MarkerPointCoordinates-',selectedWell),'markerPointCoordinates', 'filterDistance', 'nonFilterDistance');
  end
%MarkerPointCoordinates are available. Get second ring:
hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./10));
innerMask = logical(createMask(hInner));
%hOuter = impoly(handles.axes2,double(markerPointCoordinates(num2str(ringNumber * 10))./10));
%outerMask = logical(createMask(hOuter));
%if(CutInnerCircle == 1 && CutOuterCircle == 1)
%    currentRing = logical(outerMask-innerMask); 
%elseif(CutInnerCircle==1 && CutOuterCircle == 0)
    currentRing =(logical(ones(int32(sizeY/10),int32(sizeX/10)) - innerMask));
%else
%    currentRing = outerMask;
%end
currentRing = imresize(currentRing, [sizeY, sizeX]);
currentRing=sparse(currentRing);
delete(hInner);
%delete(hOuter);
innerMask=0;

EdgeCompNeurons = logical(WellNeuronDict(selectedWell));
area = optionHandler.EdgeFillNucleusAreaWithinNeurite;
[EdgeCompSizeY EdgeCompSizeX] = size(EdgeCompNeurons);
[sizeY sizeX] = size(imageHandler.NeuriteImage);
if(EdgeCompSizeX > sizeX)
      EdgeCompNeurons(:,sizeX+1:EdgeCompSizeX) = [];
      neuronHandler.NeuronPositionsEdgeComposite(selectedWell) = EdgeCompNeurons;
end
if(EdgeCompSizeY > sizeY)
      EdgeCompNeurons(sizeY+1:EdgeCompSizeY,:) = [];
      neuronHandler.NeuronPositionsEdgeComposite(selectedWell) = EdgeCompNeurons;
end


ManualNeurons = logical(neuronHandler.ManualNeuronPositionsSparse(selectedWell) .* currentRing);
filterDistance=0;
nonFilterDistance=0;
SphereArea=0;
NucleusArea64=0;
markerPointCoordinates=0;
ManualNeuronsInvert = ~ManualNeurons;
EdgeFilledNeurons = logical(neuronHandler.NeuronPositionsEdgeFill(selectedWell).* currentRing);

%EdgeFilledNeurons = EdgeFilledNeurons(:,:,2) > area;
CellomicsNeurons = neuronHandler.CellPosMatrix(selectedWell);
SkeletonNeurons = logical(neuronHandler.NeuronPositionsSkeletonization(selectedWell) .* currentRing);
currentRing=0;
NeuronsDeleted = neuronHandler.NeuronPositionsSkelDeleted(selectedWell);

SkeletonNeurons = SkeletonNeurons; %+ EdgeFilledNeurons;
SkeletonNeurons = SkeletonNeurons - (NeuronsDeleted);
ind = SkeletonNeurons<0;
SkeletonNeurons(ind) = 0;

if(isfield(handles,'NeuronCoordinates2'))
  neuronHandler2 = handles.NeuronCoordinates2;
  CellomicsNeuronsEnhanced = neuronHandler2.CellPosMatrix(selectedWell);
  bothManualCellomicsEnhanced = nnz(ManualNeurons .* CellomicsNeuronsEnhanced);
  
  CellomicsEnhancedInvert = ~CellomicsNeuronsEnhanced;
  AdditionalManualCellomicsEnhanced =  nnz(ManualNeurons .* CellomicsEnhancedInvert);
  AdditionalCellomicsEnhancedManual = nnz(CellomicsNeuronsEnhanced .* ManualNeuronsInvert);
else
    CellomicsNeuronsEnhanced = 0;
    bothManualCellomicsEnhanced = 0;
    CellomicsEnhancedInvert = 0;
    AdditionalManualCellomicsEnhanced = 0;
    AdditionalCellomicsEnhancedManual = 0;
end


disp('Neurons Manual;Neurons Cellomics;Neurons Edge Composite;Neurons Cellomics Enhanced;Neurons Skeleton;Neurons Edge Fill;Positions Manual & Edge;Additional Manual; Additional Edge Composite;;Positions Edge & Cellomics;Additional Cellomics; Additional Edge Composite;; Positions Manual & Cellomics; Additional Manual; Additional Cellomics;; Positions Manual & Enhanced Cellomics; Additional Manual; Additional Enhanced Cellomics;; Positions Manual & Edge Filled; Additioal Manual; Additional Edge Filled;; Positions Manual & Skeletonized; Additional Manual; Additional Skeletonization')
%disp('Neurons Manual;Neurons Cellomics;Neurons Edge Composite;Positions Manual & Edge;Additional Manual; Additional Edge Composite;;Positions Edge & Cellomics;Additional Cellomics; Additional Edge Composite;; Positions Manual & Cellomics; Additional Manual; Additional Cellomics;; Positions Manual & Edge Filled; Additioal Manual; Additional Edge Filled')
bothCount = nnz(EdgeCompNeurons .* ManualNeurons);
bothEdgeCellomics = nnz(EdgeCompNeurons .* CellomicsNeurons);
bothManualCellomics = nnz(ManualNeurons .* CellomicsNeurons);
bothManualEdgeFilled = nnz(EdgeFilledNeurons .* ManualNeurons);

bothManualSkeletonized = nnz(SkeletonNeurons .* ManualNeurons);






AdditionalEdge = nnz(EdgeCompNeurons .* ManualNeuronsInvert);
CellomicsNeuronsInvert = ~CellomicsNeurons;
AdditionalEdgeCellomics = nnz(EdgeCompNeurons .* CellomicsNeuronsInvert);
AdditionalManualCellomics = nnz(ManualNeurons .* CellomicsNeuronsInvert);
CellomicsNeuronsInvert=0;
AdditionalCellomicsManual = nnz(CellomicsNeurons .* ManualNeuronsInvert);
AdditionalSkeletonManual = nnz(SkeletonNeurons .* ManualNeuronsInvert);
SkeletonNeuronsInvert = ~SkeletonNeurons;
AdditionalManualSkeleton = nnz(ManualNeurons .* SkeletonNeuronsInvert);
SkeletonNeuronsInvert=0;

AdditionalEdgeFilledManual = nnz(EdgeFilledNeurons .* ManualNeuronsInvert);
ManualNeuronsInvert=0;
EdgeFilledNeuronsInvert = ~EdgeFilledNeurons;
AdditionalManualEdgeFilled = nnz(ManualNeurons .* EdgeFilledNeuronsInvert);
EdgeCompInvert = ~EdgeCompNeurons;
AdditionalCellomics = nnz(EdgeCompInvert .* CellomicsNeurons);
AdditionalManual = nnz(EdgeCompInvert .* ManualNeurons);

disp([num2str(nnz(ManualNeurons)) ';' num2str(nnz(CellomicsNeurons)) ';' num2str(nnz(EdgeCompNeurons)) ';' num2str(nnz(CellomicsNeuronsEnhanced)) ';' num2str(nnz(SkeletonNeurons)) ';' num2str(nnz(EdgeFilledNeurons)) ';' num2str(bothCount) ';' num2str(AdditionalManual) ';' num2str(AdditionalEdge) ';;' num2str(bothEdgeCellomics) ';' num2str(AdditionalCellomics) ';' num2str(AdditionalEdgeCellomics) ';;' num2str(bothManualCellomics) ';' num2str(AdditionalManualCellomics) ';' num2str(AdditionalCellomicsManual) ';;' num2str(bothManualCellomicsEnhanced) ';' num2str(AdditionalManualCellomicsEnhanced) ';' num2str(AdditionalCellomicsEnhancedManual) ';;' num2str(bothManualEdgeFilled) ';' num2str(AdditionalManualEdgeFilled) ';' num2str(AdditionalEdgeFilledManual) ';;' num2str(bothManualSkeletonized) ';' num2str(AdditionalManualSkeleton) ';' num2str(AdditionalSkeletonManual) ';;']);
%disp([num2str(nnz(ManualNeurons)) ';' num2str(nnz(CellomicsNeurons)) ';' num2str(nnz(EdgeCompNeurons)) ';' num2str(bothCount) ';' num2str(AdditionalManual) ';' num2str(AdditionalEdge) ';;' num2str(bothEdgeCellomics) ';' num2str(AdditionalCellomics) ';' num2str(AdditionalEdgeCellomics) ';;' num2str(bothManualCellomics) ';' num2str(AdditionalManualCellomics) ';' num2str(AdditionalCellomicsManual) ';;' num2str(bothManualEdgeFilled) ';' num2str(AdditionalManualEdgeFilled) ';' num2str(AdditionalEdgeFilledManual) ';;']);
handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function NeuronAlgoOptions_Callback(hObject, eventdata, handles)
% hObject    handle to NeuronAlgoOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
    optionHandler.EdgeCompositNeuriteLowerThreshold = 20;
    optionHandler.EdgeCompositMinArea = 70;
    optionHandler.EdgeCompositeDistanceNucleusWhiteArea = 3;
    optionHandler.EdgeFillNucleusAreaWithinNeurite = 0.75;
    optionHandler.EdgeFillNucleusAreaWithinNeuriteSecond = 0.55;
    optionHandler.NucleusThreshold = 15;
    
    optionHandler.SkeletonMinNeuriteLength = 40;
    %1 = Isodata for every single picture
    %2 = Isodata mean
    optionHandler.SkeletonThresholdMethod = 1;
    optionHandler.SkeletonNeuriteThresholdHardDistance = 0.4;
    SkeletonNeuriteThresholdLow = 0.75;
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
    optionHandler.EdgeCompositNeuriteLowerThreshold = 20;
    optionHandler.EdgeCompositMinArea = 70;
    optionHandler.EdgeCompositeDistanceNucleusWhiteArea = 3;
    optionHandler.EdgeFillNucleusAreaWithinNeurite = 0.75;
    optionHandler.EdgeFillNucleusAreaWithinNeurite = 0.55;
    optionHandler.NucleusThreshold = 15;
    
    optionHandler.SkeletonMinNeuriteLength = 40;
    %1 = Isodata for every single picture
    %2 = Isodata mean
    %3 = Manual Threshold
    optionHandler.SkeletonThresholdMethod = 1;
    optionHandler.SkeletonNeuriteThresholdHardDistance = 0.4;
    SkeletonNeuriteThresholdLow = 0.75;
end
prompt = {'Edge Composite: Lower Neurite Threshold:','Edge Composite: Min Neurite Area (in Pixels):','Edge Composite: Max Distance Nucleus Area:','Edge Fill: Nuc Area in Percent within Neurite','Edge Fill: Nuc Area in Percent within Neurite Lower','Skeleton: Min Neurite Length','Threshold Method (1=Isodata for single picture, 2=mean Isodata, 3=manual Threshold)','Hard Threshold Distance from Manual Threshold','Manual Low Threshold', 'Fix Threshold for Nucleus Pic'};
dlg_title = 'Options';
num_lines = 1;
def = {num2str(optionHandler.EdgeCompositNeuriteLowerThreshold),num2str(optionHandler.EdgeCompositMinArea),num2str(optionHandler.EdgeCompositeDistanceNucleusWhiteArea),num2str(optionHandler.EdgeFillNucleusAreaWithinNeurite),num2str(optionHandler.EdgeFillNucleusAreaWithinNeuriteSecond), num2str(optionHandler.SkeletonMinNeuriteLength), num2str(optionHandler.SkeletonThresholdMethod), num2str(optionHandler.SkeletonNeuriteThresholdHardDistance), num2str(optionHandler.SkeletonNeuriteThresholdLow), num2str(optionHandler.NucleusThreshold)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
optionHandler.EdgeCompositNeuriteLowerThreshold = str2num(answer{1});
optionHandler.EdgeCompositMinArea = str2num(answer{2});
optionHandler.EdgeCompositeDistanceNucleusWhiteArea = str2num(answer{3});
optionHandler.EdgeFillNucleusAreaWithinNeurite = str2num(answer{4});
optionHandler.EdgeFillNucleusAreaWithinNeuriteSecond = str2num(answer{5});
optionHandler.SkeletonMinNeuriteLength = str2num(answer{6});
optionHandler.SkeletonThresholdMethod = str2num(answer{7});
optionHandler.SkeletonNeuriteThresholdHardDistance = str2num(answer{8});
optionHandler.SkeletonNeuriteThresholdLow = str2num(answer{9});
optionHandler.NucleusThreshold = str2num(answer{10});
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles)

% --------------------------------------------------------------------
function Options_Callback(hObject, eventdata, handles)
% hObject    handle to Options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
    optionHandler.EdgeCompositNeuriteLowerThreshold = 20;
    optionHandler.EdgeCompositMinArea = 70;
    optionHandler.EdgeCompositeDistanceNucleusWhiteArea = 3;
    optionHandler.EdgeFillNucleusAreaWithinNeurite = 0.55;
     optionHandler.EdgeFillNucleusAreaWithinNeuriteSecond = 0.35;
    optionHandler.NucleusThreshold = 15;
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
    optionHandler.EdgeCompositNeuriteLowerThreshold = 20;
    optionHandler.EdgeCompositMinArea = 70;
    optionHandler.EdgeCompositeDistanceNucleusWhiteArea = 3;
    optionHandler.EdgeFillNucleusAreaWithinNeurite = 0.55;
     optionHandler.EdgeFillNucleusAreaWithinNeuriteSecond = 0.35;
    optionHandler.NucleusThreshold = 15;
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



function mapString = mapPositionToKeyString(yPos,xPos)
    mapString = [num2str(xPos) ';' num2str(yPos)];

    
function EdgeFillNeuronsAlgorithm(selectedWell, handles)
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
NucleusM = csvHandler.CellPosMatrix(selectedWell);
%Load neurite and nucleus image from file 
foldername = [imageHandler.Foldername '/ConvertedCellomics'];
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end
imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig.tif'];
imagePathNucleusSmall = [foldername '/' selectedWellLong 'NucleusSmall.tif'];
imagePathNeuriteBig = [foldername '/' selectedWellLong 'NeuriteBig.tif'];
imagePathNeuriteSmall = [foldername '/' selectedWellLong 'NeuriteSmall.tif'];

%Check if NeuronStatMatrix exists
if(~isa(neuronHandler.NeuronStatMatrix,'containers.Map'))
    neuronHandler.NeuronStatMatrix = containers.Map();
end
if(isKey(neuronHandler.NeuronStatMatrix,selectedWell))
    currentNeuronStat = neuronHandler.NeuronStatMatrix(selectedWell);
else
    currentNeuronStat = containers.Map();
end

imageHandler.ResizedNeuriteImage = imread(imagePathNeuriteSmall);
imageHandler.ResizedNucleusImage = imread(imagePathNucleusSmall);    
imageHandler.NeuriteImage = imread(imagePathNeuriteBig);
imageHandler.NucleusImage = imread(imagePathNucleusBig);

[sizeY sizeX]= size(imageHandler.NeuriteImage);

%Threshold Neurite Picture the known way
neuritePic = CutOutCircles(imageHandler.NeuriteImage,selectedWell,1,0,handles);

[neuritePic neuritePicStrong] = ThresholdPic(neuritePic,handles);
%Edge filter for Neurite Picture
%neuritePic = edge(neuritePic, 'nothinning');

SE = strel('disk', 3);
neuritePic = imdilate(neuritePic,SE);

nucleusPic = CutOutCircles(imageHandler.NucleusImage,selectedWell,1,0,handles);
%nucleusPic = zeros(7168,7168);
%[nucleusPic nucleusPicStrong] = ThresholdPic(nucleusPic,handles);
nucleusPic(nucleusPic<optionHandler.NucleusThreshold) = 0;%double(double(imageHandler.NucleusImage) ./ 255);
nucleusPic(nucleusPic>=optionHandler.NucleusThreshold) = 255;
nucleusPic = uint8(nucleusPic);
%imageHandler.BinaryImage = nucleusPic;
%imageHandler.ResizedBinaryImage = imresize(nucleusPic,0.1);
%nucleusPic = ThresholdPic(nucleusPic,0.9999999999,handles);<f
%nucleusPic = nucleusPic .* 255;
%Threshold Nucleus Picture
%nucleusPic = ThresholdPic(imageHandler.NucleusImage,handles);
fusedPic = imfuse(nucleusPic, neuritePic);
%figure(1);
%imshow(fusedPic);
neuronPointIndices = fusedPic(:,:,1) > 250 & fusedPic(:,:,2) > 250 & fusedPic(:,:,3) > 250;
%For every Nucleus:
%Check if Neurite with edges around.
%If yes: Count overlap of Nucleus and Neurites
%Count how many points of Nucleus are within the Neurite edges
EdgeFillNeuronPositionList = sparse(sizeY,sizeX);
EdgeFillNeuronPositionListSecond = sparse(sizeY,sizeX);
%stringCSVStatisticExport = sprintf('Position X;Position Y;Overlap (percent)\r\n');
[NonZeroY, NonZeroX] = find(NucleusM);
NonZeroY = uint16(NonZeroY);
NonZeroX = uint16(NonZeroX);

dbgYesCount = 0;
dbgNoCount = 0;
dbgNoStartFound = 0;

for i=1:numel(NonZeroY)    
    if(mod(i,10) == 0)
       disp([num2str(i) ' of ' num2str(numel(NonZeroY))]);
    end
    startX=0;
    startY=0;
    maxDistWhiteNuc = 0;
    for j=-maxDistWhiteNuc:maxDistWhiteNuc
        for k=-maxDistWhiteNuc:maxDistWhiteNuc
            currentX = NonZeroX(i) + k;
            currentY = NonZeroY(i) + j;
            %If point or one neighbour is white in Overlaypic.
            if(currentY < sizeY && currentX < sizeX && currentY>0 && currentX>0 && neuronPointIndices(currentY,currentX) > 0)
                startX= currentX;
                startY = currentY;
            end
        end
    end
    if(startX > 0)
        %FloodFill white Area from current point. Check also nucleus points
        stack = java.util.Stack();
        %Idea to be more efficient:
        %Cut out NucleusArea from startX - 50 until startX + 50
        %Do FloodFill Operation only on this subset
        if(startY-50>0 && startX-50>0 && startX+50 < sizeX && startY + 50 < sizeY)
            SubArea = neuronPointIndices(startY-50:startY+50,startX-50:startX+50);
            SubNucleusM = NucleusM(startY-50:startY+50,startX-50:startX+50);
            SubNucleusPic = nucleusPic(startY-50:startY+50,startX-50:startX+50);
            [maxY, maxX] = size(SubNucleusPic);
            %Recalculate startX and startY            
            stack.push([50 50]);
            %FloodFill SubArea as well as Binary NucleusPic
            %Count number of pixels of both during FloodFill operation
            [numberNeuritePixels, numberNucleusPixels] = FloodFillForEdgeFill(SubNucleusPic,SubArea,maxX,maxY,stack);
            %WhiteArea2 = imfill(SubArea,[50 50]);         
            %[yIndices xIndices] = find(WhiteArea2);
            %yIndices = yIndices + (startY-50);
            %xIndices = xIndices + (startX-50);
        else
            stack.push([startX startY]);
            WhiteArea = zeros(sizeY,sizeX);
            %WhiteArea2 = imfill(neuronPointIndices,[startY startX]);
            [numberNeuritePixels, numberNucleusPixels] = FloodFillForEdgeFill(nucleusPic,neuronPointIndices,sizeX,sizeY,stack);
             %Schwerpunkt aller Punkte berechnen.
            %[yIndices xIndices] = find(WhiteArea);
        end     
        %Check if Area has not too less and not too much pixels

        %Check relation between numberNeuritePixels and numberNucleusPixels
        
        %Add neurite independent of it's neuritePixels per
        %NucleusPixels. Save that value.
        %Data structure: Matrix with positions AND relation
        
        area = optionHandler.EdgeFillNucleusAreaWithinNeurite;
        areaSecond = optionHandler.EdgeFillNucleusAreaWithinNeuriteSecond;
        edgeFillMarked = 0;
        edgeFillMarkedSecond = 0;
        if(numberNeuritePixels/numberNucleusPixels >= area) 
            EdgeFillNeuronPositionList(NonZeroY(i), NonZeroX(i)) = 1;
            EdgeFillNeuronPositionListSecond(NonZeroY(i), NonZeroX(i)) = 1;
            edgeFillMarked=1;
            edgeFillMarkedSecond=1;
        elseif(numberNeuritePixels/numberNucleusPixels >= areaSecond) 
            EdgeFillNeuronPositionListSecond(NonZeroY(i), NonZeroX(i)) = 1;   
            edgeFillMarkedSecond=1;
        end
        %stringCSVStatisticExport = [stringCSVStatisticExport num2str(NonZeroX(i)) ';' num2str(NonZeroY(i)) ';' num2str(numberNeuritePixels/numberNucleusPixels) ';' sprintf('\r\n')];
        %EdgeFillNeuronPositionList(NonZeroY(i), NonZeroX(i),2) = (numberNeuritePixels / numberNucleusPixels)*100;
        %Check if position is marked as Manual
        NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
        if(NeuronManualM(NonZeroY(i),NonZeroX(i)))
            manualMarked=1;
        else
            manualMarked=0;
        end
        %Check if NeuronStatMatrix contains already key for selected Position
        key = mapPositionToKeyString(NonZeroY(i),NonZeroX(i));
        if(isKey(currentNeuronStat,key))
           list = currentNeuronStat(key);
        else
            list = zeros(8,1);            
        end
        list(4)=manualMarked;
        list(2) = edgeFillMarked;
        list(3) = edgeFillMarkedSecond;
        list(6) = numberNeuritePixels/numberNucleusPixels;
        currentNeuronStat(key) = list;
    else
        dbgNoStartFound = dbgNoStartFound +1;
    end    
end

%Check if Neuronhandler already has Dictionary of NeuronPositionsEdgeComposite
if(numel(neuronHandler.NeuronPositionsEdgeFill) ==0)
    neuronHandler.NeuronPositionsEdgeFill = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFillSecond) ==0)
    neuronHandler.NeuronPositionsEdgeFillSecond = containers.Map();
end
neuronHandler.NeuronStatMatrix(selectedWell) = currentNeuronStat;
%Save stringCSVStatisticExport to file
%filepath = [imageHandler.Foldername '/ConvertedCellomics/EdgeFillStat_' selectedWell '.csv'];
%fileID = fopen(filepath,'w');
%fprintf(fileID,'%s',stringCSVStatisticExport);
%fclose(fileID);
neuronHandler.NeuronPositionsEdgeFill(selectedWell) = EdgeFillNeuronPositionList;
neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = EdgeFillNeuronPositionListSecond;
handles.NeuronCoordinates = neuronHandler;
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles);    


% --------------------------------------------------------------------
function EdgeFillNeurons_Callback(hObject, eventdata, handles)
% hObject    handle to EdgeFillNeurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
EdgeFillNeuronsAlgorithm(selectedWell,handles);

function [numberNeuritePixels, numberNucleusPixels] = FloodFillForEdgeFill(NucleusBinaryPic,NeuriteOverlayPic,maxX,maxY,stack)
[sizeY, sizeX] = size(NucleusBinaryPic);
AreaResult = logical(sparse(sizeY,sizeX));
numberNeuritePixels = 0;
numberNucleusPixels = 0;
while(stack.isEmpty() == 0)
    vector = stack.pop();
    x=vector(1);
    y=vector(2);
    if(NucleusBinaryPic(y,x) >= 1 && AreaResult(y,x) == 0 && x-1 > 0 && y-1 > 0 && y+1<= maxY && x+1 <= maxX)
            AreaResult(y,x) = 1;  
            numberNucleusPixels = numberNucleusPixels+1;
            if(NeuriteOverlayPic(y,x) >= 1)
                numberNeuritePixels = numberNeuritePixels + 1;
            end
            stack.push([x y+1]);
            stack.push([x-1 y]);
            stack.push([x y-1]);
            stack.push([x+1 y]);
            %8 Neighbourhood
            stack.push([x+1 y+1]);
            stack.push([x-1 y-1]);
            stack.push([x-1 y+1]);
            stack.push([x+1 y-1]);
    elseif(NucleusBinaryPic(y,x) >= 1 && AreaResult(y,x) == 0)
            AreaResult(y,x)=1;
    end
end
result=AreaResult;


function [returnPic returnPicMin]=ThresholdPic(inputImage,handles)

    handles = guidata(handles.figure1);
    optionHandler = handles.OptionHandler;
    imageHandler = handles.ImageHandler;
    neuronHandler = handles.NeuronCoordinates;
    selectedWellNumber = get(handles.lbWell,'Value');
    wellList = get(handles.lbWell, 'string');
    selectedWell = wellList{selectedWellNumber};
    densityWidth = optionHandler.MigrationDistanceDensityImageXSize;
    densityHeight = optionHandler.MigrationDistanceDensityImageYSize;
    [sizeY sizeX]= size(imageHandler.NeuriteImage);
    ringNumber = optionHandler.DensityDistributionRingNumber;
     
%Cut out dark pixels from NeuriteImage:     
%brightNeurites(brightNeurites > 150) = 0;
%figure(1);
%imshow(brightNeurites);

%0.5 Contrast stretching
stretchlimLow = stretchlim(inputImage,[0.975 0.9999]);
stretchlimLow = neuronHandler.StretchlimResult;
brightNeurites = imadjust(inputImage,stretchlimLow,[0 1]);
%figure(2);
%imshow(brightNeurites);

%1. Median Filter
brightNeurites = medfilt2(brightNeurites);

%2. Resharp image
unsharpFilter = fspecial('unsharp');
brightNeurites = imfilter(brightNeurites,unsharpFilter);
%figure(2);
%imshow(brightNeurites);

%3. Thresholding low
%level=graythresh(brightNeurites);
brightNeurites = imcomplement(brightNeurites);
%Calculate Threshold level for every Well separately.
%level = th_intermodes(brightNeurites)./255;
%level = isodata(brightNeurites);
%level = level + 0.1
if(optionHandler.SkeletonThresholdMethod == 1)
    level = isodata(brightNeurites) + optionHandler.SkeletonNeuriteThresholdLow;
    if(level>1)
        level=isodata(brightNeurites);
    elseif(level<0)
        level=0;
    end
elseif(optionHandler.SkeletonThresholdMethod == 2)
    level=neuronHandler.IsodataResult + optionHandler.SkeletonNeuriteThresholdLow;
elseif(optionHandler.SkeletonThresholdMethod == 3)
    level = optionHandler.SkeletonNeuriteThresholdLow;
end
levelMin = level-optionHandler.SkeletonNeuriteThresholdHardDistance;
if(levelMin>1)
    levelMin=1;
elseif(levelMin<0)
    levelMin=0;
end

BW = im2bw(brightNeurites, level);
BW = imcomplement(BW);

%brightNeurites = imcomplement(brightNeurites);
BWMin = im2bw(brightNeurites,levelMin);
BWMin = imcomplement(BWMin);
%figure(3);
%imshow(BW);

%4. Remove particles
LB = 75;
UB = 4000;
%CuttedIndices = BW;
BW = xor(bwareaopen(BW,LB),  bwareaopen(BW,UB));
%figure(4);
%imshow(BW);
returnPic = BW;
returnPicMin = BWMin;

function imageMap=CutOutCircles(imageMap,selectedWell,CutInnerCircle, CutOuterCircle,handles)
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
selectedWellNumber = get(handles.lbWell,'Value');
%wellList = get(handles.lbWell, 'string');
%selectedWell = wellList{selectedWellNumber};
densityWidth = 50;
densityHeight = 50;
[sizeY sizeX]= size(imageHandler.NeuriteImage);
ringNumber = optionHandler.DensityDistributionRingNumber;
path = [imageHandler.Foldername '/ConvertedCellomics/MarkerPointCoordinates-' selectedWell '.mat'];
%Get density distribution to exclude too dense areas
%if(saveCircle)    
  if(exist(path,'file'))    
    ringNumber = optionHandler.DensityDistributionRingNumber;
    load(path);
    %hInner = impoly(handles.axes2,double(markerPointCoordinates('0')./10));
    %previousMask = createMask(hInner);
    %for i=1:ringNumber  
    %  hCurrent = impoly(handles.axes2,double(markerPointCoordinates(num2str(i*10))./10));
    %end
  else
      
     [filterDistance nonFilterDistance SphereArea markerPointCoordinates] = calculateMigrationDistance(selectedWell, densityWidth, densityHeight, handles); 
     %Save hInner to file
     foldername = imageHandler.Foldername;
     subfoldername = [foldername '/ConvertedCellomics'];
     %if(markerPointCoordinates~=0)
        save('-v7.3',strcat(subfoldername,'\MarkerPointCoordinates-',selectedWell),'markerPointCoordinates', 'filterDistance', 'nonFilterDistance');
     %end
  end
%else
%    [filterDistance nonFilterDistance SphereArea NucleusArea64 markerPointCoordinates] = calculateMigrationDistance(selectedWell, densityWidth, densityHeight, handles);   
%end
%MarkerPointCoordinates are available. Get second ring:
if(imageHandler.ZoomState)
    %Draw full pic
    imshow(imageHandler.ResizedNeuriteImage);
end
if(numel(markerPointCoordinates('10')./10) == 0)
    markerPointCoordinates=0;
end
if(markerPointCoordinates ~=0 && CutInnerCircle == 1 && CutOuterCircle == 1)
    
    hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./10));
    innerMask = logical(createMask(hInner));
    delete(hInner);
    hOuter = impoly(handles.axes2,double(markerPointCoordinates(num2str(ringNumber * 10))./10));
    outerMask = logical(createMask(hOuter));
    delete(hOuter);
    currentRing = logical(outerMask-innerMask); 
elseif(CutInnerCircle==1 && CutOuterCircle == 0 && markerPointCoordinates ~=0)
    hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./10));
    innerMask = logical(createMask(hInner));
    delete(hInner);
    currentRing =logical(ones(int32(sizeY/10),int32(sizeX/10)) - innerMask);
elseif(markerPointCoordinates ~=0);
    hOuter = impoly(handles.axes2,double(markerPointCoordinates(num2str(ringNumber * 10))./10));
    outerMask = logical(createMask(hOuter));
    delete(hOuter);
    currentRing = outerMask;
end
if(markerPointCoordinates ~=0)
  currentRing = imresize(currentRing, [sizeY, sizeX]);

if(size(imageMap,1) ~= sizeY)
  for(i=1:size(imageMap,1))
    image = imageMap(num2str(i));
    imageMap(num2str(i)) = uint8(image) .* uint8(currentRing);
  end
else
  imageMap = uint8(imageMap) .* uint8(currentRing);
end
else 
    if(size(imageMap,1) ~= sizeY)
  for(i=1:size(imageMap,1))
    image = imageMap(num2str(i));
    imageMap(num2str(i)) = uint8(image);
  end
else
  imageMap = uint8(imageMap);
end
end

function AreaResult = FloodFillNucArea(NucleusBinaryPic,yStart,xStart,NucleusM,csvHandler)
[sizeY, sizeX] = size(NucleusBinaryPic);
maxY=sizeY;
maxX=sizeX;
AreaResult = logical(zeros(sizeY,sizeX));
numberNeuritePixels = 0;
numberNucleusPixels = 0;
stack = java.util.Stack();
stack.push([xStart yStart]);
%Get from 
while(stack.isEmpty() == 0)
    vector = stack.pop();
    x=vector(1);
    y=vector(2);
    if(pdist([yStart xStart;y x]) <=5)
        add=1;
    else
        [nucCol nucRow euclidDist] = csvHandler.FindNucleusForNeuron(x, y, NucleusM ,20);
        if(nucCol == xStart && nucRow == yStart)
            add=1;
        else
            add=0;
        end
    end
    if(NucleusBinaryPic(y,x) >= 1 && AreaResult(y,x) == 0 && x-1 > 0 && y-1 > 0 && y+1<= maxY && x+1 <= maxX && add)
        
        if(add)
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
    elseif(NucleusBinaryPic(y,x) >= 1 && AreaResult(y,x) == 0 && add)
            AreaResult(y,x)=1;
        end
    end
end
result=AreaResult;

function SkeletonizationNeuronsAlgorithm(selectedWell, handles)
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
[sizeY sizeX]= size(imageHandler.NeuriteImage);
csvHandler.NeuriteLengthMatrix = sparse(sizeY, sizeX);
SkelDeletedMat = logical(sparse(sizeY,sizeX));
SkelNeurons = logical(sparse(sizeY,sizeX));
NeuronPositionsEdgeFillNeurite = logical(sparse(sizeY,sizeX));
ringNumber = optionHandler.DensityDistributionRingNumber;
selectedWellLong=selectedWell;
foldername = [imageHandler.Foldername '/ConvertedCellomics'];
neuriteAreaCount=0;

%Check if NeuronStatMatrix exists
if(~isa(neuronHandler.NeuronStatMatrix,'containers.Map'))
    neuronHandler.NeuronStatMatrix = containers.Map();
end
if(isKey(neuronHandler.NeuronStatMatrix,selectedWell))
    currentNeuronStat = neuronHandler.NeuronStatMatrix(selectedWell);
else
    currentNeuronStat = containers.Map();
end

if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end
imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig.tif'];
imagePathNucleusSmall = [foldername '/' selectedWellLong 'NucleusSmall.tif'];
imagePathNeuriteBig = [foldername '/' selectedWellLong 'NeuriteBig.tif'];
imagePathNeuriteSmall = [foldername '/' selectedWellLong 'NeuriteSmall.tif'];

imageHandler.ResizedNeuriteImage = imread(imagePathNeuriteSmall);
 
imageHandler.NeuriteImage = imread(imagePathNeuriteBig);
%composedSkeletonImage = zeros(sizeY, sizeX);
imageHandler.BinaryImage = logical(zeros(sizeY, sizeX));
imageHandler.SkeletonImage = logical(zeros(sizeY, sizeX));
%Get density distribution to exclude too dense areas
%filepath = [imageHandler.Foldername '/ConvertedCellomics/' selectedWell 'BinaryCut.tif'];
%if(~exist(filepath,'file'))
brightNeurites = CutOutCircles(imageHandler.NeuriteImage,selectedWell,1,0,handles);
nucleusPic = CutOutCircles(imageHandler.NucleusImage,selectedWell,1,0,handles);
nucleusPic(nucleusPic<optionHandler.NucleusThreshold-5) = 0;%double(double(imageHandler.NucleusImage) ./ 255);
nucleusPic(nucleusPic>=optionHandler.NucleusThreshold-5) = 1;
nucleusPic = logical(nucleusPic);
imageHandler.NucleusImage = imread(imagePathNucleusBig);
imageHandler.ResizedNucleusImage = imread(imagePathNucleusSmall);   
[BW BWStrong] = ThresholdPic(brightNeurites,handles);
%BW = imfill(BW);
%ToDo: Fill holes and imdilate
SE = strel('disk', 1);
BW = imdilate(BW,SE);


%ToDo: Save Neurite objects with single pictures to trace
%Check if already a Neuron by EdgeFillNeurons in this area available
%If not, trace Neurite and check if Neuron possible on one or
%another side  

[whitePointRows, whitePointCols] = find(BW);
alreadyVisited = uint32(ones(size(BW,1), size(BW,2)));
neuriteAreasDict = containers.Map();
i=0;
BW2= logical(zeros(size(BW,1), size(BW,2)));
secondThreshold = 20;

%FloodFill for Points, which were not visited before%As long as there are white point indices, which were not already visited
while(numel(whitePointRows) > 0)
    %Get first whitepointindex <> 0
    %FloodFill    
    stack = java.util.Stack();
    stack.push([whitePointRows(1) whitePointCols(1)]);
    currentNeuriteArea = logical(sparse(zeros(size(BW,1), size(BW,2))));
    strongCount = 0;
    visited=0;
    while(stack.isEmpty() == 0)
      vector = stack.pop();
      yPos=vector(1);
      xPos=vector(2);       
      %[visitedRows visitedCols] = ind2sub(7168,locations);
      
      if(((BW(yPos,xPos) == 1) || (brightNeurites(yPos,xPos) > secondThreshold)) && alreadyVisited(yPos, xPos) == 1 &&  yPos + 1 < size(BW,1) && xPos + 1 < size(BW,2) && xPos - 1 > 0 && yPos - 1 > 0)
            visited=1;
            if(BWStrong(yPos,xPos) == 1 )
                strongCount = strongCount+1;
            end
            BW2(yPos,xPos)= 1;
            currentNeuriteArea(yPos,xPos)= 1;
            stack.push([yPos+1 xPos]);
            stack.push([yPos xPos-1]);
            stack.push([yPos-1 xPos]);
            stack.push([yPos xPos+1]);
            %only 4 Neighbourhood!
            stack.push([yPos+1 xPos+1]);
            stack.push([yPos-1 xPos-1]);
            stack.push([yPos+1 xPos-1]);
            stack.push([yPos-1 xPos+1]);
      end
      alreadyVisited(yPos, xPos) = 0; 
    end    
    if(visited==0)
        whitePointRows(1)=[];
        whitePointCols(1)=[];
        continue;
    end
    
    newNeurite = Neurite();
    newNeuriteBig = Neurite();
    [rowIndices colIndices] = find(currentNeuriteArea);
    newNeurite.cutYPosStart = min(rowIndices) - 2;
    newNeurite.cutXPosStart = min(colIndices) - 2;    
     
    newNeurite.cutYPosEnd = max(rowIndices) + 2;
    newNeurite.cutXPosEnd = max(colIndices) + 2;
    
    if(newNeurite.cutYPosStart <= 0)
        newNeurite.cutYPosStart = 1;
    end
    if(newNeurite.cutXPosStart <= 0)
        newNeurite.cutXPosStart = 1;
    end
    if(newNeurite.cutXPosEnd > size(currentNeuriteArea,2))
        newNeurite.cutXPosEnd = size(currentNeuriteArea,2);
    end
    if(newNeurite.cutYPosEnd > size(currentNeuriteArea,1))
        newNeurite.cutYPosEnd = size(currentNeuriteArea,1);
    end
    
    newNeuriteBig.cutYPosStart = min(rowIndices) - 2;
    newNeuriteBig.cutXPosStart = min(colIndices) - 2;    
     
    newNeuriteBig.cutYPosEnd = max(rowIndices) + 2;
    newNeuriteBig.cutXPosEnd = max(colIndices) + 2;
    
    if(newNeuriteBig.cutYPosStart <= 0)
        newNeuriteBig.cutYPosStart = 1;
    end
    if(newNeuriteBig.cutXPosStart <= 0)
        newNeuriteBig.cutXPosStart = 1;
    end
    if(newNeuriteBig.cutXPosEnd > size(currentNeuriteArea,2))
        newNeuriteBig.cutXPosEnd = size(currentNeuriteArea,2);
    end
    if(newNeuriteBig.cutYPosEnd > size(currentNeuriteArea,1))
        newNeuriteBig.cutYPosEnd = size(currentNeuriteArea,1);
    end   
    
    
    currentNeuriteArea=currentNeuriteArea(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
    currentNeuriteArea = 1-currentNeuriteArea;
    
    %Copy Neurite Area to composed Skeleton image
    newNeurite.image = currentNeuriteArea;
    newNeuriteBig.image = currentNeuriteArea;
    
    
    try
        newNeurite.CalculateBaiSkeleton(6);
    catch
        while(numel(whitePointRows) > 0 && alreadyVisited(whitePointRows(1), whitePointCols(1)) == 0)
            whitePointRows(1) = [];
            whitePointCols(1) = [];
        end
        continue;
    end
    newNeuriteBig.CalculateBaiSkeleton(12);   

    newNeurite.FindBranchingCrossingPoints();
    newNeurite.KillShortBranches();
    newNeuriteBig.FindBranchingCrossingPoints();
    newNeuriteBig.KillShortBranches();
    imageHandler.SkeletonImage(newNeuriteBig.cutYPosStart:newNeuriteBig.cutYPosEnd,newNeuriteBig.cutXPosStart:newNeuriteBig.cutXPosEnd)=newNeuriteBig.skeletonImage;
    %if(numel(newNeurite.xEndpoints)>=4)
    %    newNeurite.SplitUp();
    %end
    
    %As a result there will be 2 up to 3 remaining endpoints
    %If there are 2 endpoints: There has to be one Neuron: Check if available
    %through EdgeFill-Algorithm. If yes, return.
    %If more than one, check if one has to be killed
    %If no one, check iteratively, if Neuron can be added through lowering the
    %EdgeFillThreshold.
    %If this is not possible, try to find neuron by extension of vector from
    %EndPoints.

    %Get for every Neurite the next Nucleus to the Endpoints
    %For every Neurite: Add next Nucleus to Endpoint to Resultlist

    %Next step: Check all distances from endpoint to Endpoint
    %If there are only two endpoints and distance is to short -> Kill whole
    %Neurite and if there is a Neuron by EdgeFill -> Kill it.
    newNeurite.CalculateNeuriteLength();
    newNeuriteBig.CalculateNeuriteLength();
    %Check max distance between Endpoints
    maxDistEndpoints = 0;
    for(l=1:numel(newNeuriteBig.savedXEndpoints))
        for(m=1:numel(newNeuriteBig.savedXEndpoints))
            %Check distance between both endpoints
            currentDistVec = [newNeuriteBig.savedYEndpoints(l) newNeuriteBig.savedXEndpoints(l);newNeuriteBig.savedYEndpoints(m) newNeuriteBig.savedXEndpoints(m)];
            dist = pdist(currentDistVec,'euclidean');
            if(dist>maxDistEndpoints)
                maxDistEndpoints = dist;
            end
        end
    end
   
    
    if(newNeuriteBig.neuriteLength < optionHandler.SkeletonMinNeuriteLength || maxDistEndpoints < optionHandler.SkeletonMinNeuriteLength || strongCount == 0)
        %ToDo: Delete potential Neuclei on Neurite
        %Save position of EdgeFilled Nuclei, which should be deleted
        edgeFillPositions = neuronHandler.NeuronPositionsEdgeFill(selectedWell);
        edgeFillPositionsSecond = neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
        edgeFillPositionsSecond = edgeFillPositionsSecond(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
        %Add Neurons to deleted ones
        [yPoses xPoses] = find(edgeFillPositionsSecond);
        for(f=1:numel(yPoses))
            SkelDeletedMat(newNeurite.cutYPosStart+yPoses(f)-1, newNeurite.cutXPosStart+xPoses(f)-1) = 1;
            NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
            if(NeuronManualM(newNeurite.cutYPosStart+yPoses(f)-1,newNeurite.cutXPosStart+xPoses(f)-1) == 1)
                manualMarked=1;
            else
                manualMarked=0;
            end
            key = mapPositionToKeyString(newNeurite.cutYPosStart+yPoses(f)-1,newNeurite.cutXPosStart+xPoses(f)-1);
            if(isKey(currentNeuronStat,key))
                list = currentNeuronStat(key);
            else
                list = zeros(8,1);
            end
            list(5) = 1;
            list(7) = newNeuriteBig.neuriteLength;
            list(8) = maxDistEndpoints;
            list(4) = manualMarked;
            currentNeuronStat(key) = list;
        end
        if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
            neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
        end
        neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = SkelDeletedMat;
    else
        %ToDo: Create mapping between Neurites and Neurons
        %To Kill: All Neurons without regarding Neurite
        %First, check if already marked Neuron by EdgeFill-Algorithm available.
        edgeFillPositions = neuronHandler.NeuronPositionsEdgeFill(selectedWell);
        edgeFillPositions = edgeFillPositions(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
        [yEdgeFills xEdgeFills] = find(edgeFillPositions);
        if(numel(yEdgeFills) >= 1)
            %done
            %Save Mapping between Neuron and Neurite.     
            %Check which Neuron has max. distance to next endpoint
            if(numel(yEdgeFills)==1)
                SkelNeurons(newNeurite.cutYPosStart+yEdgeFills(1)-1,newNeurite.cutXPosStart+xEdgeFills(1)-1)=1;
                NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
                %Add position of underlying nucleus to imageHandler.BinaryImage
                NucArea = nucleusPic(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                yStart=yEdgeFills(1);
                xStart=xEdgeFills(1);
                NucleusM = csvHandler.CellPosMatrix(selectedWell);
                NucleusM = NucleusM(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                %TO COMMENT IN
                 AreaResult = FloodFillNucArea(NucArea,yStart,xStart,NucleusM,csvHandler);
                imageHandler.BinaryImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd) = imageHandler.BinaryImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd) + (1-currentNeuriteArea) + AreaResult;         
                neuriteAreaCount=neuriteAreaCount+nnz(1-currentNeuriteArea);                
                if(NeuronManualM(newNeurite.cutYPosStart+yEdgeFills(1)-1,newNeurite.cutXPosStart+xEdgeFills(1)-1) == 1)
                    manualMarked=1;
                else
                    manualMarked=0;
                end
                key = mapPositionToKeyString(newNeurite.cutYPosStart+yEdgeFills(1)-1,newNeurite.cutXPosStart+xEdgeFills(1)-1);
                if(isKey(currentNeuronStat,key))
                    list = currentNeuronStat(key);
                else
                    list = zeros(8,1);
                end
                list(1) = 1;
                list(7) = newNeuriteBig.neuriteLength;
                list(4) = manualMarked;
                list(8) = maxDistEndpoints;
                currentNeuronStat(key) = list;
            else            
           %     minEndPointDist = 999;
           %     minEndPointDistInd = 0;
                markedForDeletion = ones(numel(yEdgeFills));
                for(counter=1:numel(yEdgeFills))     
                    %Check distance between different Neurons. %If distance
                    %too small, kill Neuron with farest distance to
                    %endpoint
                    for(counter2=1:numel(yEdgeFills))
                        dist = pdist([yEdgeFills(counter) xEdgeFills(counter);yEdgeFills(counter2) xEdgeFills(counter2)]);
                        if(dist<13 && dist > 0 && markedForDeletion(counter2) == 1)
                            %Mark for deletion
                            markedForDeletion(counter) = 0;
                        end
                    end
           %         endPointDist=999;
           %         endPointInd=0;
           %         for(q=1:numel(newNeurite.savedXEndpoints))
           %             if(pdist([newNeurite.cutYPosStart+newNeurite.savedYEndpoints(q)-1 newNeurite.cutXPosStart+newNeurite.savedXEndpoints(q)-1;newNeurite.cutYPosStart+yEdgeFills(counter)-1,newNeurite.cutXPosStart+xEdgeFills(counter)-1]) < endPointDist)
           %                 endPointDist = pdist([newNeurite.cutYPosStart+newNeurite.savedYEndpoints(q)-1 newNeurite.cutXPosStart+newNeurite.savedXEndpoints(q)-1;newNeurite.cutYPosStart+yEdgeFills(counter)-1,newNeurite.cutXPosStart+xEdgeFills(counter)-1]);
           %                 endPointInd = q;
           %             end
           %         end
           %         if(endPointDist < minEndPointDist)
           %             minEndPointDist = endPointDist;
           %             minEndPointDistInd=counter;
           %         end
                end
                for(counter=1:numel(yEdgeFills))
                    NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
                    if(NeuronManualM(newNeurite.cutYPosStart+yEdgeFills(counter)-1,newNeurite.cutXPosStart+xEdgeFills(counter)-1) == 1)
                        manualMarked=1;
                    else
                        manualMarked=0;
                    end
                    if(markedForDeletion(counter) == 1)
                        SkelNeurons(newNeurite.cutYPosStart+yEdgeFills(counter)-1,newNeurite.cutXPosStart+xEdgeFills(counter)-1)=1;
                        %Add position of underlying nucleus to imageHandler.BinaryImage
                        NucArea = nucleusPic(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                        yStart=yEdgeFills(counter);
                        xStart=xEdgeFills(counter);
                        NucleusM = csvHandler.CellPosMatrix(selectedWell);
                        NucleusM = NucleusM(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                        AreaResult = FloodFillNucArea(NucArea,yStart,xStart,NucleusM,csvHandler);
                        imageHandler.BinaryImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd) = imageHandler.BinaryImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd) + (1-currentNeuriteArea) + AreaResult;
                        neuriteAreaCount=neuriteAreaCount+nnz(1-currentNeuriteArea);
                        key = mapPositionToKeyString(newNeurite.cutYPosStart+yEdgeFills(counter)-1,newNeurite.cutXPosStart+xEdgeFills(counter)-1);
                        if(isKey(currentNeuronStat,key))
                            list = currentNeuronStat(key);
                        else
                            list = zeros(8,1);
                        end
                        list(1) = 2;
                        list(7) = newNeuriteBig.neuriteLength;
                        list(8) = maxDistEndpoints;
                        list(4) = manualMarked;
                        currentNeuronStat(key) = list;
                    else
                        key = mapPositionToKeyString(newNeurite.cutYPosStart+yEdgeFills(counter)-1,newNeurite.cutXPosStart+xEdgeFills(counter)-1);
                        if(isKey(currentNeuronStat,key))
                            list = currentNeuronStat(key);
                        else
                            list = zeros(8,1);
                        end
                        list(5) = 2;
                        list(7) = newNeuriteBig.neuriteLength;
                        list(8) = maxDistEndpoints;
                        list(4) = manualMarked;
                        currentNeuronStat(key) = list;
                    end
                end
            end
        else
            %Get Neuron with highest thresh in this area.
            %Check if Neuron between 0.55 and 1

            %If yes: Take this one
            %Delete all other potential neurons
            %If no: Do path extension            

            %Check for other potential Neurons in Edge Fill Neurons by 2nd
            %threshold near endpoints. Get first one with biggest thresh.
            edgeFillPositionsSecond = neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
            edgeFillPositionsSecond = edgeFillPositionsSecond(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
            %[maxA,ind] = max(A(:));
            [yInd xInd] = find(edgeFillPositionsSecond);
            if(numel(yInd) > 0)
                %Add Neuron to SkelNeurons
                SkelNeurons(newNeurite.cutYPosStart+yInd(1)-1, newNeurite.cutXPosStart+xInd(1)-1) = 1;
                NucArea = nucleusPic(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                yStart=yInd(1);
                xStart=xInd(1);
                NucleusM = csvHandler.CellPosMatrix(selectedWell);
                NucleusM = NucleusM(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                AreaResult = FloodFillNucArea(NucArea,yStart,xStart,NucleusM,csvHandler);
                imageHandler.BinaryImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd) = imageHandler.BinaryImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd) + (1-currentNeuriteArea) + AreaResult;
                neuriteAreaCount=neuriteAreaCount+nnz(1-currentNeuriteArea);
                NeuronPositionsEdgeFillNeurite(newNeurite.cutYPosStart+yInd(1)-1, newNeurite.cutXPosStart+xInd(1)-1) = 1;
                %Save Neurite length in Matrix
                csvHandler.NeuriteLengthMatrix(newNeurite.cutYPosStart+yInd(1)-1, newNeurite.cutXPosStart+xInd(1)-1) = newNeuriteBig.neuriteLength;
                NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
                if(NeuronManualM(newNeurite.cutYPosStart+yInd(1)-1,newNeurite.cutXPosStart+xInd(1)-1) == 1)
                    manualMarked=1;
                else
                    manualMarked=0;
                end
                key = mapPositionToKeyString(newNeurite.cutYPosStart+yInd(1)-1,newNeurite.cutXPosStart+xInd(1)-1);
                if(isKey(currentNeuronStat,key))
                   list = currentNeuronStat(key);
                else
                   list = zeros(8,1);
                end
                list(1) = 2;
                list(7) = newNeuriteBig.neuriteLength;
                list(8) = maxDistEndpoints;
                list(4) = manualMarked;
                currentNeuronStat(key) = list;
            else 
                minDist = 100;
                nucColSaved = 0;
                nucRowSaved = 0;
                %Centroid aller Endpunkte berechnen
                yCen=0;
                xCen=0;                
                for(i=1:numel(newNeuriteBig.xEndpoints))
                    yCen = yCen + newNeuriteBig.yEndpoints(i);
                    xCen = xCen + newNeuriteBig.xEndpoints(i);
                end
                yCen = yCen/numel(newNeuriteBig.yEndpoints);
                xCen = xCen/numel(newNeuriteBig.xEndpoints);
                for(i=1:numel(newNeurite.xEndpoints))
                    %Set startPoint ~6 px back to balance point of neurite                    
                    startPoint = [newNeurite.yEndpoints(i) newNeurite.xEndpoints(i)];
                    %stepY=yCen-startPoint(1);
                    %stepX=xCen-startPoint(2);
                    %stepY = stepY/4;
                    %stepX = stepX/4;
                    stepY = 0;
                    stepX = 0;
                    startPoint(1)=startPoint(1)+stepY;
                    startPoint(2)=startPoint(2)+stepX;
                    startPoint(1) = round(startPoint(1));
                    startPoint(2) = round(startPoint(2));
                    [nucCol nucRow euclidDist] = csvHandler.FindNucleusForNeuron(newNeurite.cutXPosStart+startPoint(2), newNeurite.cutYPosStart+startPoint(1), csvHandler.CellPosMatrix(selectedWell),20);
                    if(euclidDist<minDist)
                      minDist = euclidDist;
                      nucColSaved = nucCol;
                      nucRowSaved = nucRow;
                    end
                end
                if(nucRowSaved <= 0)
                   %Check also ot killed saved endpoints.
                   for(i=1:numel(newNeurite.savedXEndpoints))
                    %Set startPoint 6 px back to balance point of neurite
                    startPoint = [newNeurite.savedYEndpoints(i) newNeurite.savedXEndpoints(i)];
                    %stepY=yCen-startPoint(1);
                    %stepX=xCen-startPoint(2);
                    %stepY = stepY/4;
                    %stepX = stepX/4;
                    stepY = 0;
                    stepX = 0;
                    startPoint(1)=startPoint(1)+stepY;
                    startPoint(2)=startPoint(2)+stepX;
                    startPoint(1) = round(startPoint(1));
                    startPoint(2) = round(startPoint(2));
                    [nucCol nucRow euclidDist] = csvHandler.FindNucleusForNeuron(newNeurite.cutXPosStart+startPoint(2), newNeurite.cutYPosStart+startPoint(1), csvHandler.CellPosMatrix(selectedWell),20);
                    if(euclidDist<minDist)
                      minDist = euclidDist;
                      nucColSaved = nucCol;
                      nucRowSaved = nucRow;
                    end
                   end
                end
                if(nucRowSaved > 0)
                    SkelNeurons(nucRowSaved,nucColSaved)=1;
                    NeuronPositionsEdgeFillNeurite(nucRowSaved, nucColSaved) = 1;
                    neuriteAreaCount=neuriteAreaCount+nnz(1-currentNeuriteArea);
                    %NucArea = nucleusPic(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);                   
                    NucleusM = csvHandler.CellPosMatrix(selectedWell);                
                    AreaResult = FloodFillNucArea(nucleusPic,nucRowSaved,nucColSaved,NucleusM,csvHandler);
                    imageHandler.BinaryImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd) = imageHandler.BinaryImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd) + (1-currentNeuriteArea);
                    imageHandler.BinaryImage = imageHandler.BinaryImage + AreaResult;
                %Save Neurite length in Matrix
                    csvHandler.NeuriteLengthMatrix(nucRowSaved, nucColSaved) = newNeuriteBig.neuriteLength;
                    NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
                    if(NeuronManualM(nucRowSaved,nucColSaved) == 1)
                        manualMarked=1;
                    else
                        manualMarked=0;
                    end
                    key = mapPositionToKeyString(nucRowSaved,nucColSaved);
                    if(isKey(currentNeuronStat,key))
                        list = currentNeuronStat(key);
                    else
                        list = zeros(8,1);
                    end
                    list(1) = 3;
                    list(7) = newNeuriteBig.neuriteLength;
                    list(8) = maxDistEndpoints;
                    list(4) = manualMarked;
                    currentNeuronStat(key) = list;
                end
                %newNeurite.PlotSkeleton();    
                %if(newNeurite.cutYPosStart-100 > 0 && newNeurite.cutYPosEnd+100 < 7168)
                    %fusedPic = imfuse(imageHandler.NucleusImage(newNeurite.cutYPosStart-100:newNeurite.cutYPosEnd+100,newNeurite.cutXPosStart-100:newNeurite.cutXPosEnd+100),imageHandler.NeuriteImage(newNeurite.cutYPosStart-100:newNeurite.cutYPosEnd+100,newNeurite.cutXPosStart-100:newNeurite.cutXPosEnd+100));
                    %figure(2);
                    %imshow(fusedPic);
                    %if(nucRowSaved > 0)
                    %    hold on;
                    %    plot(nucColSaved-newNeurite.cutXPosStart+100, nucRowSaved-newNeurite.cutYPosStart+100, '.r');
                    %    hold off;
                    %end
                    %close(figure(2));
                %end
                %close(figure(1));
            end
            %Check if added Neuron is far enough away from Branching points 
            
        end
    end

    %newNeurite.PlotSkeleton();    
    %close(figure(1));
    %neuriteAreasDict(num2str(i)) = newNeurite;
    i=i+1;
    %if(mod(i,10) == 0)
        disp([num2str(i) ' of ' num2str(numel(whitePointRows))]);
    %end
    while(numel(whitePointRows) > 0 && alreadyVisited(whitePointRows(1), whitePointCols(1)) == 0)
        whitePointRows(1) = [];
        whitePointCols(1) = [];
    end   
end
foldername = [imageHandler.Foldername '/ConvertedCellomics'];   
imagePathBinary = [foldername '/' selectedWellLong 'Binary.tif'];
imagePathSkeleton = [foldername '/' selectedWellLong 'Skeleton.tif'];
pathNeuriteLength = [imageHandler.Foldername '/ConvertedCellomics/NeuriteLength-' selectedWell '.mat'];
neuriteLengthMatrix = csvHandler.NeuriteLengthMatrix;
save(pathNeuriteLength,'neuriteLengthMatrix');
imwrite(imageHandler.BinaryImage, imagePathBinary);
imwrite(imageHandler.SkeletonImage, imagePathSkeleton);
csvHandler.NeuriteLengthMatrix = 0;
neuriteLengthMatrix = 0;
%Add to SkelDeletedMat:
%Neurons within NeuronPositionsEdgeFill which are not in NeuronPositionsEdgeFillNeurite
%Iterate over NeuronPositonsEdgeFill
%edgeFillPositions = neuronHandler.NeuronPositionsEdgeFill(selectedWell);
%[edgeFillPosesY edgeFillPosesX] = find(edgeFillPositions);
%for(counter = 1:numel(edgeFillPosesX))
%    if(NeuronPositionsEdgeFillNeurite(edgeFillPosesY(counter),edgeFillPosesX(counter)) ~= 1 )
%        SkelDeletedMat(edgeFillPosesY(counter),edgeFillPosesX(counter)) = 1;
%    end
%end
%Check if Neuronhandler already has Dictionary of NeuronPositionsEdgeComposite

if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = SkelDeletedMat;
neuronHandler.NeuronPositionsSkeletonization(selectedWell) = SkelNeurons;
handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles);
if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
if(numel(neuronHandler.AreaDictionary) ==0)
    neuronHandler.AreaDictionary = containers.Map();
end
neuronHandler.AreaDictionary(selectedWell)=neuriteAreaCount;
neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = SkelDeletedMat;
neuronHandler.NeuronPositionsSkeletonization(selectedWell) = SkelNeurons;
neuronHandler.NeuronStatMatrix(selectedWell) = currentNeuronStat;
handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function SkeletonizationNeurons_Callback(hObject, eventdata, handles)
% hObject    handle to SkeletonizationNeurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(handles.figure1);

selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
SkeletonizationNeuronsAlgorithm(selectedWell, handles);


% --------------------------------------------------------------------
function SkeletonizeNeurites(BW, handles)
%Next steps:
%1. Cut out all Single Neurite areas for themselfs
imageHandler = handles.ImageHandler;
[whitePointRows, whitePointCols] = find(BW);
alreadyVisited = uint32(ones(size(BW,1), size(BW,2)));
neuriteAreasDict = containers.Map();
i=0;
%FloodFill for Points, which were not visited before
%As long as there are white point indices, which were not already visited
while(numel(whitePointRows) > 0)
    %Get first whitepointindex <> 0
    %FloodFill    
    stack = java.util.Stack();
    stack.push([whitePointRows(1) whitePointCols(1)]);
    currentNeuriteArea = zeros(size(BW,1), size(BW,2));
    while(stack.isEmpty() == 0)
      vector = stack.pop();
      yPos=vector(1);
      xPos=vector(2);       
      %[visitedRows visitedCols] = ind2sub(7168,locations);
      
      if(BW(yPos,xPos) == 1 && alreadyVisited(yPos, xPos) == 1 &&  yPos + 1 < size(BW,1) && xPos + 1 < size(BW,2) && xPos - 1 > 0 && yPos - 1 > 0)
            currentNeuriteArea(yPos,xPos)= 1;
            stack.push([yPos+1 xPos]);
            stack.push([yPos xPos-1]);
            stack.push([yPos-1 xPos]);
            stack.push([yPos xPos+1]);
            %only 4 Neighbourhood!
            %stack.push([yPos+1 xPos+1]);
            %stack.push([yPos-1 xPos-1]);
            %stack.push([yPos+1 xPos-1]);
            %stack.push([yPos-1 xPos+1]);
      end
      alreadyVisited(yPos, xPos) = 0; 
    end    
    %Make skeleton processing with current neurite Area
    %Cut Neurite Area on white points
    %Create Neurite object:
    %Save min Position and Max Position on big picture and cut out
    %CurrentNeuriteArea
    newNeurite = Neurite();
    [rowIndices colIndices] = find(currentNeuriteArea);
    newNeurite.cutYPosStart = min(rowIndices) - 2;
    newNeurite.cutXPosStart = min(colIndices) - 2;
     
    newNeurite.cutYPosEnd = max(rowIndices) + 2;
    newNeurite.cutXPosEnd = max(colIndices) + 2;
    currentNeuriteArea=currentNeuriteArea(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
    currentNeuriteArea = 1-currentNeuriteArea;
    newNeurite.image = currentNeuriteArea;
    newNeurite.CalculateBaiSkeleton(4);
    [yPosList, xPosList] = GetNeuronPositionsOnNeurite(newNeurite,handles);
    %newNeurite.PlotSkeleton();    
    %close(figure(1));
    neuriteAreasDict(num2str(i)) = newNeurite;
    i=i+1;
    while(numel(whitePointRows) > 0 && alreadyVisited(whitePointRows(1), whitePointCols(1)) == 0)
        whitePointRows(1) = [];
        whitePointCols(1) = [];
    end   
end

%2. Skeletonize them with Bais algorithm
    %3.1 For all Skeletons: Get Endpoints (Just 1 8 Neighbour)
    %3.2 Get Branching Points: Skeleton Points with 3 or more 8 neighbours.
    %If 4 8-neighbours: Take that as two different Neurites. Distuinguish
    %Neurites with direction
    %If 3 8-Neighbours: Check for nearest branching point and distuinguish
    %also with different directions.
    %If too many branching points: Area is too dense.

%Check if Neuronhandler already has Dictionary of NeuronPositionsEdgeComposite
if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
neuronHandler.NeuronPositionsSkeletonization(selectedWell) = neuriteAreasDict;
handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles);


function GetNeuronPositionsOnNeurite(neurite,handles)
%Distuinguish between different Neurites
%Check if a branching point is available. Two 3-neighbours or one 4-neighbour
%(Use Floyds algorithm.)
handles = guidata(handles.figure1);
%neuronHandler = handles.NeuronCoordinates;
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
[sizeY sizeX]= size(imageHandler.NeuriteImage);

SkelDeletedMat = neuronHandler.NeuronPositionsSkelDeleted(selectedWell)
SkelDeletedMat = sparse(sizeY,sizeX);
SkelNeurons = sparse(sizeY,sizeX);
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
NucleusM = csvHandler.CellPosMatrix(selectedWell);
neurite.FindBranchingCrossingPoints();
neurite.KillShortBranches();
if(numel(neurite.xEndpoints)>=4)
    neurite.SplitUp();
end

%As a result there will be 2 up to 3 remaining endpoints
%If there are 2 endpoints: There has to be one Neuron: Check if available
%through EdgeFill-Algorithm. If yes, return.
%If more than one, check if one has to be killed
%If no one, check iteratively, if Neuron can be added through lowering the
%EdgeFillThreshold.
%If this is not possible, try to find neuron by extension of vector from
%EndPoints.

%Get for every Neurite the next Nucleus to the Endpoints
%For every Neurite: Add next Nucleus to Endpoint to Resultlist

%Next step: Check all distances from endpoint to Endpoint
%If there are only two endpoints and distance is to short -> Kill whole
%Neurite and if there is a Neuron by EdgeFill -> Kill it.
neurite.CalculateNeuriteLength();
if(neurite.neuriteLength < 25)
    %ToDo: Delete potential Neuclei on Neurite
    %Save position of EdgeFilled Nuclei, which should be deleted
    edgeFillPositions = neuronHandler.NeuronPositionsEdgeFill(selectedWell);
    edgeFillThreshold = optionHandler.EdgeFillNucleusAreaWithinNeurite;
    edgeFillNeuronsFirstThresh = edgeFillPositions(neurite.cutYPosStart:neurite.cutYPosEnd,neurite.cutXPosStart:neurite.cutXPosEnd,edgeFillThreshold:1);
    %Add Neurons to deleted ones
    [yPoses xPoses] = find(edgeFillNeuronsFirstThresh);
    SkelDeletedMat(neurite.cutYPosStart+yPoses-1, neurite.cutXPosStart+xPoses-1) = 1;
    neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = SkelDeletedMat;
 
    
else
    %First, check if already marked Neuron by EdgeFill-Algorithm available.
    edgeFillPositions = neuronHandler.NeuronPositionsEdgeFill(selectedWell);
    edgeFillThreshold = optionHandler.EdgeFillNucleusAreaWithinNeurite;
    edgeFillNeuronsFirstThresh = edgeFillPositions(neurite.cutYPosStart:neurite.cutYPosEnd,neurite.cutXPosStart:neurite.cutXPosEnd,edgeFillThreshold:1);

    if(numel(find(edgeFillNeuronsFirstThresh)) >= 1)
        %done
    else
        %Get Neuron with highest thresh in this area.
        %Check if Neuron between 0.55 and 1
        
        %If yes: Take this one
        %Delete all other potential neurons
        %If no: Do path extension
        
        %Check for other potential Neurons in Edge Fill Neurons by 2nd
        %threshold near endpoints. Get first one with biggest thresh.
        edgeFillNeuronsSecondThresh = edgeFillPositions(neurite.cutYPosStart:neurite.cutYPosEnd,neurite.cutXPosStart:neurite.cutXPosEnd,:);
        %[maxA,ind] = max(A(:));
        [yInd xInd] = find(edgeFillNeuronsSecondThresh(:,:,1));
        maxIndex = 0;
        maxThresh = 0.54;
        for(i=1:numel(yInd))
            if(edgeFillNeuronsSecondThresh(yInd(i),xInd(i),2) > maxThresh)
                maxThresh = edgeFillNeuronsSecondThresh(yInd(i),xInd(i),2);
                maxIndex = i;
            end
        end
        if(maxIndex ~= 0)
            %Add Neuron to SkelNeurons
            SkelNeurons(neurite.cutYPosStart+yInd(maxIndex)-1, neurite.cutXPosStart+xInd(maxIndex)-1) = 1;
        else 
            minDist = 100;
            nucColSaved = 0;
            nucRowSaved = 0;
            for(i=1:numel(neurite.xEndpoints))
                startPoint = [neurite.yEndpoints(i) neurite.xEndpoints(i)];
                [nucCol nucRow euclidDist] = csvHandler.FindNucleusForNeuron(neurite.cutXPosStart+startPoint(2), neurite.cutYPosStart+startPoint(1), csvHandler.CellPosMatrix,50);
                if(euclidDist<minDist)
                  minDist = euclidDist;
                  nucColSaved = nucCol;
                  nucRowSaved = nucRow;
                end
            end
            if(nucRowSaved > 0)
                SkelNeurons(nucRowSaved,nucColSaved)=1;
            end
        end
    end
end
if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
    neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = SkelDeletedMat;
neuronHandler.NeuronPositionsSkeletonization = SkelNeurons;




% --------------------------------------------------------------------
function NeuritesCellomicsExport_Callback(hObject, eventdata, handles)
% hObject    handle to NeuritesCellomicsExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
NucleusM = csvHandler.CellPosMatrix(selectedWell);
densityWidth = optionHandler.MigrationDistanceDensityImageXSize;
densityHeight = optionHandler.MigrationDistanceDensityImageYSize;
[sizeY sizeX]= size(imageHandler.NeuriteImage);
ringNumber = optionHandler.DensityDistributionRingNumber;
wellList = get(handles.lbWell, 'string');
foldername = uigetdir;
for i=1:numel(wellList)
    selectedWell = wellList(i);
    selectedWell = selectedWell{1};
        %Add 0 to selected Well and image Number  if only one digit
    selectedWellLong=selectedWell;
    if(length(selectedWell) == 2)
              selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    end    
    %Get density distribution to exclude too dense areas
    filepath = [imageHandler.Foldername '/ConvertedCellomics/' selectedWellLong 'Binary.tif'];
    if(exist(filepath,'file'))
        BW = imread(filepath);
    else
        BW = zeros(sizeY,sizeX);
    end

    %Create overlay with original pic
    imagePathNeuriteBig = [imageHandler.Foldername '/ConvertedCellomics/' selectedWellLong 'NeuriteBig.tif'];

    %imageHandler.ResizedNeuriteImage = imread(imagePathNeuriteSmall); 
    %imageHandler.NeuriteImage = imread(imagePathNeuriteBig);    
    originalNeurites = imread(imagePathNeuriteBig);
    BW = originalNeurites + uint8(BW.*255); 

    %Save that pic as binary image
    %imageHandler.BinaryImage = flipdim(BW,1);
    %imageHandler.ResizedBinaryImage = imresize(imageHandler.BinaryImage,0.1);
    %handles.ImageHandler = imageHandler;
    %guidata(handles.figure1, handles);
    %figure(7);
    %imshow(imfuse(BW.*255,imadjust(originalNeurites)));
    %Next: Split picture up to 512x512 images
    %User input: Directory where pictures should be exported
    
    x=0;
    y=0;
    wellType = GetWellType(handles);
    [sizeY sizeX] = size(originalNeurites);

    %Get filename of original files
    origFolder = imageHandler.Foldername;    
    allFiles = dir(origFolder); 

   

    for(i=1:numel(allFiles))
            %Try to get orig file of Nucleus:
            nucleusRegexp = strcat('([A-Za-z].*)_(\d+)_(',selectedWellLong,')f(\d+)d0');
            tokensNucleus = regexpi(allFiles(i).name, nucleusRegexp, 'tokens');        
            if(length(tokensNucleus) > 0)
              tokensNucleus = tokensNucleus{1};
              arrayText= tokensNucleus{1};
              number = tokensNucleus{2};          
              break;
            end
    end
    while(y<(sizeY))
        currentPic=BW(y+1:y+512,x+1:x+512);
        imageNumber = csvHandler.GetPictureNumberFromPosition(x,y,wellType);
        %Write Neurite image file
        if(length(num2str(imageNumber)) == 1)
              imageNumber = ['0' num2str(imageNumber)];
        end


        filename = [arrayText '_1_' selectedWellLong 'f' num2str(imageNumber) 'd1.tif'];
        filenameNuc = [arrayText '_1_' selectedWellLong 'f' num2str(imageNumber) 'd0.tif'];
        currentPic = flipdim(currentPic,1);
        imwrite(currentPic,[foldername '/' filename]);

        %Copy file for Nuclei
        nucleusPicFile = [origFolder '/' arrayText '_' number '_' selectedWellLong 'f' num2str(imageNumber) 'd0.tif'];
        nucleusPic = imread(nucleusPicFile);


          nucleusPic = double(nucleusPic)./4095;
         % nucleusPic = imadjust(nucleusPic, [0;4095], [0;1]);
          nucleusPic = uint8(nucleusPic.*255);    


        imwrite(nucleusPic,[foldername '/' filenameNuc]);

        if(x>=(sizeX-512))
            x=0;
            y=y+512;
        else
            x=x+512;
        end    
    end
end


function wellType = GetWellType(handles)
        handles = guidata(handles.figure1);
        wellType = WellType.Undefined;
        if(get(handles.rb48,'Value'))
         wellType = WellType.Well48;
        elseif(get(handles.rb4816,'Value'))
         wellType = WellType.Well4816;         
        elseif(get(handles.rb96,'Value'))
         wellType = WellType.Well96;
         elseif(get(handles.rb96484,'Value'))
         wellType = WellType.Well96484;
        elseif(get(handles.rb9616,'Value'))
         wellType = WellType.Well9616;
        elseif(get(handles.rbWholeOt,'Value'))
         wellType = WellType.OTWhole;
        elseif(get(handles.rbWholeOt16,'Value'))
         wellType = WellType.OTWhole16;
        elseif(get(handles.rbOtSingleChamber,'Value'))
         wellType = WellType.OTSingle;
        elseif(get(handles.rbOtSingleChamber16,'Value'))
         wellType = WellType.OTSingle16;
        end


% --- Executes on button press in cbCellomicsNeurons2.
function cbCellomicsNeurons2_Callback(hObject, eventdata, handles)
% hObject    handle to cbCellomicsNeurons2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles)
% Hint: get(hObject,'Value') returns toggle state of cbCellomicsNeurons2


% --------------------------------------------------------------------
function CalcContrastStretchIndices_Callback(hObject, eventdata, handles)
% hObject    handle to CalcContrastStretchIndices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;

%Take mean values of all Pictures in experiment
isodataMean=zeros(0);
minimumMean=zeros(0);
stretchlimMean = zeros(0);
wellList = get(handles.lbWell, 'string');
waitbarHandle = waitbar(0,'Please wait. Medium Threshold is calculated.');
for i=1:numel(wellList)    
    waitbar((i/numel(wellList)),waitbarHandle);
    selectedWell = wellList{i};
    if(length(selectedWell) == 2)
      selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    elseif(length(selectedWell) == 4)
        selectedWellLong = [selectedWell(1) '0' selectedWell(2:length(selectedWell))];
    else
        selectedWellLong = selectedWell;
    end
    %Get images from image folder
    foldername = [imageHandler.Foldername '/ConvertedCellomics'];
    %imagePathNucleusBig = [foldername '/' selectedWell 'NucleusBig.tif'];
    imagePathNeuriteBig = [foldername '/' selectedWellLong 'NeuriteBig.tif'];
    neuriteImage = imread(imagePathNeuriteBig);   
    currStretchlimres = stretchlim(neuriteImage,[0.975 0.9999]);
    brightNeurites = imadjust(neuriteImage,currStretchlimres,[0 1]);
    currThresholdingLevelStrong = th_minimum(brightNeurites)./255;
    %1. Median Filter
    brightNeurites = medfilt2(brightNeurites);

    %2. Resharp image
    unsharpFilter = fspecial('unsharp');
    brightNeurites = imfilter(brightNeurites,unsharpFilter);
    %figure(2);
    %imshow(brightNeurites);

    %3. Thresholding
    %OTSU:
    %neuronHandler.ThresholdingLevel=graythresh(brightNeurites);
    %ISODATA:
    brightNeurites = imcomplement(brightNeurites);
    currThresholdingLevel = isodata(brightNeurites);
    %INTERMODES:
    %neuronHandler.ThresholdingLevel = th_intermodes(brightNeurites)./255;
 %   if(currThresholdingLevelStrong < 0.5)
 %       currThresholdingLevelStrong = 1-currThresholdingLevelStrong;
 %   end
 %   if(currThresholdingLevelStrong >= currThresholdingLevel)
        minimumMean = [minimumMean currThresholdingLevelStrong];
 %   end
    isodataMean = [isodataMean currThresholdingLevel];
    minimumMean = [minimumMean currThresholdingLevelStrong];
    currStretchlimres = currStretchlimres';
    stretchlimMean = [stretchlimMean;currStretchlimres];
end
isodataMean = mean(isodataMean);
minimumMean = mean(minimumMean);
stretchlimMean = mean(stretchlimMean);

neuronHandler.StretchlimResult = stretchlimMean;
neuronHandler.IsodataResult = isodataMean;
neuronHandler.MinimumResult = minimumMean;
handles.NeuronCoordinates = neuronHandler;

%Calculate also second Threshold for Neurite verification
%Take Minimum threshold method from ImageJ
close(waitbarHandle);
guidata(handles.figure1, handles);


% --- Executes on button press in cbBinaryPic.
function cbBinaryPic_Callback(hObject, eventdata, handles)
% hObject    handle to cbBinaryPic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbBinaryPic
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles)


function level=isodata(I)
%   Source:http://www.mathworks.com/matlabcentral/fileexchange/3195-automatic-thresholding/content/isodata.m
%   ISODATA Compute global image threshold using iterative isodata method.
%   LEVEL = ISODATA(I) computes a global threshold (LEVEL) that can be
%   used to convert an intensity image to a binary image with IM2BW. LEVEL
%   is a normalized intensity value that lies in the range [0, 1].
%   This iterative technique for choosing a threshold was developed by Ridler and Calvard .
%   The histogram is initially segmented into two parts using a starting threshold value such as 0 = 2B-1, 
%   half the maximum dynamic range. 
%   The sample mean (mf,0) of the gray values associated with the foreground pixels and the sample mean (mb,0) 
%   of the gray values associated with the background pixels are computed. A new threshold value 1 is now computed 
%   as the average of these two sample means. The process is repeated, based upon the new threshold, 
%   until the threshold value does not change any more. 
  
%
%   Class Support
%   -------------
%   The input image I can be of class uint8, uint16, or double and it
%   must be nonsparse.  LEVEL is a double scalar.
%
%   Example
%   -------
%       I = imread('blood1.tif');
%       level = graythresh(I);
%       BW = im2bw(I,level);
%       imshow(BW)
%
%   See also IM2BW.
%
% Reference :T.W. Ridler, S. Calvard, Picture thresholding using an iterative selection method, 
%            IEEE Trans. System, Man and Cybernetics, SMC-8 (1978) 630-632.

% Convert all N-D arrays into a single column.  Convert to uint8 for
% fastest histogram computation.
I = im2uint8(I(:));

% STEP 1: Compute mean intensity of image from histogram, set T=mean(I)
[counts,N]=imhist(I);
i=1;
mu=cumsum(counts);
T(i)=(sum(N.*counts))/mu(end);
T(i)=round(T(i));

% STEP 2: compute Mean above T (MAT) and Mean below T (MBT) using T from
% step 1
mu2=cumsum(counts(1:T(i)));
MBT=sum(N(1:T(i)).*counts(1:T(i)))/mu2(end);

mu3=cumsum(counts(T(i):end));
MAT=sum(N(T(i):end).*counts(T(i):end))/mu3(end);
i=i+1;
% new T = (MAT+MBT)/2
T(i)=round((MAT+MBT)/2);

% STEP 3 to n: repeat step 2 if T(i)~=T(i-1)
while abs(T(i)-T(i-1))>=1
    mu2=cumsum(counts(1:T(i)));
    MBT=sum(N(1:T(i)).*counts(1:T(i)))/mu2(end);
    
    mu3=cumsum(counts(T(i):end));
    MAT=sum(N(T(i):end).*counts(T(i):end))/mu3(end);
    
    i=i+1;
    T(i)=round((MAT+MBT)/2); 
    Threshold=T(i);
end

 % Normalize the threshold to the range [i, 1].
level = (Threshold - 1) / (N(end) - 1);


% --- Executes on button press in cbEdgeFillNeurons.
function cbEdgeFillNeurons_Callback(hObject, eventdata, handles)
% hObject    handle to cbEdgeFillNeurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
% Hint: get(hObject,'Value') returns toggle state of cbEdgeFillNeurons


% --------------------------------------------------------------------
function ExportAlgorithmStat_Callback(hObject, eventdata, handles)
% hObject    handle to ExportAlgorithmStat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
%savefoldername = uigetdir;

%Get start and stop Well

wellList = get(handles.lbWell, 'string');
exportText = sprintf('Well;xPos;yPos;Skeleton (1=Skeleton found from Edge Fill, 2=Skeleton marked from EdgeFill Second, 3=Skeleton found from endpoint of Skeleton);Edge Fill;Edge Fill Lower;Manual;Deleted by Skel (1=Deleted because too less Neurite length or Distance between Skeleton endpoints, 2=Deleted because another Neuron on same Skeleton and both Neurons are too near together);Edge Fill Overlap;Neurite length;Max dist between Skeleton Endpoints;\r\n');
for i=1:numel(wellList)    
 selectedWell = wellList{i}; 
 if(isKey(neuronHandler.NeuronStatMatrix, selectedWell))
     NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
      densityWidth = 50;
      densityHeight = 50;
    [sizeY sizeX]= size(imageHandler.NeuriteImage);
    ringNumber = optionHandler.DensityDistributionRingNumber;
    %Get density distribution to exclude too dense areas
    path = [imageHandler.Foldername '/ConvertedCellomics/MarkerPointCoordinates-' selectedWell '.mat'];
if(exist(path,'file'))    
    ringNumber = optionHandler.DensityDistributionRingNumber;
    load(path);
    filterDistance = str2double(filterDistance);
    nonFilterDistance = str2double(nonFilterDistance);
    %hInner = impoly(handles.axes2,double(markerPointCoordinates('0')./10));
    %previousMask = createMask(hInner);
    %for i=1:ringNumber  
    %  hCurrent = impoly(handles.axes2,double(markerPointCoordinates(num2str(i*10))./10));
    %end
else
      
     [filterDistance nonFilterDistance SphereArea markerPointCoordinates] = calculateMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, handles); 
     %Save hInner to file
     foldername = imageHandler.Foldername;
     subfoldername = [foldername '/ConvertedCellomics'];
     %if(markerPointCoordinates~=0)
        save('-v7.3',strcat(subfoldername,'\MarkerPointCoordinates-',selectedWell),'markerPointCoordinates', 'filterDistance', 'nonFilterDistance');
     %end
end         

    hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./10));
    innerMask = logical(createMask(hInner));

        currentRing =(logical(ones(int32(sizeY/10),int32(sizeX/10)) - innerMask));

    currentRing = imresize(currentRing, [sizeY, sizeX]);
    currentRing=sparse(currentRing);
    delete(hInner);
    innerMask=0;
     
     
     currentStatDict = neuronHandler.NeuronStatMatrix(selectedWell);
     %Get all relevant positions from matrices
     a=csvHandler.CellPosMatrix(selectedWell);
     if(size(a,2)>7168)
         a=a(:,1:7168);
         csvHandler.CellPosMatrix(selectedWell)=a;
     end
     
     [yPoses xPoses] = find(csvHandler.CellPosMatrix(selectedWell) .* currentRing);     
     %For every positive entry in Filter List: Get and Export fused pic
     for j=1:numel(yPoses)
         key=mapPositionToKeyString(yPoses(j),xPoses(j));
         if(isKey(currentStatDict,key))
           list=currentStatDict(key);
           exportText = [exportText selectedWell ';' num2str(xPoses(j)) ';' num2str(yPoses(j)) ';' num2str(list(1)) ';' num2str(list(2)) ';' num2str(list(3)) ';' num2str(list(4)) ';' num2str(list(5)) ';' num2str(list(6)) ';' num2str(list(7)) ';' num2str(list(8)) ';' sprintf('\r\n')];
         elseif(NeuronManualM(yPoses(j),xPoses(j)) == 1)
           exportText = [exportText selectedWell ';' num2str(xPoses(j)) ';' num2str(yPoses(j)) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(1) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';' sprintf('\r\n')];
         end
     end
 end
end
[FileName,PathName] = uiputfile('migdist.csv');
fileID = fopen(strcat(PathName,'\',FileName),'w');
fprintf(fileID,'%s',exportText);
fclose(fileID);


% --- Executes on button press in pbAutoThreshold.
function pbAutoThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to pbAutoThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
%optionHandler = handles.OptionHandler;
%minBrightness = optionHandler.HistogramMinNeurite;
%maxBrightness = optionHandler.HistogramMaxNeurite;
currentPicNeu = imageHandler.ResizedNeuriteImage;

stretchlimLow = stretchlim(imageHandler.ResizedNeuriteImage,[0.975 0.9999]);
currentPicNeu = double(currentPicNeu)./65536;
%currentPicNeu = imadjust(currentPicNeu, stretchlimLow, [0;1]);
currentPicNeu = uint8(currentPicNeu.*255);   
%minBrightness = optionHandler.HistogramMinNucleus;
%maxBrightness = optionHandler.HistogramMaxNucleus;
currentPicNuc = imageHandler.ResizedNucleusImage;

stretchlimLow = stretchlim(imageHandler.ResizedNucleusImage,[0.975 0.9999]);

currentPicNuc = imadjust(currentPicNuc, stretchlimLow, [0;1]);
currentPicNuc = double(currentPicNuc)./65536;
currentPicNuc = uint8(currentPicNuc.*255); 
imageHandler.ShowNucleusImage = currentPicNuc;
imageHandler.ShowNeuriteImage = currentPicNeu;
mat = imfuse(currentPicNuc,currentPicNeu);
hImage = imshow(mat);
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles)


% --------------------------------------------------------------------
function Convert_8Bit_Auto_Callback(hObject, eventdata, handles)
% hObject    handle to Convert_8Bit_Auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
foldername = [imageHandler.Foldername '/ConvertedCellomics'];

selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');

stretchLimLowNeurite = stretchlim(imageHandler.NeuriteImage,[0.975 0.9999]);
stretchLimLowNucleus = stretchlim(imageHandler.NucleusImage,[0.975 0.9999]);

allFiles = dir(foldername); 
waitbarHandle = waitbar(0,'Please wait. All images are converted to 8 Bit.');
for i=3:numel(allFiles)
    waitbar((i/numel(allFiles)),waitbarHandle);
    currentFile=allFiles(i).name;

    imageString = strcat(foldername, '\', currentFile);
    if(length(findstr(imageString, '.mat')) == 0)
      currentPic = imread(imageString);
      %deleteIndices = uint16(find(abs(currentPic) < minBrightness));
      %currentPic(deleteIndices) = 0;
    
      %currentPic = uint8(currentPic);
      if(length(findstr(currentFile, 'Neurite')) == 0)      
          stretchLimLowNucleus = stretchlim(currentPic,[0.975 0.9999]);
          currentPic = imadjust(currentPic, stretchLimLowNucleus, [0;1]);
      else
          stretchLimLowNeurite = stretchlim(currentPic,[0.975 0.9999]);
          currentPic = imadjust(currentPic, stretchLimLowNeurite, [0;1]);
      end
      %Extra
      currentPic = double(currentPic)./65536;
      currentPic = uint8(currentPic.*255);   
      
      imwrite(currentPic,imageString);
    end
end
close(waitbarHandle);


% --- Executes on button press in cbDeleted.
function cbDeleted_Callback(hObject, eventdata, handles)
% hObject    handle to cbDeleted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
% Hint: get(hObject,'Value') returns toggle state of cbDeleted


% --- Executes on button press in cbEdgeFill2.
function cbEdgeFill2_Callback(hObject, eventdata, handles)
% hObject    handle to cbEdgeFill2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbEdgeFill2
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);


% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes on button press in radiobutton26.
function radiobutton26_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton26


% --- Executes on button press in radiobutton27.
function radiobutton27_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton27


% --- Executes on button press in radiobutton28.
function radiobutton28_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton28


% --------------------------------------------------------------------
function CalculateBinaryImage_Callback(hObject, eventdata, handles)
% hObject    handle to CalculateBinaryImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
NucleusM = csvHandler.CellPosMatrix(selectedWell);
[sizeY sizeX]= size(imageHandler.NeuriteImage);

%Threshold Neurite Picture the known way
neuritePic = CutOutCircles(imageHandler.NeuriteImage,selectedWell,1,0,handles);

neuritePic = ThresholdPic(neuritePic,0,handles);
%Edge filter for Neurite Picture
%neuritePic = edge(neuritePic, 'nothinning');

SE = strel('disk', 3);
neuritePic = imdilate(neuritePic,SE);

%nucleusPic = CutOutCircles(0,1,0,handles);
%nucleusPic = zeros(7168,7168);
%nucleusPic(nucleusPic<15) = 0;%double(double(imageHandler.NucleusImage) ./ 255);
%nucleusPic(nucleusPic>=15) = 255;
%nucleusPic = uint8(nucleusPic);
imageHandler.BinaryImage = neuritePic;
imageHandler.ResizedBinaryImage = imresize(neuritePic,0.1);
handles.ImageHandler=imageHandler;
guidata(handles.figure1, handles)


% ----------------------------------------------------------------- ---
function ManualEval_Callback(hObject, eventdata, handles)
% hObject    handle to ManualEval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
manual1 = nnz(csvHandler.ManualPositions1(selectedWell));
manual2 = nnz(csvHandler.ManualPositions2(selectedWell));
manual3 = nnz(csvHandler.ManualPositions3(selectedWell));
manual4 = nnz(csvHandler.ManualPositions4(selectedWell));
disp('Manual 1; Manual 2; Manual 3; Manual4');
disp([num2str(manual1) ';' num2str(manual2) ';' num2str(manual3) ';' num2str(manual4)]);


% --- Executes on button press in cbManualEval1.
function cbManualEval1_Callback(hObject, eventdata, handles)
% hObject    handle to cbManualEval1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
% Hint: get(hObject,'Value') returns toggle state of cbManualEval1


% --------------------------------------------------------------------
function ParameterExperimentBatch_Callback(hObject, eventdata, handles)
% hObject    handle to ParameterExperimentBatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
%foldername = uigetdir;
%Call EdgeFill Neurons and Skeleton Neurons and play with parameters in
%Option Handler



foldername='C:\Users\Thoma_000\Desktop\ExperimentTest\6. Exp';

ExecuteOptionExperiment(foldername, 41,handles);

guidata(handles.figure1, handles);





function ExecuteOptionExperiment(foldername,startWell,handles)
handles = guidata(handles.figure1);
wellList = get(handles.lbWell, 'string');
%Select STARTWELL
for i=startWell:numel(wellList)     
 selectedWell = wellList{i};
 disp(['starting' startWell]);
 %if(strcmp(selectedWell,'E3'))
  %   continue;
  %end
 selectedWellLong=selectedWell;
  if(numel(selectedWell) == 2 || (~strcmp(selectedWell(2),'0') && strcmp(selectedWell(3),'_')))
      selectedWellLong = [selectedWell(1) '0' selectedWell(2:end)];
  end
  %Execute for each well:
  %ToDo: EdgeCompositeNeuronsAlgorithm
  
 % EdgeCompositAlgorithm(selectedWell, handles);
  EdgeFillNeuronsAlgorithm(selectedWell,handles);
  SkeletonizationNeuronsAlgorithm(selectedWell, handles);
  disp(['finished' selectedWell]);
  %Make Quality Check for EdgeFill Skeleton and Save value in CSV for each 
  %single well
  %Get number of Manual Neurons, EdgeFill Neurons and Skeleton Neurons
  %Get also numbers of same positions
  handles = guidata(handles.figure1);
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
WellNeuronDict = neuronHandler.NeuronPositionsEdgeComposite;


handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles);  
end


% --------------------------------------------------------------------
function CalcCellQual_Callback(hObject, eventdata, handles)
% hObject    handle to CalcCellQual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
wellList = get(handles.lbWell, 'string');
%Select STARTWELL
excelStr=['Neurons Manual;Neurons Cellomics; Positions Manual & Cellomics; Additional Manual; Additional Cellomics;' sprintf('\r\n')];
for i=1:numel(wellList)    
 selectedWell = wellList{i};

 selectedWellLong=selectedWell;
  if(numel(selectedWell) == 2 || (~strcmp(selectedWell(2),'0') && strcmp(selectedWell(3),'_')))
      selectedWellLong = [selectedWell(1) '0' selectedWell(2:end)];
  end

  handles = guidata(handles.figure1);
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
WellNeuronDict = neuronHandler.NeuronPositionsEdgeComposite;
%selectedWellNumber = get(handles.lbWell,'Value');
%wellList = get(handles.lbWell, 'string');
%selectedWell = wellList{selectedWellNumber};

densityWidth = 50;
densityHeight = 50;
[sizeY sizeX]= size(imageHandler.NeuriteImage);
%Get density distribution to exclude too dense areas
[filterDistance nonFilterDistance SphereArea NucleusArea64 markerPointCoordinates] = calculateMigrationDistance(selectedWell, densityWidth, densityHeight, handles);    
%end
%MarkerPointCoordinates are available. Get second ring:
hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./10));
innerMask = logical(createMask(hInner));
currentRing =(logical(ones(int32(sizeY/10),int32(sizeX/10)) - innerMask));
currentRing = imresize(currentRing, [sizeY, sizeX]);
currentRing=sparse(currentRing);
delete(hInner);
innerMask=0;


filterDistance=0;
nonFilterDistance=0;
SphereArea=0;
NucleusArea64=0;
markerPointCoordinates=0;
ManualNeurons = logical(neuronHandler.ManualNeuronPositionsSparse(selectedWell).* currentRing);
CellomicsNeurons = logical(neuronHandler.CellPosMatrix(selectedWell).* currentRing);
currentRing=0;
bothManualCellomics = nnz(ManualNeurons .* CellomicsNeurons);
pause(1);
ManualNeuronsInvert = ~ManualNeurons;

AdditionalCellomicsManual = nnz(CellomicsNeurons .* ManualNeuronsInvert);
ManualNeuronsInvert=0;

CellomicsNeuronsInvert = ~CellomicsNeurons;
CellomicsNeurons = nnz(CellomicsNeurons);
AdditionalManualCellomics = nnz(ManualNeurons .* CellomicsNeuronsInvert);
CellomicsNeuronsInvert=0;



excelStr=[excelStr selectedWell ';' num2str(nnz(ManualNeurons)) ';' num2str(CellomicsNeurons) ';'  num2str(bothManualCellomics) ';' num2str(AdditionalManualCellomics) ';' num2str(AdditionalCellomicsManual) ';;' sprintf('\r\n')];
%Write excelStr to excelFile
disp(excelStr);

%Delete saved matrices after saving to excel for memory reasons.

handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles);  
end
disp('finished');


% --- Executes on button press in cbSkeletonPic.
function cbSkeletonPic_Callback(hObject, eventdata, handles)
% hObject    handle to cbSkeletonPic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
% Hint: get(hObject,'Value') returns toggle state of cbSkeletonPic


% --- Executes on button press in cbManualEval2.
function cbManualEval2_Callback(hObject, eventdata, handles)
% hObject    handle to cbManualEval2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
% Hint: get(hObject,'Value') returns toggle state of cbManualEval2


% --- Executes on button press in cbManualEval3.
function cbManualEval3_Callback(hObject, eventdata, handles)
% hObject    handle to cbManualEval3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
% Hint: get(hObject,'Value') returns toggle state of cbManualEval3


% --- Executes on button press in cbManualEval4.
function cbManualEval4_Callback(hObject, eventdata, handles)
% hObject    handle to cbManualEval4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
% Hint: get(hObject,'Value') returns toggle state of cbManualEval4


% --- Executes on selection change in popupPlus1.
function popupPlus1_Callback(hObject, eventdata, handles)
% hObject    handle to popupPlus1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupPlus1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupPlus1


% --- Executes during object creation, after setting all properties.
function popupPlus1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupPlus1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupMinus1.
function popupMinus1_Callback(hObject, eventdata, handles)
% hObject    handle to popupMinus1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupMinus1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupMinus1


% --- Executes during object creation, after setting all properties.
function popupMinus1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupMinus1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupPlus2.
function popupPlus2_Callback(hObject, eventdata, handles)
% hObject    handle to popupPlus2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupPlus2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupPlus2


% --- Executes during object creation, after setting all properties.
function popupPlus2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupPlus2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupMinus2.
function popupMinus2_Callback(hObject, eventdata, handles)
% hObject    handle to popupMinus2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupMinus2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupMinus2


% --- Executes during object creation, after setting all properties.
function popupMinus2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupMinus2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cbRemoveCore.
function cbRemoveCore_Callback(hObject, eventdata, handles)
% hObject    handle to cbRemoveCore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
% Hint: get(hObject,'Value') returns toggle state of cbRemoveCore


% --- Executes on selection change in popupMinus3.
function popupMinus3_Callback(hObject, eventdata, handles)
% hObject    handle to popupMinus3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupMinus3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupMinus3


% --- Executes during object creation, after setting all properties.
function popupMinus3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupMinus3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function ExportFullQualityCheck_Callback(hObject, eventdata, handles)
% hObject    handle to ExportFullQualityCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
%savefoldername = uigetdir;
prompt = {'First Well', 'Last Well'};
dlg_title = 'Quality Check Export';
answer = inputdlg(prompt,dlg_title,1);
wellList = get(handles.lbWell, 'string');
exportText = sprintf('Well;Ring Number;Nuclei;Nuclei without cut;Neurons Manual with cut;Neurons Cellomics with cut;Neurons Manual without cut; Neurons Cellomics without cut;Neurons Cellomics 2;Neurons Edge Fill;Neurons Skeleton;;Neurons Cellomics & Manual with cut;Additional Manual;Additional Cellomics;;Neurons Cellomics & Manual without cut;Additional Manual;Additional Cellomics;;Neurons Cellomics 2 & Manual;Additional Manual;Additional Cellomics 2;;Neurons Edge Fill & Manual;Additional Manual;Additional Edge Fill;;Neurons Skeleton & Manual;Additional Manual;Additional Skeleton;\r\n');
startindex = find(ismember(wellList,answer(1)));
stopindex=find(ismember(wellList,answer(2)));
for i=startindex:stopindex
 selectedWell = wellList{i}; 
 densityWidth = 50;
densityHeight = 50;
[sizeY sizeX]= size(imageHandler.NeuriteImage);
ringNumber = optionHandler.DensityDistributionRingNumber;
%Get density distribution to exclude too dense areas
path = [imageHandler.Foldername '/ConvertedCellomics/MarkerPointCoordinates-' selectedWell '.mat'];
%Get density distribution to exclude too dense areas
%if(saveCircle)    
  if(exist(path,'file'))    
    ringNumber = optionHandler.DensityDistributionRingNumber;
    load(path);
    %hInner = impoly(handles.axes2,double(markerPointCoordinates('0')./10));
    %previousMask = createMask(hInner);
    %for i=1:ringNumber  
    %  hCurrent = impoly(handles.axes2,double(markerPointCoordinates(num2str(i*10))./10));
    %end
  else
     [filterDistance nonFilterDistance SphereArea markerPointCoordinates] = calculateMigrationDistance(selectedWell, densityWidth, densityHeight, handles); 
     %Save hInner to file
     foldername = imageHandler.Foldername;
     subfoldername = [foldername '/ConvertedCellomics'];
     save('-v7.3',strcat(subfoldername,'\MarkerPointCoordinates-',selectedWell),'markerPointCoordinates', 'filterDistance', 'nonFilterDistance');
  end
  

  if(markerPointCoordinates ~=0 && numel(markerPointCoordinates('10') > 0))
    hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./10));
    innerMask = logical(createMask(hInner));
%hOuter = impoly(handles.axes2,double(markerPointCoordinates(num2str(ringNumber * 10))./10));
%outerMask = logical(createMask(hOuter));
%if(CutInnerCircle == 1 && CutOuterCircle == 1)
%    currentRing = logical(outerMask-innerMask); 
%elseif(CutInnerCircle==1 && CutOuterCircle == 0)
    currentRing =(logical(ones(int32(sizeY/10),int32(sizeX/10)) - innerMask));
%else
%    currentRing = outerMask;
%end
    currentRing = imresize(currentRing, [sizeY, sizeX]);
    currentRing=sparse(currentRing);
    delete(hInner);
%delete(hOuter);
    innerMask=0;
  else
      currentRing=ones(sizeY,sizeX);
  end
if(isfield(handles,'NeuronCoordinates2'))
    neuronHandler2 = handles.NeuronCoordinates2;
    EdgeCompNeurons = logical(neuronHandler2.CellPosMatrix(selectedWell));
    [EdgeCompSizeY EdgeCompSizeX] = size(EdgeCompNeurons);
end
area = optionHandler.EdgeFillNucleusAreaWithinNeurite;

[sizeY sizeX] = size(imageHandler.NeuriteImage);
NucleusM = logical(csvHandler.CellPosMatrix(selectedWell));
[NucleusMSizeY NucleusMSizeX]= size(NucleusM);
if(NucleusMSizeX > sizeX)
     NucleusM(:,sizeX+1:NucleusMSizeX) = [];
     csvHandler.CellPosMatrix(selectedWell) = NucleusM;
end

if(NucleusMSizeY > sizeY)
     NucleusM(sizeY+1:NucleusMSizeY,:) = [];
     csvHandler.CellPosMatrix(selectedWell) = NucleusM;
end
NucleusMWithCut = logical(NucleusM.*currentRing); 
ManualNeurons = logical(neuronHandler.ManualNeuronPositionsSparse(selectedWell) .* currentRing);

filterDistance=0;
nonFilterDistance=0;
SphereArea=0;
NucleusArea64=0;
markerPointCoordinates=0;
ManualNeuronsInvert = ~ManualNeurons;
if(isKey(neuronHandler.NeuronPositionsEdgeFill,selectedWell))

    EdgeFilledNeurons = logical(neuronHandler.NeuronPositionsEdgeFill(selectedWell).* currentRing);
    SkeletonNeurons = logical(neuronHandler.NeuronPositionsSkeletonization(selectedWell) .* currentRing);
    NeuronsDeleted = neuronHandler.NeuronPositionsSkelDeleted(selectedWell);

    SkeletonNeurons = SkeletonNeurons - (NeuronsDeleted);
    ind = SkeletonNeurons<0;
    SkeletonNeurons(ind) = 0;
    bothManualEdgeFilled = nnz(EdgeFilledNeurons .* ManualNeurons);
    bothManualSkeletonized = nnz(SkeletonNeurons .* ManualNeurons);
        
    AdditionalSkeletonManual = nnz(SkeletonNeurons .* ManualNeuronsInvert);

    AdditionalEdgeFilledManual = nnz(EdgeFilledNeurons .* ManualNeuronsInvert);
   % ManualNeuronsInvert=0;
    EdgeFilledNeuronsInvert = ~EdgeFilledNeurons;
    AdditionalManualEdgeFilled = nnz(ManualNeurons .* EdgeFilledNeuronsInvert);
    EdgeFilledNeuronsInvert=0;
    
        SkeletonNeuronsInvert = ~SkeletonNeurons;
      AdditionalManualSkeleton = nnz(ManualNeurons .* SkeletonNeuronsInvert);
    SkeletonNeuronsInvert=0;
else
    EdgeFilledNeurons=0;
    SkeletonNeurons=0;
    bothManualEdgeFilled=0;
    bothManualSkeletonized=0;
    AdditionalManualEdgeFilled=0;
    AdditionalEdgeFilledManual=0;
    AdditionalSkeletonManual=0;
    AdditionalManualEdgeFilled=0;
end
%EdgeFilledNeurons = EdgeFilledNeurons(:,:,2) > area;
CellomicsNeuronsWithoutCut = neuronHandler.CellPosMatrix(selectedWell);

CellomicsNeuronsWithCut = logical(CellomicsNeuronsWithoutCut.* currentRing);
currentRing=0;

ManualNeuronsWithoutCut = logical(neuronHandler.ManualNeuronPositionsSparse(selectedWell));

%disp('Neurons Manual;Neurons Cellomics with cut;Neurons Edge Composite;Positions Manual & Edge;Additional Manual; Additional Edge Composite;;Positions Edge & Cellomics;Additional Cellomics; Additional Edge Composite;; Positions Manual & Cellomics; Additional Manual; Additional Cellomics;; Positions Manual & Edge Filled; Additioal Manual; Additional Edge Filled')
if(isfield(handles,'NeuronCoordinates2'))
    bothCount = nnz(EdgeCompNeurons .* ManualNeurons);
else
    bothCount=0;
end
bothManualCellomics = nnz(ManualNeurons .* CellomicsNeuronsWithCut);
bothManualCellomicsWithoutCut = nnz(ManualNeuronsWithoutCut .* CellomicsNeuronsWithoutCut);


if(isfield(handles,'NeuronCoordinates2'))
    AdditionalEdge = nnz(EdgeCompNeurons .* ManualNeuronsInvert);   
else
    AdditionalEdge=0;
end
CellomicsNeuronsInvertWithCut = ~CellomicsNeuronsWithCut;

AdditionalManualCellomics = nnz(ManualNeurons .* CellomicsNeuronsInvertWithCut);
CellomicsNeuronsInvertWithCut=0;
CellomicsNeuronsInvertWithoutCut = ~CellomicsNeuronsWithoutCut;
ManualNeuronsInvertWithoutCut = ~ManualNeuronsWithoutCut;
AdditionalManualCellomicsWithoutCut = nnz(ManualNeuronsWithoutCut .* CellomicsNeuronsInvertWithoutCut);
CellomicsNeuronsInvertWithoutCut=0;
CellomicsNeuronsInvert=0;
AdditionalCellomicsManual = nnz(CellomicsNeuronsWithCut .* ManualNeuronsInvert);
AdditionalCellomicsManualWithoutCut = nnz(CellomicsNeuronsWithoutCut .* ManualNeuronsInvertWithoutCut);
ManualNeuronsInvertWithoutCut=0;

if(isfield(handles,'NeuronCoordinates2'))
EdgeCompInvert = ~EdgeCompNeurons;
AdditionalManual = nnz(EdgeCompInvert .* ManualNeurons);
EdgeCompInvert=0;
else
    AdditionalManual=0;
    EdgeCompNeurons=0;
end


 %Export for whole ring
 exportText = [exportText selectedWell ';' 'whole Well' ';' num2str(nnz(NucleusMWithCut)) ';' num2str(nnz(NucleusM)) ';' num2str(nnz(ManualNeurons)) ';' num2str(nnz(CellomicsNeuronsWithCut)) ';' num2str(nnz(ManualNeuronsWithoutCut)) ';' num2str(nnz(CellomicsNeuronsWithoutCut)) ';'  num2str(nnz(EdgeCompNeurons)) ';' num2str(nnz(EdgeFilledNeurons)) ';' num2str(nnz(SkeletonNeurons)) ';;' num2str(bothManualCellomics) ';' num2str(AdditionalManualCellomics) ';' num2str(AdditionalCellomicsManual) ';;'  num2str(bothManualCellomicsWithoutCut) ';' num2str(AdditionalManualCellomicsWithoutCut) ';' num2str(AdditionalCellomicsManualWithoutCut) ';;' num2str(bothCount) ';' num2str(AdditionalManual) ';' num2str(AdditionalEdge) ';;' num2str(bothManualEdgeFilled) ';' num2str(AdditionalManualEdgeFilled) ';' num2str(AdditionalEdgeFilledManual) ';;' num2str(bothManualSkeletonized) ';' num2str(AdditionalManualSkeleton) ';' num2str(AdditionalSkeletonManual) ';;' sprintf('\r\n')]
end




[FileName,PathName] = uiputfile('NeuronStatistics.csv');
fileID = fopen(strcat(PathName,'\',FileName),'w');
fprintf(fileID,'%s',exportText);
fclose(fileID);
