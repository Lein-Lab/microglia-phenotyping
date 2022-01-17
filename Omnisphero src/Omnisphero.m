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

    function varargout = Omnisphero(varargin)
% OMNISPHERO MATLAB code for Omnisphero.fig
%      OMNISPHERO, by itself, creates a new OMNISPHERO or raises the existing
%      singleton*.
%
%      H = OMNISPHERO returns the handle to a new OMNISPHERO or the handle to
%      the existing singleton*.
%
%      OMNISPHERO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OMNISPHERO.M with the given input arguments.
%
%      OMNISPHERO('Property','Value',...) creates a new OMNISPHERO or raises the
%      existing singleton*.  Starting from the  left, property value pairs are
%      applied to the GUI before Omnisphero_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Omnisphero_OpeningFcn via varargin.
%
%      *See GUI Opt on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Omnisphero

% Last Modified by GUIDE v2.5 09-Mar-2018 12:03:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Omnisphero_OpeningFcn, ...
                   'gui_OutputFcn',  @Omnisphero_OutputFcn, ...
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


% --- Executes just before Omnisphero is made visible.
function Omnisphero_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Omnisphero (see VARARGIN)

% Choose default command line output for Omnisphero
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes Omnisphero wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Omnisphero_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%Open GUI in Full Screen Mode
addpath(genpath('./lib'));
%maxfig(handles.figure1,1);


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
try
    [filename, pathname] = uigetfile('*.csv', 'Select CSV File with Cellomics Data');
    LoadCSV(filename,pathname,0,0,handles);
catch errorObj
% If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
% --------------------------------------------------------------------

function LoadCSV(filename,pathname,wellType,readMDB,handles)
    handles = guidata(handles.figure1);
if (~isempty(filename) && ~isempty(pathname))
    csvHandler = CSVCoordinates();
    filepath = strcat(pathname,filename);
    if(readMDB==0)
        csvHandler.ReadCSVFile(filepath, wellType);
    else
        csvHandler.ReadMDBFile(filepath, wellType,0);    
    end
    handles.CSVCoordinates = csvHandler;    
    guidata(handles.figure1, handles)
end

% --------------------------------------------------------------------
function LoadNeuronCSV2_Callback(hObject, eventdata, handles)
% hObject    handle to LoadNeuronCSV2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)
try
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
        wellType = GetWellType(handles);
        %%if(readMDB==1)
        %%    csvHandler.ReadMDBFile(filepath, wellType,1);
        %%else
            neuronHandler.ReadCSVFile(filepath, wellType, saveValuesArray, filterArray);       
        %%end
        neuronHandler.ManualNeuronPositions = containers.Map();
        neuronHandler.CellPosMatrix = WellNeuronDict;
        handles.NeuronCoordinates2 = neuronHandler;
        guidata(handles.figure1, handles);
        WellNeuronDict = CalcNucleusNeuronDist(handles, neuronHandler.CellPosMatrix);
        neuronHandler.CellPosMatrix = WellNeuronDict;
        handles.NeuronCoordinates2 = neuronHandler;
        guidata(handles.figure1, handles);
end
catch errorObj
% If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function LoadNeuronCSV(filename,pathname,wellType,readMDB ,handles)
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
    if(readMDB==1)
        neuronHandler.ReadMDBFile(filepath, wellType,1);
    else
        neuronHandler.ReadCSVFile(filepath, wellType, saveValuesArray, filterArray);
    end
    neuronHandler.ManualNeuronPositions = containers.Map();
        handles.NeuronCoordinates = neuronHandler;
    guidata(handles.figure1, handles)
    
   % WellNeuronDict = CalcNucleusNeuronDist(handles, neuronHandler.CellPosMatrix);
   % neuronHandler.CellPosMatrix = WellNeuronDict;
    handles.NeuronCoordinates = neuronHandler;
    guidata(handles.figure1, handles)
end


function LoadNeuronCSV_Callback(hObject, eventdata, handles)
% hObject    handle to LoadNeuronCSV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    [filename, pathname] = uigetfile('*.csv', 'Select CSV File with Cellomics Data');
    LoadNeuronCSV(filename,pathname,handles,0);
catch errorObj
% If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
% --------------------------------------------------------------------
function FixNeuronPositions_Callback(hObject, eventdata, handles)
% hObject    handle to FixNeuronPositions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    handles = guidata(handles.figure1);
    FixNeuronPositions(hObject, eventdata, handles);
    guidata(handles.figure1, handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
    

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
    waitbarHandle = waitbar(0,'Please wait. I am currently mapping Neurons to Nuclei.', 'WindowStyle', 'modal');
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

          [xMapped,yMapped,distance] = neuronHandler.FindNucleusForNeuron(col,row,NucleusM, 100,0,0);
          distsum=distsum+distance;
          if(yMapped > 0)
            NeuronM(row,col) = 0;
            NeuronM(yMapped,xMapped) = 1;
          end
       end
       WellNeuronDict(selectedWell) = NeuronM;
       neuronHandler.CellPosMatrix = WellNeuronDict;
       distsum=distsum/numel(nucleusRows);
       disp('Medium distance at well ');
       disp(selectedWell)
       disp(distsum);        
      end
    end
    close(waitbarHandle);
else    
    errordlg('ERROR. Neuron and Nucleus CSV need both to be loaded.');
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
    waitbarHandle = waitbar(0,'Please wait. I am currently mapping Neurons to Nuclei.', 'WindowStyle', 'modal');
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
          xMapped=col;
          yMapped=row;
          distance=0;
          %[xMapped,yMapped,distance] = csvHandler.FindNucleusForNeuron(col,row,NucleusM,100,0,0);
          distsum=distsum+distance;
          if(yMapped ~= -1)
            NeuronM(row,col) = 0;
            NeuronM(yMapped,xMapped) = 1;
          else
              break;
          end
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


function LoadImageFolder(foldername,wellType,channel2,channel3,channel4,v2,v3,v4,handles)
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
  if(isfield(handles,'OptionHandler'))
        optionHandler = handles.OptionHandler;
  else
        CreateOptionHandler(handles);
        handles = guidata(handles.figure1);
        optionHandler = handles.OptionHandler;
  end
  handles.CSVCoordinates = csvHandler;
  guidata(handles.figure1, handles)
  imageHandler.Foldername = foldername;
  maxValueWholeImage=0;
  maxFilterList = zeros(0);
  
    %Check if Subfolder with files already exists.
    subfoldername = strcat(foldername, '/ConvertedCellomics');
    if(~exist(subfoldername, 'dir') || ~numel(dir(subfoldername)) > 10)        
        %Select channel -> Cell Type relation
        str={'Neurons','Oligodendrocytes','Astrocytes','Channel Not Available'}
%         [channel2,v2] = listdlg('PromptString','Which cell Type is Channel 2?',...
%                 'SelectionMode','single',...
%                 'ListString',str);
%             
%         [channel3,v3] = listdlg('PromptString','Which cell Type is Channel 3?',...
%                 'SelectionMode','single',...
%                 'ListString',str);    
%             
%         [channel4,v4] = listdlg('PromptString','Which cell Type is Channel 4?',...
%                 'SelectionMode','single',...
%                 'ListString',str); 
%         if(channel2 == 4)
%             v2=0;
%         end
%         if(channel3==4)
%             v3=0;
%         end
%         if(channel4==4)
%             v4=0;
%         end
        waitbarHandle = waitbar(0,'Please wait. All Cellomics images will be preprocessed. This could take a while. Go and have a coffee!', 'WindowStyle', 'modal');
        selectedWellNumber = get(handles.lbWell,'Value');
        wellList = get(handles.lbWell, 'string');
        selectedWell = wellList{selectedWellNumber};
    
        % Load all images in folder.
        allFiles = dir(foldername); 
        %Check if Subfolder for own files exists. If not: Create
        mkdir(foldername,'ConvertedCellomics');
        subfoldername = [foldername '/ConvertedCellomics'];
    
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
            c1Regexp = strcat('ARRAY.*\d*_(',currentWell,')f(\d+)d0');
            c2Regexp = strcat('ARRAY.*\d*_(',currentWell,')f(\d+)d1');
            c3Regexp = strcat('ARRAY.*\d*_(',currentWell,')f(\d+)d2');
            c4Regexp = strcat('ARRAY.*\d*_(',currentWell,')f(\d+)d3');
            tokensc1 = regexpi(allFiles(i).name, c1Regexp, 'tokens');
            tokensc2 = regexpi(allFiles(i).name, c2Regexp, 'tokens');
            tokensc3 = regexpi(allFiles(i).name, c3Regexp, 'tokens');
            tokensc4 = regexpi(allFiles(i).name, c4Regexp, 'tokens');
        
            if(length(tokensc1) > 0)
                tokensc1 = tokensc1{1};
                fileWell = tokensc1{1};
                oneFileFound = 1;
            elseif(length(tokensc2) > 0)
                tokensc2 = tokensc2{1};
                fileWell = tokensc2{1};
                oneFileFound = 1;
            elseif(length(tokensc3) > 0)
                tokensc3 = tokensc3{1};
                fileWell = tokensc3{1};
                oneFileFound = 1;
            elseif(length(tokensc4) > 0)
                tokensc4 = tokensc4{1};
                fileWell = tokensc4{1};
                oneFileFound = 1;
            end
            if((length(tokensc1) > 1) && strcmp(fileWell, currentWell))
                %Nucleus Pic. Get index number and add to nucleus Dictionary
                index = tokensc1{2};
                imageHandler.NucleusPicArray{str2num(index) + 1} = allFiles(i).name;
            elseif((length(tokensc2) > 1) && strcmp(fileWell, currentWell))
                %Neurite Pic. Get index number and add to neurite Dictionary
                index = tokensc2{2};
                imageHandler.NeuritePicArray{str2num(index) + 1} = allFiles(i).name;
            elseif((length(tokensc3) > 1) && strcmp(fileWell, currentWell))
                %Neurite Pic. Get index number and add to neurite Dictionary
                index = tokensc3{2};
                imageHandler.C3PicArray{str2num(index) + 1} = allFiles(i).name;
           elseif((length(tokensc4) > 1) && strcmp(fileWell, currentWell))
                %Neurite Pic. Get index number and add to neurite Dictionary
                index = tokensc4{2};
                imageHandler.C4PicArray{str2num(index) + 1} = allFiles(i).name;
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
       
        
        
        if(wellType == WellType.Well48)
           maxIndex = numel(csvHandler.Well48Offsets);
        elseif(wellType == WellType.Well4816)
            maxIndex = numel(csvHandler.Well48Offsets16Bit);
        elseif(wellType == WellType.Well96)
            maxIndex = numel(csvHandler.Well96Offsets);
        elseif(wellType == WellType.Well96484)
            maxIndex = numel(csvHandler.Well96Offsets484);
        elseif(wellType == WellType.Well9616)
            maxIndex = numel(csvHandler.Well96Offsets16Bit);
        elseif(wellType == WellType.Well9616255)
            maxIndex = numel(csvHandler.Well96Offsets16Bit255);
        elseif(wellType == WellType.OTWhole)
            maxIndex = numel(csvHandler.WellWholeOTOffsets);
        elseif(wellType == WellType.OTWhole16)
            maxIndex = numel(csvHandler.WellWholeOTOffsets16Bit);
        elseif(wellType == WellType.OTSingle)
            maxIndex = numel(csvHandler.WellSingleChamberOffsets);
        elseif(wellType == WellType.OTSingle16)
            maxIndex = numel(csvHandler.WellSingleChamberOffsets16Bit);
        end
        KeyArray = zeros(maxIndex,3);  
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
        c3Image = zeros(0,0);
        c4Image = zeros(0,0);
        sixteenBit = 0;
    
        %Fill Image Row by Row, with SortedKeyArray
        breaked = 0;
        highestNucleusValue=0;
        highestNeuriteValue=0;
        maxValueWhole=0;
        
        while(i<maxIndex)
            neuriteImageLine = zeros(0,0);
            nucleusImageLine = zeros(0,0);
            c3ImageLine = zeros(0,0);
            c4ImageLine = zeros(0,0);
            while(y==previousY && i<=maxIndex)
              x = SortedKeyArray(i,2);
              y = SortedKeyArray(i,3);
              imageNumber = SortedKeyArray(i,1);
              filenamenuc = imageHandler.NucleusPicArray(imageNumber);
              filenamenuc = filenamenuc{1};
              if(v2)
                filenameneu = imageHandler.NeuritePicArray(imageNumber);
                filenameneu = filenameneu{1};
                imageString = strcat(foldername, '/', filenameneu);
                if(exist(imageString,'file') && ~ isempty(filenameneu) && v2)
                    currentNeuritePic = uint8(imread(imageString));
                else(v2)
                    currentNeuritePic = uint8(zeros(512,512));
                end
              end
              if(v3)
                filenamec3 = imageHandler.C3PicArray(imageNumber);
                filenamec3 = filenamec3{1};
                imageString = strcat(foldername, '/', filenamec3);
                if(exist(imageString,'file') && ~isempty(filenamec3) && v3)
                  currentC3Pic = imread(imageString);              
                else
                  currentC3Pic = zeros(512,512);
                end
              end
              if(v4)
                filenamec4 = imageHandler.C4PicArray(imageNumber);
                filenamec4 = filenamec4{1};
                imageString = strcat(foldername, '/', filenamec4);
                if(exist(imageString,'file') && ~isempty(filenamec4))
                    currentC4Pic = imread(imageString);              
                else
                    currentC4Pic = zeros(512,512);
                end
              end
             
             
              imageString = strcat(foldername, '/', filenamenuc);
              if(exist(imageString,'file') && ~isempty(filenamenuc))
                currentNucleusPic = imread(imageString);
              else
                currentNucleusPic = zeros(512,512);
              end
              maxValueNucleus = uint8(max(currentNucleusPic(:)));
              
              maxValueNeurite=0;
              maxValueC3=0;
              maxValueC4=0;
              
              if(v2)
                maxValueNeurite = uint8(max(currentNeuritePic(:)));
                neuriteImageLine = uint8(cat(2,neuriteImageLine,currentNeuritePic));
              end
              if(v3)
                maxValueC3 = uint8(max(currentC3Pic(:)));
                c3ImageLine = uint8(cat(2,c3ImageLine,currentC3Pic));
              end
              if(v4)
                maxValueC4 = uint8(max(currentC4Pic(:)));
                c4ImageLine = uint8(cat(2,c4ImageLine,currentC4Pic));
              end
              maxValueWhole = uint8(max([maxValueNucleus maxValueNeurite maxValueC3 maxValueC4]));
              maxValueWholeImage = uint8(max([uint8(maxValueWholeImage) maxValueWhole]));
              nucleusImageLine = uint8(cat(2,nucleusImageLine,currentNucleusPic));
              
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
            
            if(wellType == WellType.OTSingle16)            
                if(v2)
                    neuriteImageLine = uint8(flipdim(neuriteImageLine,1));
                end
                nucleusImageLine = uint8(flipdim(nucleusImageLine,1));
                if(v3)
                    c3ImageLine = uint8(flipdim(c3ImageLine,1));
                end
                if(v4)
                    c4ImageLine = uint8(flipdim(c4ImageLine,1));
                end
            end
            if(v2)
                neuriteImage = uint8(cat(1,neuriteImage,neuriteImageLine));
            end
            nucleusImage = uint8(cat(1,nucleusImage,nucleusImageLine));
            if(v3)
                c3Image = uint8(cat(1,c3Image,c3ImageLine));
            end
            if(v4)
                c4Image = uint8(cat(1,c4Image,c4ImageLine));
            end
            end
            if(breaked==1)
                continue;
            end
            
            %If 16 Bit image, invert y Achses on image
            if(wellType == WellType.Well4816 || wellType == WellType.Well9616 || wellType == WellType.OTWhole16) % || wellType == WellType.Well9616255) %|| wellType == WellType.OTSingle16)    
                if(v2)
                   neuriteImage = flipdim(neuriteImage,1);                    
                end
                if(v3)
                    c3Image = flipdim(c3Image,1);                    
                end
                if(v4)
                    c4Image = flipdim(c4Image,1);                    
                end
                nucleusImage = flipdim(nucleusImage,1);
            end
            if(v2)
                neuriteImageSmall = imresize(neuriteImage,optionHandler.ScalingFactor);
                SaveStitchedChannelImage(neuriteImage,neuriteImageSmall,channel2,currentWell,subfoldername, optionHandler.FileEnding)  
            end
            if(v3)
                c3ImageSmall = imresize(c3Image,optionHandler.ScalingFactor);
                SaveStitchedChannelImage(c3Image,c3ImageSmall,channel3,currentWell,subfoldername, optionHandler.FileEnding)  
            end
            if(v4)
                c4ImageSmall = imresize(c4Image,optionHandler.ScalingFactor);
                SaveStitchedChannelImage(c4Image,c4ImageSmall,channel4,currentWell,subfoldername, optionHandler.FileEnding)  
            end            
                        
            imageHandler.ResizedNucleusImage = imresize(nucleusImage,optionHandler.ScalingFactor);                
            imageHandler.NucleusImage = nucleusImage;
    
            %For every well: Save Original Pics in Subfolder ConvertedCellomics
            
            filepath = [subfoldername '/' currentWell 'NucleusBig' optionHandler.FileEnding];
            imwrite(nucleusImage,filepath); 
            
            filepath = [subfoldername '/' currentWell 'NucleusSmall' optionHandler.FileEnding];
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

function SaveStitchedChannelImage(bigImage,smallImage,channelIndex,currentWell,subfoldername, fileEnding)
if(channelIndex==1)
    channelName='Neurite';
elseif(channelIndex==2)
    channelName='Oligo';
else
    channelName='Astro';
end
filepath = [subfoldername '/' currentWell channelName 'Big' fileEnding];
imwrite(bigImage,filepath);
filepath = [subfoldername '/' currentWell channelName 'Small' fileEnding];
imwrite(smallImage,filepath);
    
% --------------------------------------------------------------------
function LoadImageFolder_Callback(hObject, eventdata, handles)
% hObject    handle to LoadImageFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Restructure this method.
% Currently: Image Folder is saved. Only current Well is saved in memory.
% New: Save Image Folder. Create Pics for all Wells and save them in own
% SubFolder of Folder.
try
foldername = uigetdir;
wellType = WellType.Undefined;
subfoldername = strcat(foldername, '/ConvertedCellomics');
channel2=0;
channel3=0;
channel4=0;
v2=0;
v3=0;
v4=0;
if(numel(foldername) > 1 && ~exist(subfoldername, 'dir') || ~numel(dir(subfoldername)) > 10)        
    %Select channel -> Cell Type relation
    str={'Neurons','Oligodendrocytes','Astrocytes','Channel Not Available'}
    [channel2,v2] = listdlg('PromptString','Which cell Type is Channel 2?',...
            'SelectionMode','single',...
            'ListString',str);

    [channel3,v3] = listdlg('PromptString','Which cell Type is Channel 3?',...
            'SelectionMode','single',...
            'ListString',str);    

    [channel4,v4] = listdlg('PromptString','Which cell Type is Channel 4?',...
            'SelectionMode','single',...
            'ListString',str); 
    if(channel2 == 4)
        v2=0;
    end
    if(channel3==4)
        v3=0;
    end
    if(channel4==4)
        v4=0;
    end
end
if(foldername ~=0)
    LoadImageFolder(foldername,wellType,channel2,channel3,channel4,v2,v3,v4,handles);
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
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
optionHandler = handles.OptionHandler;
ScalingFactor = optionHandler.ScalingFactor;
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
    h = imrect(gca,[xStart/(1/ScalingFactor),yStart/(1/ScalingFactor),512/(1/ScalingFactor), 512/(1/ScalingFactor)]);
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
    [xNuc,yNuc,distance] = csvHandler.FindNucleusForNeuron(xMapped,yMapped,NucleusM,30,0,0);
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
    [xNuc,yNuc,distance] = neuronHandler.FindNucleusForNeuron(xMapped,yMapped,NucleusM,0,0);
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
    [sizeY sizeX]=size(imageHandler.NeuriteImage);
    %Get state of 'ManualButtonGroup'
    if(numel(neuronHandler.ManualNeuronPositionsSparse) ==0) 
        neuronHandler.ManualNeuronPositionsSparse = containers.Map()
    end
    if(numel(csvHandler.ManualPositions1) ==0)        
        csvHandler.ManualPositions1 = containers.Map();
    end 
    if(numel(csvHandler.ManualPositions2) ==0)        
        csvHandler.ManualPositions2 = containers.Map();
    end
    if(numel(csvHandler.ManualPositions3) ==0)        
        csvHandler.ManualPositions3 = containers.Map();
    end
    if(numel(csvHandler.ManualPositions4) ==0)        
        csvHandler.ManualPositions4 = containers.Map();
    end
    if(isKey(csvHandler.ManualPositions1,selectedWell))
            currentNuclei1 = csvHandler.ManualPositions1(selectedWell);           
    else     
            currentNuclei1 = sparse(sizeY,sizeX);            
            csvHandler.ManualPositions1(selectedWell) = currentNuclei1;
    end      
    if(isKey(csvHandler.ManualPositions2,selectedWell))
         currentNuclei2 = csvHandler.ManualPositions2(selectedWell);
    else
        currentNuclei2 = sparse(sizeY,sizeX);
        csvHandler.ManualPositions2(selectedWell) = currentNuclei2;
    end
    if(isKey(csvHandler.ManualPositions3,selectedWell))
         currentNuclei3 = csvHandler.ManualPositions3(selectedWell);
    else
        currentNuclei3 = sparse(sizeY,sizeX);
        csvHandler.ManualPositions3(selectedWell) = currentNuclei1;
    end
    if(isKey(csvHandler.ManualPositions4,selectedWell))
         currentNuclei4 = csvHandler.ManualPositions4(selectedWell);
    else
        currentNuclei4 = sparse(sizeY,sizeX);
        csvHandler.ManualPositions4(selectedWell) = currentNuclei4;
    end
    if(isKey(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
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
    [xNuc,yNuc,distance] = csvHandler.FindNucleusForNeuron(xMapped,yMapped,NucleusM,30,0,0);
    
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
    elseif(ScalingFactor ==1)
    %Left click OK
    %Check if csvHandler.ManualPositions 1-4 exist
    clickType = get(gcf,'SelectionType');
    [sizeY sizeX]=size(imageHandler.NeuriteImage);
    %Get state of 'ManualButtonGroup'
    if(numel(neuronHandler.ManualNeuronPositionsSparse) ==0) 
        neuronHandler.ManualNeuronPositionsSparse = containers.Map()
    end
    if(numel(csvHandler.ManualPositions1) ==0)        
        csvHandler.ManualPositions1 = containers.Map();
    end 
    if(numel(csvHandler.ManualPositions2) ==0)        
        csvHandler.ManualPositions2 = containers.Map();
    end
    if(numel(csvHandler.ManualPositions3) ==0)        
        csvHandler.ManualPositions3 = containers.Map();
    end
    if(numel(csvHandler.ManualPositions4) ==0)        
        csvHandler.ManualPositions4 = containers.Map();
    end
    if(isKey(csvHandler.ManualPositions1,selectedWell))
            currentNuclei1 = csvHandler.ManualPositions1(selectedWell);           
    else     
            currentNuclei1 = sparse(sizeY,sizeX);            
            csvHandler.ManualPositions1(selectedWell) = currentNuclei1;
    end      
    if(isKey(csvHandler.ManualPositions2,selectedWell))
         currentNuclei2 = csvHandler.ManualPositions2(selectedWell);
    else
        currentNuclei2 = sparse(sizeY,sizeX);
        csvHandler.ManualPositions2(selectedWell) = currentNuclei2;
    end
    if(isKey(csvHandler.ManualPositions3,selectedWell))
         currentNuclei3 = csvHandler.ManualPositions3(selectedWell);
    else
        currentNuclei3 = sparse(sizeY,sizeX);
        csvHandler.ManualPositions3(selectedWell) = currentNuclei1;
    end
    if(isKey(csvHandler.ManualPositions4,selectedWell))
         currentNuclei4 = csvHandler.ManualPositions4(selectedWell);
    else
        currentNuclei4 = sparse(sizeY,sizeX);
        csvHandler.ManualPositions4(selectedWell) = currentNuclei4;
    end
    if(isKey(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
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
    
    xMapped = int32(x);
    yMapped = int32(y);
    WellNucleusDict = csvHandler.CellPosMatrix;
    NucleusM = WellNucleusDict(selectedWell);
    [xNuc,yNuc,distance] = csvHandler.FindNucleusForNeuron(xMapped,yMapped,NucleusM,30,0,0);
    
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
    markerX = xNuc;
    markerY = yNuc;    
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
try
if(isfield(handles,'ImageHandler'))
    handles = guidata(handles.figure1);
    imageHandler = handles.ImageHandler;
    [x, y, BW, xi, yi] = roipoly;
    handles.ImageHandler = imageHandler;
    guidata(handles.figure1, handles) 
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
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
if(numel(neuronHandler.ManualNeuronPositions) > 0 && isKey(neuronHandler.ManualNeuronPositions,selectedWell) && numel(imageHandler.ZoomState) > 1)
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
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
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
    
    %[subNucImage subNeuImage] = adjustBrightness(subNucImage, subNeuImage, optionHandler);
    
    %if(~numel(imageHandler.NucleusImage) == 0)              
    %        mat = imfuse(subNucImage,subNeuImage);
    %        hImage = imshow(mat);        
    %        set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
    %        hPixelInfo = impixelinfo;        
    %        imageHandler.MousePosition = hPixelInfo;
    %end
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
    imageHandler.ZoomState = [zoomOnZoomedImage(1)*(1/optionHandler.ScalingFactor) zoomOnZoomedImage(2)*(1/optionHandler.ScalingFactor) zoomOnZoomedImage(3)*(1/optionHandler.ScalingFactor) zoomOnZoomedImage(4)*(1/optionHandler.ScalingFactor)];
%     currentPicNuc = imageHandler.NucleusImage(zoomOnZoomedImage(2)*10:zoomOnZoomedImage(2)*10+zoomOnZoomedImage(4)*10, zoomOnZoomedImage(1)*10:zoomOnZoomedImage(1)*10+zoomOnZoomedImage(3)*10);
%     currentPicNeu = imageHandler.NeuriteImage(zoomOnZoomedImage(2)*10:zoomOnZoomedImage(2)*10+zoomOnZoomedImage(4)*10, zoomOnZoomedImage(1)*10:zoomOnZoomedImage(1)*10+zoomOnZoomedImage(3)*10);
%     %minBrightness = optionHandler.HistogramMinNeurite;
%     %maxBrightness = optionHandler.HistogramMaxNeurite;
%     %currentPicNeu = double(currentPicNeu)./3996;
%     %currentPicNeu = imadjust(currentPicNeu, [double(minBrightness)/3996;double(maxBrightness)/3996], [0;1]);
%     %currentPicNeu = uint8(currentPicNeu.*255);   
%     %minBrightness = optionHandler.HistogramMinNucleus;
%     %maxBrightness = optionHandler.HistogramMaxNucleus;
%     %currentPicNuc = double(currentPicNuc)./3996;
%     %currentPicNuc = imadjust(currentPicNuc, [double(minBrightness)/3996;double(maxBrightness)/3996], [0;1]);
%     %currentPicNuc = uint8(currentPicNuc.*255); 
%     [currentPicNuc currentPicNeu] = adjustBrightness(currentPicNuc, currentPicNeu, optionHandler);
%     
%     mat = imfuse(currentPicNuc,currentPicNeu);    
%     hImage = imshow(mat);        
%     set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
%     hPixelInfo = impixelinfo;        
%     imageHandler.MousePosition = hPixelInfo;
    %Load selected Matrices to GUI
    RefreshZoomedImage(handles);
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
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
optionHandler = handles.OptionHandler;
[subNucImage subNeuImage] = adjustBrightness(imageHandler.ResizedNucleusImage, imageHandler.ResizedNeuriteImage, optionHandler);
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'rbNeuritePic'
        if(isfield(handles,'ImageHandler'))
            imageHandler = handles.ImageHandler;
            imageHandler.ZoomState = zeros(1);            
            if(~numel(imageHandler.NucleusImage) == 0 && ~numel(imageHandler.ResizedNeuriteImage) == 0)                      
                hImage = imshow(subNeuImage);
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
                hImage = imshow(subNucImage);  
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
                mat = imfuse(subNucImage,subNeuImage);
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
            plot((nucleusCols./(1/ScalingFactor)),(nucleusRows./(1/ScalingFactor)),'Linestyle','none','Marker','.','Markersize',2,'Color',[1 .5 0]);       
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
            plot((neuronCols./(1/ScalingFactor)),(neuronRows./(1/ScalingFactor)),'Linestyle','none','Marker','.','Markersize',10,'Color','red');       
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
ChangeWell(selectedWellNumber,handles);


function ChangeWell(selectedWellNumber,handles)
optionHandler = handles.OptionHandler;
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
    if(numel(optionHandler.FileEnding) == 0)
        optionHandler.FileEnding = '.tif';
    end
    foldername = [imageHandler.Foldername '/ConvertedCellomics'];
    imagePathNucleusBig = [foldername '/' selectedWell 'NucleusBig' optionHandler.FileEnding];
    imagePathNucleusSmall = [foldername '/' selectedWell 'NucleusSmall' optionHandler.FileEnding];
    imagePathNeuriteBig = [foldername '/' selectedWell 'NeuriteBig' optionHandler.FileEnding];
    imagePathNeuriteSmall = [foldername '/' selectedWell 'NeuriteSmall' optionHandler.FileEnding];
    imagePathOligoBig = [foldername '/' selectedWell 'OligoBig' optionHandler.FileEnding];
    imagePathOligoSmall = [foldername '/' selectedWell 'OligoSmall' optionHandler.FileEnding];    
    imagePathAstroBig = [foldername '/' selectedWell 'AstroBig' optionHandler.FileEnding];
    imagePathAstroSmall = [foldername '/' selectedWell 'AstroSmall' optionHandler.FileEnding];    
    imagePathSkeleton = [foldername '/' selectedWell 'Skeleton' optionHandler.FileEnding];
    imagePathBinary = [foldername '/' selectedWell 'Binary' optionHandler.FileEnding];
    if(exist(imagePathNeuriteSmall,'file'))
        imageHandler.NeuriteImage = imread(imagePathNeuriteBig);
        imageHandler.ResizedNeuriteImage = imread(imagePathNeuriteSmall);
    end
    if(exist(imagePathOligoSmall,'file'))
        imageHandler.OligoImage = imread(imagePathOligoBig);
        imageHandler.ResizedOligoImage = imread(imagePathOligoSmall);
    end
    if(exist(imagePathAstroSmall,'file'))
        imageHandler.AstroImage = imread(imagePathAstroBig);
        imageHandler.ResizedAstroImage = imread(imagePathAstroSmall);
    end
    
    imageHandler.ResizedNucleusImage = imread(imagePathNucleusSmall);
    imageHandler.NucleusImage = imread(imagePathNucleusBig);
    
    if(exist(imagePathBinary,'file'))
        imageHandler.BinaryImage = imread(imagePathBinary);
        imageHandler.ResizedBinaryImage = imresize(imageHandler.BinaryImage,optionHandler.ScalingFactor);
    end
    if(exist(imagePathSkeleton,'file'))
        imageHandler.SkeletonImage = imread(imagePathSkeleton);
        imageHandler.ResizedSkeletonImage = imresize(imageHandler.SkeletonImage,optionHandler.ScalingFactor);
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
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
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
    if(rowNumber == 0)
    [rowNumber, colNumber]=size(imageHandler.NucleusImage);
    end
    A = optionHandler.RectangleSizeA;
    B = optionHandler.RectangleSizeB;
    if(imageHandler.ZoomState)
        h = imrect(gca,[1,1,(A), B]);
    else
        h = imrect(gca,[1,1,(A * optionHandler.ScalingFactor),(B * optionHandler.ScalingFactor)]);
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
        %Wandle alle Positionen in Position f�r das gro�e Bild um        
        currentFilter{filterEntries+1} = poly2mask(polyPos(:,1)+xmin,polyPos(:,2)+ymin,colNumber, rowNumber);
        %currentFilter{filterEntries+2} = poly2mask(polyPos2(:,1)+xmin,polyPos2(:,2)+ymin,colNumber, rowNumber);
    else
        currentFilter{filterEntries+1} = createMask(h);
        %currentFilter{filterEntries+2} = createMask(h2);
        currentFilter{filterEntries+1} = sparse(imresize(currentFilter{filterEntries+1},[rowNumber, colNumber]));
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
    [rowNumber, colNumber]=size(imageHandler.NucleusImage);
    if(imageHandler.ZoomState)
        xmin = imageHandler.ZoomState(1);
        ymin = imageHandler.ZoomState(2);  
        polyPos = getPosition(h);
        %Wandle alle Positionen in Position f�r das gro�e Bild um        
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
        waitbarHandle = waitbar(0,'Set Filter for Sphere automatically.', 'WindowStyle', 'modal');
        selectedWellNumber = get(handles.lbWell,'Value');
        wellList = get(handles.lbWell, 'string');
        selectedWell = wellList{selectedWellNumber};
        waitbar((20/100),waitbarHandle);
        SphereAreaSizeX = optionHandler.MigrationDistanceDensityImageXSize;
        SphereAreaSizeY = optionHandler.MigrationDistanceDensityImageYSize;
        [sizeY sizeX] = size(imageHandler.NeuriteImage);
        NucleusM = csvHandler.CellPosMatrix(selectedWell);
        [filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, optionHandler, NucleusM, imageHandler.Foldername, sizeY, sizeX);         
        
        waitbar((60/100),waitbarHandle);
        [rowNumber, colNumber]=size(imageHandler.NeuriteImage);
        sizingFactorX2 = colNumber/sizingFactorX;
        sizingFactorY2 = rowNumber/sizingFactorY;     
        %Calculate Factor for preview Pic
        sizingFactorX2 = sizingFactorX2 * optionHandler.ScalingFactor;
        sizingFactorY2 = sizingFactorY2 * optionHandler.ScalingFactor;
        allowedDistanceNucleus = 3*((sizingFactorX2 + sizingFactorY2)/2);
        allowedDistanceSphere = 10*((sizingFactorX2 + sizingFactorY2)/2);
        rowColMatSorted = sortMatrixByNearestNeighbour(SphereArea,allowedDistanceSphere,sizingFactorX,sizingFactorY,handles);
        waitbar((70/100),waitbarHandle);
        %negativeRowColMatSorted = sortMatrixByNearestNeighbour(NucleusArea,allowedDistanceNucleus,sizingFactorX,sizingFactorY,handles);
      %  newNucleusMask = imresize(NucleusArea,[rowNumber,colNumber]);
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
    if(numel(lbFilterString) == 0)
        lbFilterString=cell(0);
        lbFilterString{1} = num2str(filterEntries);
    else
        lbFilterString = [lbFilterString;num2str(filterEntries)];
    end
    set(handles.lbFilter, 'String', lbFilterString);

    %Save to file
    foldername = imageHandler.Foldername;
    subfoldername = [foldername '/ConvertedCellomics'];
    FMask = FilterMask();
    FMask.PositiveFilters = imageHandler.PositiveFilters;
    FMask.NegativeFilters = imageHandler.NegativeFilters;
    save('-v7.3',strcat(subfoldername,'/',selectedWell),'FMask');
    FMask=0;
    imageHandler.PositiveFilters = 0;
    imageHandler.NegativeFilters = 0;
end
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles)
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


function rowColMatSorted = sortMatrixByNearestNeighbour(SphereArea,allowedMaxDistance,AreaSizeX, AreaSizeY,handles)
imageHandler=handles.ImageHandler;
[NonZeroRows NonZeroCols] = find(SphereArea);
[rowNumber, colNumber]=size(imageHandler.NeuriteImage);
sizingFactorX = colNumber/AreaSizeX;
sizingFactorY = rowNumber/AreaSizeY;
%Divide by 10 for Preview Pic
sizingFactorX = sizingFactorX *  optionHandler.ScalingFactor;
sizingFactorY = sizingFactorY *  optionHandler.ScalingFactor;
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
try
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
        %Wandle alle Positionen in Position f�r das gro�e Bild um        
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
save('-v7.3',strcat(subfoldername,'/',selectedWell),'FMask');
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles)
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


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
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
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
filterMaskSmall = imresize(filterMask,optionHandler.ScalingFactor);    
     nucImageResized = double(imageHandler.ResizedNucleusImage).*double(filterMaskSmall);
    
     nucResizedTransposed = transpose(nucImageResized);
     %nucImageTransposed = transpose(nucImage);
     %Cut Matrix and set new Zoom State
     [maxPossibleRow, maxPossibleCol] = (size(nucImageResized));
     maxPossibleRow = maxPossibleRow /optionHandler.ScalingFactor;
     maxPossibleCol = maxPossibleCol /optionHandler.ScalingFactor;
     [minRow, minCol] = (find(nucImageResized,1,'first'));
     minRow = minRow /optionHandler.ScalingFactor;
     minCol = minCol /optionHandler.ScalingFactor;
     [minColTrans, minRowTrans] = find(nucResizedTransposed,1,'first');
     minRowTrans = minRowTrans / optionHandler.ScalingFactor;
     minColTrans = minColTrans / optionHandler.ScalingFactor;
     minRow = min(minRow, minRowTrans);
     minCol = min(minCol, minColTrans);
     [maxRow, maxCol] = (find(nucImageResized,1,'last'));
     maxRow = maxRow / optionHandler.ScalingFactor;
     maxCol = maxCol / optionHandler.ScalingFactor;
     [maxColTrans, maxRowTrans] = find(nucResizedTransposed,1,'last');    
     maxColTrans = maxColTrans / optionHandler.ScalingFactor;    
     maxRowTrans = maxRowTrans / optionHandler.ScalingFactor;
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
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --------------------------------------------------------------------
function DeleteNeuronsManual_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteNeuronsManual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
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
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --------------------------------------------------------------------
function CountNeuronsManual_Callback(hObject, eventdata, handles)
% hObject    handle to CountNeuronsManual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
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
guidata(handles.figure1, handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on selection change in lbFilter.
function lbFilter_Callback(hObject, ~, handles)
% hObject    handle to lbFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% If available: Apply Filter mask on combined image.
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
filterNumber = get(hObject,'Value');
if(isfield(handles,'OptionHandler'))
	optionHandler = handles.OptionHandler;
else
	optionHandler = Option();
end

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
    filterMask = (filterHandle);
    %Remove excluded parts of picture
    if(iscell(imageHandler.NegativeFilters) > 0)
        negativeMask = imageHandler.NegativeFilters;
        for i=1:numel(negativeMask)
          workMask = negativeMask{i};
          workMask = int8(workMask -1);
          workMask = logical(workMask .* (-1));
          workMask = sparse(logical(workMask));
          filterMask = (logical((filterMask .* workMask)));
        end
    end
    workMask=0;
    %resize filter mask for pic
    filterMaskSmall = sparse(imresize(full(filterMask),optionHandler.ScalingFactor));
    %Multiply mask with image
    nucImage = uint8(imageHandler.NucleusImage).*uint8(full(filterMask));
    %We need to also do this for the other channels:
    if(size(imageHandler.NeuriteImage)>0)
    neuImage = uint8(imageHandler.NeuriteImage).*uint8(full(filterMask));
    neuImagePoly = neuImage;
    else
    neuImagePoly = 0;
    neuImage = 0;
    end
    if(size(imageHandler.OligoImage)>0)
    oliImage = uint8(imageHandler.OligoImage).*uint8(full(filterMask)); 
    oliImagePoly = oliImage;
    else
    oliImagePoly = 0;
    oliImage = 0;
    end
    if(size(imageHandler.AstroImage)>0)
    astImage = uint8(imageHandler.AstroImage).*uint8(full(filterMask));
    astImagePoly = astImage;
    else
    astImagePoly = 0; 
    astImage = 0;
    end
    if(size(imageHandler.BinaryImage)>0)
    binImage = uint8(imageHandler.BinaryImage).*uint8(full(filterMask)); 
    binImagePoly = binImage;
    else
    binImagePoly = 0;   
    binImage = 0;
    end
    if(size(imageHandler.SkeletonImage)>0)
    skelImage = uint8(imageHandler.SkeletonImage).*uint8(full(filterMask)); 
    skelImagePoly = skelImage;
    else
    skelImagePoly = 0;  
    skelImage = 0;
    end
    % Here calculate the area of the Mask:
    AreaMask = sum(sum(filterMask));
    nucImagePoly = nucImage;
    %neuImage = uint8(imageHandler.NeuriteImage).*uint8(filterMask);
    % Here we have to discriminate between rectangular filters and
    % polygonal filters!
    % Therefore we will store the filterMask under a new name!
    filterMaskPoly = filterMask;
    filterMask=0;    
    %Transpose only resized image, because Transposing is very cost
    %intensive
    nucImageResized = logical(logical(imageHandler.ResizedNucleusImage).*logical(filterMaskSmall));
    nucResizedTransposed = transpose(nucImageResized);
    %nucImageTransposed = transpose(nucImage);
    %Cut Matrix and set new Zoom State
    [maxPossibleRow, maxPossibleCol] = size(nucImage);
    [minRow, minCol] = find(nucImage,1,'first');
    [minColTrans, minRowTrans] = find(nucResizedTransposed,1,'first');
    minRowTrans = minRowTrans  / optionHandler.ScalingFactor;
    minColTrans = minColTrans  / optionHandler.ScalingFactor;
    minRow = min(minRow, minRowTrans);
    minCol = min(minCol, minColTrans);
    [maxRow, maxCol] = find(nucImage,1,'last');
    [maxColTrans, maxRowTrans] = find(nucResizedTransposed,1,'last');
    
    maxColTrans = maxColTrans / optionHandler.ScalingFactor;
    
    maxRowTrans = maxRowTrans / optionHandler.ScalingFactor;
    
    maxRow = max(maxRow,maxRowTrans);
    maxCol = max(maxCol,maxColTrans);
    maxRow = min(maxRow, maxPossibleRow);
    maxCol = min(maxCol, maxPossibleCol);
    %Set correct zoom State
    %imageHandler.ShowNucleusImage = nucImage(minRow:maxRow,minCol:maxCol);
    %imageHandler.ShowNeuriteImage = neuImage(minRow:maxRow,minCol:maxCol);
    
%     if(~numel(imageHandler.NucleusImage) == 0 && ~numel(imageHandler.ResizedNeuriteImage) == 0)              
%             [subNucImage subNeuImage] = adjustBrightness(imageHandler.ShowNucleusImage, imageHandler.ShowNeuriteImage, handles.OptionHandler);    
%             mat = imfuse(subNucImage,subNeuImage);
%             hImage = imshow(mat);        
%             set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});
%             hPixelInfo = impixelinfo;        
%             imageHandler.MousePosition = hPixelInfo;
%     end
    
    imageHandler.ZoomState = [minCol minRow maxCol-minCol maxRow-minRow];
    % Here we calculate the area of the Zoom to sea whether we got a
    % rectengular shape or a polygon:
    AreaBox = ((maxCol-minCol)+1)*((maxRow-minRow)+1);
    if((AreaBox/AreaMask)==1)
        D= 0
    else
        D =1
    end    
    imageHandler.PositiveFilters = 0;
    imageHandler.NegativeFilters = 0;
    FMask = 0;
    if((AreaBox/AreaMask)==1)
        imageHandler.NucleusImage = nucImage;
        imageHandler.NeuriteImage = neuImage;
        imageHandler.OligoImage = oliImage;
        imageHandler.AstroImage = astImage;
        imageHandler.BinaryImage = binImage;
        imageHandler.SkeletonImage = skelImage;
        handles.ImageHandler = imageHandler;
        guidata(handles.figure1, handles);
    else    
    imageHandler.NucleusImage = nucImagePoly;
    if(size(neuImagePoly)>0);
        imageHandler.NeuriteImage = neuImagePoly;
    end  
    if(size(oliImagePoly)>0);
        imageHandler.OligoImage = oliImagePoly;
    end  
    if(size(astImagePoly)>0);
        imageHandler.AstroImage = astImagePoly;
    end  
    if(size(binImagePoly)>0);
        imageHandler.BinaryImage = binImagePoly;
    end  
    if(size(skelImagePoly)>0);
        imageHandler.SkeletonImage = skelImagePoly;
    end  
    handles.ImageHandler = imageHandler;
    guidata(handles.figure1, handles);
    end
    RefreshZoomedImage(handles);
    
    handles.ImageHandler = imageHandler;
    guidata(handles.figure1, handles);
    LoadManualNeuronPositionsToGUI(handles);
    guidata(handles.figure1, handles);
    
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

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

function SaveFixedNeuronPositions(filename,pathname,handles)
handles = guidata(handles.figure1);
imageHandler = 0;
neuronHandler = 0;
csvHandler = 0;
optionHandler = 0;
oligoHandler = 0;
if(isfield(handles,'ImageHandler'))
    imageHandler = handles.ImageHandler;
    
    %Set not needed Fields of imageHandler to 0, to reduce amount of needed
    %memory and prevent Out of Memory exception
    neuriteImage = imageHandler.NeuriteImage;
    nucleusImage = imageHandler.NucleusImage;
    oligoImage = imageHandler.OligoImage;
    astroImage = imageHandler.AstroImage;
    resizedNeuriteImage = imageHandler.ResizedNeuriteImage;
    resizedNucleusImage = imageHandler.ResizedNucleusImage;
    resizedOligoImage = imageHandler.ResizedOligoImage;
    resizedAstroImage = imageHandler.ResizedAstroImage;
    if(isprop(imageHandler,'ShowNucleusImage'))
        showNucleusImage = imageHandler.ShowNucleusImage;
        showNeuriteImage = imageHandler.ShowNeuriteImage;
    else
        showNucleusImage=0;
        showNeuriteImage=0;
    end
    imageHandler.NeuriteImage=0;
    imageHandler.NucleusImage=0;
    imageHandler.OligoImage=0;
    imageHandler.AstroImage=0;
    imageHandler.ResizedOligoImage=0;
    imageHandler.ResizedAstroImage=0;
    imageHandler.ResizedNeuriteImage=0;
    imageHandler.ResizedNucleusImage=0;
    imageHandler.ShowNucleusImage=0;
    imageHandler.ShowNeuriteImage=0;
    imageHandler.NeuritePicArray = 0;
    imageHandler.NucleusPicArray = 0;
    imageHandler.MousePosition=0;
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
if(isfield(handles,'OligoCoordinates'))
    oligoHandler = handles.OligoCoordinates;
end
save('-v7.3',strcat(pathname,filename),'imageHandler','neuronHandler','csvHandler','optionHandler','oligoHandler');
%Zip file
zip(strcat(pathname,filename,'.zip'),strcat(pathname,filename));
delete(strcat(pathname,filename));
if(isfield(handles,'ImageHandler'))
    imageHandler.NeuriteImage = neuriteImage;
    imageHandler.NucleusImage = nucleusImage;
    imageHandler.ResizedNeuriteImage = resizedNeuriteImage;
    imageHandler.ResizedNucleusImage = resizedNucleusImage;
    imageHandler.OligoImage=oligoImage;
    imageHandler.AstroImage=astroImage;
    imageHandler.ResizedOligoImage=resizedOligoImage;
    imageHandler.ResizedAstroImage=resizedAstroImage;
    imageHandler.ShowNucleusImage = showNucleusImage;
    imageHandler.ShowNeuriteImage = showNeuriteImage;
end
disp('Saved');

% --------------------------------------------------------------------
function SaveFixedNeuronPositions_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFixedNeuronPositions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Iterate over Neuron Coordinates and Save mapped data to CSV
try
[FileName,PathName] = uiputfile;
SaveFixedNeuronPositions(FileName,PathName,handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


function Load(pathname,filename,handles)
handles = guidata(handles.figure1);
if(any(regexp(filename,'.zip$')))
    unzip(strcat(pathname,filename));
    filename=filename(1:length(filename)-4)
    load(filename);
    delete(filename);
else
    load(strcat(pathname,filename));
end


title(filename);
%New: One file per well for all Filters
%Load and save on every lbWellClick
set(gcf,'name',filename);
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
        save('-v7.3',strcat(subfoldername,'/',currentWell),'FMask');
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
if(exist('oligoHandler'))
    oligoHandler.ManualCountingMode=0;
    handles.OligoCoordinates = oligoHandler;
end
if(exist('csvHandler'))
    handles.CSVCoordinates = csvHandler;
    guidata(handles.figure1, handles)
    ReloadWellList(handles,0);
end
if(exist('optionHandler'))
    handles.OptionHandler = optionHandler;
end
guidata(handles.figure1, handles)

% --------------------------------------------------------------------
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
[filename, pathname] = uigetfile({'*.mat;*.zip'}, 'Select MAT File or zipped MAT File');
if(~strcmp(class(filename),'double'))
    Load(pathname,filename,handles);
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


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
try
ConvertPicturesTo8Bit(handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
function ConvertPicturesTo8Bit(handles,foldername)
handles = guidata(handles.figure1);
allFiles = dir(foldername); 
waitbarHandle = waitbar(0,'Please wait. All images are converted to 8 Bit.', 'WindowStyle', 'modal');
for i=3:numel(allFiles)
    waitbar((i/numel(allFiles)),waitbarHandle);
    currentFile=allFiles(i).name;
    if(strfind(currentFile, 'Neurite'))
        minBrightness = 0;
        maxBrightness = 4095;
    else
        minBrightness = 0;
        maxBrightness = 4095;
    end
    imageString = strcat(foldername, '/', currentFile);
    if(length(findstr(imageString, '.mat')) == 0 && (length(findstr(imageString, '.tif')) > 0 || length(findstr(imageString, '.TIF')) > 0))
      currentPic = imread(imageString);
      deleteIndices = uint16(find(abs(currentPic) < minBrightness));
      currentPic(deleteIndices) = 0;
    
      %currentPic = uint8(currentPic);
      currentPic = uint8(uint32(uint32(currentPic).*255)./4096);    
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


function [returnPic] = adjustBrightness(inputPic,picIndex, optionHandler)
if(numel(picIndex) == 4 & picIndex == 'anuc')
    minBrightness = optionHandler.HistogramMinNucleus;
    maxBrightness = optionHandler.HistogramMaxNucleus;    
elseif(numel(picIndex) == 3 & picIndex == 'neu')
    minBrightness = optionHandler.HistogramMinNeurite;
    maxBrightness = optionHandler.HistogramMaxNeurite;    
elseif(numel(picIndex) == 5 & picIndex == 'oligo')
    minBrightness = optionHandler.HistogramMinOligo;
    maxBrightness = optionHandler.HistogramMaxOligo;    
elseif(numel(picIndex) == 5 & picIndex == 'astro')
    minBrightness = optionHandler.HistogramMinAstro;
    maxBrightness = optionHandler.HistogramMaxAstro;    
end
if(exist('maxBrightness','var') && maxBrightness > 255)
        returnPic = imadjust(inputPic);
elseif(exist('maxBrightness','var'))
        returnPic = imadjust((inputPic), [double(minBrightness)/255;double(maxBrightness)/255], [0;1]);
else
    returnPic = inputPic;
end    

% --- Executes on button press in pb8bitPreview.
function pb8bitPreview_Callback(hObject, eventdata, handles)
% hObject    handle to pb8bitPreview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
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
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

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
optionHandler = handles.OptionHandler;
SphereAreaSizeX = optionHandler.MigrationDistanceDensityImageXSize;
SphereAreaSizeY = optionHandler.MigrationDistanceDensityImageYSize;
NucleusM = csvHandler.CellPosMatrix(selectedWell);
[sizeY sizeX] = size(imageHandler.NeuriteImage);
[filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, optionHandler, NucleusM, imageHandler.Foldername, sizeY, sizeX); 

hPixelInfo = impixelinfo;        
imageHandler.MousePosition = hPixelInfo;
disp(['Migration distance is ' num2str(filterDistance) ' pixel. Without filter it would be ' num2str(nonFilterDistance) ' pixel.']); 


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
if(isfield(handles,'OligoCoordinates'))
    oligoHandler = handles.OligoCoordinates;
end
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
      selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
elseif(length(selectedWell) == 4)
        selectedWellLong = [selectedWell(1) '0' selectedWell(2:length(selectedWell))];
end
%Check which checkboxes are checked.
checkedPicturesMap = containers.Map();
if(get(handles.cbNucleusPicture ,'Value'))
    if(imageHandler.ZoomState)
        checkedPicturesMap('anuc') =  imageHandler.NucleusImage;
    else    
        checkedPicturesMap('anuc') =  imageHandler.ResizedNucleusImage;
    end
end
if(get(handles.cbNeuritePicture ,'Value'))
    if(imageHandler.ZoomState)
        checkedPicturesMap('neu') = imageHandler.NeuriteImage;
    else
        checkedPicturesMap('neu') = imageHandler.ResizedNeuriteImage;
    end
end
if(get(handles.cbOligoPicture, 'Value'))
    if(imageHandler.ZoomState)
        checkedPicturesMap('oligo') = imageHandler.OligoImage;
    else
        checkedPicturesMap('oligo') = imageHandler.ResizedOligoImage;
    end
end
if(get(handles.cbAstroPicture, 'Value'))
    if(imageHandler.ZoomState)
        checkedPicturesMap('astro') = imageHandler.AstroImage;
    else
        checkedPicturesMap('astro') = imageHandler.ResizedAstroImage;
    end
end
if(get(handles.cbSkeletonPic, 'Value'))
    if(imageHandler.ZoomState)
        checkedPicturesMap('skel') = imageHandler.SkeletonImage;
    else
        checkedPicturesMap('skel') = imageHandler.ResizedSkeletonImage;
    end
end
if(get(handles.cbBinaryPic , 'Value'))
    if(imageHandler.ZoomState)
        checkedPicturesMap('bin') = imageHandler.BinaryImage;
    else
        checkedPicturesMap('bin') = imageHandler.ResizedBinaryImage;
    end
end

cbRemoveCoreChecked = get(handles.cbRemoveCore, 'Value');
cbPreprocessedImagesChecked = get(handles.cbPreprocessedImages, 'Value');

%Analog to RefreshUnzoomed image but care about cutting Big Pictures before
zoomOnZoomedImage = imageHandler.ZoomState;


if(cbPreprocessedImagesChecked)
    for(i=1:numel(checkedPicturesMap.keys))
        currentKey = checkedPicturesMap.keys;
        currentKey = currentKey(i);
        currentKey=currentKey{1};
        currentImage=checkedPicturesMap(currentKey);

        [sizeYorig sizeXorig] = size(imageHandler.NucleusImage);

        NucleusM = csvHandler.CellPosMatrix(selectedWell);
        foldername = imageHandler.Foldername;
    
        if(optionHandler.CutOutSphereCore == 1)
            currentImage = CutOutCircles(currentImage,selectedWell,1,1,0,optionHandler, foldername, sizeYorig, sizeXorig, NucleusM, 0);
        end
        
       if(strcmp(currentKey, 'anuc'))    
           imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
           if(exist(imagePathWatershed,'file'))
               currentImage = imread(imagePathWatershed);
           else
               currentImage(currentImage<optionHandler.NucleusThreshold) = 0;
               currentImage(currentImage>=optionHandler.NucleusThreshold) = 255;
               currentImage = uint8(currentImage);
           end
        %currentImage(find(currentImage<optionHandler.NucleusThreshold)) = 0;%double(double(imageHandler.NucleusImage) ./ 255);
       end
        stretchlimLow = neuronHandler.StretchlimResult;
       if(strcmp(currentKey,'neu'))
        [currentImage neuritePicStrong lMin] = ThresholdPic(currentImage,optionHandler, stretchlimLow, sizeYorig, sizeXorig, 0, 'Neurite');
        SE = strel('disk', 3);
        currentImage = imdilate(currentImage,SE);
        %currentImage = imadjust(currentImage,stretchlimLow,[0 1]);
        %currentImage = medfilt2(currentImage);
        %unsharpFilter = fspecial('unsharp');
        %currentImage = imfilter(currentImage,unsharpFilter);    
        %currentPicNeu = imcomplement(currentPicNeu);
        %level = optionHandler.SkeletonNeuriteThresholdLow;
        %currentImage(find(currentImage<(255-(level*255))))=0;
       elseif(strcmp(currentKey,'oligo'))
        [currentImage neuritePicStrong lMin] = ThresholdPic(currentImage,optionHandler, stretchlimLow, sizeYorig, sizeXorig, 0, 'Oligo');
        SE = strel('disk', 3);
        currentImage = imdilate(currentImage,SE);
        %currentImage = imadjust(currentImage,stretchlimLow,[0 1]);
        %currentImage = medfilt2(currentImage);
        %unsharpFilter = fspecial('unsharp');
        %currentImage = imfilter(currentImage,unsharpFilter);    
        %currentPicNeu = imcomplement(currentPicNeu);
        %level = optionHandler.SkeletonNeuriteThresholdLow;
        %currentImage(find(currentImage<(255-(level*255))))=0;
       end       
       %checkedPicturesMap(currentKey) = adjustBrightness(currentImage,currentKey, optionHandler);
       checkedPicturesMap(currentKey) = currentImage;
    end
elseif(cbRemoveCoreChecked)    
    NucleusM = csvHandler.CellPosMatrix(selectedWell);
    foldername = imageHandler.Foldername;
    [sizeY sizeX] = size(imageHandler.NeuriteImage);
    for(i=1:numel(checkedPicturesMap.keys))
        currentKey = checkedPicturesMap.keys;
        currentKey = currentKey(i);
        currentKey=currentKey{1};       
        if(optionHandler.CutOutSphereCore == 1)
            checkedPicturesMap(currentKey) = CutOutCircles(checkedPicturesMap(currentKey),selectedWell,1, 1,0,optionHandler, foldername, sizeY, sizeX, NucleusM, imageHandler);    
        end
    end
end
for(i=1:numel(checkedPicturesMap.keys))
        currentKey = checkedPicturesMap.keys;
        currentKey = currentKey(i);
        currentKey=currentKey{1};
        currentValue = checkedPicturesMap(currentKey);
        if(imageHandler.ZoomState)
            checkedPicturesMap(currentKey) = currentValue(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));     
        end
        if(cbPreprocessedImagesChecked==0)
            checkedPicturesMap(currentKey) = adjustBrightness(checkedPicturesMap(currentKey),currentKey,handles.OptionHandler);
        end
        %checkedPicturesMap(currentKey) = checkedPicturesMap(currentKey);
end

WellNucleusDict = csvHandler.CellPosMatrix;
NucleusM = WellNucleusDict(selectedWell);
if(imageHandler.ZoomState)
    NucleusM = NucleusM(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
end
WellNeuronDict = neuronHandler.CellPosMatrix;
%NeuronM = WellNeuronDict(selectedWell);
% if(imageHandler.ZoomState)
%     NeuronM = NeuronM(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
% end


WellSkeletonDict = neuronHandler.NeuronPositionsSkeletonization;

if(numel(WellSkeletonDict) > 0 && isKey(WellSkeletonDict,selectedWell))
    SkeletonM = WellSkeletonDict(selectedWell);
    if(imageHandler.ZoomState)
        SkeletonM = SkeletonM(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
    end
end
if(isfield(handles,'NeuronCoordinates2'))
 WellNeuronDict2 = neuronHandler2.CellPosMatrix;
 NeuronM2 = WellNeuronDict2(selectedWell);
 if(imageHandler.ZoomState)
    NeuronM2 = NeuronM2(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
 end
end

% WellManualDict = neuronHandler.ManualNeuronPositionsSparse;
% ManualM = WellManualDict(selectedWell);
% if(imageHandler.ZoomState)
%     ManualM = ManualM(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
% end
%WellEdgeCompDict = neuronHandler.NeuronPositionsEdgeComposite;
WellEdgeFilledDict = neuronHandler.NeuronPositionsEdgeFill;
%EdgeCompM = WellEdgeCompDict(selectedWell);
%EdgeCompM = EdgeCompM(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
if(isa(WellEdgeFilledDict,'containers.Map') && isKey(WellEdgeFilledDict,selectedWell))
    EdgeFilledM = WellEdgeFilledDict(selectedWell);
    if(imageHandler.ZoomState)
        EdgeFilledM = EdgeFilledM(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
    end
    WellDeletedDict = neuronHandler.NeuronPositionsSkelDeleted;

    if(numel(WellDeletedDict) > 0 && isKey(WellDeletedDict,selectedWell))
        DeletedM = WellDeletedDict(selectedWell);
        if(imageHandler.ZoomState)
            DeletedM = DeletedM(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
        end
    end
end


if(numel(checkedPicturesMap.keys) == 4)
    keyOne = checkedPicturesMap.keys;
    keyOne = keyOne(1);
    keyOne=keyOne{1};
    keyTwo = checkedPicturesMap.keys;
    keyTwo = keyTwo(2);
    keyTwo=keyTwo{1};
    keyThree = checkedPicturesMap.keys;
    keyThree = keyThree(3);
    keyThree=keyThree{1};
    keyFour = checkedPicturesMap.keys;
    keyFour = keyFour(4);
    keyFour=keyFour{1};
    image = cat(3,checkedPicturesMap(keyThree)+checkedPicturesMap(keyFour),checkedPicturesMap(keyTwo)+checkedPicturesMap(keyFour),checkedPicturesMap(keyOne));   
elseif(numel(checkedPicturesMap.keys) == 3)
    keyOne = checkedPicturesMap.keys;
    keyOne = keyOne(1);
    keyOne=keyOne{1};
    keyTwo = checkedPicturesMap.keys;
    keyTwo = keyTwo(2);
    keyTwo=keyTwo{1};
    keyThree = checkedPicturesMap.keys;
    keyThree = keyThree(3);
    keyThree=keyThree{1};
    
    image = cat(3,checkedPicturesMap(keyThree),checkedPicturesMap(keyTwo),checkedPicturesMap(keyOne));          
else
    if(numel(checkedPicturesMap.keys) == 0)
        [sizeYorig sizeXorig] = size(imageHandler.NucleusImage);
        sizeYorig = 717;
        sizeXorig = 717;
        image = logical(ones(sizeYorig, sizeXorig));
    end
    for(i=1:numel(checkedPicturesMap.keys))
        currentKey = checkedPicturesMap.keys;
        currentKey = currentKey(i);
        currentKey=currentKey{1};
        currentValue = checkedPicturesMap(currentKey);
        if(i==1)
            image=currentValue;
        else
            image=imfuse(image,currentValue);
        end
    end
end
if(exist('image','var'))
    hImage = imshow(image);
    set(hImage,'ButtonDownFcn',{@executeOnImageClick, handles});         
    hPixelInfo = impixelinfo;        
    imageHandler.MousePosition = hPixelInfo;
end
         


%New:
listboxstrings = get(handles.popupPlus1);
listboxindex = get(handles.popupPlus1, 'value');
listboxstring = listboxstrings.String(listboxindex); 
MatrixPlus1 = GetMatrixFromListboxString(listboxstring,1,selectedWell,handles);
listboxindex = get(handles.popupPlus2, 'value');
listboxstring = listboxstrings.String(listboxindex); 
MatrixPlus2 = GetMatrixFromListboxString(listboxstring,1,selectedWell,handles);
listboxindex = get(handles.popupMinus1, 'value');
listboxstring = listboxstrings.String(listboxindex); 
MatrixMinus1 = GetMatrixFromListboxString(listboxstring,1,selectedWell,handles);
listboxindex = get(handles.popupMinus2, 'value');
listboxstring = listboxstrings.String(listboxindex); 
MatrixMinus2 = GetMatrixFromListboxString(listboxstring,1,selectedWell,handles);
listboxindex = get(handles.popupMinus3, 'value');
listboxstring = listboxstrings.String(listboxindex); 
MatrixMinus3 = GetMatrixFromListboxString(listboxstring,1,selectedWell,handles);
foldername = imageHandler.Foldername;
[sizeY sizeX] = size(imageHandler.NucleusImage);
if(MatrixPlus2 ~= -1)
        %MatrixPlus2 = CutOutCircles(MatrixPlus2,selectedWell,1, 0,0,optionHandler, foldername, sizeY, sizeX, NucleusM, NeuronM, imageHandler);    
        if(imageHandler.ZoomState)
            MatrixPlus2 = MatrixPlus2(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));        
        end
end
if(MatrixPlus1 ~= -1)
    %MatrixPlus1 = CutOutCircles(MatrixPlus1,selectedWell,1, 0,0,optionHandler, foldername, sizeY, sizeX, NucleusM, NeuronM, imageHandler);    
    if(imageHandler.ZoomState)
        MatrixPlus1 = MatrixPlus1(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));    
    end
    if(MatrixMinus1 ~= -1)
        if(imageHandler.ZoomState)
            MatrixMinus1 = MatrixMinus1(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
        end
        MatrixPlus1 = MatrixPlus1 - MatrixMinus1;
        ind = MatrixPlus1<0;
        MatrixPlus1(ind) = 0;
    end
    if(MatrixMinus2 ~= -1)
        if(imageHandler.ZoomState)
            MatrixMinus2 = MatrixMinus2(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
        end
        MatrixPlus1 = MatrixPlus1 - MatrixMinus2;
        ind = MatrixPlus1<0;
        MatrixPlus1(ind) = 0;
    end
    if(MatrixMinus3 ~= -1)
        if(imageHandler.ZoomState)
            MatrixMinus3 = MatrixMinus3(zoomOnZoomedImage(2):zoomOnZoomedImage(2)+zoomOnZoomedImage(4), zoomOnZoomedImage(1):zoomOnZoomedImage(1)+zoomOnZoomedImage(3));
        end
        MatrixPlus2 = MatrixPlus2 - MatrixMinus3;
        ind = MatrixPlus2<0;
        MatrixPlus2(ind) = 0;
    end
    [nucleusRows, nucleusCols] = find(MatrixPlus1);   
     hold on;
     if(imageHandler.ZoomState)
         plot((nucleusCols),(nucleusRows),'Linestyle','none','Marker','.','Markersize',15,'Color','red');       
     else
        plot((nucleusCols./(1/optionHandler.ScalingFactor)),(nucleusRows./(1/optionHandler.ScalingFactor)),'Linestyle','none','Marker','.','Markersize',10,'Color','red');
     end     
     hold off;
end
if(MatrixPlus2 ~= -1)    
    [nucleusRows, nucleusCols] = find(MatrixPlus2);   
     hold on;
     if(imageHandler.ZoomState)
        plot((nucleusCols),(nucleusRows),'Linestyle','none','Marker','.','Markersize',15,'Color','blue');       
     else
        plot((nucleusCols./(1/optionHandler.ScalingFactor)),(nucleusRows./(1/optionHandler.ScalingFactor)),'Linestyle','none','Marker','.','Markersize',10,'Color','blue');
		%plot((nucleusCols./10),(nucleusRows./10),'Linestyle','none','Marker','.','Markersize',10,'Color','blue');
     end
     hold off;
end


function mat=GetMatrixFromListboxString(listboxstring,convert,selectedWell,handles)
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
if(convert==1)
    listboxstring=listboxstring{1};
end
if(isfield(handles,'NeuronCoordinates'))
    neuronHandler = handles.NeuronCoordinates;
end
if(isfield(handles,'NeuronCoordinates2'))
    neuronHandler2 = handles.NeuronCoordinates2;
end
if(isfield(handles,'OligoCoordinates'))
    oligoHandler = handles.OligoCoordinates;
end
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
    case 'Composite Fill Neurons'
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
    case 'Neuronal Tracing'
        del = neuronHandler.NeuronPositionsSkelDeleted(selectedWell);
        mat = neuronHandler.NeuronPositionsSkeletonization(selectedWell);
        mat=mat-del;
        ind = mat<0;
        mat(ind) = 0;
    case 'Composite Fill Oligos'
         mat = oligoHandler.NeuronPositionsEdgeFill;
        if(mat ~= -1)
            mat = mat(selectedWell);
        end
    case 'Oligo Tracing'
        del = oligoHandler.NeuronPositionsSkelDeleted(selectedWell);
        mat = oligoHandler.NeuronPositionsSkeletonization(selectedWell);
        mat=mat-del;
        ind = mat<0;
        mat(ind) = 0;
    otherwise
        mat=-1;        
end

function SetMatrixFromListboxString(listboxstring,srcmat,selectedWell,handles)
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
if(isfield(handles,'NeuronCoordinates'))
    neuronHandler = handles.NeuronCoordinates;
end
if(isfield(handles,'NeuronCoordinates2'))
    neuronHandler2 = handles.NeuronCoordinates2;
end
if(isfield(handles,'OligoCoordinates'))
    oligoHandler = handles.OligoCoordinates;
end
listboxstring=listboxstring{1};
if(numel(csvHandler.ManualPositions1) ==0)
     csvHandler.ManualPositions1 = containers.Map();
     csvHandler.ManualPositions2 = containers.Map();
     csvHandler.ManualPositions3 = containers.Map();
     csvHandler.ManualPositions4 = containers.Map();
end


switch listboxstring
    case 'Nucleus Matrix'
        csvHandler.CellPosMatrix(selectedWell) = srcmat;         
    case 'Cellomics Neurons'
        neuronHandler.CellPosMatrix(selectedWell) = srcmat;
    case 'Cellomics Neurons 2'
        neuronHandler2.CellPosMatrix(selectedWell) = srcmat;
    case 'Manual Neurons Main'
        neuronHandler.ManualNeuronPositionsSparse(selectedWell) = srcmat;
    case 'Manual Neurons 1'
        csvHandler.ManualPositions1(selectedWell)=srcmat;        
    case 'Manual Neurons 2'
        csvHandler.ManualPositions2(selectedWell)=srcmat;        
    case 'Manual Neurons 3'
        csvHandler.ManualPositions3(selectedWell)=srcmat;        
    case 'Manual Neurons 4'
        csvHandler.ManualPositions4(selectedWell)=srcmat;        
    case 'Edge Composit Neurons'
        neuronHandler.NeuronPositionsEdgeComposite(selectedWell)=srcmat;        
    case 'Composite Fill Neurons'
        neuronHandler.NeuronPositionsEdgeFill(selectedWell)=srcmat;        
    case 'Skeleton Neurons'        
        neuronHandler.NeuronPositionsSkeletonization(selectedWell)=srcmat;        
    case 'Neuronal Tracing'
        neuronHandler.NeuronPositionsSkeletonization(selectedWell)=srcmat;
    case 'Composite Fill Oligos'
        oligoHandler.NeuronPositionsEdgeFill(selectedWell)=srcmat;
    case 'Oligo Tracing'
        oligoHandler.NeuronPositionsSkeletonization(selectedWell)=srcmat;
    otherwise
        srcmat=-1;        
end
guidata(handles.figure1, handles);
    
function RefreshUnzoomedImage(handles)
optionHandler = Option();
ScalingFactor = optionHandler.ScalingFactor;
RefreshZoomedImage(handles);
imageHandler = handles.ImageHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
path = [imageHandler.Foldername '/ConvertedCellomics/' selectedWell '.mat'];
%Load filter for preview
if(exist(path,'file'))
    load(path);
    imageHandler.PositiveFilters = FMask.PositiveFilters;
    imageHandler.NegativeFilters = FMask.NegativeFilters;
    filterEntries = numel(imageHandler.PositiveFilters);
    currentFilterList = imageHandler.PositiveFilters;    
    for(i=1:numel(currentFilterList))
        currentFilter=currentFilterList{i};
        currentFilterY = transpose(currentFilter);
        [yOfX xStart] = find(currentFilter,1);
        [yOfx xStop] = find(currentFilter,1,'last');
        [xOfY yStart ] = find(currentFilterY,1);
        [xOfy yStop ] = find(currentFilterY,1,'last');
        if(xStop<xStart)
            xDummy = xStart;
            xStart = xStop;
            xStop = xDummy;
        end
        if(yStop<yStart)
            yDummy = yStart;
            yStart = yStop;
            yStop = yDummy;
        end
        width = xStop - xStart;
        height = yStop - yStart;
        
        
        h = imrect(gca,[xStart/(1/ScalingFactor),yStart/(1/ScalingFactor),width/(1/ScalingFactor), height/(1/ScalingFactor)]); 
        f= 1;
    end
end


function [DensityM] = getDensityDistributionsFromCSV(NucleusM,sizeXTarget,sizeYTarget)
    [sizeY, sizeX] = size(NucleusM);
    [row col] = find(NucleusM);
    row = [row;sizeY];
    col = [col;sizeX];
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


function [filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates result resultCompFill AnalyzeHandler] = calculateDensityDistribution(selectedWell, SphereAreaSizeX, SphereAreaSizeY, AnalyzeHandler, handles)
markerPointCoordinates=-1;
filterDistance = -1;
filterDistanceCompFill = -1;
nonFilterDistance=-1;
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
neuronHandler = AnalyzeHandler;
[sizeY sizeX] = size(imageHandler.NucleusImage);
SphereArea=-1;
optionHandler = handles.OptionHandler;
ringNumber = optionHandler.DensityDistributionRingNumber;
NucleusM = csvHandler.CellPosMatrix(selectedWell);
[filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, optionHandler, NucleusM, imageHandler.Foldername, sizeY, sizeX); 

%Migration Distance calculation gives distance from Sphere and Nucleus
%We need: LinePic
if(filterDistance == -1)
   result = zeros(ringNumber,9);
   resultCompFill = zeros(ringNumber,9);
elseif (markerPointCoordinates==0)
   result = zeros(ringNumber,9);
   resultCompFill = zeros(ringNumber,9);
elseif(numel(markerPointCoordinates('0')) == 0)
    result = zeros(ringNumber,9);
    resultCompFill = zeros(ringNumber,9);
else

%NucleusM = csvHandler.CellPosMatrix;
%NeuronM = neuronHandler.CellPosMatrix;

if(~isa(neuronHandler.CellPosMatrix,'containers.Map'))
       neuronHandler.CellPosMatrix = containers.Map();
end


if(~isKey(neuronHandler.CellPosMatrix,selectedWell))
    NeuronM = sparse(logical(zeros(sizeY,sizeX)));
    neuronHandler.CellPosMatrix(selectedWell) = NeuronM;
else
    NeuronM = neuronHandler.CellPosMatrix(selectedWell);
end
if(numel(neuronHandler.ManualNeuronPositionsSparse) ==0) 
        neuronHandler.ManualNeuronPositionsSparse = containers.Map()
end
if(~isKey(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = sparse(sizeY, sizeX);
end
NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);

if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
if(~isKey(neuronHandler.NeuronPositionsSkeletonization,selectedWell))
    neuronHandler.NeuronPositionsSkeletonization(selectedWell) = sparse(sizeY, sizeX);
end
NeuronSkelM = neuronHandler.NeuronPositionsSkeletonization(selectedWell);
NeuronCompFillM = neuronHandler.NeuronPositionsEdgeFill(selectedWell);
if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
    neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
end
if(~isKey(neuronHandler.NeuronPositionsSkelDeleted,selectedWell))
    neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = sparse(sizeY, sizeX);
end
NeuronsDeleted = neuronHandler.NeuronPositionsSkelDeleted(selectedWell);

    NeuronSkelM = NeuronSkelM - (NeuronsDeleted);
    ind = NeuronSkelM<0;
    NeuronSkelM(ind) = 0;



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
hInner = impoly(handles.axes2,double(markerPointCoordinates('0')./(1/ScalingFactor)));
previousMask = createMask(hInner);
cbRings = get(handles.cbPolys ,'Value');
%if(~cbRings)
    delete(hInner);
%end
    for i=1:ringNumber
    %Create masks from LinePic and cut them from CSVs
    %hTest = impoly(handles.axes2, [1 1; 500 500]);   
        hCurrent = impoly(handles.axes2,double(markerPointCoordinates(num2str(i*10))./(1/ScalingFactor)));    
        currentMask = logical(createMask(hCurrent));
        currentRing = logical(currentMask-previousMask); 
        currentRing = imresize(currentRing, [sizeY, sizeX]);
        previousMask = currentMask;
        %Save hCurrent to file
     %   if(~cbRings)
            delete(hCurrent);
     %   end
    %Filter each ring on CSVs to calculate Density            
        Nucleus2 = logical(logical(NucleusM) .* currentRing);
        NeuronCell2 = logical(logical(NeuronM) .* currentRing);
        NeuronManual2 = logical(logical(NeuronManualM) .* currentRing);
        NeuronSkel2 = logical(logical(NeuronSkelM) .* currentRing);
        NeuronCompFill2 = logical(logical(NeuronCompFillM) .* currentRing);
        result(i,1) = nnz(currentRing);
        result(i,2) = nnz(Nucleus2);
        result(i,3) = nnz(NeuronCell2);
        result(i,4) = nnz(NeuronManual2);
        result(i,5) = nnz(NeuronSkel2);
        result(i,6) = result(i,3) / result(i,2);
        result(i,7) = result(i,4) / result(i,2);
        result(i,8) = result(i,5) / result(i,2);
        result(i,9) = result(i,2) / result(i,1);
        
        resultCompFill(i,1) = nnz(currentRing);
        resultCompFill(i,2) = nnz(Nucleus2);
        resultCompFill(i,3) = nnz(NeuronCell2);
        resultCompFill(i,4) = nnz(NeuronManual2);
        resultCompFill(i,5) = nnz(NeuronCompFill2);
        resultCompFill(i,6) = resultCompFill(i,3) / resultCompFill(i,2);
        resultCompFill(i,7) = resultCompFill(i,4) / resultCompFill(i,2);
        resultCompFill(i,8) = resultCompFill(i,5) / resultCompFill(i,2);
        resultCompFill(i,9) = resultCompFill(i,2) / resultCompFill(i,1);
    end
end


function Imx = keepMaxObj(X)
    %Error checking:
    assert(islogical(X)||isstruct(X),'The first input argument is expected to be a struct or a logical');
    if isstruct(X)
        CC = X;
        parts = {'PixelIdxList','ImageSize'};
        assert(all(ismember(parts,fieldnames(CC))),'CC is expected to be the output from bwconncomp');
    else
        CC = bwconncomp(X);
    end  
    clear X;
    %Preallocate and find number of voxels/object
    Nvox = zeros(CC.NumObjects,1);
    for ii = 1:CC.NumObjects
        Nvox(ii) = numel(CC.PixelIdxList{ii});
    end
    %Find the biggest object's index, warn and save all if there are multiples
    [mx,midx] = max(Nvox);
    more_than1_max = sum(mx==Nvox);
    if more_than1_max > 1
        midx = find(mx == Nvox);
        warning('Multiple:Maxima', 'There were %i objects with the maximum size.\n  They are all left on!',more_than1_max);
    end    
    %Create the final image
    Imx = false(CC.ImageSize);
    Imx([CC.PixelIdxList{midx}]) = true;

function [filterDistance nonFilterDistance sphereRadiusFilter SphereArea markerPointCoordinates] = calculateMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, optionHandler, NucleusM, foldername, sizeY, sizeX)
ringNumber = optionHandler.DensityDistributionRingNumber;
SphereArea=0;
NucleusArea=0;
filterDistance=-1;
nonFilterDistance=-1;
sphereRadiusFilter=-1;
markerPointCoordinates = 0;
NucleusArea64=0;
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
      selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
elseif(length(selectedWell) == 4)
        selectedWellLong = [selectedWell(1) '0' selectedWell(2:length(selectedWell))];
end

    
    %Get images from image folder
imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
imagePathNucleusSmall = [foldername '/' selectedWellLong 'NucleusSmall' optionHandler.FileEnding];
imagePathNeuriteBig = [foldername '/' selectedWellLong 'NeuriteBig' optionHandler.FileEnding];
imagePathNeuriteSmall = [foldername '/' selectedWellLong 'NeuriteSmall' optionHandler.FileEnding];
if(exist(imagePathNucleusSmall,'file'))
    ResizedNeuriteImage = imread(imagePathNeuriteSmall);
    ResizedNucleusImage = imread(imagePathNucleusSmall);    
    NeuriteImage = imread(imagePathNeuriteBig);
    NucleusImage = imread(imagePathNucleusBig);
    selectedWellLong=selectedWell;
    if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    elseif(length(selectedWell) == 4)
        selectedWellLong = [selectedWell(1) '0' selectedWell(2:length(selectedWell))];
    end   

    [DensityM] = getDensityDistributionsFromCSV(NucleusM,32,32);
    %[DensityNeuron] = getDensityDistributionsFromCSV(NeuronM,32,32);
    
    %Find circle in Density Edge
    %Get point of Maximum Density
    [C,rowmaxarray]=max(DensityM);
    [maxvalue,xMaxSphere]=max(C);
    yMaxSphere=rowmaxarray(xMaxSphere);
    %densityNeurons = DensityNeuron(uint8(yMaxSphere/2), uint8(xMaxSphere/2));
    i=0;    
     thresholdNucImage = optionHandler.MigDistLowerNucleusThreshold;
     %New Calculation method for Nucleus Area
     [thresholdedRows thresholdedCols] = find(NucleusImage > thresholdNucImage);
     thresholdedRows = [thresholdedRows;sizeY];
     thresholdedCols = [thresholdedCols;sizeX];
     thresholdedRows = [thresholdedRows;0];
     thresholdedCols = [thresholdedCols;0];
     
          
     %1. Cut out border
     %Only corners and not so much!     
     NucleusArea = NucleusImage(401:sizeY-400, 401:sizeX-400);
     NucleusArea = logical(im2bw(NucleusArea, thresholdNucImage/255));
     NucleusArea = logical(keepMaxObj(NucleusArea));
     vertical = logical(zeros(sizeY-(2*400),400));
     horizontal = logical(zeros((400),sizeX));
     NucleusArea = logical([vertical NucleusArea vertical]);
     NucleusArea = logical([horizontal; NucleusArea; horizontal]);
     %2. Set corners black!
     %Upper left
     NucleusArea(1:1024,1:1024)=0;
     %Upper right
     NucleusArea(1:1024,sizeX-1024:sizeX)=0;
     %Lower left
     NucleusArea(sizeY-1024:sizeY,1:1024)=0;
     %Lower right
     NucleusArea(sizeY-1024:sizeY,sizeX-1024:sizeX)=0;
     

    CuttedIndicesNuc=0;
    NucleusAreaTemp=0;
    
     [row,col] = find(NucleusArea);
     if(length(row) > 0 && length(col) > 0)
       
         
       %sx and sy have to be on area.
       %If not get next area point and set it there.
       sYMapped = sum(row)/length(row);
       sXMapped = sum(col)/length(col);
       sYMappedorig = sum(row)/length(row);
       sXMappedorig = sum(col)/length(col);
       [sXorig sYorig] = MapPoint(sXMapped,sYMapped,SphereAreaSizeX,SphereAreaSizeY,sizeX,sizeY,1);
       if(NucleusArea(uint16(sYMapped),uint16(sXMapped)) == 0)
            [D IDX] = bwdist(NucleusArea);
            [sYMapped sXMapped] = ind2sub([sizeY sizeX],IDX(uint16(sYMapped),uint16(sXMapped)));
       end
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
       [DensityMSized] = getDensityDistributionsFromCSV(NucleusMSorted,SphereAreaSizeX, SphereAreaSizeY);
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
        xMaxSphereSmall = (xMaxSphere/SphereAreaSizeX)*32;
        yMaxSphereSmall = (yMaxSphere/SphereAreaSizeX)*32;
        yMaxSphereSmall=uint8(yMaxSphereSmall);
        xMaxSphereSmall=uint8(xMaxSphereSmall);
        if(yMaxSphereSmall < 1)
            yMaxSphereSmall=1;
        end
        if(xMaxSphereSmall < 1)
            xMaxSphereSmall=1;
        end
        %densityNeurons = DensityNeuron(uint8(yMaxSphereSmall), uint8(xMaxSphereSmall));
        i=i+1;
        stack.push([xMaxSphere yMaxSphere]);
      end
      stack.push([sX sY]);
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
%F�r weitere Verarbeitung:
% Berechne den Schwerpunkt aller wei�en Pixel in NucleusArea
     
      minimizingFactorX = 32/SphereAreaSizeX;
      minimizingFactorY = 32/SphereAreaSizeY;
      mediumDistanceList = zeros(0,0);
      mediumDistanceListInInterval = zeros(0);
      mediumDistanceListRadius = zeros(0,0);
      mediumDistanceListRadiusInInterval = zeros(0);
     
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
        x=sXorig;
        y=sYorig;
        xBig=sXMappedorig;
        yBig=sYMappedorig;
        if(uint8(x)<1 | uint8(y)<1)
            continue;
        end
        SphereArea(uint8(y),uint8(x)) =0;
        steps=0;
        while((1==1 && (sXorig==sX && sYorig == sY) || (sXorig~=sX || sYorig~=sY)) && uint8(x) <= SphereAreaSizeX && uint8(y) <= SphereAreaSizeY && uint8(x) > 0 && uint8(y) > 0 && uint8(xBig) > 0 && uint8(yBig) >0 && uint16(xBig) < sizeX && uint16(yBig) < sizeY)
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
          mediumDistanceListRadius = [mediumDistanceListRadius; pdist([xOuterSphere yOuterSphere;sXMapped sYMapped]) successCount];
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
            mediumDistanceListRadius = [mediumDistanceList; pdist([xOuterSphere yOuterSphere;sXMapped sYMapped]) successCount];  
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
      currentDistanceRadius = mediumDistanceListRadius(i,1);
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
      threshold = (maxDistance/3);
      %threshold = (maxDistance-(2*distancestddev));
    %Check if current distance is in Interval. Check also if one of both
    %neighbourpoints is also in interval
      %if(currentDistance >= 150 && currentDistance >= (maxDistance-(2.0*distancestddev)))
      if(currentDistance >= threshold && (leftNeighbourDistance > threshold || rightNeighbourDistance > threshold))
        mediumDistanceListInInterval = [mediumDistanceListInInterval currentDistance];
        mediumDistanceListRadiusInInterval = [mediumDistanceListRadiusInInterval currentDistanceRadius];
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
     sphereRadiusFilter = num2str(mean(mediumDistanceListRadiusInInterval), '%f');
     disp('Medium distance total (in Pixel): ');
     disp(mediumDistance);
     disp('Medium distance selected (in Pixel): ');
     disp(mediumDistanceSelected);
     disp('Sphere Radius (in Pixel): ');
     disp(sphereRadiusFilter);
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
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
% Hint: get(hObject,'Value') returns toggle state of cbNucleusPicture


% --- Executes on button press in cbNeuritePicture.
function cbNeuritePicture_Callback(hObject, eventdata, handles)
% hObject    handle to cbNeuritePicture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
% Hint: get(hObject,'Value') returns toggle state of cbNeuritePicture


% --- Executes on button press in cbNucleusMatrix.
function cbNucleusMatrix_Callback(hObject, eventdata, handles)
% hObject    handle to cbNucleusMatrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
% Hint: get(hObject,'Value') returns toggle state of cbNucleusMatrix


% --- Executes on button press in cbCellomicsNeuronMatrix.
function cbCellomicsNeuronMatrix_Callback(hObject, eventdata, handles)
% hObject    handle to cbCellomicsNeuronMatrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
% Hint: get(hObject,'Value') returns toggle state of cbCellomicsNeuronMatrix


% --- Executes on button press in cbManualNeuronMatrix.
function cbManualNeuronMatrix_Callback(hObject, eventdata, handles)
% hObject    handle to cbManualNeuronMatrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
% Hint: get(hObject,'Value') returns toggle state of cbManualNeuronMatrix


function DensDistAll(FileName,PathName,handles)
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;

fileID = fopen(strcat(PathName,'/',FileName),'w');
handles = guidata(handles.figure1);
ringCount = optionHandler.DensityDistributionRingNumber;
densityWidth = optionHandler.MigrationDistanceDensityImageXSize;
densityHeight = optionHandler.MigrationDistanceDensityImageYSize
wellList = get(handles.lbWell, 'string');
waitbarHandle = waitbar(0,'Calculating and exporting Migration Distances for all wells', 'WindowStyle', 'modal');
selectedWellNumber = get(handles.lbWell,'Value');
exportText = sprintf('Well;Ring Number;Migration Distance(px);Manual measured Migration Distance (px);Migration Distance without filter (px);Number of Pixels;Number Nuclei;Cellomics Neurons;Manual Neurons;Skeleton Neurons;Cellomics Neurons per Nuclei;Manual Neurons per Nuclei;Skeleton Neurons per Nuclei;Nuclei per Pixels;\r\n');
for currentWellIndex=optionHandler.StartWell:numel(wellList)
    currentWell = wellList(currentWellIndex);
    currentWell = currentWell{1};
    disp(currentWell);    
    [distFilter, distWithoutFilter, SphereRadius, SphereArea, markerPointCoordinates densDist densDistCompFill resultHandler]=calculateDensityDistribution(currentWell,densityWidth,densityHeight,handles.NeuronCoordinates,handles);
    handles.NeuronCoordinates = resultHandler;
    if(numel(densDist) == 1)
        %Fill up densDist with dummy values
        densDist = zeros(ringCount,9);        
    end
    if(distFilter >-1)
        distFilter = num2str(distFilter,'%f');
        distWithoutFilter = num2str(distWithoutFilter,'%f');
        for i=1:optionHandler.DensityDistributionRingNumber        
          if(isKey(imageHandler.MigDistWellMapping,currentWell))
            manualDist = num2str(imageHandler.MigDistWellMapping(currentWell), '%f')
            manualDist = strrep(manualDist,'.',',');
            exportText = [exportText currentWell ';' num2str(i) ';' distFilter ';' manualDist ';' distWithoutFilter ';' num2str(densDist(i,1)),'.',',' ';' num2str(densDist(i,2)) ';' num2str(densDist(i,3)) ';' num2str(densDist(i,4)) ';' num2str(densDist(i,5)) ';' num2str(densDist(i,6)) ';' num2str(densDist(i,7)) ';' num2str(densDist(i,8)) ';' num2str(densDist(i,9)) sprintf('\r\n')];
          else
            exportText = [exportText currentWell ';' num2str(i) ';' distFilter ';;' distWithoutFilter ';' num2str(densDist(i,1)) ';' num2str(densDist(i,2)) ';' num2str(densDist(i,3)) ';' num2str(densDist(i,4)) ';' num2str(densDist(i,5)) ';' num2str(densDist(i,6)) ';' num2str(densDist(i,7)) ';' num2str(densDist(i,8)) ';' num2str(densDist(i,9)) sprintf('\r\n')];
          end
        end
    else
        for i=1:optionHandler.DensityDistributionRingNumber
            exportText = [exportText currentWell ';' num2str(i) ';' num2str(0) ';;' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(0) sprintf('\r\n')];
        end
    end
    waitbar((currentWellIndex/numel(wellList)),waitbarHandle);
end
%Write exportText to File
fprintf(fileID,'%s',exportText);
close(waitbarHandle);

% --------------------------------------------------------------------
function MigrationDistanceAll_Callback(hObject, eventdata, handles)
% hObject    handle to MigrationDistanceAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Outer Loop: Iterate over all Wells.
try
[FileName,PathName] = uiputfile('migdist.csv');
DensDistAll(FileName,PathName,handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

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
imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
imagePathNucleusSmall = [foldername '/' selectedWellLong 'NucleusSmall' optionHandler.FileEnding];
imagePathNeuriteBig = [foldername '/' selectedWellLong 'NeuriteBig' optionHandler.FileEnding];
imagePathNeuriteSmall = [foldername '/' selectedWellLong 'NeuriteSmall' optionHandler.FileEnding];

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


[filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, densityWidth, densityHeight, optionHandler, NucleusM, imageHandler.Foldername, sizeY, sizeX); 
%MarkerPointCoordinates are available. Get second ring:
if(markerPointCoordinates ~= 0)
    hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./(1/ScalingFactor)));
    innerMask = logical(createMask(hInner));
    hOuter = impoly(handles.axes2,double(markerPointCoordinates(num2str(ringNumber * 10))./(1/ScalingFactor)));
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
        if(startY-50>0 && startX-50>0 && startX+50 < sizeXsmall && startY + 50 < sizeYsmall)
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
            [xNuc,yNuc,distance] = csvHandler.FindNucleusForNeuron(xFocus,yFocus,NucleusM,75,0,0);
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
try
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
EdgeCompositAlgorithm(selectedWell,handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --------------------------------------------------------------------
function StatisticsCSVExport_Callback(hObject, eventdata, handles)
    try
% hObject    handle to StatisticsCSVExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Extract Rectangles from ImageHandler. 
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
[FileName,PathName] = uiputfile('WellStatistics.csv');
fileID = fopen(strcat(PathName,'/',FileName),'w');
waitbarHandle = waitbar(0,'Calculating and exporting Statistics for all wells', 'WindowStyle', 'modal');
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
[sizeY sizeX]=size(imageHandler.NucleusImage);
%Preprocessing: Change existing filters to Sparse filters
%Iterate over all wells

exportText = sprintf('Well;Filter;#Nuclei;#Neurons Cellomics;#Manual Main;#Manual 1;#Manual 2;#Manual 3;#Manual 4;%%Neurons by Cellomics;%%Neurons by Clicks;%%Neurons by Neuronal Tracing;#Cellomics & Manuell;#Additional Manual;#Additional Cellomics;#Composite Fill & Manual;#Additional Manual;#Additional Composite Fill;#Neuronal Tracing & Manual;#Additional Manual;#Additional Neuronal Tracing\r\n');
for i=1:numel(wellList)    
    selectedWell = wellList{i};
    disp(strcat('Current well: ',selectedWell)); 
    manualDist = '';
    if(isKey(imageHandler.MigDistWellMapping,selectedWell))
            manualDist = num2str(imageHandler.MigDistWellMapping(currentWell), '%f')
            manualDist = strrep(manualDist,'.',',');
            exportText = [exportText currentWell ';' distFilter ';' manualDist ';' distWithoutFilter sprintf('\r\n')];
    end
     NucleusM = csvHandler.CellPosMatrix(selectedWell);
     foldername = imageHandler.Foldername;
     [sizeY sizeX] = size(imageHandler.NeuriteImage);
    [distFilter, distWithoutFilter, sphereRadius, SphereArea, NucleusArea]=calculateMigrationDistance(selectedWell,64,64,optionHandler, NucleusM, foldername, sizeY, sizeX);
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
    filterList = imageHandler.PositiveFilters;
    for j=1:numel(filterList)        
        if(iscell(filterList))
            filterMask = filterList(j);
            filterMask = sparse(filterMask{1});
            filterList{j}=filterMask;
        end        
    end
    imageHandler.PositiveFilters=filterList;
    %Iterate over all filters
    filterList = imageHandler.PositiveFilters;    
    for j=1:numel(filterList)
        if(iscell(filterList))
            filterMask = filterList(j);
            filterMask = sparse(filterMask{1});
            
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
           	if(numel(csvHandler.ManualPositions1) ==0)        
                csvHandler.ManualPositions1 = containers.Map();
                csvHandler.ManualPositions2 = containers.Map();
                csvHandler.ManualPositions3 = containers.Map();
                csvHandler.ManualPositions4 = containers.Map();        
            end    
            if(~isKey(csvHandler.ManualPositions1,selectedWell))                         
                currentNuclei1 = logical(sparse(sizeY,sizeX));               
                csvHandler.ManualPositions1(selectedWell) = currentNuclei1;
            end  
            if(~isKey(csvHandler.ManualPositions2,selectedWell))                         
                currentNuclei2 = logical(sparse(sizeY,sizeX));               
                csvHandler.ManualPositions2(selectedWell) = currentNuclei2;
            end  
            if(~isKey(csvHandler.ManualPositions3,selectedWell))                         
                currentNuclei3 = logical(sparse(sizeY,sizeX));               
                csvHandler.ManualPositions3(selectedWell) = currentNuclei3;
            end  
            if(~isKey(csvHandler.ManualPositions4,selectedWell))                         
                currentNuclei4 = logical(sparse(sizeY,sizeX));               
                csvHandler.ManualPositions4(selectedWell) = currentNuclei4;
            end  
            if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
                neuronHandler.NeuronPositionsSkeletonization = containers.Map();
            end
            if(~isKey(neuronHandler.NeuronPositionsSkeletonization,selectedWell))  
                neuronHandler.NeuronPositionsSkeletonization(selectedWell) = logical(sparse(sizeY,sizeX));
            end
            NucleusM = csvHandler.CellPosMatrix(selectedWell);
            
            if(~isa(neuronHandler.CellPosMatrix,'containers.Map'))
                neuronHandler.CellPosMatrix = containers.Map();
            end
            if(isKey(neuronHandler.CellPosMatrix,selectedWell))
                NeuronM = neuronHandler.CellPosMatrix(selectedWell);
            else
                NeuronM = logical(zeros(sizeY,sizeX)); 
            end
            if(~isa(neuronHandler.NeuronPositionsEdgeFill,'containers.Map'))
                neuronHandler.NeuronPositionsEdgeFill = containers.Map();
            end
            %NeuronM = neuronHandler.CellPosMatrix(selectedWell);
            if(isKey(neuronHandler.NeuronPositionsEdgeFill,selectedWell))
                CompositeFillM = neuronHandler.NeuronPositionsEdgeFill(selectedWell);
            else
                CompositeFillM = logical(zeros(sizeY,sizeX)); 
            end
            if(~isa(neuronHandler.ManualNeuronPositionsSparse,'containers.Map'))
                neuronHandler.ManualNeuronPositionsSparse = containers.Map();
            end
            if(isKey(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
                NeuronManualM = logical(neuronHandler.ManualNeuronPositionsSparse(selectedWell));
                NeuronManualM1 = logical(csvHandler.ManualPositions1(selectedWell));
                NeuronManualM2 = logical(csvHandler.ManualPositions2(selectedWell));
                NeuronManualM3 = logical(csvHandler.ManualPositions3(selectedWell));
                NeuronManualM4 = logical(csvHandler.ManualPositions4(selectedWell));
            else
                NeuronManualM = sparse(logical(zeros(sizeY,sizeX))); 
                NeuronManualM1 = sparse(logical(zeros(sizeY,sizeX))); 
                NeuronManualM2 = sparse(logical(zeros(sizeY,sizeX))); 
                NeuronManualM3 = sparse(logical(zeros(sizeY,sizeX))); 
                NeuronManualM4 = sparse(logical(zeros(sizeY,sizeX))); 
            end
            
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
            NeuronM = logical(logical(NeuronM) .* filterMask);               
            numberNeurons = nnz(NeuronM);
            NeuronManualM = logical(logical(full(NeuronManualM)) .* filterMask);
            NeuronManualM1 = logical(logical(full(NeuronManualM1)) .* filterMask);
            NeuronManualM2 = logical(logical(full(NeuronManualM2)) .* filterMask);
            NeuronManualM3 = logical(logical(full(NeuronManualM3)) .* filterMask);
            NeuronManualM4 = logical(logical(full(NeuronManualM4)) .* filterMask);
            CompositeFillM = logical(logical(full(CompositeFillM)) .* filterMask); 
            numberManual = nnz(NeuronManualM);
            numberManual1 = nnz(NeuronManualM1);
            numberManual2 = nnz(NeuronManualM2);
            numberManual3 = nnz(NeuronManualM3);
            numberManual4 = nnz(NeuronManualM4);
            numberCompositeFill = nnz(CompositeFillM);
            SkeletonNeurons = logical(logical(neuronHandler.NeuronPositionsSkeletonization(selectedWell)) .* filterMask);
            filterMask=0;
           
%disp('Neurons Manual;Neurons Cellomics;Neurons Edge Composite;Positions Manual & Edge;Additional Manual; Additional Edge Composite;;Positions Edge & Cellomics;Additional Cellomics; Additional Edge Composite;; Positions Manual & Cellomics; Additional Manual; Additional Cellomics;; Positions Manual & Edge Filled; Additioal Manual; Additional Edge Filled')

bothManualCellomics = nnz(NeuronManualM .* NeuronM);
bothManualSkeletonized = nnz(SkeletonNeurons .* NeuronManualM);
bothManualComposite = nnz(CompositeFillM .* NeuronManualM);
CellomicsNeuronsInvert =  logical(full(~NeuronM));



AdditionalManualCellomics = nnz(NeuronManualM .* CellomicsNeuronsInvert);
CellomicsNeuronsInvert=0;
ManualNeuronsInvert = logical(full(~NeuronManualM));
AdditionalCellomicsManual = nnz(NeuronM .* ManualNeuronsInvert);
AdditionalSkeletonManual = nnz(SkeletonNeurons .* ManualNeuronsInvert);
AdditionalCompositeManual = nnz(CompositeFillM .* ManualNeuronsInvert);
ManualNeuronsInvert=0;
SkeletonNeuronsInvert = logical(full(~SkeletonNeurons));
CompositeFillNeuronsInvert = logical(full(~CompositeFillM));
AdditionalManualSkeleton = nnz(NeuronManualM .* SkeletonNeuronsInvert);
AdditionalManualCompositeFill = nnz(NeuronManualM .* CompositeFillNeuronsInvert);
SkeletonNeuronsInvert=0;
EdgeCompNeurons=0;

AdditionalCellomics = 0;%;nnz(EdgeCompInvert .* CellomicsNeurons);
AdditionalManual = 0;%nnz(EdgeCompInvert .* ManualNeurons);
            
            
            percentageCellomics = (numberNeurons/numberNuclei) * 100;
            percentageManual = (numberManual/numberNuclei) * 100;
            percentageTracing = (nnz(SkeletonNeurons)/numberNuclei) * 100;
            numberNuclei = num2str(numberNuclei,'%f');
            %numberNuclei = strrep(numberNuclei, '.',',');
            numberNeurons = num2str(numberNeurons,'%f');
            %numberNeurons = strrep(numberNeurons, '.',',');
            numberManual = num2str(numberManual,'%f');
            %numberManual = strrep(numberManual, '.',',');
            numberManual1 = num2str(numberManual1,'%f');
            %numberManual1 = strrep(numberManual1, '.',',');
            numberManual2 = num2str(numberManual2,'%f');
            %numberManual2 = strrep(numberManual2, '.',',');
            numberManual3 = num2str(numberManual3,'%f');
            %numberManual3 = strrep(numberManual3, '.',',');
            numberManual4 = num2str(numberManual4,'%f');
            %numberManual4 = strrep(numberManual4, '.',',');
            percentageCellomics = num2str(percentageCellomics,'%f');
            %percentageCellomics = strrep(percentageCellomics, '.',',');
            percentageManual = num2str(percentageManual,'%f');
            %percentageManual = strrep(percentageManual, '.',',');
            percentageTracing = num2str(percentageTracing,'%f');
            %percentageTracing = strrep(percentageTracing, '.',',');
            bothManualCellomics = num2str(bothManualCellomics,'%f');
            %bothManualCellomics = strrep(bothManualCellomics, '.',',');
            AdditionalManualCellomics = num2str(AdditionalManualCellomics,'%f');
            %AdditionalManualCellomics = strrep(AdditionalManualCellomics, '.',',');
            AdditionalCellomicsManual = num2str(AdditionalCellomicsManual,'%f');
            %AdditionalCellomicsManual = strrep(AdditionalCellomicsManual, '.',',');
            bothManualSkeletonized = num2str(bothManualSkeletonized,'%f');
            %bothManualSkeletonized = strrep(bothManualSkeletonized, '.',',');
            AdditionalManualSkeleton = num2str(AdditionalManualSkeleton,'%f');
            %AdditionalManualSkeleton = strrep(AdditionalManualSkeleton, '.',',');
            AdditionalSkeletonManual = num2str(AdditionalSkeletonManual,'%f');
            %AdditionalSkeletonManual = strrep(AdditionalSkeletonManual, '.',',');
            
            %Add entry to CSV Exort string
            exportText = [exportText selectedWell ';' num2str(j) ';' numberNuclei ';' numberNeurons ';' numberManual ';' numberManual1 ';' numberManual2 ';' numberManual3 ';' numberManual4 ';' percentageCellomics ';' percentageManual ';' percentageTracing ';' bothManualCellomics ';' AdditionalManualCellomics ';' AdditionalCellomicsManual ';' num2str(bothManualComposite) ';' num2str(AdditionalManualCompositeFill) ';' num2str(AdditionalCompositeManual) ';' bothManualSkeletonized ';' AdditionalManualSkeleton ';' AdditionalSkeletonManual sprintf('\r\n')];
        end
    end
    waitbar((i/numel(wellList)),waitbarHandle);
end
%Write exportText to File
fprintf(fileID,'%s',exportText);
close(waitbarHandle);
handles.NeuronCoordinates = neuronHandler;
handles.CSVCoordinates = csvHandler;
guidata(handles.figure1, handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

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
try
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
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


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
try
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
densityWidth = optionHandler.MigrationDistanceDensityImageXSize;
densityHeight = optionHandler.MigrationDistanceDensityImageYSize;
[filterDistance nonFilterDistance SphereRadius SphereArea markerPointCoordinates result resultCompFill resultHandler] = calculateDensityDistribution(selectedWell, densityWidth, densityHeight, handles.NeuronCoordinates, handles);
handles.NeuronCoordinates = resultHandler;
if(numel(result) > 1)
    disp('Ring Number;Pixels in Ring;Nuclei in Ring;Neurons by Cellomics in Ring;Manual counted Neurons in Ring;Skeletonized Neurons in Ring; Cellomics Neurons per Nuclei;Manual Neurons per Nuclei;Skeleton Neurons per Nuclei; Nuclei per Pixels');
    for i=1:optionHandler.DensityDistributionRingNumber
        disp([num2str(i) ';' strrep(num2str(result(i,1)),'.',',') ';' strrep(num2str(result(i,2)),'.',',') ';' strrep(num2str(result(i,3)),'.',',') ';' strrep(num2str(result(i,4)),'.',',') ';' strrep(num2str(result(i,5)),'.',',') ';' strrep(num2str(result(i,6)),'.',',') ';' strrep(num2str(result(i,7)),'.',',') ';' strrep(num2str(result(i,8)),'.',',') ';' strrep(num2str(result(i,9)),'.',',')]);
    end
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
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

try
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
          hInner = impoly(handles.axes2,double(markerPointCoordinates('0')./(1/0.1)));
          previousMask = createMask(hInner);
          for i=10:ringNumber  
            hCurrent = impoly(handles.axes2,double(markerPointCoordinates(num2str(i*10))./(1/0.1)));
          end
        end
    end
else
    if(imageHandler.ZoomState)
        RefreshZoomedImage(handles);
    else
        RefreshUnzoomedImage(handles);
    end
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --------------------------------------------------------------------
function ConvertSingleChambTo96er_Callback(hObject, eventdata, handles)
% hObject    handle to ConvertSingleChambTo96er (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
optionHandler = handles.OptionHandler;
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
imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
imagePathNucleusSmall = [foldername '/' selectedWellLong 'NucleusSmall' optionHandler.FileEnding];
imagePathNeuriteBig = [foldername '/' selectedWellLong 'NeuriteBig' optionHandler.FileEnding];
imagePathNeuriteSmall = [foldername '/' selectedWellLong 'NeuriteSmall' optionHandler.FileEnding];   

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
    filterMaskSmall = imresize(filterMask,optionHandler.ScalingFactor);    
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
    nucImageSmall = imresize(nucImage,optionHandler.ScalingFactor);
    neuImageSmall = imresize(neuImage,optionHandler.ScalingFactor);
    %Save new pics to ConvertedCellomics Folder
    imagePathNucleusBigNew = [foldername '/' selectedWellLong '_' num2str(i) 'NucleusBig' optionHandler.FileEnding];
    imagePathNucleusSmallNew = [foldername '/' selectedWellLong '_' num2str(i) 'NucleusSmall' optionHandler.FileEnding];
    imagePathNeuriteBigNew = [foldername '/' selectedWellLong '_' num2str(i) 'NeuriteBig' optionHandler.FileEnding];
    imagePathNeuriteSmallNew = [foldername '/' selectedWellLong '_' num2str(i) 'NeuriteSmall' optionHandler.FileEnding]; 
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
ReloadWellList(handles,0);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --- Executes when selected object is changed in uipanel2.
function ReloadWellList(handles, method)
% hObject    handle to the selected object in uipanel2 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

%If csvHandler is available: Get all available Wells from CSVHandler in
%CellPosMatrix

handles = guidata(handles.figure1);
if(method==0)
    if(isfield(handles,'CSVCoordinates'))
               %Map Neuron Positions to Nucleus Positions
        wellList = get(handles.lbWell, 'string');
        csvHandler = handles.CSVCoordinates;
        %neuronHandler = handles.NeuronCoordinates;
        if(isprop(csvHandler, 'CellPosMatrix'))
            WellNeuronDict = csvHandler.CellPosMatrix;
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
else
    %Get latest well from ConvertedCellomicsFolder
    if(isfield(handles,'ImageHandler'))        
        imageHandler = handles.ImageHandler;
        foldername = [imageHandler.Foldername '/ConvertedCellomics/'];
        allFiles = dir(foldername); 
        fulltext = '';
        wellList = get(handles.lbWell, 'string');

        for j=1:numel(wellList)
            selectedWell = wellList(j);
            selectedWell = selectedWell{1};
            selectedWellLong=selectedWell;
            if(length(selectedWell) == 2)
                      selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
            end
            for(i=1:numel(allFiles))
                currentName=allFiles(i).name;
                c1Regexp = strcat(selectedWellLong,'NucleusSmall');
                %c1Regexp = strcat('ARRAY.*\d*_(',currentWell,')f(\d+)d0');
                tokensc1 = regexpi(currentName, c1Regexp, 'tokens');
                if(length(tokensc1) > 0)
                    fileWell = tokensc1{1};
                    fulltext = [fulltext;wellList(j)];
                end
            end
        end
        set(handles.lbWell, 'String', fulltext);
    end
end


% --- Executes on button press in cbEdgeCompositeNeurons.
function cbEdgeCompositeNeurons_Callback(hObject, eventdata, handles)
% hObject    handle to cbEdgeCompositeNeurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on button press in cbSkeletonNeurons.
function cbSkeletonNeurons_Callback(hObject, eventdata, handles)
% hObject    handle to cbSkeletonNeurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end



% --------------------------------------------------------------------
function QualityCheckOld_Callback(hObject, eventdata, handles)
% hObject    handle to QualityCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
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
     NucleusM = csvHandler.CellPosMatrix(selectedWell);
     foldername = imageHandler.Foldername;
     [sizeY sizeX] = size(imageHandler.NeuriteImage);
    [filterDistance nonFilterDistance sphereRadius SphereArea NucleusArea64 markerPointCoordinates] = calculateMigrationDistance(selectedWell, densityWidth, densityHeight, optionHandler, NucleusM, foldername, sizeY, sizeX);    
%end
%MarkerPointCoordinates are available. Get second ring:
hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./(1/ScalingFactor)));
innerMask = logical(createMask(hInner));
%hOuter = impoly(handles.axes2,double(markerPointCoordinates(num2str(ringNumber * 10))./10));
%outerMask = logical(createMask(hOuter));
%if(CutInnerCircle == 1 && CutOuterCircle == 1)
%    currentRing = logical(outerMask-innerMask); 
%elseif(CutInnerCircle==1 && CutOuterCircle == 0)
    currentRing =(logical(ones(int32(sizeY/(1/ScalingFactor)),int32(sizeX/(1/ScalingFactor))) - innerMask));
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
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function [DP FP TN] = CalculateQuality(foldername,neuronHandler,selectedWell,sizeY,sizeX,yStartPos,yEndPos,xStartPos,xEndPos,optionHandler,NucleusM)
majorMask = logical(zeros(sizeY,sizeX));
majorMask(yStartPos:yEndPos,xStartPos:xEndPos)=1;
SphereAreaSizeX = optionHandler.MigrationDistanceDensityImageXSize;
SphereAreaSizeY = optionHandler.MigrationDistanceDensityImageYSize;
ringNumber = optionHandler.DensityDistributionRingNumber;
%Get density distribution to exclude too dense areas
%Get density distribution to exclude too dense areas
%if(saveCircle)
if(optionHandler.CutOutSphereCore == 1)
  [filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, optionHandler, NucleusM, foldername, sizeY, sizeX); 
else
    markerPointCoordinates=0;
end
%MarkerPointCoordinates are available. Get second ring:
innerMask = ones(uint16(sizeY/(1/ScalingFactor)),uint16(sizeX/(1/ScalingFactor)));
if(markerPointCoordinates ~=0 && numel(markerPointCoordinates('10') > 0))
    mk = markerPointCoordinates('10')./(1/ScalingFactor);
    innerMask = roipoly(innerMask,mk(:,1),mk(:,2));
    outerMask = ones(uint16(sizeY/(1/ScalingFactor)),uint16(sizeX/(1/ScalingFactor)));
    mk = markerPointCoordinates(num2str(ringNumber * 10))./(1/ScalingFactor);
    outerMask = roipoly(innerMask,mk(:,1),mk(:,2));
    SE = strel('disk',100);
    outerMask=imdilate(outerMask,SE);
    currentRing =logical(outerMask - innerMask);    
    currentRing = imresize(currentRing, [sizeY, sizeX]);
else
    currentRing = logical(ones(sizeY,sizeX));
end
innerMask=0;
if(~size(currentRing,1)) == size(neuronHandler.ManualNeuronPositionsSparse(selectedWell),1) && size(currentRing,2) == size(neuronHandler.ManualNeuronPositionsSparse(selectedWell),2)
    currentRing = logical(ones(sizeY,sizeX));
end
ManualNeurons = logical(neuronHandler.ManualNeuronPositionsSparse(selectedWell) .* currentRing);
ManualNeurons = logical(ManualNeurons .* majorMask);
filterDistance=0;
nonFilterDistance=0;
SphereArea=0;
NucleusArea64=0;
markerPointCoordinates=0;
ManualNeuronsInvert = logical(full(logical(~ManualNeurons)));

%EdgeFilledNeurons = EdgeFilledNeurons(:,:,2) > area;
if(~isa(neuronHandler.CellPosMatrix,'containers.Map'))
    neuronHandler.CellPosMatrix = containers.Map();
end
if(isKey(neuronHandler.CellPosMatrix, selectedWell))
    CellomicsNeurons = logical(neuronHandler.CellPosMatrix(selectedWell));
else
    CellomicsNeurons = logical(zeros(sizeY,sizeX));
end
SkeletonNeurons = logical(full(neuronHandler.NeuronPositionsSkeletonization(selectedWell)) .* currentRing);

currentRing=0;
NeuronsDeleted = neuronHandler.NeuronPositionsSkelDeleted(selectedWell);

SkeletonNeurons = SkeletonNeurons - (NeuronsDeleted);
SkeletonNeurons = logical(SkeletonNeurons .* majorMask);
ind = SkeletonNeurons<0;
SkeletonNeurons(ind) = 0;


bothManualSkeletonized = nnz(SkeletonNeurons .* ManualNeurons);

AdditionalEdge = 0;% nnz(EdgeCompNeurons .* ManualNeuronsInvert);
CellomicsNeuronsInvert = full(logical(~CellomicsNeurons));
AdditionalManualCellomics = nnz(ManualNeurons .* CellomicsNeuronsInvert);
CellomicsNeuronsInvert=0;
AdditionalSkeletonManual = nnz(SkeletonNeurons .* ManualNeuronsInvert);
SkeletonNeuronsInvert = ~SkeletonNeurons;
AdditionalManualSkeleton = nnz(ManualNeurons .* SkeletonNeuronsInvert);
SkeletonNeuronsInvert=0;
EdgeCompNeurons=0;
ManualNeuronsInvert=0;
%EdgeCompInvert = ~EdgeCompNeurons;
AdditionalCellomics = 0;%;nnz(EdgeCompInvert .* CellomicsNeurons);
AdditionalManual = 0;%nnz(EdgeCompInvert .* ManualNeurons);
DP = bothManualSkeletonized;
FP = AdditionalSkeletonManual;
TN = AdditionalManualSkeleton;


% --------------------------------------------------------------------
function QualityCheck_Callback(hObject, eventdata, handles)
% hObject    handle to QualityCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
optionHandler = handles.OptionHandler;
WellNeuronDict = neuronHandler.NeuronPositionsEdgeComposite;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
foldername = [imageHandler.Foldername '/ConvertedCellomics'];


SphereAreaSizeX = optionHandler.MigrationDistanceDensityImageXSize;
SphereAreaSizeY = optionHandler.MigrationDistanceDensityImageYSize;
[sizeY sizeX]= size(imageHandler.NeuriteImage);
NucleusM = csvHandler.CellPosMatrix(selectedWell);

ringNumber = optionHandler.DensityDistributionRingNumber;
%Get density distribution to exclude too dense areas
[filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, optionHandler, NucleusM, imageHandler.Foldername, sizeY, sizeX); 
%MarkerPointCoordinates are available. Get second ring:
hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./(1/ScalingFactor)));
innerMask = logical(createMask(hInner));
outerMask = ones(uint16(sizeY/(1/ScalingFactor)),uint16(sizeX/(1/ScalingFactor)));
mk = markerPointCoordinates(num2str(ringNumber * 10))./(1/ScalingFactor);
outerMask = roipoly(innerMask,mk(:,1),mk(:,2));
SE = strel('disk',100);
outerMask=imdilate(outerMask,SE);
currentRing =logical(outerMask - innerMask);
currentRing = imresize(currentRing, [sizeY, sizeX]);
%currentRing=sparse(currentRing);
delete(hInner);
%delete(hOuter);
innerMask=0;

%EdgeCompNeurons = logical(WellNeuronDict(selectedWell));
area = optionHandler.EdgeFillNucleusAreaWithinNeurite;
%[EdgeCompSizeY EdgeCompSizeX] = size(EdgeCompNeurons);
[sizeY sizeX] = size(imageHandler.NeuriteImage);
if(~isKey(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = sparse(sizeY, sizeX);
end
ManualNeurons = logical(neuronHandler.ManualNeuronPositionsSparse(selectedWell) .* currentRing);
filterDistance=0;
nonFilterDistance=0;
SphereArea=0;
NucleusArea64=0;
markerPointCoordinates=0;
ManualNeuronsInvert = logical(full(logical(~ManualNeurons)));
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
bothCount = 0;%nnz(EdgeCompNeurons .* ManualNeurons);
bothEdgeCellomics = 0;%nnz(EdgeCompNeurons .* CellomicsNeurons);
bothManualCellomics = nnz(ManualNeurons .* CellomicsNeurons);
bothManualEdgeFilled = nnz(EdgeFilledNeurons .* ManualNeurons);
bothManualSkeletonized = nnz(SkeletonNeurons .* ManualNeurons);


AdditionalEdge = 0;% nnz(EdgeCompNeurons .* ManualNeuronsInvert);
CellomicsNeuronsInvert = ~CellomicsNeurons;
AdditionalEdgeCellomics = 0;% nnz(EdgeCompNeurons .* CellomicsNeuronsInvert);
AdditionalManualCellomics = nnz(ManualNeurons .* CellomicsNeuronsInvert);
CellomicsNeuronsInvert=0;
AdditionalCellomicsManual = nnz(CellomicsNeurons .* ManualNeuronsInvert);
AdditionalSkeletonManual = nnz(SkeletonNeurons .* ManualNeuronsInvert);
SkeletonNeuronsInvert = ~SkeletonNeurons;
AdditionalManualSkeleton = nnz(ManualNeurons .* SkeletonNeuronsInvert);
SkeletonNeuronsInvert=0;
EdgeCompNeurons=0;
AdditionalEdgeFilledManual = nnz(EdgeFilledNeurons .* ManualNeuronsInvert);
ManualNeuronsInvert=0;
EdgeFilledNeuronsInvert = ~EdgeFilledNeurons;
AdditionalManualEdgeFilled = nnz(ManualNeurons .* EdgeFilledNeuronsInvert);
%EdgeCompInvert = ~EdgeCompNeurons;
AdditionalCellomics = 0;%;nnz(EdgeCompInvert .* CellomicsNeurons);
AdditionalManual = 0;%nnz(EdgeCompInvert .* ManualNeurons);

disp([num2str(nnz(ManualNeurons)) ';' num2str(nnz(CellomicsNeurons)) ';' num2str(nnz(EdgeCompNeurons)) ';' num2str(nnz(CellomicsNeuronsEnhanced)) ';' num2str(nnz(SkeletonNeurons)) ';' num2str(nnz(EdgeFilledNeurons)) ';' num2str(bothCount) ';' num2str(AdditionalManual) ';' num2str(AdditionalEdge) ';;' num2str(bothEdgeCellomics) ';' num2str(AdditionalCellomics) ';' num2str(AdditionalEdgeCellomics) ';;' num2str(bothManualCellomics) ';' num2str(AdditionalManualCellomics) ';' num2str(AdditionalCellomicsManual) ';;' num2str(bothManualCellomicsEnhanced) ';' num2str(AdditionalManualCellomicsEnhanced) ';' num2str(AdditionalCellomicsEnhancedManual) ';;' num2str(bothManualEdgeFilled) ';' num2str(AdditionalManualEdgeFilled) ';' num2str(AdditionalEdgeFilledManual) ';;' num2str(bothManualSkeletonized) ';' num2str(AdditionalManualSkeleton) ';' num2str(AdditionalSkeletonManual) ';;']);
%disp([num2str(nnz(ManualNeurons)) ';' num2str(nnz(CellomicsNeurons)) ';' num2str(nnz(EdgeCompNeurons)) ';' num2str(bothCount) ';' num2str(AdditionalManual) ';' num2str(AdditionalEdge) ';;' num2str(bothEdgeCellomics) ';' num2str(AdditionalCellomics) ';' num2str(AdditionalEdgeCellomics) ';;' num2str(bothManualCellomics) ';' num2str(AdditionalManualCellomics) ';' num2str(AdditionalCellomicsManual) ';;' num2str(bothManualEdgeFilled) ';' num2str(AdditionalManualEdgeFilled) ';' num2str(AdditionalEdgeFilledManual) ';;']);
handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end



% --------------------------------------------------------------------
function AdvancedNeuronParameter_Callback(hObject, eventdata, handles)
% hObject    handle to AdvancedNeuronParameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
if(isfield(handles,'OptionHandler'))
    optionHandler = handles.OptionHandler;
    if(optionHandler == 0)
        CreateOptionHandler(handles);
        handles = guidata(handles.figure1);
        optionHandler = handles.OptionHandler;
    end
else
    CreateOptionHandler(handles);
    handles = guidata(handles.figure1);
    optionHandler = handles.OptionHandler;
end
prompt = {'Skeleton: Min Neurite Length', 'Tolerance Angle for Neuron from Endpoints', 'Max Distance of Neuron from Endpoints', 'Min pixels strong thresholded', 'Min distance between Neurons', 'Min size Neurite area', 'Fuzz Filter activated', 'Second Fuzz Filter activated', 'Skeleton: Max Overlay Size'};
dlg_title = 'Options';
num_lines = 1;
if(isnumeric(optionHandler.SkeletonThresholdMethod))
    optionHandler.SkeletonThresholdMethod = num2str(optionHandler.SkeletonThresholdMethod);
end
def = {num2str(optionHandler.SkeletonMinNeuriteLength), num2str(optionHandler.ToleranceAngleFromEndpoint), num2str(optionHandler.MaxDistanceFromEndpoint), num2str(optionHandler.MinNumberStrongThresholdPixels), num2str(optionHandler.MinDistanceBetweenNeurons), num2str(optionHandler.MinSizeNeuriteArea), num2str(optionHandler.FuzzFilterActivated), num2str(optionHandler.FuzzFilterSecondActivated), num2str(optionHandler.MaxOverlaySize)};
answer = inputdlg(prompt,dlg_title,num_lines,def);

optionHandler.SkeletonMinNeuriteLength = str2num(answer{1});
optionHandler.ToleranceAngleFromEndpoint = str2num(answer{2});
optionHandler.MaxDistanceFromEndpoint = str2num(answer{3});
optionHandler.MinNumberStrongThresholdPixels = str2num(answer{4});
optionHandler.MinDistanceBetweenNeurons = str2num(answer{5});
optionHandler.MinSizeNeuriteArea = str2num(answer{6});
optionHandler.FuzzFilterActivated = str2num(answer{7});
optionHandler.FuzzFilterSecondActivated = str2num(answer{8});
optionHandler.MaxOverlaySize = str2num(answer{9});
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles)

catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --------------------------------------------------------------------
function AdvancedOligoParam_Callback(hObject, eventdata, handles)
% hObject    handle to AdvancedOligoParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
if(isfield(handles,'OptionHandler'))
    optionHandler = handles.OptionHandler;
    if(optionHandler == 0)
        CreateOptionHandler(handles);
        handles = guidata(handles.figure1);
        optionHandler = handles.OptionHandler;
    end
else
    CreateOptionHandler(handles);
    handles = guidata(handles.figure1);
    optionHandler = handles.OptionHandler;
end
prompt = {'Composite Fill: Nuc Area in Percent within Neurite','Composite Fill: Nuc Area in Percent within Neurite Lower','Composite Fill: Neurite Area in Percent within Nucleus Min','Composite Fill: Neurite Area in Percent within Nucleus Max','Composite Fill: Min Size Neurite Area','Skeleton: Min Neurite Length','Threshold Method (1=Isodata for single picture, 2=mean Isodata, 3=manual Threshold)','Hard Threshold Distance from Manual Threshold','Manual Low Threshold Neurite', 'Manual Low Threshold Oligo', 'Fix Threshold for Nucleus Pic', 'Tolerance Angle for Neuron from Endpoints', 'Max Distance of Neuron from Endpoints', 'Min pixels strong thresholded', 'Min distance between Neurons', 'Min size Neurite area', 'Fuzz Filter activated', 'Second Fuzz Filter activated', 'Composite Fill: Look around', 'Skeleton: Max Overlay Size'};
dlg_title = 'Options';
num_lines = 1;
if(isnumeric(optionHandler.SkeletonThresholdMethod))
    optionHandler.SkeletonThresholdMethod = num2str(optionHandler.SkeletonThresholdMethod);
end
def = {num2str(optionHandler.EdgeFillNucleusAreaWithinNeurite),num2str(optionHandler.EdgeFillNucleusAreaWithinNeuriteSecond), num2str(optionHandler.CompositeFillNeuriteAreaWithinNucleus),num2str(optionHandler.CompositeFillNeuriteAreaWithinNucleusMax), num2str(optionHandler.MinSizeNeuriteAreaCompFill),  num2str(optionHandler.SkeletonMinNeuriteLength), optionHandler.SkeletonThresholdMethod, num2str(optionHandler.SkeletonNeuriteThresholdHardDistance),  num2str(optionHandler.SkeletonNeuriteThresholdLow), num2str(optionHandler.SkeletonOligoThresholdLow), num2str(optionHandler.NucleusThreshold), num2str(optionHandler.ToleranceAngleFromEndpoint), num2str(optionHandler.MaxDistanceFromEndpoint), num2str(optionHandler.MinNumberStrongThresholdPixels), num2str(optionHandler.MinDistanceBetweenNeurons), num2str(optionHandler.MinSizeNeuriteArea), num2str(optionHandler.FuzzFilterActivated), num2str(optionHandler.FuzzFilterSecondActivated), num2str(optionHandler.EdgeFillLookAround), num2str(optionHandler.MaxOverlaySize)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
optionHandler.EdgeFillNucleusAreaWithinNeurite = str2num(answer{1});
optionHandler.EdgeFillNucleusAreaWithinNeuriteSecond = str2num(answer{2});
optionHandler.CompositeFillNeuriteAreaWithinNucleus  = str2num(answer{3});
optionHandler.CompositeFillNeuriteAreaWithinNucleusMax = str2num(answer{4});
optionHandler.MinSizeNeuriteAreaCompFill = str2num(answer{5});
optionHandler.SkeletonMinNeuriteLength = str2num(answer{6});
optionHandler.SkeletonThresholdMethod = (answer{7});
optionHandler.SkeletonNeuriteThresholdHardDistance = str2num(answer{8});
optionHandler.SkeletonNeuriteThresholdLow = str2num(answer{9});
optionHandler.SkeletonOligoThresholdLow = str2num(answer{10});
optionHandler.NucleusThreshold = str2num(answer{11});
optionHandler.ToleranceAngleFromEndpoint = str2num(answer{12});
optionHandler.MaxDistanceFromEndpoint = str2num(answer{13});
optionHandler.MinNumberStrongThresholdPixels = str2num(answer{14});
optionHandler.MinDistanceBetweenNeurons = str2num(answer{15});
optionHandler.MinSizeNeuriteArea = str2num(answer{16});
optionHandler.FuzzFilterActivated = str2num(answer{17});
optionHandler.FuzzFilterSecondActivated = str2num(answer{18});
optionHandler.EdgeFillLookAround = str2num(answer{19});
optionHandler.MaxOverlaySize = str2num(answer{20});
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles)

catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --------------------------------------------------------------------
function NeuronAlgoOptions_Callback(hObject, eventdata, handles)
% hObject    handle to NeuronAlgoOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
if(isfield(handles,'OptionHandler'))
    optionHandler = handles.OptionHandler;
    if(optionHandler == 0)
        CreateOptionHandler(handles);
        handles = guidata(handles.figure1);
        optionHandler = handles.OptionHandler;
    end
else
    CreateOptionHandler(handles);
    handles = guidata(handles.figure1);
    optionHandler = handles.OptionHandler;
end
prompt = {'Threshold Method (1=Isodata for single picture, 2=mean Isodata, 3=manual Threshold)','Manual Low Threshold Neurite', 'Hard Threshold Distance from Man. Thresh.', 'Manual Low Threshold Oligo', 'Fix Threshold for Nucleus Pic'};
dlg_title = 'Options';
num_lines = 1;
if(isnumeric(optionHandler.SkeletonThresholdMethod))
    optionHandler.SkeletonThresholdMethod = num2str(optionHandler.SkeletonThresholdMethod);
end
def = {optionHandler.SkeletonThresholdMethod, num2str(optionHandler.SkeletonNeuriteThresholdLow), num2str(optionHandler.SkeletonNeuriteThresholdHardDistance), num2str(optionHandler.SkeletonOligoThresholdLow), num2str(optionHandler.NucleusThreshold)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
optionHandler.SkeletonThresholdMethod = (answer{1});
optionHandler.SkeletonNeuriteThresholdLow = str2num(answer{2});
optionHandler.SkeletonNeuriteThresholdHardDistance = str2num(answer{3});
optionHandler.SkeletonOligoThresholdLow = str2num(answer{4});
optionHandler.NucleusThreshold = str2num(answer{5});
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles)

catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --------------------------------------------------------------------
function Options_Callback(hObject, eventdata, handles)
% hObject    handle to Options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
if(isfield(handles,'OptionHandler'))
    optionHandler = handles.OptionHandler;
    if(optionHandler == 0)
        CreateOptionHandler(handles);
        handles = guidata(handles.figure1);
        optionHandler = handles.OptionHandler;
    end
else
    CreateOptionHandler(handles);
end
prompt = {'Migration Distance: Lower Threshold for Nucleus Pixels:', 'Density Distribution: Number of Rings:','Density Distribution: Width of Density Picture','Density Distribution: Height of Density Picture','Image File Format','Startwell','Excluded Wells','Cut Out Spherecore', 'ScalingFactor', 'Neuron Type'};
dlg_title = 'Options';
num_lines = 1;
if(numel(optionHandler.FileEnding) == 0)
    optionHandler.FileEnding = '.tif';
end
def = {num2str(optionHandler.MigDistLowerNucleusThreshold), num2str(optionHandler.DensityDistributionRingNumber),num2str(optionHandler.MigrationDistanceDensityImageXSize), num2str(optionHandler.MigrationDistanceDensityImageYSize),optionHandler.FileEnding,num2str(optionHandler.StartWell),optionHandler.ExcludedWells,num2str(optionHandler.CutOutSphereCore),num2str(optionHandler.ScalingFactor),num2str(optionHandler.NeuronType)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
optionHandler.MigDistLowerNucleusThreshold = str2num(answer{1});
optionHandler.DensityDistributionRingNumber = str2num(answer{2});
optionHandler.MigrationDistanceDensityImageXSize = str2num(answer{3});
optionHandler.MigrationDistanceDensityImageYSize = str2num(answer{4});
optionHandler.FileEnding = answer{5};

optionHandler.StartWell = str2num(answer{6});
optionHandler.ExcludedWells = answer{7};
optionHandler.CutOutSphereCore = str2num(answer{8});
optionHandler.ScalingFactor = str2num(answer{9});
optionHandler.NeuronType = answer{10};
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --------------------------------------------------------------------
function MultiParamAnalysisOptions_Callback(hObject, eventdata, handles)
% hObject    handle to MultiParamAnalysisOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
if(isfield(handles,'OptionHandler'))
    optionHandler = handles.OptionHandler;
    if(optionHandler == 0)
        CreateOptionHandler(handles);
        handles = guidata(handles.figure1);
        optionHandler = handles.OptionHandler;
    end
else
    CreateOptionHandler(handles);
    handles = guidata(handles.figure1);
    optionHandler = handles.OptionHandler;
end
prompt = {'First FP Threshold (in %)', 'Second FP Threshold (in %)', 'Use Filters (1=yes, First filter on each well used), 0=no (Calculation on whole Well),', 'Filter Size X', 'Filter Size Y', 'Excluded Wells (separated by ;)', 'Number of used CPU cores'};
dlg_title = 'Options';
num_lines = 1;
def = {num2str(optionHandler.MaxAllowedFPFirst), num2str(optionHandler.MaxAllowedFPSecond), num2str(optionHandler.FilterMulPa), num2str(optionHandler.RectangleSizeA), num2str(optionHandler.RectangleSizeB), optionHandler.ExcludedWellsMulPa,num2str(optionHandler.NumberOfCores)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
optionHandler.MaxAllowedFPFirst = str2num(answer{1});
optionHandler.MaxAllowedFPSecond = str2num(answer{2});
optionHandler.FilterMulPa = str2num(answer{3});
optionHandler.RectangleSizeA = str2num(answer{4});
optionHandler.RectangleSizeB = str2num(answer{5});
optionHandler.ExcludedWellsMulPa = answer{6};
optionHandler.NumberOfCores = str2num(answer{7});
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles)
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function CreateOptionHandler(handles,controlWells)
handles = guidata(handles.figure1);
optionHandler = Option();
    
     if(exist('controlWells'))
         optionHandler.SkeletonThresholdMethod = controlWells;
     else
         optionHandler.SkeletonThresholdMethod = 3;
     end    
    handles.OptionHandler = optionHandler;
    guidata(handles.figure1, handles);



% --------------------------------------------------------------------
function AdvancedCompFillParam_Callback(hObject, eventdata, handles)
% hObject    handle to AdvancedCompFillParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
if(isfield(handles,'OptionHandler'))
    optionHandler = handles.OptionHandler;
    if(optionHandler == 0)
        CreateOptionHandler(handles);
        handles = guidata(handles.figure1);
        optionHandler = handles.OptionHandler;
    end
else
    CreateOptionHandler(handles);
    handles = guidata(handles.figure1);
    optionHandler = handles.OptionHandler;
end
prompt = {'Composite Fill: Nuc Area in Percent within Neurite','Composite Fill: Nuc Area in Percent within Neurite Lower','Composite Fill: Neurite Area in Percent within Nucleus Min','Composite Fill: Neurite Area in Percent within Nucleus Max','Composite Fill: Min Size Neurite Area', 'Composite Fill: Look around'};
dlg_title = 'Options';
num_lines = 1;
if(isnumeric(optionHandler.SkeletonThresholdMethod))
    optionHandler.SkeletonThresholdMethod = num2str(optionHandler.SkeletonThresholdMethod);
end
def = {num2str(optionHandler.EdgeFillNucleusAreaWithinNeurite),num2str(optionHandler.EdgeFillNucleusAreaWithinNeuriteSecond), num2str(optionHandler.CompositeFillNeuriteAreaWithinNucleus),num2str(optionHandler.CompositeFillNeuriteAreaWithinNucleusMax), num2str(optionHandler.MinSizeNeuriteAreaCompFill), num2str(optionHandler.EdgeFillLookAround)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
optionHandler.EdgeFillNucleusAreaWithinNeurite = str2num(answer{1});
optionHandler.EdgeFillNucleusAreaWithinNeuriteSecond = str2num(answer{2});
optionHandler.CompositeFillNeuriteAreaWithinNucleus  = str2num(answer{3});
optionHandler.CompositeFillNeuriteAreaWithinNucleusMax = str2num(answer{4});
optionHandler.MinSizeNeuriteAreaCompFill = str2num(answer{5});
optionHandler.EdgeFillLookAround = str2num(answer{6});
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles)

catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --------------------------------------------------------------------
function HistogramOptions_Callback(hObject, eventdata, handles)
% hObject    handle to HistogramOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  try
handles = guidata(handles.figure1);
if(isfield(handles,'OptionHandler'))
    optionHandler = handles.OptionHandler;
    if(optionHandler == 0)
        CreateOptionHandler(handles);
        handles = guidata(handles.figure1);
        optionHandler = handles.OptionHandler;
    end
else
    CreateOptionHandler(handles);
    handles = guidata(handles.figure1);
    optionHandler = handles.OptionHandler;
end
prompt = {'Histogramm Correction: Max Intensity Nucleus:','Histogramm Correction: Max Intensity Neurite:','Histogramm Correction: Max Intensity Oligo:','Histogramm Correction: Max Intensity Astro:','Histogramm Correction: Min Intensity Nucleus:','Histogramm Correction: Min Intensity Neurite:','Histogramm Correction: Min Intensity Oligo:','Histogramm Correction: Min Intensity Astro:'};
dlg_title = 'Options';
num_lines = 1;
if(numel(optionHandler.FileEnding) == 0)
    optionHandler.FileEnding = '.tif';
end
def = {num2str(optionHandler.HistogramMaxNucleus),num2str(optionHandler.HistogramMaxNeurite),num2str(optionHandler.HistogramMaxOligo),num2str(optionHandler.HistogramMaxAstro),num2str(optionHandler.HistogramMinNucleus),num2str(optionHandler.HistogramMinNeurite),num2str(optionHandler.HistogramMinOligo),num2str(optionHandler.HistogramMinAstro)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
optionHandler.HistogramMaxNucleus = str2num(answer{1});
optionHandler.HistogramMaxNeurite = str2num(answer{2});
optionHandler.HistogramMaxOligo = str2num(answer{3});
optionHandler.HistogramMaxAstro = str2num(answer{4});
optionHandler.HistogramMinNucleus = str2num(answer{5});
optionHandler.HistogramMinNeurite = str2num(answer{6});
optionHandler.HistogramMinOligo = str2num(answer{7});
optionHandler.HistogramMinAstro = str2num(answer{8});
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
  
    
function mapString = mapPositionToKeyString(yPos,xPos)
    mapString = [num2str(xPos) ';' num2str(yPos)];


    
function [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithm(selectedWell, yStartPos, yEndPos, xStartPos, xEndPos, save, optionHandler, neuronHandler, NucleusM, foldername, channelName)
[sizeYorig sizeXorig] = size(NucleusM);
NucleusM = NucleusM(yStartPos:yEndPos,xStartPos:xEndPos);

%[D IDX] = bwdist(full(NucleusM));
TN=0;
if(isfield(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
    NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
    NeuronManualM = NeuronManualM(yStartPos:yEndPos,xStartPos:xEndPos);
    TN = nnz(NeuronManualM);
else
    NeuronManualM = sparse(logical(zeros(sizeYorig,sizeXorig)));
    if(isa(neuronHandler.ManualNeuronPositionsSparse,'containers.Map') == 0)
        neuronHandler.ManualNeuronPositionsSparse = containers.Map();
    end
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
end
%Load neurite and nucleus image from file 
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end
imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
imagePathNeuriteBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NeuriteBig' optionHandler.FileEnding];

%Check if NeuronStatMatrix exists
if(~isa(neuronHandler.NeuronStatMatrix,'containers.Map'))
    neuronHandler.NeuronStatMatrix = containers.Map();
end
if(~isa(neuronHandler.OmniNucleusPositions,'containers.Map'))
    neuronHandler.OmniNucleusPositions = containers.Map();
end
currentNeuronStat = containers.Map();
neuritePic = imread(imagePathNeuriteBig);
nucleusImage = imread(imagePathNucleusBig);
sizeY = yEndPos - yStartPos+1;
sizeX = xEndPos - xStartPos+1;
[sizeYOrig sizeXOrig]= size(neuritePic);

%Threshold Neurite Picture the known way
    
[sizeY sizeX] = size(neuritePic);
if(optionHandler.FilterMulPa~=1 && optionHandler.CutOutSphereCore == 1)
    neuritePic = CutOutCircles(neuritePic,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM, 0);
end
[neuritePic neuritePicStrong lMin] = ThresholdPic(neuritePic,optionHandler, neuronHandler.StretchlimResult, sizeY, sizeX,0,channelName);

SE = strel('disk', 3);
neuritePic = imdilate(neuritePic,SE);
neuritePic = neuritePic(yStartPos:yEndPos,xStartPos:xEndPos);
neuritePicStrong = neuritePicStrong(yStartPos:yEndPos,xStartPos:xEndPos);
nucleusPic=nucleusImage;
if(optionHandler.FilterMulPa~=1 && optionHandler.CutOutSphereCore == 1)
    nucleusPic = CutOutCircles(nucleusPic,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM,0, 0);
end
nucleusPic(nucleusPic<optionHandler.NucleusThreshold) = 0;
nucleusPic(nucleusPic>=optionHandler.NucleusThreshold) = 255;
nucleusPic=imcomplement(nucleusPic);
nucleusPic = uint8(nucleusPic);

%Do Watershed
%Check if there is already a watershed image. If not, create it.
imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
if(exist(imagePathWatershed,'file'))
    nucleusPic = imread(imagePathWatershed);
else
    javaaddpath('lib/ij.jar');
    javaaddpath('lib/mij.jar');
    MIJ.start(false);
    MIJ.createImage(nucleusPic);
    MIJ.run('Options...', 'iterations=1 count=1 edm=Overwrite do=Nothing');
    MIJ.run('Watershed');
    %nucleusPicBin = logical(MIJ.getCurrentImage);
    MIJ.run('Save', 'save=temp.tif');
    MIJ.closeAllWindows;
    MIJ.exit;
    nucleusPic = imread('temp.tif');
    delete temp.tif;
    nucleusPic = logical(1-nucleusPic);
    imwrite(nucleusPic,imagePathWatershed);
end
nucleusPic = nucleusPic(yStartPos:yEndPos,xStartPos:xEndPos);
[sizeYsmall sizeXsmall] = size(nucleusPic);
fusedPic = imfuse(nucleusPic, neuritePic);
%figure(1);
%imshow(fusedPic);
neuronPointIndices = fusedPic(:,:,1) > 250 & fusedPic(:,:,2) > 250 & fusedPic(:,:,3) > 250;
[sizeYneuronPoints sizeXneuronPoints] = size(neuronPointIndices);
%For every Nucleus:
%Check if Neurite with edges around.
%If yes: Count overlap of Nucleus and Neurites
%Count how many points of Nucleus are within the Neurite edges
EdgeFillNeuronPositionList = sparse(double(sizeY),double(sizeX));
EdgeFillNeuronPositionListSecond = sparse(double(sizeY),double(sizeX));
EdgeFillNeuronPositionList = EdgeFillNeuronPositionList(yStartPos:yEndPos,xStartPos:xEndPos);
EdgeFillNeuronPositionListSecond = EdgeFillNeuronPositionListSecond(yStartPos:yEndPos,xStartPos:xEndPos);
%stringCSVStatisticExport = sprintf('Position X;Position Y;Overlap (percent)\r\n');
[NonZeroY, NonZeroX] = find(NucleusM);
NonZeroY = uint16(NonZeroY);
NonZeroX = uint16(NonZeroX);

dbgYesCount = 0;
dbgNoCount = 0;
dbgNoStartFound = 0;
TP = 0;
FP = 0;
meanSize=zeros(0);
meanBrightness=zeros(0);
count=0;
nucleusPicBin = logical(nucleusPic);
% Get all connected components and label them
LabelComps=bwlabel(neuronPointIndices,8);
LabelCompsNeu = bwlabel(neuritePic);
%Check for each connected component, if also a place on NucleusPic is given
LabelCompsNuc = bwlabel(nucleusPicBin);
for i=1:numel(NonZeroY)    
    if(mod(i,10) == 0)
       disp([num2str(i) ' of ' num2str(numel(NonZeroY))]);
    end
    startX=0;
    startY=0;
    maxDistWhiteNuc = optionHandler.EdgeFillLookAround;
    maxDistWhiteNucIndex=-9;
    maxDistWhiteNucIndexSaved=-9;
    for j=-maxDistWhiteNuc:maxDistWhiteNuc
        for k=-maxDistWhiteNuc:maxDistWhiteNuc
            maxDistWhiteNucIndex = abs(abs(j)+abs(k));
            currentX = NonZeroX(i) + k;
            currentY = NonZeroY(i) + j;
            %If point or one neighbour is white in Overlaypic.            
            if(currentY < sizeY && currentX < sizeX && currentY>0 && currentX>0 && currentX<sizeXneuronPoints && currentY < sizeYneuronPoints && neuronPointIndices(currentY,currentX) > 0 && (startX ==0 || maxDistWhiteNucIndex < maxDistWhiteNucIndexSaved))
                startX= currentX;
                startY = currentY;                
            end
        end
    end
    if(startX > 0)
        %FloodFill white Area from current point. Check also nucleus points
        %stack = java.util.Stack();
        %Idea to be more efficient:
        %Cut out NucleusArea from startX - 50 until startX + 50
        %Do FloodFill Operation only on this subset
        if(startY-150>0 && startX-150>0 && startX+150 < sizeXsmall && startY + 150 < sizeYsmall)
            SubArea = neuronPointIndices(startY-150:startY+150,startX-150:startX+150);
            SubNucleusM = NucleusM(startY-150:startY+150,startX-150:startX+150);
            SubNucleusPic = nucleusPic(startY-150:startY+150,startX-150:startX+150);
            SubNucleusImage = nucleusImage(startY-150:startY+150,startX-150:startX+150);
            
            SubLabelComps = LabelComps(startY-150:startY+150,startX-150:startX+150);
            SubLabelCompsNuc = LabelCompsNuc(startY-150:startY+150,startX-150:startX+150);
            SubLabelCompsNeu = LabelCompsNeu(startY-150:startY+150,startX-150:startX+150);
            conCompIndex = SubLabelComps(150, 150);
            conCompIndexNuc = SubLabelCompsNuc(150, 150);
            conCompIndexNeu = SubLabelCompsNeu(150, 150);
            
            if(conCompIndex == 0)
               conCompIndex = SubLabelComps(151, 151); 
               if(conCompIndex == 0)
                   conCompIndex = SubLabelComps(149, 149); 
                   if(conCompIndex == 0)
                     conCompIndex = SubLabelComps(151, 149); 
                        if(conCompIndex == 0)
                            conCompIndex = SubLabelComps(149, 151); 
                        end
                   end
               end
            end
            if(conCompIndexNuc == 0)
               conCompIndexNuc = SubLabelCompsNuc(151, 151); 
               if(conCompIndexNuc == 0)
                   conCompIndexNuc = SubLabelCompsNuc(149, 149); 
                   if(conCompIndexNuc == 0)
                     conCompIndexNuc = SubLabelCompsNuc(151, 149); 
                        if(conCompIndexNuc == 0)
                            conCompIndexNuc = SubLabelCompsNuc(149, 151); 
                        end
                   end
               end
            end
            if(conCompIndexNeu == 0)
               conCompIndexNeu = SubLabelCompsNeu(151, 151); 
               if(conCompIndexNeu == 0)
                   conCompIndexNeu = SubLabelCompsNeu(149, 149); 
                   if(conCompIndexNeu == 0)
                     conCompIndexNeu = SubLabelCompsNeu(151, 149); 
                        if(conCompIndexNeu == 0)
                            conCompIndexNeu = SubLabelCompsNeu(149, 151); 
                        end
                   end
               end
            end
            
            
            %Extract area
            numberNeuritePixels = bwarea(SubArea == 1 & SubLabelCompsNuc == conCompIndexNuc);
            numberNucleusPixels = bwarea(SubLabelCompsNuc == conCompIndexNuc);            
            numberNeuritePixelsTotal = bwarea(conCompIndexNeu == SubLabelCompsNeu);
            
            count=count+1;
            meanSize = [meanSize numberNucleusPixels];
            meanBrightness = [meanBrightness CalculateNucleusBrightness(SubNucleusPic, SubNucleusImage)];
            [maxY, maxX] = size(SubNucleusPic);
            %Recalculate startX and startY            
         %   stack.push([50 50]);
            %FloodFill SubArea as well as Binary NucleusPic
            %Count number of pixels of both during FloodFill operation
      %      [numberNeuritePixels, numberNucleusPixels] = FloodFillForEdgeFill(SubNucleusPic,SubArea,maxX,maxY,stack);
            %WhiteArea2 = imfill(SubArea,[50 50]);         
            %[yIndices xIndices] = find(WhiteArea2);
            %yIndices = yIndices + (startY-50);
            %xIndices = xIndices + (startX-50);
        else
            %Extract connected component at given point
            conCompIndex = LabelComps(startY, startX);
            conCompIndexNuc = LabelCompsNuc(startY, startX);
            conCompIndexNeu = LabelCompsNeu(startY, startX);
            %Extract area
            numberNeuritePixels = bwarea(LabelCompsNuc == conCompIndexNuc & neuronPointIndices == 1);
            numberNucleusPixels = bwarea(LabelCompsNuc == conCompIndexNuc);
            numberNeuritePixelsTotal = bwarea(LabelCompsNeu == conCompIndexNeu);
%            stack.push([startX startY]);
            WhiteArea = zeros(sizeY,sizeX);
    %        [numberNeuritePixels, numberNucleusPixels] = FloodFillForEdgeFill(nucleusPic,neuronPointIndices,sizeX,sizeY,stack);
             %Schwerpunkt aller Punkte berechnen.
            %[yIndices xIndices] = find(WhiteArea);
        end     
        %Check if Area has not too less and not too much pixels

        %Check relation between numberNeuritePixels and numberNucleusPixels
        
        %Add neurite independent of it's neuritePixels per
        %NucleusPixels. Save that value.
        %Data structure: Matrix with positions AND relation
        
        %Mean calculation of brightness and NucleusArea
        
        area = optionHandler.EdgeFillNucleusAreaWithinNeurite;
        areaSecond = optionHandler.EdgeFillNucleusAreaWithinNeuriteSecond;
        areaAround = optionHandler.CompositeFillNeuriteAreaWithinNucleus;
        areaAroundMax = optionHandler.CompositeFillNeuriteAreaWithinNucleusMax;
        edgeFillMarked = 0;
        edgeFillMarkedSecond = 0;
        
        %Find next position of NucleusM from current NucleusOM
        %[curY curX] = ind2sub([sizeY,sizeX],IDX(NonZeroY(i),NonZeroX(i)));
        
        if(numberNeuritePixels/numberNucleusPixels >= area && numberNucleusPixels/numberNeuritePixelsTotal >= areaAround && numberNucleusPixels/numberNeuritePixelsTotal <= areaAroundMax && numberNeuritePixelsTotal > optionHandler.MinSizeNeuriteAreaCompFill) 
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
        
        if(NeuronManualM(uint16(yStartPos)-1+uint16(NonZeroY(i)),uint16(xStartPos)-1+uint16(NonZeroX(i))))
            manualMarked=1;            
            if(edgeFillMarked)
                TP = TP + 1;
                TN = TN - 1;
            end
        else
            manualMarked=0;
            if(edgeFillMarked)
                FP = FP + 1;
            end
        end
        if(save)
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
            eFillFirst = 0;
            eFillSecond=0;
        else
            eFillFirst = EdgeFillNeuronPositionList;
            eFillSecond=EdgeFillNeuronPositionListSecond;
        end
    else
        dbgNoStartFound = dbgNoStartFound +1;
    end    
end

meanSize = median(meanSize);
meanBrightness = median(meanBrightness);

if(numel(neuronHandler.NucleusBrightnessDict) ==0)% || ~isKey(neuronHandler.NucleusBrightnessDict, selectedWell)
    neuronHandler.NucleusBrightnessDict = containers.Map();
    neuronHandler.NucleusAreaDict = containers.Map();
    neuronHandler.NucleusBrightnessDict(selectedWell) = meanBrightness;
    neuronHandler.NucleusAreaDict(selectedWell) = meanSize;
else
    neuronHandler.NucleusAreaDict(selectedWell) = meanBrightness;
    neuronHandler.NucleusAreaDict(selectedWell) = meanSize;
end
if(numel(neuronHandler.NeuronPositionsEdgeFill) ==0)
    neuronHandler.NeuronPositionsEdgeFill = containers.Map();
end
if(~isKey(neuronHandler.NeuronPositionsEdgeFill, selectedWell)) 
    neuronHandler.NeuronPositionsEdgeFill(selectedWell) = sparse(logical(zeros(sizeYOrig,sizeXOrig)));
    fullMatrix = neuronHandler.NeuronPositionsEdgeFill(selectedWell);
    fullMatrix(yStartPos:yEndPos,xStartPos:xEndPos) = EdgeFillNeuronPositionList;
    neuronHandler.NeuronPositionsEdgeFill(selectedWell) = fullMatrix;
    %neuronHandler.NeuronPositionsEdgeFill(selectedWell) = EdgeFillNeuronPositionList;
else
    fullMatrix = neuronHandler.NeuronPositionsEdgeFill(selectedWell);
    fullMatrix(yStartPos:yEndPos,xStartPos:xEndPos) = EdgeFillNeuronPositionList;
    neuronHandler.NeuronPositionsEdgeFill(selectedWell) = fullMatrix;
end
if(numel(neuronHandler.NeuronPositionsEdgeFillSecond) ==0)
    neuronHandler.NeuronPositionsEdgeFillSecond = containers.Map();
end
if(~isKey(neuronHandler.NeuronPositionsEdgeFillSecond, selectedWell))
    neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = sparse(logical(zeros(sizeYOrig,sizeXOrig)));
    fullMatrix = neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
    fullMatrix(yStartPos:yEndPos,xStartPos:xEndPos) = EdgeFillNeuronPositionListSecond;
    neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = fullMatrix;
    %neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = EdgeFillNeuronPositionListSecond;
else
    fullMatrix = neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
    fullMatrix(yStartPos:yEndPos,xStartPos:xEndPos) = EdgeFillNeuronPositionListSecond;
    neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = fullMatrix;
end
neuronHandler.NeuronStatMatrix(selectedWell) = currentNeuronStat;
%handles.NeuronCoordinates = neuronHandler;
%handles.ImageHandler = imageHandler;
%guidata(handles.figure1, handles);    
%Check if Neuronhandler already has Dictionary of NeuronPositionsEdgeComposite


function watershaded = SplitUpNucleis(watershaded)
props = regionprops(watershaded,'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Centroid');
[sizeY sizeX] = size(watershaded);
%If major axis length too high, split up the area.
%Get orientation and from 
for(i=1:numel(props))      
    prop = props(i);
    if(prop.MajorAxisLength > 27)
        %Split up against the orientation of the major axis.
        %Vektor: Starte vom Mittelpunkt der Ellipse
        %In Richtung der k�rzeren Achse, also 90� und -90� -Winkel aus der
        %Hauptachsenrichtung.
        %Wende Sinus und Kosinussatz an
        deltaX = cosd(prop.Orientation+90);
        deltaY = sind(prop.Orientation+90);
        currentX = prop.Centroid(1);
        currentY = prop.Centroid(2);
        pathX = zeros(0);
        pathY = zeros(0);
        while(uint16(currentX) > 2 && uint16(currentY) > 2 && uint16(currentX) < 1022 && uint16(currentY) < 1022 && (watershaded(uint16(currentY),uint16(currentX)) == 1 || (watershaded(uint16(currentY+deltaY),uint16(currentX+deltaX)) == 1) || (watershaded(uint16(currentY+deltaY*2),uint16(currentX+deltaX*2)) == 1)))
            pathX = [pathX uint16(currentX)];
            pathY = [pathY uint16(currentY)];
            currentX = currentX + deltaX;
            currentY = currentY + deltaY;
        end
        deltaX = cosd(prop.Orientation-90);
        deltaY = sind(prop.Orientation-90);
        currentX = prop.Centroid(1);
        currentY = prop.Centroid(2);
        while(uint16(currentX) > 2 && uint16(currentY) > 2 && uint16(currentX) < sizeX-2 && uint16(currentY) < sizeY-2 && (watershaded(uint16(currentY),uint16(currentX)) == 1 || (watershaded(uint16(currentY+deltaY),uint16(currentX+deltaX)) == 1) || (watershaded(uint16(currentY+deltaY*2),uint16(currentX+deltaX*2)) == 1)))
            pathX = [pathX uint16(currentX)];
            pathY = [pathY uint16(currentY)];
            currentX = currentX + deltaX;
            currentY = currentY + deltaY;
        end
        for(i=1:numel(pathX))
            watershaded(pathY(i),pathX(i)) = 0;
            watershaded(pathY(i)-1,pathX(i)) = 0;
            watershaded(pathY(i),pathX(i)-1) = 0;
        end
    end
end

% --------------------------------------------------------------------
function EdgeFillNeurons_Callback(hObject, eventdata, handles)
% hObject    handle to EdgeFillNeurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
waitbarHandle = waitbar(0.5,'Executing Composite Fill', 'WindowStyle', 'modal');
%flappybird;
try
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
neuronHandler = handles.NeuronCoordinates;
csvHandler = handles.CSVCoordinates;
imageHandler = handles.ImageHandler;
foldername = imageHandler.Foldername;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
NucleusM = csvHandler.CellPosMatrix(selectedWell);
imageHandler = handles.ImageHandler;
[sizeY sizeX] = size(imageHandler.NucleusImage);
if(imageHandler.ZoomState)
   xMin = uint32(imageHandler.ZoomState(1));
    yMin = uint32(imageHandler.ZoomState(2));
    xMax = uint32(imageHandler.ZoomState(3) + xMin);
    yMax = uint32(imageHandler.ZoomState(4) + yMin);    
    [neuronHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell, yMin, yMax, xMin, xMax, 0, optionHandler, neuronHandler, NucleusM,foldername, 'NeuriteBig');
    RefreshZoomedImage(handles);
else
    firstY=1;
    firstX=1;
    if(sizeY - firstY > 10000 || sizeX - firstX > 10000)        
        sizeXsmall = (sizeX-firstX)/2;
        sizeYsmall = (sizeY-firstY)/2;
        [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithm(selectedWell,firstY,sizeYsmall,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,'Neurite');
        [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithm(selectedWell,firstY,sizeYsmall,sizeXsmall+1,sizeX,1,optionHandler,neuronHandler,NucleusM,foldername,'Neurite');
        [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithm(selectedWell,sizeYsmall+1,sizeY,sizeXsmall+1,sizeX,1,optionHandler,neuronHandler,NucleusM,foldername,'Neurite');
        [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithm(selectedWell,sizeYsmall+1,sizeY,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,'Neurite');
    else
        [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithm(selectedWell,firstY,sizeY, firstX, sizeX,1,optionHandler,neuronHandler,NucleusM,foldername,'Neurite');  
    end
    %[neuronHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell,1,sizeY,1,sizeX,1,optionHandler,neuronHandler,NucleusM,foldername);
end
handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles);   
%system('taskkill /fi "WINDOWTITLE eq FBird');
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
close(waitbarHandle);

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


function [returnPic returnPicMin levelMin]=ThresholdPic(inputImage,optionHandler,stretchlimLow,sizeY,sizeX,type,channelName)
    densityWidth = optionHandler.MigrationDistanceDensityImageXSize;
    densityHeight = optionHandler.MigrationDistanceDensityImageYSize;
    ringNumber = optionHandler.DensityDistributionRingNumber;
if(type <1)     
%0.5 Contrast stretching
%stretchlimLow = stretchlim(inputImage,[0.975 0.9999]);
%stretchlimLow = neuronHandler.StretchlimResult;
brightNeurites=inputImage;
if(strcmp(channelName,'Neurite'))
    brightNeurites = imadjust(inputImage,stretchlimLow,[0 1]);
    %1. Median Filter
    brightNeurites = medfilt2(brightNeurites);

    %2. Resharp image
    unsharpFilter = fspecial('unsharp');
    brightNeurites = imfilter(brightNeurites,unsharpFilter);
end

%3. Thresholding low
%level=graythresh(brightNeurites);
brightNeurites = imcomplement(brightNeurites);
%Calculate Threshold level for every Well separately.
%level = th_intermodes(brightNeurites)./255;
%level = isodata(brightNeurites);
%level = level + 0.1
if(strcmp(channelName, 'Oligo'))
    thresholdLevel = optionHandler.SkeletonOligoThresholdLow;
else
    thresholdLevel = optionHandler.SkeletonNeuriteThresholdLow;
end

if(optionHandler.SkeletonThresholdMethod == '1' || optionHandler.SkeletonThresholdMethod == 1)
    level = isodata(brightNeurites) + thresholdLevel;
    if(level>1)
        level=isodata(brightNeurites);
    elseif(level<0)
        level=0;
    end
elseif(optionHandler.SkeletonThresholdMethod == '2' || optionHandler.SkeletonThresholdMethod == 2)
    level=neuronHandler.IsodataResult + thresholdLevel;
elseif(optionHandler.SkeletonThresholdMethod == '3' || optionHandler.SkeletonThresholdMethod == 3)
    level = thresholdLevel;
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
%UB = 4000;
%CuttedIndices = BW;
%BW = xor(bwareaopen(BW,LB),  bwareaopen(BW,UB));
BW = bwareaopen(BW,LB);
%figure(4);
%imshow(BW);
returnPic = BW;
returnPicMin = BWMin;

else
        backgroundball =  imopen(inputImage,strel('disk',50));
        imgballbackground = inputImage - backgroundball;
        %Calculates the histogram of img1
        [lehisto x] = imhist(imgballbackground);
        %Calculates new background for img1 using the triangle method
        level = triangle_th(lehisto,256);
        %Creates a binary image of img1 using level
        BW = im2bw(imgballbackground,level);
        % Excludes all areas in 'imgnucleusbinary' smaller than N
        BW = bwareaopen(binaryimg,4000);
        % Neurites are eroded to obtain only cell bodies
        BW = imerode(binaryimg,SE);
        returnPic = BW;
        returnPicMin = BW;
end        
function selectedWellLong = GetSelectedWellLong(selectedWell)
selectedWellLong = selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end
    

function imageMap=CutOutCircles(imageMap,selectedWell,CutInnerCircle, CutOuterCircle,refreshIM,optionHandler,foldername,sizeY,sizeX, NucleusM, imageHandler)
%wellList = get(handles.lbWell, 'string');
%selectedWell = wellList{selectedWellNumber};
densityWidth = optionHandler.MigrationDistanceDensityImageXSize;
densityHeight = optionHandler.MigrationDistanceDensityImageYSize;
[sizeYorig sizeXorig] = size(imageMap);
ringNumber = optionHandler.DensityDistributionRingNumber;
% if(numel(strfind(foldername, 'ConvertedCellomics')) <= 0)
%     foldername = [foldername '/ConvertedCellomics'];
% end


[filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, densityWidth, densityHeight, optionHandler, NucleusM, foldername, sizeY, sizeX); 
%MarkerPointCoordinates are available. Get second ring:
selectedWellLong = GetSelectedWellLong(selectedWell);
%figure;
% if(~(imageHandler == 0) & imageHandler.ZoomState)
%     %Draw full pic
%     imshow(imageHandler.ResizedNeuriteImage);
% else
%     imshow(imread(imagePathNeuriteSmall));
% end
ScalingFactor = 0.1;
if(markerPointCoordinates ~= 0 && numel(markerPointCoordinates('10')./(1/ScalingFactor)) == 0)
    markerPointCoordinates=0;
end


if(markerPointCoordinates ~=0 && CutInnerCircle == 1 && CutOuterCircle == 0)    
    innerMask = ones(uint16(sizeY/(1/ScalingFactor)),uint16(sizeX/(1/ScalingFactor)));
    mk = markerPointCoordinates('10')./(1/ScalingFactor);
    innerMask = roipoly(innerMask,mk(:,1),mk(:,2));    
    currentRing =logical(ones(int32(sizeY/(1/ScalingFactor)),int32(sizeX/(1/ScalingFactor))) - innerMask);    
elseif(CutInnerCircle==1 && CutOuterCircle == 1 && markerPointCoordinates ~=0)
    innerMask = ones(uint16(sizeY/(1/ScalingFactor)),uint16(sizeX(1/ScalingFactor)));
    mk = markerPointCoordinates('10')./(1/ScalingFactor);
    innerMask = roipoly(innerMask,mk(:,1),mk(:,2));
    outerMask = ones(uint16(sizeY/(1/ScalingFactor)),uint16(sizeX/(1/ScalingFactor)));
    mk = markerPointCoordinates(num2str(ringNumber * 10))./(1/ScalingFactor);
    outerMask = roipoly(innerMask,mk(:,1),mk(:,2));
    SE = strel('disk',100);
    outerMask=imdilate(outerMask,SE);
    %currentRing =logical(ones(int32(sizeY/10),int32(sizeX/10)) - innerMask);
    currentRing =logical(outerMask - innerMask);
elseif(markerPointCoordinates ~=0);
    outerMask = ones(uint16(sizeY/(1/ScalingFactor)),uint16(sizeX/(1/ScalingFactor)));
    mk = markerPointCoordinates(num2str(ringNumber * 10))./(1/ScalingFactor);
    outerMask = roipoly(innerMask,mk(:,1),mk(:,2));
    currentRing = outerMask;
end
if(markerPointCoordinates ~=0)
  currentRing = imresize(currentRing, [sizeYorig, sizeXorig]);
    if(size(imageMap,1) ~= sizeY)
      %for(i=1:size(imageMap,1))
        image = imageMap;
        imageMap = uint8(image) .* uint8(currentRing);
      %end
    else
      if(issparse(imageMap))
        imageMap = (imageMap) .* (currentRing);
      else
        imageMap = uint8(imageMap) .* uint8(currentRing);
      end
    end
    else 
        if(size(imageMap,1) ~= sizeY)
      %for(i=1:size(imageMap,1))
        image = imageMap;
        imageMap = uint8(image);
      %end
    else
      imageMap = uint8(imageMap);
    end
end
if(~(imageHandler == 0) & imageHandler.ZoomState)
    if(refreshIM)
        RefreshZoomedImage(handles);
    end
end

function AreaResult = FloodFillNucArea(NucleusBinaryPic,yStart,xStart,NucleusM,D,Ind,origSizeY,origSizeX)
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
    if(x>=maxX || y >= maxY)
        continue;
    end
    if(pdist([yStart xStart;y x]) <=5)
        add=1;
    else        
        %[nucCol nucRow euclidDist] = csvHandler.FindNucleusForNeuron(x, y, NucleusM,20,0,0);        
        [nucCol nucRow euclidDist] = FindNucleusForNeuronWithoutDirection(x, y, D,Ind,origSizeY,origSizeX);
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


function [nucCol nucRow euclidDist] = FindNucleusForNeuronWithoutDirection(neuCol, neuRow, D, Ind,sizeY,sizeX)
       nucCol = -1;
       nucRow = -1;
       found = 0;
       iterator=0;
       euclidDist=99999;
       x=neuCol;
       y=neuRow;           
       %Distance Transformation for NucleusM
       [maxSmallY maxSmallX] = size(D);
       if(neuCol <= maxSmallX && neuRow <= maxSmallY)
           [nucRow nucCol] = ind2sub([sizeY sizeX],Ind(neuRow, neuCol));
           euclidDist = D(neuRow, neuCol);           
       end

function fuzzy = CheckForFuzzOnPosition(yMid,xMid,yFactor,xFactor,optionHandler, LowDensityNeurite, BrightDensity)
   fuzzy=0;
   xDens = uint8((xMid)/xFactor);
   yDens = uint8((yMid)/yFactor);   
   if(xDens<1)
       xDens=1;
   end
   if(xDens>=128)
       xDens=127;
   end
   if(yDens<1)
       yDens=1;
   end
   if(yDens >= 128)
       yDens = 127;
   end
       
    fuzzCounter=0;
    if(LowDensityNeurite(yDens,xDens) < optionHandler.FuzzFilterSecondActivated)
        fuzzCounter=fuzzCounter+4;
    end
    if(LowDensityNeurite(yDens+1,xDens) < optionHandler.FuzzFilterSecondActivated)
        fuzzCounter=fuzzCounter+1;
    end
    if(LowDensityNeurite(yDens,xDens+1) < optionHandler.FuzzFilterSecondActivated)
        fuzzCounter=fuzzCounter+1;
    end
    if(LowDensityNeurite(yDens+1,xDens+1) < optionHandler.FuzzFilterSecondActivated)
        fuzzCounter=fuzzCounter+1;
    end
    
    
    if(optionHandler.FuzzFilterActivated <= 11 && BrightDensity(yDens,xDens) < 3)
        fuzzCounter=fuzzCounter+1;
    end
    if(optionHandler.FuzzFilterActivated <= 11 && BrightDensity(yDens+1,xDens) < 3)
        fuzzCounter=fuzzCounter+1;
    end
    if(optionHandler.FuzzFilterActivated <= 11 && BrightDensity(yDens+1,xDens+1) < 3)
        fuzzCounter=fuzzCounter+1;
    end
    if(optionHandler.FuzzFilterActivated <= 11 && BrightDensity(yDens,xDens+1) < 3)
        fuzzCounter=fuzzCounter+1;
    end

    
    if (optionHandler.FuzzFilterActivated <= 11 && (fuzzCounter >= optionHandler.FuzzFilterActivated))
        fuzzy=1;  
    elseif(optionHandler.FuzzFilterActivated > 11 && fuzzCounter >= 3)
        fuzzy=1;
    end
       
function [neuronHandler TP FP TN] = SkeletonizationNeuronsAlgorithm(selectedWell,yStartPos, yEndPos, xStartPos, xEndPos, sav, optionHandler, sizeYorig, sizeXorig, foldername, neuronHandler, NucleusM, AnalyzeNeurites, channelName, competingHandler)
AnalyzeNeurites=0;
NucleusM = logical(NucleusM);
sizeY = yEndPos - yStartPos+1;
sizeX = xEndPos - xStartPos+1;
neuriteLengthMatrix = containers.Map();
SkelDeletedMat = logical(sparse(double(sizeY), double(sizeX)));
SkelNeurons = logical(sparse(double(sizeY), double(sizeX)));
NeuronPositionsEdgeFillNeurite = logical(sparse(double(sizeY), double(sizeX)));
ringNumber = optionHandler.DensityDistributionRingNumber;
selectedWellLong=selectedWell;
neuriteAreaCount=0;
neuriteCounter=0;
TP = 0;
FP = 0;
TN=0;
%NucleusM = NucleusM(yStartPos:yEndPos,xStartPos:xEndPos);
TruePositives = zeros(0,0);
%Just for testing purposes.
%foldernam = [foldername '/ConvertedCellomics'];
path = [foldername '/MarkerPointCoordinates-' selectedWell '.mat'];
mkdir(foldername,selectedWell);

%Get density distribution to exclude too dense areas
%if(saveCircle)     
ringNumber = optionHandler.DensityDistributionRingNumber;
% if(optionHandler.CutOutSphereCore==1)
%     load(path);
%     innerMask = ones(uint16(sizeY/10),uint16(sizeX/10));
%     if(markerPointCoordinates ~= 0)
%     mk = markerPointCoordinates('10')./10;
%     innerMask = roipoly(innerMask,mk(:,1),mk(:,2));
%     currentRing =logical(ones(int32(sizeY/10),int32(sizeX/10)) - innerMask);
%     currentRing =imresize(currentRing,[sizeY,sizeX]);
%     end
% end

if(~isa(neuronHandler.CellPosMatrix,'containers.Map'))
    neuronHandler.CellPosMatrix = containers.Map();
end
if(isKey(neuronHandler.CellPosMatrix,selectedWell))
    NeuronM = neuronHandler.CellPosMatrix(selectedWell);
else
    NeuronM = logical(zeros(sizeYorig,sizeXorig)); 
end
NeuronM = NeuronM(yStartPos:yEndPos,xStartPos:xEndPos);
if(~isa(neuronHandler.ManualNeuronPositionsSparse,'containers.Map'))
    neuronHandler.ManualNeuronPositionsSparse = containers.Map();
end
if(isKey(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
    ManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
else
    ManualM = sparse(logical(zeros(sizeYorig,sizeXorig))); 
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) =ManualM;
end

if(competingHandler==0)
    CompetingM = logical(zeros(sizeYorig,sizeXorig)); 
elseif(competingHandler ~=0)
    if(~isa(competingHandler.CellPosMatrix,'containers.Map'))
        competingHandler.CellPosMatrix = containers.Map();
    end
    if(isKey(neuronHandler.CellPosMatrix,selectedWell))
        CompetingM = competingHandler.CellPosMatrix(selectedWell);
    else
        CompetingM = logical(zeros(sizeYorig,sizeXorig)); 
    end
end
CompetingM = CompetingM(yStartPos:yEndPos,xStartPos:xEndPos);

ManualM = ManualM(yStartPos:yEndPos,xStartPos:xEndPos);
TN = nnz(ManualM);
[D Ind]= bwdist(full(logical(NucleusM)));
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
imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
imagePathNucleusSmall = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusSmall' optionHandler.FileEnding];
imagePathNeuriteBig = [foldername '/ConvertedCellomics/' selectedWellLong channelName 'Big' optionHandler.FileEnding];
imagePathNeuriteSmall = [foldername '/ConvertedCellomics/' selectedWellLong channelName 'Small' optionHandler.FileEnding];
imagePathSkeleton = [foldername '/ConvertedCellomics/' selectedWellLong 'Skeleton' optionHandler.FileEnding];
imagePath = [foldername '/ConvertedCellomics/'];
NeuriteImage = imread(imagePathNeuriteBig);
BinaryImage = sparse(logical(zeros(sizeYorig, sizeXorig)));
allfiles = dir(imagePath);
SkeletonImage = sparse(logical(zeros(sizeYorig, sizeXorig)));
Skelimg = 0;
% for(p=1:numel(allfiles))
% ind_Sk = strfind([foldername '/' allfiles(p).name],'Skeleton');
% if(numel(ind_Sk) > 0)
%     Skelimg = 1;
% end
% end
if(Skelimg == 1)    
SkeletonImageOld = imread(imagePathSkeleton);
end

if(optionHandler.CutOutSphereCore == 1)
    brightNeurites = CutOutCircles(NeuriteImage,selectedWell,1,0,0,optionHandler, foldername, sizeYorig, sizeXorig, NucleusM, 0);
else
    brightNeurites=NeuriteImage;
end

%Implement Top hat filter with rolling ball
%SE = strel('ball',8,7);
%brightNeurites = imtophat(brightNeurites,SE);

stretchlimLow = neuronHandler.StretchlimResult;
[BW BWStrong lMin] = ThresholdPic(brightNeurites,optionHandler,stretchlimLow,sizeYorig, sizeXorig,0,channelName);
%BW = imfill(BW);
%ToDo: Fill holes and imdilate
SE = strel('disk', 1);
BW = imdilate(BW,SE);
BW=BW(yStartPos:yEndPos,xStartPos:xEndPos);
BWStrong=BWStrong(yStartPos:yEndPos,xStartPos:xEndPos);
NucleusImage = imread(imagePathNucleusBig);
if(optionHandler.CutOutSphereCore == 1)
    nucleusPic = CutOutCircles(NucleusImage,selectedWell,1,0,1,optionHandler, foldername, sizeY, sizeX, NucleusM, 0);
else
    nucleusPic=NucleusImage;
end
nucleusPic(nucleusPic<optionHandler.NucleusThreshold) = 0;%double(double(imageHandler.NucleusImage) ./ 255);
nucleusPic(nucleusPic>=optionHandler.NucleusThreshold) = 1;
nucleusPic = uint8(nucleusPic);
nucleusPic2 = nucleusPic(yStartPos:yEndPos,xStartPos:xEndPos);
overlayAreaIndices = find(nucleusPic2 > 0 & BW > 0);
nucleusPic2=0;
overlayArea=logical(zeros(sizeY,sizeX));
overlayArea(overlayAreaIndices)=1;



%Do Watershed
%Check if there is already a watershed image. If not, create it.
imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
if(exist(imagePathWatershed,'file'))
    nucleusPic = imread(imagePathWatershed);
else
    javaaddpath('lib/ij.jar');
    javaaddpath('lib/mij.jar');
    MIJ.start(false);
    MIJ.createImage(nucleusPic);
    MIJ.run('Options...', 'iterations=1 count=1 edm=Overwrite do=Nothing');
    MIJ.run('Watershed');
    %nucleusPicBin = logical(MIJ.getCurrentImage);
    MIJ.run('Save', 'save=temp.tif');
    MIJ.closeAllWindows;
    MIJ.exit;
    nucleusPic = imread('temp.tif');
    delete temp.tif;
    nucleusPic = logical(1-nucleusPic);
    imwrite(nucleusPic,imagePathWatershed);
end
nucleusPic = nucleusPic(yStartPos:yEndPos,xStartPos:xEndPos);


[thresholdedRowsBig thresholdedColsBig] = find(NeuriteImage > 25);
[thresholdedRowsSmall thresholdedColsSmall] = find(NeuriteImage > 11);
thresholdedRowsBig = uint16(thresholdedRowsBig);
thresholdedColsBig = uint16(thresholdedColsBig);
thresholdedRowsSmall = uint16(thresholdedRowsSmall);
thresholdedColsSmall = uint16(thresholdedColsSmall);

thresholdedRowsBig = [thresholdedRowsBig;sizeYorig];
thresholdedColsBig = [thresholdedColsBig;sizeXorig];
thresholdedRowsBig = [thresholdedRowsBig;0];
thresholdedColsBig = [thresholdedColsBig;0];
thresholdedRowsSmall = [thresholdedRowsSmall;sizeYorig];
thresholdedColsSmall = [thresholdedColsSmall;sizeXorig];
thresholdedRowsSmall = [thresholdedRowsSmall;0];
thresholdedColsSmall = [thresholdedColsSmall;0];
BrightDensity = (hist3([double(thresholdedRowsBig),double(thresholdedColsBig)],[128,128]));
LowDensity = (hist3([double(thresholdedRowsSmall),double(thresholdedColsSmall)],[128,128]));


%[thresholdedRowsBigNuclei thresholdedColsBigNuclei] = find(NucleusImage > optionHandler.NucleusThreshold);
%[thresholdedRowsSmallNuclei thresholdedColsSmallNuclei] = find(NucleusImage < 20);
%thresholdedRowsBigNuclei = uint16(thresholdedRowsBigNuclei);
%thresholdedColsBigNuclei = uint16(thresholdedColsBigNuclei);
%thresholdedRowsSmallNuclei = uint16(thresholdedRowsSmallNuclei);
%thresholdedColsSmallNuclei = uint16(thresholdedColsSmallNuclei);

%thresholdedRowsBigNuclei = [thresholdedRowsBigNuclei;sizeYorig];
%thresholdedColsBigNuclei = [thresholdedColsBigNuclei;sizeXorig];
%thresholdedRowsBigNuclei = [thresholdedRowsBigNuclei;0];
%thresholdedColsBigNuclei = [thresholdedColsBigNuclei;0];
%thresholdedRowsSmallNuclei = [thresholdedRowsSmallNuclei;sizeYorig];
%thresholdedColsSmallNuclei = [thresholdedColsSmallNuclei;sizeXorig];
%thresholdedRowsSmallNuclei = [thresholdedRowsSmallNuclei;0];
%thresholdedColsSmallNuclei = [thresholdedColsSmallNuclei;0];
%BrightDensityNuclei = (hist3([double(thresholdedRowsBigNuclei),double(thresholdedColsBigNuclei)],[128,128]));
%LowDensityNuclei = (hist3([double(thresholdedRowsSmallNuclei),double(thresholdedColsSmallNuclei)],[128,128]));


%[thresholdedRowsBigNeurite thresholdedColsBigNeurite] = find(NeuriteImage > optionHandler.NucleusThreshold);
[thresholdedRowsSmallNeurite thresholdedColsSmallNeurite] = find(NeuriteImage < 5);
%thresholdedRowsBigNeurite = uint16(thresholdedRowsBigNeurite);
%thresholdedColsBigNeurite = uint16(thresholdedColsBigNeurite);
thresholdedRowsSmallNeurite = uint16(thresholdedRowsSmallNeurite);
thresholdedColsSmallNeurite = uint16(thresholdedColsSmallNeurite);

%thresholdedRowsBigNeurite = [thresholdedRowsBigNeurite;sizeYorig];
%thresholdedColsBigNeurite = [thresholdedColsBigNeurite;sizeXorig];
%thresholdedRowsBigNeurite = [thresholdedRowsBigNeurite;0];
%thresholdedColsBigNeurite = [thresholdedColsBigNeurite;0];
thresholdedRowsSmallNeurite = [thresholdedRowsSmallNeurite;sizeYorig];
thresholdedColsSmallNeurite = [thresholdedColsSmallNeurite;sizeXorig];
thresholdedRowsSmallNeurite = [thresholdedRowsSmallNeurite;0];
thresholdedColsSmallNeurite = [thresholdedColsSmallNeurite;0];
%BrightDensityNeurite = (hist3([double(thresholdedRowsBigNeurite),double(thresholdedColsBigNeurite)],[128,128]));
%Catch out of memory exception and deactivate Fuzz, because it is not so
%important and leads to out of memory exceptions on OT
try
    LowDensityNeurite = (hist3([double(thresholdedRowsSmallNeurite),double(thresholdedColsSmallNeurite)],[128,128]));
catch
    LowDensityNeurite = logical(zeros(128,128));
end
thresholdedRowsSmallNeurite=0;
thresholdedRowsSmallNeurite=0;
thresholdedColsSmallNeurite=0;
%LowDensity = imcomplement(LowDensity);
%ToDo: Save Neurite objects with single pictures to trace
%Check if already a Neuron by EdgeFillNeurons in this area available
%If not, trace Neurite and check if Neuron possible on one or
%another side  

%The same maybe faster:
L = uint16(bwlabel(BW,8));
NotList = [0];

%Comment in this part, if you want to use SubThresholding for single big
%areas of the Neurite image
% while 1
%     %Get label with highest occurance in L:
%     mostOccuredLabel = mode(L(~ismember(L,NotList)));
%     if(mostOccuredLabel==0)
%         break;
%     end
%     NotList = [NotList mostOccuredLabel];
%     currentNeuriteArea = logical(sparse(zeros(size(BW,1), size(BW,2))));
%     
%     if(numel(find(L==mostOccuredLabel)) < 1000000)
%     currentNeuriteArea(find(L == mostOccuredLabel))=1;  
%     
%     areaSize = bwarea(full(currentNeuriteArea));
%     if(areaSize > optionHandler.MinSizeNeuriteArea)
%         %Apply harder threshold in relevant part of BW
%         newNeurite = Neurite();
%         newNeuriteBig = Neurite();
%         [rowIndices colIndices] = find(currentNeuriteArea);
%         
%         newNeurite.cutYPosStart = min(rowIndices) - 2;
%         newNeurite.cutXPosStart = min(colIndices) - 2;    
%         newNeurite.cutYPosEnd = max(rowIndices) + 2;
%         newNeurite.cutXPosEnd = max(colIndices) + 2;
%         if(newNeurite.cutYPosStart <= 0)
%             newNeurite.cutYPosStart = 1;
%         end
%         if(newNeurite.cutXPosStart <= 0)
%             newNeurite.cutXPosStart = 1;
%         end
%         if(newNeurite.cutXPosEnd > size(currentNeuriteArea,2))
%             newNeurite.cutXPosEnd = size(currentNeuriteArea,2);
%         end
%         if(newNeurite.cutYPosEnd > size(currentNeuriteArea,1))
%             newNeurite.cutYPosEnd = size(currentNeuriteArea,1);
%         end
%         %Try to split up neurite
%         NeuriteSubRegion = NeuriteImage(newNeurite.cutYPosStart+yStartPos-1:newNeurite.cutYPosEnd+yStartPos-1,newNeurite.cutXPosStart+xStartPos-1:newNeurite.cutXPosEnd+xStartPos-1);
%         NeuriteSubRegion = imcomplement(NeuriteSubRegion);
%         subHist = imhist(NeuriteSubRegion);
%         thresholdSub = IJIsoData(subHist);
%         currentNeuriteArea = im2bw(NeuriteSubRegion,thresholdSub/255);
%         BW(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd) = 1-currentNeuriteArea;
%     else
%         break;
%     end
%     else
%         break;
%     end
% end
%Again: Remove to small areas of BW and Cut out Circles again
LB = 75;
%UB = 4000;
%CuttedIndices = BW;
%BW = xor(bwareaopen(BW,LB),  bwareaopen(BW,UB));

BW = bwareaopen(BW,LB);
if(optionHandler.CutOutSphereCore == 1)
    BW = logical(CutOutCircles(BW,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM, 0));
else
    BW=logical(BW);
end
CC = bwconncomp(BW);

%[whitePointRows, whitePointCols] = find(BW);
alreadyVisited = logical(ones(size(BW,1), size(BW,2)));
neuriteAreasDict = containers.Map();
i=0;
%BW2= logical(zeros(size(BW,1), size(BW,2)));
secondThreshold = 20;

%Get different areas were nucleus Pic and BW are both > 0

%currentNeuriteAreas = logical(sparse(zeros(size(BW,1),size(BW,2),CC.NumObjects)));
for cou=1:CC.NumObjects
%while(numel(whitePointRows) > 0)
    %Get first whitepointindex <> 0
    %FloodFill    
    %stack = java.util.Stack();
    %stack.push([whitePointRows(1) whitePointCols(1)]);   
    %currentNeuriteArea = currentNeuriteAreas(:,:,cou);
    %if(numel(CC.PixelIdxList{cou}) < 1000000)
    currentNeuriteArea = sparse(false(sizeY, sizeX));    
    strongCount = 0;

    
    currentNeuriteArea(CC.PixelIdxList{cou}) = 1;
    strongCount = BWStrong(CC.PixelIdxList{cou});
    strongCount = bwarea(strongCount);
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
    BrightDensityNuclei = imhist(NucleusImage(newNeuriteBig.cutYPosStart:newNeuriteBig.cutYPosEnd,newNeuriteBig.cutXPosStart:newNeuriteBig.cutXPosEnd));
    nucVariance = var(BrightDensityNuclei);
    
    %ToDo: Evtl. differentiate between Zoomed and Non Zoomed Mode!
    %currentNeuriteArea=currentNeuriteArea(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
    currentNeuriteArea = logical(1-currentNeuriteArea);
    currentNeuriteAreaSmall=currentNeuriteArea(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
    
   xMid = uint16(((newNeurite.cutXPosStart + newNeurite.cutXPosEnd)/2) +xStartPos-1);
   xMid2 = uint16((newNeurite.cutXPosStart + 5) +xStartPos-1);
   yMid = uint16(((newNeurite.cutYPosStart + newNeurite.cutYPosEnd)/2) +yStartPos-1);
   yMid2 = uint16((newNeurite.cutYPosStart + 5) +yStartPos-1);
   xFactor = sizeXorig/128;
   yFactor = sizeYorig/128;
   xDens = uint8((xMid)/xFactor);
   yDens = uint8((yMid)/yFactor);
   xDens2 = uint8((xMid2)/xFactor);
   yDens2 = uint8((yMid2)/yFactor);
   currentX =xMid;
   currentY = yMid;
   if(xDens<1)
       xDens=1;
   end
   if(xDens>=128)
       xDens=127;
   end
   if(yDens<1)
       yDens=1;
   end
   if(yDens >= 128)
       yDens = 127;
   end
   if(xDens2<1)
       xDens2=1;
   end
   if(xDens2>=128)
       xDens2=127;
   end
   if(yDens2<1)
       yDens2=1;
   end
   if(yDens2 >= 128)
       yDens2 = 127;
   end
    
    %Copy Neurite Area to composed Skeleton image
    newNeurite.image = currentNeuriteAreaSmall;
    newNeuriteBig.image = currentNeuriteAreaSmall;    
    fuzzy=0;

    if((optionHandler.FuzzFilterSecondActivated==1 || optionHandler.FuzzFilterSecondActivated==2) && (LowDensityNeurite(yDens2,xDens2) < 2 || LowDensityNeurite(yDens,xDens) < 2 || LowDensityNeurite(yDens2,xDens2) < 2 || LowDensityNeurite(yDens2,xDens2) < 2 || LowDensityNeurite(yDens2,xDens2) < 2))   
        fuzzy=1;
    end
    if (optionHandler.FuzzFilterActivated > 0.5 && optionHandler.FuzzFilterActivated < 2.5 && ((BrightDensity(yDens2,xDens2) < 3 && BrightDensity(yDens2+1,xDens2) < 3 && BrightDensity(yDens2+1,xDens2+1) < 3 && BrightDensity(yDens2,xDens2+1) < 3 && BrightDensity(yDens,xDens) < 3) || LowDensity(yDens,xDens) > 2000 || LowDensity(yDens+1,xDens+1) > 2000 || LowDensity(yDens+1,xDens) > 2000 || LowDensity(yDens,xDens+1) > 2000))
        fuzzy=1;
    end
    aSize = bwarea(1-currentNeuriteAreaSmall);
	NucAreaSmall = nucleusPic(newNeuriteBig.cutYPosStart:newNeuriteBig.cutYPosEnd,newNeuriteBig.cutXPosStart:newNeuriteBig.cutXPosEnd);  
	SE = strel('disk', 3);
	[sizeYskel sizeXskel] = size(NucAreaSmall);
	newNeuriteBig.nucleusImage = NucAreaSmall;
    if(aSize > optionHandler.MinSizeNeuriteArea)
          noNeurite=1;
%          %Do separate Thresholding on neurite image for current region
%          NeuriteSubRegion = imageHandler.NeuriteImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
%          NeuriteSubRegion = imcomplement(NeuriteSubRegion);
%          thresoldSub = isodata(NeuriteSubRegion)
%          subHist = imhist(NeuriteSubRegion);
%          thresholdSub = IJIsoData(subHist);
%          currentNeuriteArea = im2bw(NeuriteSubRegion./255,thresholdSub/255);
%         %FuzzTest Cut off evtl.
         
     else
        noNeurite=0;
		if(optionHandler.NeuronType == 0)
			try        

				disp([num2str(cou) ' of ' num2str(CC.NumObjects)]);
				newNeurite.CalculateBaiSkeleton(6); 
				newNeuriteBig.CalculateBaiSkeleton(16);            
				newNeurite.FindBranchingCrossingPoints();
				newNeurite.KillShortBranches();
				newNeuriteBig.FindBranchingCrossingPoints();



				newNeuriteBig.KillShortBranches();
				newNeurite.CalculateNeuriteLength();
				newNeuriteBig.CalculateNeuriteLength();
				if(newNeurite.neuriteLength > newNeuriteBig.neuriteLength)
					newNeuriteBig.neuriteLength = newNeurite.neuriteLength;
				end

			SkeletonImage(newNeuriteBig.cutYPosStart+yStartPos-1:newNeuriteBig.cutYPosEnd+yStartPos-1,newNeuriteBig.cutXPosStart+xStartPos-1:newNeuriteBig.cutXPosEnd+xStartPos-1)=logical(SkeletonImage(newNeuriteBig.cutYPosStart+yStartPos-1:newNeuriteBig.cutYPosEnd+yStartPos-1,newNeuriteBig.cutXPosStart+xStartPos-1:newNeuriteBig.cutXPosEnd+xStartPos-1)+newNeuriteBig.skeletonImage);


			catch
				disp('Skeletonization error');
				noNeurite=1;
			end        
		else
			if(Skelimg == 0)
             disp([num2str(cou) ' of ' num2str(CC.NumObjects)]);
            try
              newNeurite.CalculateBaiSkeleton(6); 
              newNeuriteBig.CalculateBaiSkeleton(40); 
            catch
              disp('Skeletonization error');
              noNeurite=1;
            end
            
           else
             %newNeuriteBig.skeletonImage = SkeletonImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
             %newNeuriteBig.skeletonImageOrig = SkeletonImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
             %neneurwNeurite.skeletonImage = newNeuriteBig.skeletonImage;
             newNeuriteBig.skeletonImage = zeros(sizeYskel,sizeXskel);
             newNeuriteBig.skeletonImage(SkeletonImageOld(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd) == 1 & (currentNeuriteAreaSmall == 0)) = 1;
             newNeurite.skeletonImage = newNeuriteBig.skeletonImage;
             newNeuriteBig.skeletonImageOrig = newNeuriteBig.skeletonImage;
             
         end
           try     
            disp([num2str(cou) ' of ' num2str(CC.NumObjects)]);
            newNeurite.CalculateBaiSkeleton(6); 
            %newNeuriteBig.CalculateBaiSkeleton(40);            
            newNeurite.FindBranchingCrossingPoints();
            newNeurite.KillShortBranches();
            newNeuriteBig.FindBranchingCrossingPoints();
              newNeuriteBig.KillShortBranches();
              newNeuriteBig.CalculateNeuriteLength();
              %newNeuriteBig.CompleteNeuritePostprocessing();
            %newNeuriteBig.KillShortBranches();
            newNeurite.CalculateNeuriteLength();
            newNeuriteBig.CalculateNeuriteLength();
            if(newNeurite.neuriteLength > newNeuriteBig.neuriteLength)
                newNeuriteBig.neuriteLength = newNeurite.neuriteLength;
            end
            
            SkeletonImage(newNeuriteBig.cutYPosStart+yStartPos-1:newNeuriteBig.cutYPosEnd+yStartPos-1,newNeuriteBig.cutXPosStart+xStartPos-1:newNeuriteBig.cutXPosEnd+xStartPos-1)=logical(SkeletonImage(newNeuriteBig.cutYPosStart+yStartPos-1:newNeuriteBig.cutYPosEnd+yStartPos-1,newNeuriteBig.cutXPosStart+xStartPos-1:newNeuriteBig.cutXPosEnd+xStartPos-1)+newNeuriteBig.skeletonImage);
            
            
			catch
				disp('Skeletonization error');
				noNeurite=1;
			end        
		end
        
    end    
    
    
    %%For debugging purposes:
%     if(imageHandler.ZoomState)
%      xMin = uint32(imageHandler.ZoomState(1));
%      yMin = uint32(imageHandler.ZoomState(2));
%      xMax = uint32(imageHandler.ZoomState(3) + xMin);
%      yMax = uint32(imageHandler.ZoomState(4) + yMin);
%      
%      
%      neuriteChannelAround=imageHandler.NeuriteImage(yMin-50:yMax+50,xMin-50:xMax+50);
%      nucleusChannelAround=imageHandler.NucleusImage(yMin-50:yMax+50,xMin-50:xMax+50);
%      figure(2);
%      imshow(imfuse(neuriteChannelAround,nucleusChannelAround));
%     end
     
    % newNeurite.PlotSkeleton();
    
    
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

   
   %Check also size of connected Nuc Area around
   %NucArea = nucleusPic;   
   %nucSize = size(find(NucArea),1);
   %neuSize = numel(find(currentNeuriteArea));
   %If relation between skeleton length and neurite area is wrong (big
   %neurite area and small skeleton, then go ahead)
   %ToDo: It would be better to use the full CompositeFillResults instead!
    if(optionHandler.CutOutSphereCore == 0 && aSize > 10000)
%       %Use CompositeFill Results for this area
       edgeFillPositions = neuronHandler.NeuronPositionsEdgeFill(selectedWell);
%       edgeFillPositions = edgeFillPositions(yStartPos:yEndPos,xStartPos:xEndPos);
%       edgeFillPositions = edgeFillPositions(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
%       [yEdgeFills xEdgeFills] = find(edgeFillPositions);
%       
    end
    if(( noNeurite == 0 && (newNeuriteBig.neuriteLength < optionHandler.SkeletonMinNeuriteLength || maxDistEndpoints < optionHandler.SkeletonMinNeuriteLength) && bwarea(1-currentNeuriteAreaSmall) < 10000) || strongCount < optionHandler.MinNumberStrongThresholdPixels || fuzzy == 1)% || numel(find(currentNeuriteArea)) > 60000) 
        %ToDo: Delete potential Neuclei on Neurite
        %Save position of EdgeFilled Nuclei, which should be deleted
        %edgeFillPositions = neuronHandler.NeuronPositionsEdgeFill(selectedWell);
        edgeFillPositionsSecond = neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
        edgeFillPositionsSecond = edgeFillPositionsSecond(yStartPos:yEndPos,xStartPos:xEndPos);
        edgeFillPositionsSecond = edgeFillPositionsSecond(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
        %Add Neurons to deleted ones
        [yPoses xPoses] = find(edgeFillPositionsSecond);
        if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
            neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
        end        
        
        %Check if Matrix already exists
        if(isKey(neuronHandler.NeuronPositionsSkelDeleted,selectedWell))
            fullMatrix = neuronHandler.NeuronPositionsSkelDeleted(selectedWell);
            fullMatrix(yStartPos:yEndPos,xStartPos:xEndPos) = SkelDeletedMat;
            neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = fullMatrix;
        else               
            neuronHandler.NeuronPositionsSkelDeleted(selectedWell)= sparse(logical(zeros(sizeYorig,sizeXorig)));
            fullMatrix = neuronHandler.NeuronPositionsSkelDeleted(selectedWell);
            fullMatrix(yStartPos:yEndPos,xStartPos:xEndPos) = SkelDeletedMat;
            neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = fullMatrix;
        end                
    else
        %ToDo: Create mapping between Neurites and Neurons
        %To Kill: All Neurons without regarding Neurite
        %First, check if already marked Neuron by EdgeFill-Algorithm available.
        overlayAreaSmall = overlayArea(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
        overlayAreaSmall = bwlabel(overlayAreaSmall);
        
        edgeFillPositions = neuronHandler.NeuronPositionsEdgeFill(selectedWell);
        edgeFillPositions = edgeFillPositions(yStartPos:yEndPos,xStartPos:xEndPos);
        edgeFillPositions = edgeFillPositions(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
        [yEdgeFills xEdgeFills] = find(edgeFillPositions);
        if(numel(yEdgeFills) >= 1)
            %done
            %Save Mapping between Neuron and Neurite.     
            %Check which Neuron has max. distance to next endpoint
            if(numel(yEdgeFills)==1)
                NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
                %Add position of underlying nucleus to imageHandler.BinaryImage
                NucAreaSmall = nucleusPic(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                
                
                if(overlayAreaSmall(yEdgeFills(1),xEdgeFills(1)) == 0)
                    overlaySize=0;
                else
                    overlayLabel = overlayAreaSmall(yEdgeFills(1),xEdgeFills(1));
                    overlayAreaSmall(overlayAreaSmall~=overlayLabel)=0;
                    overlayAreaSmall(overlayAreaSmall==overlayLabel)=1;
                    overlaySize = bwarea(overlayAreaSmall);
                end
                
                yStart=yEdgeFills(1);
                xStart=xEdgeFills(1);                
                NucAreaSmall = nucleusPic(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd); 
                SE = strel('disk', 3);
                NucAreaSmallDil = imdilate(NucAreaSmall,SE);
                
                NucleusMSub = NucleusM(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                DSub = D(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                IndSub = Ind(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);                
                AreaResultSmall = bwlabel(NucAreaSmallDil);
                               
                nucLabel = AreaResultSmall(yStart,xStart);
                if(nucLabel==0)
                    nucleusSize=-1;
                    AreaResultSmall(:)=0;
                    brightness=-1;
                    NucleusImageSmall=NucleusImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                else
                    AreaResultSmall(AreaResultSmall~=nucLabel)=0;
                    AreaResultSmall(AreaResultSmall==nucLabel)=1;      
                    NucleusImageSmall=NucleusImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                    brightness = CalculateNucleusBrightness(AreaResultSmall,NucleusImageSmall);
                    nucleusSize = bwarea(AreaResultSmall);
                end  
                
                %AreaResultSmall = FloodFillNucArea(NucAreaSmall,yStart,xStart,NucleusMSub,DSub,IndSub,sizeY,sizeX);                
                NucleusImageSmall=NucleusImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                if(optionHandler.FuzzFilterActivated > 2.5)
                    fuzzy = CheckForFuzzOnPosition(newNeurite.cutYPosStart+yEdgeFills(1)-1+yStartPos,newNeurite.cutXPosStart+xEdgeFills(1)-1+xStartPos,yFactor,xFactor,optionHandler,LowDensityNeurite, BrightDensity);
                end
                %AreaResultSmall = AreaResult(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                if((nucleusSize > 1 || nucleusSize == -1) && nucLabel ~=0 && overlaySize < optionHandler.MaxOverlaySize && fuzzy < 0.5 && CompetingM(newNeurite.cutYPosStart+yEdgeFills(1)-1,newNeurite.cutXPosStart+xEdgeFills(1)-1)==0)
                    AreaResultSmall = bwlabel(NucAreaSmall);
                    nucLabel = AreaResultSmall(yStart,xStart);
                    if(nucLabel==0)
                        nucleusSize=-1;
                        AreaResultSmall(:)=0;
                        brightness=-1;
                        NucleusImageSmall=NucleusImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                    else
                        AreaResultSmall(AreaResultSmall~=nucLabel)=0;
                        AreaResultSmall(AreaResultSmall==nucLabel)=1;      
                        NucleusImageSmall=NucleusImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                        brightness = CalculateNucleusBrightness(AreaResultSmall,NucleusImageSmall);
                        nucleusSize = bwarea(AreaResultSmall);
                    end        
                    if(SkelNeurons(newNeurite.cutYPosStart+yEdgeFills(1)-1,newNeurite.cutXPosStart+xEdgeFills(1)-1) == 0)
                        newNeurite.numberNuclei=newNeurite.numberNuclei+1;
                        newNeuriteBig.numberNuclei=newNeuriteBig.numberNuclei+1;
                        newNeuriteBig.nucleusImage(AreaResultSmall==1)=1;
                    end
                    SkelNeurons(newNeurite.cutYPosStart+yEdgeFills(1)-1,newNeurite.cutXPosStart+xEdgeFills(1)-1)=1;                    
                    BinaryImage(newNeurite.cutYPosStart+yStartPos-1:newNeurite.cutYPosEnd+yStartPos-1,newNeurite.cutXPosStart+xStartPos-1:newNeurite.cutXPosEnd+xStartPos-1) = logical(logical(BinaryImage(newNeurite.cutYPosStart+yStartPos-1:newNeurite.cutYPosEnd+yStartPos-1,newNeurite.cutXPosStart+xStartPos-1:newNeurite.cutXPosEnd+xStartPos-1)) + (1-currentNeuriteAreaSmall) + AreaResultSmall);         
                    neuriteAreaCount=neuriteAreaCount+nnz(1-currentNeuriteArea);   
                    %newNeuriteBig.nucleusImage = logical(newNeuriteBig.nucleusImage + logical(AreaResultSmall));
                    %newNeurite.nucleusImage = logical(newNeurite.nucleusImage + logical(AreaResultSmall));

                    if(NeuronManualM(newNeurite.cutYPosStart+yEdgeFills(1)+yStartPos-2,newNeurite.cutXPosStart+xEdgeFills(1)+xStartPos-2) == 1)
                        manualMarked=1;
                        TP = TP + 1;
                        TN = TN - 1;
                        TruePositives = [TruePositives;newNeurite.cutYPosStart+yEdgeFills(1)-1 newNeurite.cutXPosStart+xEdgeFills(1)-1];
                    else
                        manualMarked=0;
                        FP = FP + 1;
                    end
                    if(sav)
                        key = mapPositionToKeyString(newNeurite.cutYPosStart+yEdgeFills(1)-1,newNeurite.cutXPosStart+xEdgeFills(1)-1);
                        if(isKey(currentNeuronStat,key))
                            list = currentNeuronStat(key);
                        else
                            list = zeros(10,1);
                        end
                        list(1) = 1;
                        if(noNeurite == 0)
                            list(7) = newNeuriteBig.neuriteLength;
                        end
                        list(4) = manualMarked;
                        list(8) = maxDistEndpoints;                        
                        list(9) = brightness;
                        list(10) = nucleusSize;
                        list(11) = nucVariance;
                        currentNeuronStat(key) = list;
                    end
                end
            else
                markedForDeletion = ones(numel(yEdgeFills));
                for(counter=1:numel(yEdgeFills))     
                    %Check distance between different Neurons. %If distance
                    %too small, kill Neuron with farest distance to
                    %endpoint
                    
                    
                    %Check size of Neurite. If area > 1000, take all
                    %Neurons. Else only the one with highest overlap.
                    for(counter2=1:numel(yEdgeFills))
                        dist = pdist([yEdgeFills(counter) xEdgeFills(counter);yEdgeFills(counter2) xEdgeFills(counter2)]);
                       % key = mapPositionToKeyString(newNeurite.cutYPosStart+yEdgeFills(counter)-1+yStartPos-1,newNeurite.cutXPosStart+xEdgeFills(counter)-1+xStartPos-1);
                      %  list1 = currentNeuronStat(key);
                       % key = mapPositionToKeyString(newNeurite.cutYPosStart+yEdgeFills(counter2)-1+yStartPos-1,newNeurite.cutXPosStart+xEdgeFills(counter2)-1+xStartPos-1);
                      %  list2 = currentNeuronStat(key);
                        
                        if(dist<optionHandler.MinDistanceBetweenNeurons && dist > 0 && markedForDeletion(counter2) == 1)                            
                            %Check first, which has most overlap.
                           % if(list2(6) > list1(6))                                
                            %Mark for deletion
                                markedForDeletion(counter) = 0;
                           % else
                                %markedForDeletion(counter) = 0;
                           % end
                        end
                    end
                end                
               
                    startCou=1;
                    endeCou=numel(yEdgeFills);
                  
                      
                for(counter=startCou:endeCou)
                    NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
                    if(NeuronManualM(newNeurite.cutYPosStart+yEdgeFills(counter)+yStartPos-2,newNeurite.cutXPosStart+xEdgeFills(counter)+xStartPos-2) == 1)
                        manualMarked=1;
                    else
                        manualMarked=0;
                    end
                    yStart=yEdgeFills(counter);
                    xStart=xEdgeFills(counter);
                    NucleusMSub = NucleusM(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                    DSub = D(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                    IndSub = Ind(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                    NucAreaSmall = nucleusPic(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                    
                    SE = strel('disk', 3);
                    NucAreaSmallDil = imdilate(NucAreaSmall,SE);
                    AreaResultSmall = bwlabel(NucAreaSmallDil);
                    nucLabel = AreaResultSmall(yStart,xStart);
                    if(nucLabel==0)
                        nucleusSize=-1;
                        AreaResultSmall(:)=0;
                        brightness=-1;
                        NucleusImageSmall=NucleusImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                    else
                        AreaResultSmall(AreaResultSmall~=nucLabel)=0;
                        AreaResultSmall(AreaResultSmall==nucLabel)=1;      
                        NucleusImageSmall=NucleusImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                        brightness = CalculateNucleusBrightness(AreaResultSmall,NucleusImageSmall);
                        nucleusSize = bwarea(AreaResultSmall);
                    end
                    overlayAreaSmall = overlayArea(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                    overlayAreaSmall = bwlabel(overlayAreaSmall);
                    if(overlayAreaSmall(yEdgeFills(counter),xEdgeFills(counter)) == 0)
                        overlaySize=0;
                    else
                        overlayLabel = overlayAreaSmall(yEdgeFills(counter),xEdgeFills(counter));
                        overlayAreaSmall(overlayAreaSmall~=overlayLabel)=0;
                        overlayAreaSmall(overlayAreaSmall==overlayLabel)=1;
                        overlaySize = bwarea(overlayAreaSmall);
                    end
                    if(optionHandler.FuzzFilterActivated > 2.5)
                        fuzzy = CheckForFuzzOnPosition(newNeurite.cutYPosStart+yEdgeFills(counter)-1+yStartPos,newNeurite.cutXPosStart+xEdgeFills(counter)-1+xStartPos,yFactor,xFactor,optionHandler,LowDensityNeurite, BrightDensity);
                    end
                    if(markedForDeletion(counter) == 1 && fuzzy < 0.5)
                        %AreaResultSmall = AreaResult(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                        if((nucleusSize > 1 || nucleusSize==-1) && nucLabel ~=0 && overlaySize < optionHandler.MaxOverlaySize && CompetingM(newNeurite.cutYPosStart+yEdgeFills(counter)-1,newNeurite.cutXPosStart+xEdgeFills(counter)-1)==0)
                            AreaResultSmall = bwlabel(NucAreaSmall);
                            if(nucLabel==0)
                                nucleusSize=-1;
                                AreaResultSmall(:)=0;
                                brightness=-1;
                                NucleusImageSmall=NucleusImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                            else
                                AreaResultSmall(AreaResultSmall~=nucLabel)=0;
                                AreaResultSmall(AreaResultSmall==nucLabel)=1;      
                                NucleusImageSmall=NucleusImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                                brightness = CalculateNucleusBrightness(AreaResultSmall,NucleusImageSmall);
                                nucleusSize = bwarea(AreaResultSmall);
                            end
                            if(SkelNeurons(newNeurite.cutYPosStart+yEdgeFills(counter)-1,newNeurite.cutXPosStart+xEdgeFills(counter)-1)==0)
                               newNeurite.numberNuclei=newNeurite.numberNuclei+1;
                               newNeuriteBig.numberNuclei=newNeuriteBig.numberNuclei+1;  
                               newNeuriteBig.nucleusImage(AreaResultSmall==1)=1;
                            end
                            SkelNeurons(newNeurite.cutYPosStart+yEdgeFills(counter)-1,newNeurite.cutXPosStart+xEdgeFills(counter)-1)=1;    
%                             newNeuriteBig.nucleusImage = logical(newNeuriteBig.nucleusImage + AreaResultSmall);
%                             newNeurite.nucleusImage = logical(newNeurite.nucleusImage + AreaResultSmall);

                            if(manualMarked == 1)
                                TP=TP+1;
                                TN = TN - 1;
                                TruePositives = [TruePositives;newNeurite.cutYPosStart+yEdgeFills(counter)-1 newNeurite.cutXPosStart+xEdgeFills(counter)-1];
                            else
                                FP=FP+1;
                            end
                            %Add position of underlying nucleus to imageHandler.BinaryImage

                            BinaryImage(newNeurite.cutYPosStart+yStartPos-1:newNeurite.cutYPosEnd+yStartPos-1,newNeurite.cutXPosStart+xStartPos-1:newNeurite.cutXPosEnd+xStartPos-1) = logical(logical(BinaryImage(newNeurite.cutYPosStart+yStartPos-1:newNeurite.cutYPosEnd+yStartPos-1,newNeurite.cutXPosStart+xStartPos-1:newNeurite.cutXPosEnd+xStartPos-1)) + (1-currentNeuriteAreaSmall) + AreaResultSmall);
                            neuriteAreaCount=neuriteAreaCount+nnz(1-currentNeuriteArea);
                            if(sav)
                                key = mapPositionToKeyString(newNeurite.cutYPosStart+yEdgeFills(counter)-1,newNeurite.cutXPosStart+xEdgeFills(counter)-1);
                                if(isKey(currentNeuronStat,key))
                                    list = currentNeuronStat(key);
                                else
                                    list = zeros(10,1);
                                end
                                list(1) = 2;
                                if(noNeurite == 0)
                                    list(7) = newNeuriteBig.neuriteLength;
                                end
                                list(8) = maxDistEndpoints;
                                list(4) = manualMarked;
                                list(9) = brightness;
                                list(10) = nucleusSize;
                                list(11) = nucVariance;
                                currentNeuronStat(key) = list;
                            end
                        end
                    elseif(sav)
                        key = mapPositionToKeyString(newNeurite.cutYPosStart+yEdgeFills(counter)-1,newNeurite.cutXPosStart+xEdgeFills(counter)-1);
                        if(isKey(currentNeuronStat,key))
                            list = currentNeuronStat(key);
                        else
                            list = zeros(10,1);
                        end
                        list(5) = 2;
                        if(noNeurite == 0)
                            list(7) = newNeuriteBig.neuriteLength;
                        end
                        list(8) = maxDistEndpoints;
                        list(4) = manualMarked;
                        list(9) = brightness;
                        list(10) = nucleusSize;
                        list(11) = nucVariance;
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
            edgeFillPositionsSecond = edgeFillPositionsSecond(newNeurite.cutYPosStart+yStartPos-1:newNeurite.cutYPosEnd+yStartPos-1,newNeurite.cutXPosStart+xStartPos-1:newNeurite.cutXPosEnd+xStartPos-1);
            %[maxA,ind] = max(A(:));
            [yInd xInd] = find(edgeFillPositionsSecond);
            maxOverlap = 0;
            maxInd = 0;
            for(i=1:numel(yInd))
                key = mapPositionToKeyString(newNeurite.cutYPosStart+yInd(i)-1+yStartPos-1,newNeurite.cutXPosStart+xInd(i)-1+xStartPos-1);
                if(isKey(currentNeuronStat,key))
                   list = currentNeuronStat(key);
                   if(list(6) > maxOverlap)
                       maxOverlap = list(6);
                       maxInd = i;
                   end
                end
            end          
            
            if(maxInd > 0 && maxOverlap > optionHandler.EdgeFillNucleusAreaWithinNeuriteSecond)
                %Add Neuron to SkelNeurons
                yStart=yInd(maxInd);
                xStart=xInd(maxInd);
                NucleusMSub = NucleusM(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                DSub = D(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                IndSub = Ind(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                NucAreaSmall = nucleusPic(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd); 
                SE = strel('disk', 3);
                NucAreaSmallDil = imdilate(NucAreaSmall,SE);
                NucleusImageSmall=NucleusImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                AreaResultSmall = bwlabel(NucAreaSmallDil);
                nucLabel=AreaResultSmall(yStart,xStart);
                
                if(overlayAreaSmall(yInd(1),xInd(1)) == 0)
                    overlaySize=0;
                else
                    overlayLabel = overlayAreaSmall(yInd(1),xInd(1));
                    overlayAreaSmall(overlayAreaSmall~=overlayLabel)=0;
                    overlayAreaSmall(overlayAreaSmall==overlayLabel)=1;
                    overlaySize = bwarea(overlayAreaSmall);
                end
                
                
                if(nucLabel==0)
                    nucleusSize=-1;
                    AreaResultSmall(:)=0;
                    brightness=-1;
                    NucleusImageSmall=NucleusImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                else
                    AreaResultSmall(AreaResultSmall~=nucLabel)=0;
                    AreaResultSmall(AreaResultSmall==nucLabel)=1;      
                    NucleusImageSmall=NucleusImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                    brightness = CalculateNucleusBrightness(AreaResultSmall,NucleusImageSmall);
                    nucleusSize = bwarea(AreaResultSmall);
                end
                if(optionHandler.FuzzFilterActivated > 2.5)
                    fuzzy = CheckForFuzzOnPosition(newNeurite.cutYPosStart+yInd(1)-1+yStartPos,newNeurite.cutXPosStart+xInd(1)-1+xStartPos,yFactor,xFactor,optionHandler,LowDensityNeurite, BrightDensity);
                end
                if((nucleusSize > 0 || nucleusSize==-1) && nucLabel ~=0 && overlaySize < optionHandler.MaxOverlaySize && fuzzy < 0.5 && CompetingM(newNeurite.cutYPosStart+yInd(1)-1, newNeurite.cutXPosStart+xInd(1)-1)==0)
                    AreaResultSmall = bwlabel(NucAreaSmall);
                    nucLabel = AreaResultSmall(yStart,xStart);
                    if(nucLabel==0)
                        nucleusSize=-1;
                        AreaResultSmall(:)=0;
                        brightness=-1;
                        NucleusImageSmall=NucleusImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                    else
                        AreaResultSmall(AreaResultSmall~=nucLabel)=0;
                        AreaResultSmall(AreaResultSmall==nucLabel)=1;      
                        NucleusImageSmall=NucleusImage(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                        brightness = CalculateNucleusBrightness(AreaResultSmall,NucleusImageSmall);
                        nucleusSize = bwarea(AreaResultSmall);
                    end
                    if(SkelNeurons(newNeurite.cutYPosStart+yInd(1)-1, newNeurite.cutXPosStart+xInd(1)-1) ==0)
                        newNeurite.numberNuclei=newNeurite.numberNuclei+1;
                        newNeuriteBig.numberNuclei=newNeuriteBig.numberNuclei+1;
                        newNeuriteBig.nucleusImage(AreaResultSmall==1)=1;
                    end
                    
                    SkelNeurons(newNeurite.cutYPosStart+yInd(1)-1, newNeurite.cutXPosStart+xInd(1)-1) = 1;
%                     newNeuriteBig.nucleusImage = logical(newNeuriteBig.nucleusImage + AreaResultSmall);
%                     newNeurite.nucleusImage = logical(newNeurite.nucleusImage + AreaResultSmall);
                    %NucArea = nucleusPic(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                    BinaryImage(newNeurite.cutYPosStart+yStartPos-1:newNeurite.cutYPosEnd+yStartPos-1,newNeurite.cutXPosStart+xStartPos-1:newNeurite.cutXPosEnd+xStartPos-1) = logical(logical(BinaryImage(newNeurite.cutYPosStart+yStartPos-1:newNeurite.cutYPosEnd+yStartPos-1,newNeurite.cutXPosStart+xStartPos-1:newNeurite.cutXPosEnd+xStartPos-1)) + (1-currentNeuriteAreaSmall) + AreaResultSmall);
                    neuriteAreaCount=neuriteAreaCount+nnz(1-currentNeuriteArea);
                    NeuronPositionsEdgeFillNeurite(newNeurite.cutYPosStart+yInd(1)-1, newNeurite.cutXPosStart+xInd(1)-1) = 1;
                    brightness = CalculateNucleusBrightness(AreaResultSmall,NucleusImageSmall);
                    %Save Neurite length in Matrix
                   % if(noNeurite==0)
                   %     neuronHandler.NeuriteLengthMatrix(newNeurite.cutYPosStart+yInd(1)-1, newNeurite.cutXPosStart+xInd(1)-1) = newNeuriteBig.neuriteLength;
                   % end
                    NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
                    if(NeuronManualM(newNeurite.cutYPosStart+yInd(1)+yStartPos-2,newNeurite.cutXPosStart+xInd(1)+xStartPos-2) == 1)
                        manualMarked=1;
                        TP = TP + 1;
                        TN = TN - 1;
                        TruePositives = [TruePositives;newNeurite.cutYPosStart+yInd(1)-1 newNeurite.cutXPosStart+xInd(1)-1];
                    else
                        manualMarked=0;
                        FP = FP + 1;
                    end
                    if(sav)
                        key = mapPositionToKeyString(newNeurite.cutYPosStart+yInd(maxInd)-1,newNeurite.cutXPosStart+xInd(maxInd)-1);
                        if(isKey(currentNeuronStat,key))
                           list = currentNeuronStat(key);
                        else
                           list = zeros(8,1);
                        end
                        list(1) = 2;
                        if(noNeurite==0)
                            list(7) = newNeuriteBig.neuriteLength;
                        end
                        list(8) = maxDistEndpoints;
                        list(4) = manualMarked;
                        list(9) = brightness;
                        list(10) = nucleusSize;
                        list(11) = nucVariance;
                        currentNeuronStat(key) = list;
                    end
                end
            elseif(noNeurite==0) 
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
                for(i=1:numel(newNeurite.xShortedEndpoints))
                    %ToDo: Add criterias direction and overlap!
                    %Berechne Geradendefinition zwischen Centroid und Endpunkt
                    u = [yCen xCen];
                    v = [newNeurite.yShortedEndpoints(i) newNeurite.xShortedEndpoints(i)];
                    rVec = v-u;
                    
                    %CosTheta = dot(u,v)/(norm(u)*norm(v));
                    %ThetaInDegrees = acos(CosTheta)*180/pi
                    
                    %Set startPoint ~6 px back to balance point of neurite                    
                    startPoint = [newNeurite.yShortedEndpoints(i) newNeurite.xShortedEndpoints(i)];
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
                    %[nucCol nucRow euclidDist] = FindNucleusForNeuron(newNeurite.cutXPosStart+startPoint(2)+double(xStartPos-1), newNeurite.cutYPosStart+startPoint(1)+double(yStartPos-1), NucleusM, optionHandler.MaxDistanceFromEndpoint, rVec, optionHandler.ToleranceAngleFromEndpoint);
                    [nucCol nucRow euclidDist] = FindNucleusForNeuronNew(newNeurite.cutXPosStart+startPoint(2)+double(xStartPos-1), newNeurite.cutYPosStart+startPoint(1)+double(yStartPos-1), NucleusM, optionHandler.MaxDistanceFromEndpoint, rVec, optionHandler.ToleranceAngleFromEndpoint);
                    nucCol = nucCol - (xStartPos-1);
                    nucRow = nucRow - (yStartPos-1);
                    if(euclidDist<minDist)
                      minDist = euclidDist;
                      nucColSaved = nucCol;
                      nucRowSaved = nucRow;
                    end
                end
                if(nucRowSaved <= 0)
                   %Check also ot killed saved endpoints.
                   for(i=1:numel(newNeuriteBig.yShortedEndpoints))
                    %Set startPoint 6 px back to balance point of neurite
                    startPoint = [newNeuriteBig.yShortedEndpoints(i) newNeuriteBig.xShortedEndpoints(i)];
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
                    u = [yCen xCen];
                    v = [newNeuriteBig.yShortedEndpoints(i) newNeuriteBig.xShortedEndpoints(i)];
                    rVec = v-u;
                    %[nucCol nucRow euclidDist] = FindNucleusForNeuron(newNeuriteBig.cutXPosStart+startPoint(2)+double(xStartPos-1), newNeuriteBig.cutYPosStart+startPoint(1)+double(yStartPos-1), NucleusM, optionHandler.MaxDistanceFromEndpoint, rVec, optionHandler.ToleranceAngleFromEndpoint);
                    [nucCol nucRow euclidDist] = FindNucleusForNeuronNew(newNeuriteBig.cutXPosStart+startPoint(2)+double(xStartPos-1), newNeuriteBig.cutYPosStart+startPoint(1)+double(yStartPos-1), NucleusM, optionHandler.MaxDistanceFromEndpoint, rVec, optionHandler.ToleranceAngleFromEndpoint);
                    nucCol = nucCol - (xStartPos-1);
                    nucRow = nucRow - (yStartPos-1);
                    if(euclidDist<minDist)
                      minDist = euclidDist;
                      nucColSaved = nucCol;
                      nucRowSaved = nucRow;
                    end
                   end
                end
                if(nucRowSaved > 0 && nucColSaved > 0 && nucRowSaved < size(NeuronPositionsEdgeFillNeurite,1) && nucColSaved < size(NeuronPositionsEdgeFillNeurite,2))
                    NeuronPositionsEdgeFillNeurite(nucRowSaved, nucColSaved) = 1;
                    neuriteAreaCount=neuriteAreaCount+nnz(1-currentNeuriteArea);
                    %NucArea = nucleusPic(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);  
                     SE = strel('disk', 3);
                     nucleusPicDil = imdilate(nucleusPic,SE);
                    AreaResult = FloodFillNucArea(nucleusPicDil,nucRowSaved,nucColSaved,NucleusM,D,Ind,sizeY,sizeX);
                    AreaResultSmall = AreaResult;%(yStartPos:yEndPos,xStartPos:xEndPos);
                    NucleusImageSmall=NucleusImage(yStartPos:yEndPos,xStartPos:xEndPos);
                    nucleusSize = bwarea(AreaResultSmall);
                    overlayAreaSmall=overlayArea;
                    overlayAreaSmall = bwlabel(overlayAreaSmall);
                    if(overlayArea(nucRowSaved,nucColSaved) == 0)
                        overlaySize=0;
                    else
                        overlayLabel = overlayAreaSmall(nucRowSaved,nucColSaved);
                        overlayAreaSmall(overlayArea~=overlayLabel)=0;
                        overlayAreaSmall(overlayArea==overlayLabel)=1;
                        overlaySize = bwarea(overlayAreaSmall);
                    end
                    if(optionHandler.FuzzFilterActivated > 2.5)
                        fuzzy = CheckForFuzzOnPosition(nucRowSaved+yStartPos,nucColSaved+xStartPos,yFactor,xFactor,optionHandler,LowDensityNeurite, BrightDensity);
                    end
                    if(nucleusSize > 1 && overlaySize < optionHandler.MaxOverlaySize && fuzzy < 0.5 && CompetingM(nucRowSaved,nucColSaved)==0)
                        AreaResult = FloodFillNucArea(nucleusPic,nucRowSaved,nucColSaved,NucleusM,D,Ind,sizeY,sizeX);
                        AreaResultSmall = AreaResult;%(yStartPos:yEndPos,xStartPos:xEndPos);
                        NucleusImageSmall=NucleusImage(yStartPos:yEndPos,xStartPos:xEndPos);
                        nucleusSize = bwarea(AreaResultSmall);
                        brightness = CalculateNucleusBrightness(AreaResultSmall,NucleusImageSmall);
                        if(SkelNeurons(nucRowSaved,nucColSaved)==0)
                            newNeurite.numberNuclei=newNeurite.numberNuclei+1;
                            newNeuriteBig.numberNuclei=newNeuriteBig.numberNuclei+1;
                            %newNeurite.nucleusImage(AreaResultSmall==1)=1;
                            if(yStartPos==1 && xStartPos == 1 && yEndPos == sizeYorig && xEndPos == sizeXorig)
                                AreaResultSmallArea = AreaResultSmall(newNeurite.cutYPosStart+yStartPos-1:newNeurite.cutYPosEnd+yStartPos-1,newNeurite.cutXPosStart+xStartPos-1:newNeurite.cutXPosEnd+xStartPos-1);
                            else
                                AreaResultSmallArea = AreaResultSmall(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd);
                            end                            
                            newNeuriteBig.nucleusImage = logical(newNeuriteBig.nucleusImage + AreaResultSmallArea);
                        end
                        SkelNeurons(nucRowSaved,nucColSaved)=1;                          
                        BinaryImage(newNeurite.cutYPosStart+yStartPos-1:newNeurite.cutYPosEnd+yStartPos-1,newNeurite.cutXPosStart+xStartPos-1:newNeurite.cutXPosEnd+xStartPos-1) = logical(logical(BinaryImage(newNeurite.cutYPosStart+yStartPos-1:newNeurite.cutYPosEnd+yStartPos-1,newNeurite.cutXPosStart+xStartPos-1:newNeurite.cutXPosEnd+xStartPos-1)) + (1-currentNeuriteAreaSmall));
                        BinaryImage(yStartPos:yEndPos,xStartPos:xEndPos) = logical(logical(BinaryImage(yStartPos:yEndPos,xStartPos:xEndPos)) +AreaResultSmall);
                        
%                         newNeuriteBig.nucleusImage = logical(newNeuriteBig.nucleusImage + AreaResultSmallArea);
%                         
%                         newNeurite.numberNuclei = newNeurite.numberNuclei+1;
                    %Save Neurite length in Matrix
                    %    if(noNeurite==0)
                    %        neuronHandler.NeuriteLengthMatrix(nucRowSaved, nucColSaved) = newNeuriteBig.neuriteLength;                    
                    %    end
                        NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
                        if(NeuronManualM(nucRowSaved+yStartPos-1,nucColSaved+xStartPos-1) == 1)
                            TruePositives = [TruePositives;nucRowSaved nucColSaved];
                            manualMarked=1;
                            TP = TP + 1;
                            TN = TN - 1;
                        else
                            manualMarked=0;
                            FP = FP + 1;
                        end
                        if(sav)
                            key = mapPositionToKeyString(nucRowSaved,nucColSaved);
                            if(isKey(currentNeuronStat,key))
                                list = currentNeuronStat(key);
                            else
                                list = zeros(10,1);
                            end
                            list(1) = 3;
                            if(noNeurite==0)
                                list(7) = newNeuriteBig.neuriteLength;
                            end
                            list(8) = maxDistEndpoints;
                            list(4) = manualMarked;
                            list(9) = brightness;
                            list(10) = nucleusSize;
                            list(11) = nucVariance;
                            currentNeuronStat(key) = list;
                        end
                    end
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
        
        if(AnalyzeNeurites > 0 && newNeuriteBig.numberNuclei > 0 && noNeurite == 0)
        
        %Evtl: Do this all in Neurite class.
        %ToDo: Neurite Length calculation
        %1. Get number of nuclei
        %if(newNeurite.numberNuclei==1)
            %Just start without special cases and deal only with cases with one
            %nucleus.
            %1. Check if Nucleus is in mid or at border of Neurite
            %For that: Subtract Nucleus from Neurite and check if Neurite
            %is split up in two different areas. If yes, check if both
            %areas have at least (MinNeuriteLength).
            oldSkel = newNeuriteBig.skeletonImage;
            
            %figure(2);
            %imshow(oldSkel);
            
            bothIndices = find(logical(newNeuriteBig.skeletonImage) & logical(newNeuriteBig.nucleusImage));
            
            %Doppelt gemoppelt h�lt besser :-D
            newNeuriteBig.FindBranchingCrossingPoints2();           
            newNeuriteBig.FindBranchingCrossingPoints2(); 
            %Finde �berdeckungsl�nge zwischen Nucleus und Neuritenbild
            
            %minNeuriteLength = optionHandler.SkeletonMinNeuriteLength - numel(bothIndices);
            %if(minNeuriteLength < 11)
                minNeuriteLength=11;
            %end
            [sizeYskel sizeXskel] = size(newNeuriteBig.skeletonImage);
            subtractImage = newNeuriteBig.skeletonImage - newNeuriteBig.nucleusImage;
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
            newNeuriteBig.skeletonImage(labelSubtract==maxIndex)=1;
            subtractImage(labelSubtract==maxIndex) = 1;
            %ToDo: Ensure that at least one subtractImageLabel remains and 
            %ensure that bothIndices are counted before finding new
            %BranchingCrossingPoints.
            
            
            %Do this only if there is at least no Nucleus Point, touching a
            %Neurite Point
            
            [DSubtract,IDX] = bwdist(subtractImage);
            [DNuc,IDXNuc] = bwdist(newNeuriteBig.nucleusImage);
            neuriteCounter=neuriteCounter+1;
            
            %Get minimum of DSubtract on Indices where
            %newNeuriteBig.nucleusImage is 1
            minDist = min(DSubtract(newNeuriteBig.nucleusImage==1));
            if(minDist >= 2)
            
                %1. Get mid of Nucleus and find shortest path to next Neurite
                [Ay Ax]= find(newNeuriteBig.nucleusImage>=1);
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
                end
            end
            newNeuriteBig.skeletonImage = subtractImage;
            %newNeuriteBig.FindBranchingCrossingPoints2();           
            
            
            
            %Create figure with overlayed Nucleus, Skeleton and original
            %Neurite in /SkeletonExport/selectedWell/ID.png
            
            currentObject = zeros(9,20);
            newNeuriteBig.skeletonImage = parsiSkel(newNeuriteBig.skeletonImage);
            branchImage = logical(bwmorph(newNeuriteBig.skeletonImage, 'branchpoints'));
            endpointImage = logical(bwmorph(newNeuriteBig.skeletonImage, 'endpoints'));
            %2 Endpoints: No branches, maybe just circles
            newNeuriteBig.xBranchingPoints = zeros(0);
            newNeuriteBig.yBranchingPoints = zeros(0);
            [yBranchingPoints xBranchingPoints] = find(branchImage==1);

            newNeuriteBig.skeletonImage = parsiSkel(newNeuriteBig.skeletonImage);
            
            [yBranchingPoints xBranchingPoints] = KillCloseBranches(yBranchingPoints,xBranchingPoints);
   
            newNeuriteBig.yBranchingPoints = yBranchingPoints;
            newNeuriteBig.xBranchingPoints = xBranchingPoints; 
            
            %figure(3);
            %imshow(imfuse(newNeuriteBig.skeletonImage,newNeuriteBig.nucleusImage));
            exportRGB = zeros(sizeYskel,sizeXskel,3);
            exportRGB(:,:,3) = 1-newNeuriteBig.image;
            exportRGB(:,:,1) = newNeuriteBig.nucleusImage;
            exportRGB(:,:,2) = newNeuriteBig.skeletonImage;
            
            
            [rowIndices colIndices] = find(currentNeuriteArea);
            yIndex = min(rowIndices);
            xIndex = min(colIndices);
            cutYPosStartRef = min(rowIndices);
            cutXPosStartRef = min(colIndices);
            positionID = sub2ind([sizeY sizeX],cutYPosStartRef, cutXPosStartRef);
            
            
            %figure(1);
            %imshow(exportRGB);
            imwrite(exportRGB,[foldername '/' selectedWell '/' num2str(neuriteCounter) '.png']);
            connections=SplitUpNeuronsAndNuclei(SubNucM,newNeuriteBig.nucleusImage,newNeuriteBig.skeletonImage);
            labelSubtract = bwlabel(newNeuriteBig.skeletonImage);   
            for(t=1:size(connections,1))
                %Save to Excel:
                %Values for Sub Neurite
                %0. Neurite ID
                currentObject(1,1)=connections(t,1);
                %1. Subneurite Length
                currentObject(2,1)=connections(t,2);
                %2. Subeurite Branching Points
                currentObject(3,1)=connections(t,3);
                %3. Comma separated Nucleus IDs
                currentObject(4,1:end)=connections(t,1:end);
                
                %3. Values for whole Neurite
                %1. Number Nuclei
                currentObject(5,1) = newNeuriteBig.numberNuclei;
                %2. Number total Branching Points
                currentObject(6,1) = numel(newNeuriteBig.yBranchingPoints);
                %3. Neurite length
                currentObject(7,1) = nnz(newNeuriteBig.skeletonImage);
                %4. Number of distinct Neurons per Skeleton
                currentObject(8,1) = max(labelSubtract(:));
                currentObject(9,1) = positionID;
                
                 %5. Export skeleton picture with ID
                 neuriteLengthMatrix([selectedWell ';' num2str(neuriteCounter) ';' num2str(t)]) = currentObject;
            end
        end
       % elseif(newNeurite.numberNuclei>1)
       % end        
        %2. Get neurite length minus nucleus positions
        %3. Neurite length medium = neurite length / number nuclei
        
        %Write result into dictionary
    end

    %newNeurite.PlotSkeleton();    
    %close(figure(1));
    %neuriteAreasDict(num2str(i)) = newNeurite;
    i=i+1;
    %if(mod(i,10) == 0)
        
    %end
%     while(numel(whitePointRows) > 0 && alreadyVisited(whitePointRows(1), whitePointCols(1)) == 0)
%         whitePointRows(1) = [];
%         whitePointCols(1) = [];
%     end   
    %end
end

imagePathBinary = [foldername '/' selectedWellLong 'Binary' optionHandler.FileEnding];
imagePathSkeleton = [foldername '/ConvertedCellomics/' selectedWellLong 'Skeleton' optionHandler.FileEnding];
pathNeuriteLength = [foldername '/ConvertedCellomics/NeuriteLength-' selectedWell '.mat'];
%neuriteLengthMatrix = 0;
%neuriteLengthMatrix = neuronHandler.NeuriteLengthMatrix;
%save(pathNeuriteLength,'neuriteLengthMatrix');
imwrite(full(BinaryImage), imagePathBinary);
imwrite(full(SkeletonImage), imagePathSkeleton);
%neuronHandler.NeuriteLengthMatrix = 0;

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

%Check if Neuronhandler already has Dictionary of NeuronPositionsEdgeComposite
if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
if(~isKey(neuronHandler.NeuronPositionsSkeletonization, selectedWell))
    neuronHandler.NeuronPositionsSkeletonization(selectedWell)= sparse(logical(zeros(sizeYorig,sizeXorig)));
    fullMatrix = neuronHandler.NeuronPositionsSkeletonization(selectedWell);
    fullMatrix(yStartPos:yEndPos,xStartPos:xEndPos) = SkelNeurons(1:yEndPos-yStartPos+1,1:xEndPos-xStartPos+1);
    neuronHandler.NeuronPositionsSkeletonization(selectedWell) = fullMatrix;
else
    fullMatrix = neuronHandler.NeuronPositionsSkeletonization(selectedWell);
    fullMatrix(yStartPos:yEndPos,xStartPos:xEndPos) = SkelNeurons(1:yEndPos-yStartPos+1,1:xEndPos-xStartPos+1);
    neuronHandler.NeuronPositionsSkeletonization(selectedWell) = fullMatrix;
end

%if(numel(neuronHandler.NeuriteLengthMatrix) ==0 || ~isKey(neuronHandler.NeuriteLengthMatrix, selectedWell))
    neuronHandler.NeuriteLengthMatrix = containers.Map();
    neuronHandler.NeuriteLengthMatrix(selectedWell) = neuriteLengthMatrix;
%else
%    neuronHandler.NeuriteLengthMatrix(selectedWell) = neuriteLengthMatrix;
%end


%handles.NeuronCoordinates = neuronHandler;
%guidata(handles.figure1, handles);
if(numel(neuronHandler.AreaDictionary) ==0)
    neuronHandler.AreaDictionary = containers.Map();
end
neuronHandler.AreaDictionary(selectedWell)=neuriteAreaCount;

if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
    neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
end
if(isKey(neuronHandler.NeuronPositionsSkelDeleted, selectedWell))
    fullMatrix = neuronHandler.NeuronPositionsSkelDeleted(selectedWell);
else
    fullMatrix = zeros(sizeYorig,sizeXorig);
end
fullMatrix(yStartPos:yEndPos,xStartPos:xEndPos) = SkelDeletedMat;
neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = fullMatrix;
%neuronHandler.NeuronPositionsSkeletonization(selectedWell) = SkelNeurons;
neuronHandler.NeuronStatMatrix(selectedWell) = currentNeuronStat;

%if(sizeY == sizeYorig && sizeX==sizeXorig)
    [TP FP TN] = CalculateQuality(foldername,neuronHandler,selectedWell,sizeYorig,sizeXorig,yStartPos,yEndPos,xStartPos,xEndPos,optionHandler,NucleusM);    
%end


function [yBranchingPoints xBranchingPoints] = KillCloseBranches(yBranchingPoints,xBranchingPoints)
%Again: Kill branches which are close together.
xMid = zeros(0);
yMid = zeros(0);
killIndicesList = zeros(0);
newPointsListY = zeros(0);
newPointsListX = zeros(0);

%Solange es zwei Branching Points mit pdist < 30 gibt, f�hre folgenden
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


function connections=SplitUpNeuronsAndNuclei(subNucM,nucleusImage,neuriteImage,newNeuriteBig,sizeYorig,sizeXorig)
labeledNuclei = bwlabel(nucleusImage);
labeledNeurites = bwlabel(neuriteImage);
[sizeY sizeX] = size(nucleusImage);
[D IDX2] = bwdist(nucleusImage);
[D2 IDX] = bwdist(full(subNucM));
%Now check for connections between Neurites and Nuclei
connections=zeros(max(labeledNeurites(:)),20);
for(i=1:max(labeledNeurites(:)))
    currentNeurite = zeros(sizeY,sizeX);
    currentNeurite(labeledNeurites==i)=1;
    connections(i,1)=i;
    %Calculate Subneurite length
    connections(i,2) = numel(labeledNeurites(labeledNeurites==i));
    %Calculate Subneurite number of Branching Points
    branchImage = logical(bwmorph(currentNeurite, 'branchpoints'));
    [yBranchingPoints xBranchingPoints] = find(branchImage>0);
    [yBranchingPoints xBranchingPoints] = KillCloseBranches(yBranchingPoints,xBranchingPoints);
    connections(i,3) = numel(yBranchingPoints);
    %For every pixel in current Nucleus, get all Neurite pixels with
    %distance 1
    nucleusCount=4;
    distanceOneIndices = find(currentNeurite & D <= 1.5);
    for(k=1:numel(distanceOneIndices))
        nucleusIndex = IDX(distanceOneIndices(k));
        nucleusIndexGlobal = IDX(distanceOneIndices(k));
        [idxGlobalY idxGlobalX] = ind2sub([sizeY sizeX],nucleusIndexGlobal);
        idxGlobalY = idxGlobalY + newNeuriteBig.cutYPosStart;
        idxGlobalX = idxGlobalX + newNeuriteBig.cutXPosStart;
        nucleusIndexGlobal = sub2ind([sizeYorig sizeXorig], idxGlobalY, idxGlobalX);
        nucleusID = nucleusIndexGlobal;
        %nucleusID = labeledNuclei(nucleusIndex);
        if(~any(connections(i,4:end)==nucleusID))
            connections(i,nucleusCount) = nucleusID;
            nucleusCount=nucleusCount+1;
        end
    end    
end

function  variance=CheckForFuzzAround(subNucImage,neuImage,BrightDensity,LowDensity,x,y,sizeY,sizeX)
%1. Split up Nucleus image in different areas
%2. Calculate for each area the histogram
%3. Calculate for each histogram the variance


    

function [medBrightness] = CalculateNucleusBrightness(NucAreaBinary, NucAreaFull)
    indices = find(NucAreaBinary);    
    medBrightness = mean(NucAreaFull(indices));
    
        

%Get For NeuronPosition col, row nearest Nucleus in NucleusM
function [nucCol nucRow euclidDist] = FindNucleusForNeuron(neuCol, neuRow, NucleusM, maxDist, rVec, tolerance)
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
               %Check if point is in bounds, if applicable
               if(tolerance ~=0)
                   %Calculate rVec                       
                   u=[neuRow neuCol];
                   v=[i j];
                   rVec2 = v-u;
                   %Calculate angle
                   CosTheta = dot(rVec,rVec2)/(norm(rVec)*norm(rVec2));
                   ThetaInDegrees = acos(CosTheta)*180/pi;
               end
               if(((tolerance == 0) || (ThetaInDegrees > -tolerance && ThetaInDegrees < tolerance)) && i>0 && j<=sizeX && i<=sizeY && j>0)
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

%Get For NeuronPosition col, row nearest Nucleus in NucleusM
function [nucCol nucRow euclidDist] = FindNucleusForNeuronNew(neuCol, neuRow, NucleusM, maxDist, rVec, tolerance)
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
               %Check if point is in bounds, if applicable
               if(i>0 && j<=sizeX && i<=sizeY && j>0 && NucleusM(i,j) == 1)
                   if(tolerance ~=0)
                       %Calculate rVec                       
                       u=[neuRow neuCol];
                       v=[i j];
                       rVec2 = v-u;
                       %Calculate angle
                       CosTheta = dot(rVec,rVec2)/(norm(rVec)*norm(rVec2));
                       ThetaInDegrees = acos(CosTheta)*180/pi;
                   end
                   if(((tolerance == 0) || (ThetaInDegrees > -tolerance && ThetaInDegrees < tolerance)) && i>0 && j<=sizeX && i<=sizeY && j>0)                            
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


% --------------------------------------------------------------------
function SkeletonizationNeurons_Callback(hObject, eventdata, handles)
% hObject    handle to SkeletonizationNeurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
waitbarHandle = waitbar(0.5,'Executing Skeletonization. This could take up to 45 minutes.', 'WindowStyle', 'modal');
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
csvHandler = handles.CSVCoordinates;
foldername = imageHandler.Foldername;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
[sizeY sizeX] = size(imageHandler.NeuriteImage);

NucleusM = csvHandler.CellPosMatrix(selectedWell);
if(imageHandler.ZoomState)
    xMin = uint32(imageHandler.ZoomState(1));
    yMin = uint32(imageHandler.ZoomState(2));
    xMax = uint32(imageHandler.ZoomState(3) + xMin);
    yMax = uint32(imageHandler.ZoomState(4) + yMin);
     imageHandler.NucleusImage = 0;
     imageHandler.NeuriteImage = 0;
    [neuronHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell, yMin, yMax, xMin, xMax, 0, optionHandler, neuronHandler, NucleusM,foldername,'Neurite');
    [neuronHandler a b c] = SkeletonizationNeuronsAlgorithm(selectedWell,yMin, yMax, xMin, xMax, 0,optionHandler, sizeY,sizeX,foldername,neuronHandler,NucleusM,0,'Neurite',0);
     handles.NeuronCoordinates = neuronHandler;
    guidata(handles.figure1, handles);
    ChangeWell(selectedWellNumber,handles);
    RefreshZoomedImage(handles);
else
    %Check size of image handler. If image is bigger than 10.000 pixels,
    %split it up into four pics
    if(size(imageHandler.NucleusImage,1) > 10000 || size(imageHandler.NucleusImage,2) > 10000)        
        sizeXsmall = sizeX/2;
        sizeYsmall = sizeY/2;
        imageHandler.NucleusImage = 0;
        imageHandler.NeuriteImage = 0;
        [neuronHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell,1,sizeYsmall,1,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,'Neurite');
        [neuronHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell,1,sizeYsmall,sizeXsmall+1,sizeX,1,optionHandler,neuronHandler,NucleusM,foldername,'Neurite');
        [neuronHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell,sizeYsmall+1,sizeY,sizeXsmall+1,sizeX,1,optionHandler,neuronHandler,NucleusM,foldername,'Neurite');
        [neuronHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell,sizeYsmall+1,sizeY,1,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,'Neurite');
    else
         imageHandler.NucleusImage = 0;
         imageHandler.NeuriteImage = 0;
        [neuronHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell,1,sizeY,1,sizeX,1,optionHandler,neuronHandler,NucleusM,foldername,'Neurite');
    end
    [neuronHandler a b c] = SkeletonizationNeuronsAlgorithm(selectedWell,1, sizeY, 1, sizeX, 1,optionHandler, sizeY,sizeX,foldername,neuronHandler,NucleusM,1,'Neurite',0);    
     handles.NeuronCoordinates = neuronHandler;
     guidata(handles.figure1, handles);
    
    ChangeWell(selectedWellNumber,handles);
    RefreshUnzoomedImage(handles);
end
handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles);

catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
close(waitbarHandle);

% --------------------------------------------------------------------
function SkeletonizeNeurites(BW, handles)
%Next steps:
%1. Cut out all Single Neurite areas for themselfs
imageHandler = handles.ImageHandler;
[whitePointRows, whitePointCols] = find(BW);
alreadyVisited = logical(ones(size(BW,1), size(BW,2)));
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
                [nucCol nucRow euclidDist] = csvHandler.FindNucleusForNeuron(neurite.cutXPosStart+startPoint(2), neurite.cutYPosStart+startPoint(1), csvHandler.CellPosMatrix,50,0,0);
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
try
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
    filepath = [imageHandler.Foldername '/ConvertedCellomics/' selectedWellLong 'Binary' optionHandler.FileEnding];
    if(exist(filepath,'file'))
        BW = imread(filepath);
    else
        BW = zeros(sizeY,sizeX);
    end

    %Create overlay with original pic
    imagePathNeuriteBig = [imageHandler.Foldername '/ConvertedCellomics/' selectedWellLong 'NeuriteBig' optionHandler.FileEnding];

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
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function wellType = GetWellType(handles)
handles = guidata(handles.figure1);
wellType = WellType.Undefined;

%Reihenfolge:
%1: 48 - 8
%2: 48 - 16
%3: 96 - 484
%4: 96 - 8
%5: 96 - 16
%6: OT - 8
%7: OT - 16
%8: Whole - 8
%9: Whole - 16

[wellType,ok] = listdlg('PromptString','Select a Well Type:','SelectionMode','single','ListString',{'48 Well Plate 8 Bit         ';'48 Well Plate 16 Bit        ';'96 Well Plate 8 Bit 484 Pics';'96 Well Plate 8 Bit         ';'96 Well Plate 16 Bit        ';'OT Single Chamber 8 Bit     ';'OT Single Chamber 16 Bit    ';'Whole OT 8 Bit              ';'Whole OT 16 Bit         ';'96 Well Plate 16255 Bit     '});
if(ok>0)
    if(wellType==1)
     wellType = WellType.Well48;
    elseif(wellType==2)
     wellType = WellType.Well4816;         
    elseif(wellType==4)
     wellType = WellType.Well96;
    elseif(wellType==3)
     wellType = WellType.Well96484;
    elseif(wellType==5)
     wellType = WellType.Well9616;
    elseif(wellType==9)
     wellType = WellType.OTWhole16;
    elseif(wellType==8)
     wellType = WellType.OTWhole;
    elseif(wellType==6)
     wellType = WellType.OTSingle;
    elseif(wellType==7)
     wellType = WellType.OTSingle16;
    elseif(wellType==10)
     wellType = WellType.Well9616255;
    end
end


% --- Executes on button press in cbCellomicsNeurons2.
function cbCellomicsNeurons2_Callback(hObject, eventdata, handles)
% hObject    handle to cbCellomicsNeurons2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles)
% Hint: get(hObject,'Value') returns toggle state of cbCellomicsNeurons2
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function CalcContrastStretchIndices(handles)
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
if(~isfield(handles,'NeuronCoordinates'))
    handles.NeuronCoordinates=CSVCoordinates();
end
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
skeletonThresholdMethod = optionHandler.SkeletonThresholdMethod;

if(strfind(skeletonThresholdMethod,';'))
    c = textscan(skeletonThresholdMethod, '%s', 'delimiter', ';');
    %Mean isodata of all wells given in C have to be calculated.
    %Set values of optionHandler
    optionHandler.SkeletonThresholdMethod = '4';
    a=c{1};
    %optionHandler.SkeletonNeuriteThresholdLow = meanOfIsodataResult
else
    a=cell(1);
    a{1}=skeletonThresholdMethod;
    optionHandler.SkeletonThresholdMethod = '4';
end
%Take mean values of all Pictures in experiment
isodataMean=zeros(0);
isodataMeanNuc=zeros(0);
%minimumMean=zeros(0);
stretchlimMean = zeros(0);
wellList = get(handles.lbWell, 'string');
waitbarHandle = waitbar(0,'Please wait. Medium Threshold is calculated.', 'WindowStyle', 'modal');
for i=optionHandler.StartWell:numel(wellList)    
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
    imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
    imagePathNeuriteBig = [foldername '/' selectedWellLong 'NeuriteBig' optionHandler.FileEnding];
    neuriteImage = imread(imagePathNeuriteBig); 
    nucleusImage = imread(imagePathNucleusBig);
    currStretchlimres = stretchlim(neuriteImage,[0.975 0.9999]);
    brightNeurites = imadjust(neuriteImage,currStretchlimres,[0 1]);
    %currThresholdingLevelStrong = th_minimum(brightNeurites)./255;
    %1. Median Filter
    brightNeurites = medfilt2(brightNeurites);

    %2. Resharp image
    unsharpFilter = fspecial('unsharp');
    brightNeurites = imfilter(brightNeurites,unsharpFilter);

    %ISODATA:
    
    if(optionHandler.SkeletonThresholdMethod == '2' | optionHandler.SkeletonThresholdMethod == 2 | any(ismember(a,selectedWell)))
        brightNeurites = imcomplement(brightNeurites);
        brightNucleus = imcomplement(nucleusImage);
        brightNucleus = imadjust(brightNucleus);
        brightNucleus=im2bw(brightNucleus);        
        %Get darkest 
        indices=find(brightNucleus==0);
        currThresholdingLevelNuc = min(nucleusImage(indices));
        %currThresholdingLevelNuc = isodata(brightNucleus);
        currThresholdingLevel = isodata(brightNeurites);
        isodataMean = [isodataMean currThresholdingLevel];
        isodataMeanNuc = [isodataMeanNuc currThresholdingLevelNuc];        
    end
   % minimumMean = [minimumMean currThresholdingLevelStrong];
    
    currStretchlimres = currStretchlimres';
    stretchlimMean = [stretchlimMean;currStretchlimres];
end
if(optionHandler.SkeletonThresholdMethod =='4');
    optionHandler.SkeletonThresholdMethod = '3';
end
isodataMean = mean(isodataMean);
isodataMeanNuc = mean(isodataMeanNuc);
%minimumMean = mean(minimumMean);
stretchlimMean = mean(stretchlimMean,1);

%Check for each pixel max regarding value in original image


neuronHandler.StretchlimResult = stretchlimMean;
neuronHandler.IsodataResult = isodataMean;
%neuronHandler.MinimumResult = minimumMean;
handles.NeuronCoordinates = neuronHandler;

%Calculate also second Threshold for Neurite verification
%Take Minimum threshold method from ImageJ
close(waitbarHandle);
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function CalcContrastStretchIndices_Callback(hObject, eventdata, handles)
% hObject    handle to CalcContrastStretchIndices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
CalcContrastStretchIndices(handles);
guidata(handles.figure1, handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on button press in cbBinaryPic.
function cbBinaryPic_Callback(hObject, eventdata, handles)
% hObject    handle to cbBinaryPic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbBinaryPic
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

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
Threshold=T(i);
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



function level=IJIsoData(data)
maxValue = numel(data)-1;
count0 = data(1);
data(1) = 0; %set to zero so erased areas aren't included
countMax = data(maxValue);
data(maxValue) = 0;
min = 1;
while ((data(min)==0) && (min<maxValue))
	min = min + 1;
end
maximum = maxValue;
while ((data(maximum)==0) && (maximum>0))
        maximum = maximum -1;
    if (min>=maximum) 
        data(1)= count0;
        data(maxValue)=countMax;
        level = (numel(data))/2;
        break;
    end
    movingIndex = min;
    inc = max(maximum/40, 1);
    while 1
        sum1=0.0;
        sum2=0.0;
        sum3=0.0;
        sum4=0.0;
        for (i=min:movingIndex) 
            sum1 = sum1 + i*data(i);
            sum2 = sum2 + data(i);
        end
        for (i=(movingIndex+1):maximum)
            sum3 = sum3 + i*data(i);
            sum4 = sum4 + data(i);
        end
        result = (sum1/sum2 + sum3/sum4)/2.0;
        movingIndex = movingIndex + 1;
        if(~(((movingIndex+1)<=result && movingIndex<maximum-1)))
            break;
        end
    end% while ((movingIndex+1)<=result && movingIndex<max-1);
    data(1)= count0;
    data(maxValue)=countMax;
    level = round(result);
    break;
end

% --- Executes on button press in cbEdgeFillNeurons.
function cbEdgeFillNeurons_Callback(hObject, eventdata, handles)
% hObject    handle to cbEdgeFillNeurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
% Hint: get(hObject,'Value') returns toggle state of cbEdgeFillNeurons
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --------------------------------------------------------------------
function ExportAlgorithmStat_Callback(hObject, eventdata, handles)
% hObject    handle to ExportAlgorithmStat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
foldername = imageHandler.Foldername;
[sizeY sizeX]=size(imageHandler.NucleusImage);
%Get start and stop Well
wellList = get(handles.lbWell, 'string');
exportText = sprintf('Well;xPos;yPos;Skeleton (1=Skeleton found from Edge Fill, 2=Skeleton marked from EdgeFill Second, 3=Skeleton found from endpoint of Skeleton);Edge Fill;Edge Fill Lower;Manual;Deleted by Skel (1=Deleted because too less Neurite length or Distance between Skeleton endpoints, 2=Deleted because another Neuron on same Skeleton and both Neurons are too near together);Edge Fill Overlap;Neurite length;Max dist between Skeleton Endpoints;Nucleus Brightness;Nucleus Size;Nucleus Variance;\r\n');
for i=optionHandler.StartWell:numel(wellList)    
     selectedWell = wellList{i}; 
     NucleusM = csvHandler.CellPosMatrix(selectedWell);
     NeuronM = neuronHandler.CellPosMatrix(selectedWell);
     if(length(selectedWell) == 2)
      selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    elseif(length(selectedWell) == 4)
        selectedWellLong = [selectedWell(1) '0' selectedWell(2:length(selectedWell))];
    else
        selectedWellLong = selectedWell;
     end
    
    imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
    nucleusImage = imread(imagePathNucleusBig);
    if(optionHandler.CutOutSphereCore == 1)
        nucleusPic = CutOutCircles(nucleusImage,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM, 0);
    else
        nucleusPic=nucleusImage;
    end
    nucleusPic(nucleusPic<optionHandler.NucleusThreshold) = 0;
    nucleusPic(nucleusPic>=optionHandler.NucleusThreshold) = 255;
    nucleusPic = uint8(nucleusPic);
    %nucleusPic = nucleusPic(yStartPos:yEndPos,xStartPos:xEndPos);
    nucleusPicBin = logical(nucleusPic);
    nucleusPicBin = logical(xor(bwareaopen(nucleusPicBin,1),  bwareaopen(nucleusPicBin,15000)));
    nucleusPicBinBigNuclei = xor(bwareaopen(nucleusPicBin,250),  bwareaopen(nucleusPicBin,15000));
    D = (~nucleusPicBinBigNuclei);
    nucleusPicBinBigNuclei=0;
    D = -D;
    D(~nucleusPicBin) = -Inf;
    L = watershed(D);
    nucleusPicBin(find(~logical(L))) = 0;
    L=0;
    nucleusPicBin = xor(bwareaopen(nucleusPicBin,35),  bwareaopen(nucleusPicBin,15000));
 
 
    if(isKey(neuronHandler.NeuronStatMatrix, selectedWell))
         NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
        densityWidth = optionHandler.MigrationDistanceDensityImageXSize;
        densityHeight = optionHandler.MigrationDistanceDensityImageYSize;
        [sizeY sizeX]= size(imageHandler.NeuriteImage);
        ringNumber = optionHandler.DensityDistributionRingNumber;
        %Get density distribution to exclude too dense areas
        [filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, densityWidth, densityHeight, optionHandler, NucleusM, imageHandler.Foldername, sizeY, sizeX); 
    
    hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./(1/ScalingFactor)));
    innerMask = logical(createMask(hInner));

     currentRing =(logical(ones(int32(sizeY/(1/ScalingFactor)),int32(sizeX/(1/ScalingFactor))) - innerMask));

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
           if(numel(list) <= 8)
               [siz brightn]=CalculateNucleusSizeAndBrightness(yPoses(j),xPoses(j),nucleusPicBin,nucleusImage,handles)
               list(9) = brightn;
               list(10) = siz;
           end
           if(numel(list) <= 10)
               list(11)=-1;
           end
           exportText = [exportText selectedWell ';' num2str(xPoses(j)) ';' num2str(yPoses(j)) ';' num2str(list(1)) ';' num2str(list(2)) ';' num2str(list(3)) ';' num2str(list(4)) ';' num2str(list(5)) ';' num2str(list(6)) ';' num2str(list(7)) ';' num2str(list(8)) ';' strrep(num2str(list(9)),'.',',') ';' strrep(num2str(list(10)),'.',',') ';' strrep(num2str(list(11)),'.',',') ';' sprintf('\r\n')];
         elseif(NeuronManualM(yPoses(j),xPoses(j)) == 1)
                  [siz brightn]=CalculateNucleusSizeAndBrightness(yPoses(j),xPoses(j),nucleusPicBin,nucleusImage,handles)  
                  exportText = [exportText selectedWell ';' num2str(xPoses(j)) ';' num2str(yPoses(j)) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(1) ';' num2str(0) ';' num2str(0) ';;;' strrep(num2str(brightn),'.',',') ';' strrep(num2str(siz),'.',',') ';' sprintf('\r\n')];
         end
     end
 end
end

[FileName,PathName] = uiputfile('migdist.csv');
fileID = fopen(strcat(PathName,'/',FileName),'w');
fprintf(fileID,'%s',exportText);
fclose(fileID);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function [siz brightness]=CalculateNucleusSizeAndBrightness(y,x,nucleusBin,nucleusOriginal,handles)
%Mark label on part of nucleusBinaryPic
%For efficiency: Cut it out
[sizeY sizeX]=size(nucleusBin);
nucleusBinary = nucleusBin;
nucleusOrig = nucleusOriginal;
if(x-50 > 1 && y - 50 > 1 && x+50 < sizeX && y+50 < sizeY) 
    nucleusBinary = nucleusBinary(y-50:y+50,x-50:x+50);
    nucleusOrig=nucleusOriginal(y-50:y+50,x-50:x+50);
    %Set x and y to mid
    x=50;
    y=50;
end
nucleusBinary = bwlabel(nucleusBinary);
nucLabel = nucleusBinary(y,x);
nucleusBinary(nucleusBinary~=nucLabel)=0;
nucleusBinary(nucleusBinary==nucLabel)=1;   
siz = bwarea(nucleusBinary);
ind = find(nucleusBinary==1);
brightness = mean(nucleusOrig(ind));




% --- Executes on button press in pbAutoThreshold.
function pbAutoThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to pbAutoThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
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
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --------------------------------------------------------------------
function Convert_8Bit_Auto_Callback(hObject, eventdata, handles)
% hObject    handle to Convert_8Bit_Auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
foldername = [imageHandler.Foldername '/ConvertedCellomics'];

selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');

stretchLimLowNeurite = stretchlim(imageHandler.NeuriteImage,[0.975 0.9999]);
stretchLimLowNucleus = stretchlim(imageHandler.NucleusImage,[0.975 0.9999]);

allFiles = dir(foldername); 
waitbarHandle = waitbar(0,'Please wait. All images are converted to 8 Bit.', 'WindowStyle', 'modal');
for i=3:numel(allFiles)
    waitbar((i/numel(allFiles)),waitbarHandle);
    currentFile=allFiles(i).name;

    imageString = strcat(foldername, '/', currentFile);
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
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --- Executes on button press in cbDeleted.
function cbDeleted_Callback(hObject, eventdata, handles)
% hObject    handle to cbDeleted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
% Hint: get(hObject,'Value') returns toggle state of cbDeleted


% --- Executes on button press in cbEdgeFill2.
function cbEdgeFill2_Callback(hObject, eventdata, handles)
% hObject    handle to cbEdgeFill2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbEdgeFill2
try
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


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
try
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

    NucleusM = csvHandler.CellPosMatrix(selectedWell);
    NeuronM = neuronHandler.CellPosMatrix(selectedWell);
    foldername = imageHandler.Foldername;
    [sizeY sizeX] = size(imageHandler.NeuriteImage);
%Threshold Neurite Picture the known way
if(optionHandler.CutOutSphereCore == 1)
    neuritePic = CutOutCircles(imageHandler.NeuriteImage,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM, NeuronM,0);
else
    neuritePic = imageHandler.NeuriteImage;
end

neuritePic = ThresholdPic(neuritePic,0,optionHandler, neuronHandler.StretchlimResult, sizeY, sizeX,0,'Neurite');
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
imageHandler.ResizedBinaryImage = imresize(neuritePic,optionHandler.ScalingFactor);
handles.ImageHandler=imageHandler;
guidata(handles.figure1, handles)
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% ----------------------------------------------------------------- ---
function ManualEval_Callback(hObject, eventdata, handles)
% hObject    handle to ManualEval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
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
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on button press in cbManualEval1.
function cbManualEval1_Callback(hObject, eventdata, handles)
% hObject    handle to cbManualEval1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
RefreshUnzoomedImage(handles);
% Hint: get(hObject,'Value') returns toggle state of cbManualEval1
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --------------------------------------------------------------------
function ParameterExperimentBatch_Callback(hObject, eventdata, handles)
% hObject    handle to ParameterExperimentBatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
%foldername = uigetdir;
%Call EdgeFill Neurons and Skeleton Neurons and play with parameters in
%Option Handler
ExecuteOptionExperiment(handles,0,'Neurite');
guidata(handles.figure1, handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


function ExecuteOptionExperiment(handles,onlyCompFill,channelName)
h = waitbar(0.5,'Neuronal quantification in progress. This could take a few hours. Please be patient.', 'WindowStyle', 'modal');
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
optionHandler = handles.OptionHandler;

if(strcmp(channelName,'Oligo'))
    competingHandler = handles.NeuronCoordinates;
    if(isfield(handles,'OligoCoordinates'))
       neuronHandler = handles.OligoCoordinates;
    else
       neuronHandler = CSVCoordinates();
    end
else
    neuronHandler = handles.NeuronCoordinates;
    competingHandler = 0;
end
wellList = get(handles.lbWell, 'string');
foldername = imageHandler.Foldername;
[sizeY sizeX] = size(imageHandler.NucleusImage);
%Select STARTWELL

if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
if(numel(neuronHandler.AreaDictionary) ==0)
    neuronHandler.AreaDictionary = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
    neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFill) ==0)
    neuronHandler.NeuronPositionsEdgeFill = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFillSecond) ==0)
    neuronHandler.NeuronPositionsEdgeFillSecond = containers.Map();
end
if(numel(neuronHandler.NeuronStatMatrix) ==0)
    neuronHandler.NeuronStatMatrix = containers.Map();
end

%Create Matlab job
%Leave one Core out to ensure that GUI isn't freezing
maxCores =optionHandler.NumberOfCores;
if(maxCores > 1)
    c = parcluster();
    c.NumWorkers = maxCores;
    job = createJob(c,'JobData',optionHandler);
end

for(j=optionHandler.StartWell:numel(wellList))
    selectedWell=wellList(j);
    selectedWell=selectedWell{1};
    selectedWellLong=selectedWell;
    if(length(selectedWell) == 2)
              selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    end
    NucleusM = csvHandler.CellPosMatrix(selectedWell);    
    if(size(NucleusM,2) > sizeX || size(NucleusM,1) > sizeY)
        NucleusM = NucleusM(1:sizeY,1:sizeX);
        csvHandler.CellPosMatrix(selectedWell) = NucleusM;
        handles.CSVCoordinates = csvHandler;
        guidata(handles.figure1, handles);
    end
    if(maxCores > 1)
        %Check if watershed images are available. If not create them before
        %processing as ImageJ operation can't be done in parallel.
        %Check if there is already a watershed image. If not, create it.
        imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
        if(~exist(imagePathWatershed,'file'))
            WatershedNuclei(handles);
        end
        createTask(job, @ExecuteCompositeFillAndSkeleton, 4, {selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler});
    else
       [neuronHandler, TP, FP, TN] = ExecuteCompositeFillAndSkeleton(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler);
    end
end
if(maxCores > 1)
    submit(job);
    wait(job);
    %job = out(2);
    out=fetchOutputs(job);
    %ToDo: Merge NeuronHandler
    TP = sum(cell2mat(out(:,2)));
    FP = sum(cell2mat(out(:,3)));
    TN = sum(cell2mat(out(:,4)));
end
% TPPerWell = cell2mat(out(:,2));
% FPPerWell = cell2mat(out(:,3));
% TNPerWell = cell2mat(out(:,4));
% TPPerWell = TPPerWell ./ (TPPerWell+FPPerWell);
% FPPerWell = FPPerWell ./ (TPPerWell+FPPerWell);
i=0;

if(maxCores > 1)
    for(j=optionHandler.StartWell:numel(wellList))
     %Get current NeuronHandler result and write it in main Neuronhandler
        i=i+1;
        selectedWell=wellList(j);
        selectedWell=selectedWell{1};    
    
        currentNeuronHandler = out(i,1);
        currentNeuronHandler=currentNeuronHandler{1};
        %Merge Neuron positions as well as neurite length matrix.
        neuronHandler.NeuronPositionsSkeletonization(selectedWell) = currentNeuronHandler.NeuronPositionsSkeletonization(selectedWell);
        neuronHandler.AreaDictionary(selectedWell) = currentNeuronHandler.AreaDictionary(selectedWell);
        neuronHandler.NeuronPositionsEdgeFill(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFill(selectedWell);
        neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
        neuronHandler.NeuronStatMatrix(selectedWell) = currentNeuronHandler.NeuronStatMatrix(selectedWell);
        neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = currentNeuronHandler.NeuronPositionsSkelDeleted(selectedWell);
        %neuronHandler.NeuriteLengthMatrix(selectedWell) = currentNeuronHandler.NeuriteLengthMatrix;
    end
end
close(h);
if(strcmp(channelName,'Oligo'))
     handles.OligoCoordinates = neuronHandler;
else
    handles.NeuronCoordinates = neuronHandler;
end
guidata(handles.figure1, handles);  



% --------------------------------------------------------------------
function CalcCellQual_Callback(hObject, eventdata, handles)
% hObject    handle to CalcCellQual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
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
imageHandler = handles.ImageHandler;us
optionHandler = handles.OptionHandler;
WellNeuronDict = neuronHandler.NeuronPositionsEdgeComposite;
%selectedWellNumber = get(handles.lbWell,'Value');
%wellList = get(handles.lbWell, 'string');
%selectedWell = wellList{selectedWellNumber};

densityWidth = 50;
densityHeight = 50;
[sizeY sizeX]= size(imageHandler.NeuriteImage);
%Get density distribution to exclude too dense areas
 NucleusM = csvHandler.CellPosMatrix(selectedWell);
 foldername = imageHandler.Foldername;
 [sizeY sizeX] = size(imageHandler.NeuriteImage);
[filterDistance nonFilterDistance sphereRadius SphereArea NucleusArea64 markerPointCoordinates] = calculateMigrationDistance(selectedWell, densityWidth, densityHeight, optionHandler, NucleusM,  foldername, sizeY, sizeX);    
%end
%MarkerPointCoordinates are available. Get second ring:
hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./(1/ScalingFactor)));
innerMask = logical(createMask(hInner));
currentRing =(logical(ones(int32(sizeY/(1/ScalingFactor)),int32(sizeX/(1/ScalingFactor))) - innerMask));
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
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on button press in cbSkeletonPic.
function cbSkeletonPic_Callback(hObject, eventdata, handles)
% hObject    handle to cbSkeletonPic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
% Hint: get(hObject,'Value') returns toggle state of cbSkeletonPic
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on button press in cbManualEval2.
function cbManualEval2_Callback(hObject, eventdata, handles)
% hObject    handle to cbManualEval2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
% Hint: get(hObject,'Value') returns toggle state of cbManualEval2
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on button press in cbManualEval3.
function cbManualEval3_Callback(hObject, eventdata, handles)
% hObject    handle to cbManualEval3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
% Hint: get(hObject,'Value') returns toggle state of cbManualEval3
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on button press in cbManualEval4.
function cbManualEval4_Callback(hObject, eventdata, handles)
% hObject    handle to cbManualEval4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
% Hint: get(hObject,'Value') returns toggle state of cbManualEval4
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on selection change in popupPlus1.
function popupPlus1_Callback(hObject, eventdata, handles)
% hObject    handle to popupPlus1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
% Hints: contents = cellstr(get(hObject,'String')) returns popupPlus1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupPlus1
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --- Executes during object creation, after setting all properties.
function popupPlus1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupPlus1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
try
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on selection change in popupMinus1.
function popupMinus1_Callback(hObject, eventdata, handles)
% hObject    handle to popupMinus1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
% Hints: contents = cellstr(get(hObject,'String')) returns popupMinus1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupMinus1


% --- Executes during object creation, after setting all properties.
function popupMinus1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupMinus1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
try
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on selection change in popupPlus2.
function popupPlus2_Callback(hObject, eventdata, handles)
% hObject    handle to popupPlus2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
% Hints: contents = cellstr(get(hObject,'String')) returns popupPlus2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupPlus2
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes during object creation, after setting all properties.
function popupPlus2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupPlus2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
try
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on selection change in popupMinus2.
function popupMinus2_Callback(hObject, eventdata, handles)
% hObject    handle to popupMinus2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
% Hints: contents = cellstr(get(hObject,'String')) returns popupMinus2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupMinus2
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes during object creation, after setting all properties.
function popupMinus2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupMinus2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
try
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on button press in cbRemoveCore.
function cbRemoveCore_Callback(hObject, eventdata, handles)
% hObject    handle to cbRemoveCore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
% Hint: get(hObject,'Value') returns toggle state of cbRemoveCore
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on selection change in popupMinus3.
function popupMinus3_Callback(hObject, eventdata, handles)
% hObject    handle to popupMinus3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
% Hints: contents = cellstr(get(hObject,'String')) returns popupMinus3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupMinus3
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

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
try
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
if(isfield(handles,'OligoCoordinates') && isprop(handles.OligoCoordinates,'NeuronPositionsSkeletonization'))
    oligoHandler = handles.OligoCoordinates;
end
optionHandler = handles.OptionHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
%savefoldername = uigetdir;
%prompt = {'First Well', 'Last Well'};
%dlg_title = 'Quality Check Export';
%answer = inputdlg(prompt,dlg_title,1);
wellList = get(handles.lbWell, 'string');
exportText = sprintf('Well;Ring Number;Nuclei;Nuclei without cut;Neurons Manual with cut;Neurons Cellomics with cut;Neurons Manual without cut; Neurons Cellomics without cut;Neurons Cellomics 2;Neurons Edge Fill;Neurons Skeleton;;Neurons Cellomics & Manual with cut;Additional Manual;Additional Cellomics;;Neurons Cellomics & Manual without cut;Additional Manual;Additional Cellomics;;Neurons Cellomics 2 & Manual;Additional Manual;Additional Cellomics 2;;Neurons Edge Fill & Manual;Additional Manual;Additional Edge Fill;;Neurons Skeleton & Manual;Additional Manual;Additional Skeleton;\r\n');
startindex = optionHandler.StartWell;
stopindex=numel(wellList);
waitbarHandle = waitbar(0,'Exporting Quality Check', 'WindowStyle', 'modal');
for i=startindex:stopindex
    waitbar(i/stopindex);
 selectedWell = wellList{i}; 
densityWidth = optionHandler.MigrationDistanceDensityImageXSize;
densityHeight = optionHandler.MigrationDistanceDensityImageYSize;
[sizeY sizeX]= size(imageHandler.NeuriteImage);
ringNumber = optionHandler.DensityDistributionRingNumber;
%Get density distribution to exclude too dense areas   
NucleusM = logical(csvHandler.CellPosMatrix(selectedWell));
[filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, densityWidth, densityHeight, optionHandler, NucleusM, imageHandler.Foldername, sizeY, sizeX); 
  if(markerPointCoordinates ~=0 && numel(markerPointCoordinates('10') > 0))
    hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./(1/ScalingFactor)));
    innerMask = logical(createMask(hInner));
    outerMask = ones(uint16(sizeY/(1/ScalingFactor)),uint16(sizeX/(1/ScalingFactor)));
    mk = markerPointCoordinates(num2str(ringNumber * 10))./(1/ScalingFactor);
    outerMask = roipoly(innerMask,mk(:,1),mk(:,2));
    SE = strel('disk',100);
    outerMask=imdilate(outerMask,SE);
    %currentRing =logical(ones(int32(sizeY/10),int32(sizeX/10)) - innerMask);
    currentRing =logical(outerMask - innerMask);    
    currentRing = imresize(currentRing, [sizeY, sizeX]);
    delete(hInner);
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
if(~isa(neuronHandler.ManualNeuronPositionsSparse,'containers.Map'))
    neuronHandler.ManualNeuronPositionsSparse = containers.Map();
end
if(isKey(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
    ManualNeurons = logical(neuronHandler.ManualNeuronPositionsSparse(selectedWell) .* currentRing);
else
    ManualNeurons = sparse(logical(zeros(sizeY,sizeX)));
end


filterDistance=0;
nonFilterDistance=0;
SphereArea=0;
NucleusArea64=0;
markerPointCoordinates=0;
ManualNeuronsInvert = logical(full(logical(~ManualNeurons)));

if(~isa(neuronHandler.NeuronPositionsEdgeFill,'containers.Map'))
    neuronHandler.NeuronPositionsEdgeFill = containers.Map();
end

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
if(isKey(neuronHandler.CellPosMatrix,selectedWell))
    CellomicsNeuronsWithoutCut = neuronHandler.CellPosMatrix(selectedWell);
    CellomicsNeuronsWithCut = logical(CellomicsNeuronsWithoutCut.* currentRing);
else
    CellomicsNeuronsWithoutCut = 0;
    CellomicsNeuronsWithCut = 0;
end
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
CellomicsNeuronsInvertWithoutCut = logical(full(~CellomicsNeuronsWithoutCut));
ManualNeuronsInvertWithoutCut = logical(full(~ManualNeuronsWithoutCut));
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
fileID = fopen(strcat(PathName,'/',FileName),'w');
fprintf(fileID,'%s',exportText);
fclose(fileID);
close(waitbarHandle);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --------------------------------------------------------------------
function CopyMat_Callback(hObject, eventdata, handles)
% hObject    handle to CopyMat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
answer = questdlg('Copy Matrices on all Wells or only on current Well?','On all Wells?', 'all Wells', 'current Well', 'current Well');
str={'Empty','Nucleus Matrix','Cellomics Neurons','Cellomics Neurons 2','Manual Neurons Main','Manual Neurons 1','Manual Neurons 2','Manual Neurons 3','Manual Neurons 4','Composite Fill Neurons','Neuronal Tracing', 'Composite Fill Oligos','Oligo Tracing'}
[src,v] = listdlg('PromptString','Select matrix you want to copy','SelectionMode','single','ListString',str);
[trg,v] = listdlg('PromptString','Select target matrix','SelectionMode','single','ListString',str);
listboxstrings = get(handles.popupPlus1);
listboxstring = listboxstrings.String(src);
trg = listboxstrings.String(trg);
if(strcmp(answer,'current Well'))
    selectedWell = wellList{selectedWellNumber};
    srcmat=GetMatrixFromListboxString(listboxstring,1,selectedWell,handles);
    SetMatrixFromListboxString(trg,srcmat,selectedWell,handles);
else
    for(j=optionHandler.StartWell:numel(wellList))
        selectedWell=wellList(j);
        selectedWell=selectedWell{1};
        srcmat=GetMatrixFromListboxString(listboxstring,1,selectedWell,handles);
        SetMatrixFromListboxString(trg,srcmat,selectedWell,handles);
    end
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --------------------------------------------------------------------
function CopyMatAnotherSave_Callback(hObject, eventdata, handles)
% hObject    handle to CopyMatAnotherSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
imageHandler = handles.ImageHandler;
if(isfield(handles,'OligoCoordinates'))
    oligoHandler = handles.OligoCoordinates;
end
optionHandler = handles.OptionHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
str={'Empty','Nucleus Matrix','Cellomics Neurons','Cellomics Neurons 2','Manual Neurons Main','Manual Neurons 1','Manual Neurons 2','Manual Neurons 3','Manual Neurons 4','Composite Fill Neurons','Neuronal Tracing', 'Composite Fill Oligos','Oligo Tracing'}
[src,v] = listdlg('PromptString','Select matrix you want to copy','SelectionMode','single','ListString',str);
[trg,v] = listdlg('PromptString','Select target matrix','SelectionMode','single','ListString',str);
listboxstrings = get(handles.popupPlus1);
listboxstring = listboxstrings.String(src);
trg = listboxstrings.String(trg);
trg = trg{1};
listboxstring=listboxstring{1};
%Copy regarding Matrix in Temp object
if(strcmp(listboxstring,'Composite Fill Neurons'))
    tempMatrices = neuronHandler.NeuronPositionsEdgeFill;
elseif(strcmp(listboxstring,'Neuronal Tracing'))
    tempMatrices = neuronHandler.NeuronPositionsSkeletonization;
elseif(strcmp(listboxstring,'Manual Neurons Main'))
    tempMatrices = neuronHandler.ManualNeuronPositionsSparse;
elseif(strcmp(listboxstring,'Nucleus Matrix'))
    tempMatrices = csvHandler.CellPosMatrix;
elseif(strcmp(listboxstring,'Manual Neurons 1'))
     tempMatrices = csvHandler.ManualPositions1;
elseif(strcmp(listboxstring,'Manual Neurons 2'))
    tempMatrices = csvHandler.ManualPositions2;
elseif(strcmp(listboxstring,'Manual Neurons 3'))
    tempMatrices = csvHandler.ManualPositions3;
elseif(strcmp(listboxstring,'Manual Neurons 4'))
    tempMatrices = csvHandler.ManualPositions4;
elseif(strcmp(listboxstring,'Manual Neurons 4'))
    tempMatrices = csvHandler.ManualPositions4;
elseif(strcmp(listboxstring,'Composite Fill Oligos'))
    tempMatrices = oligoHandler.NeuronPositionsEdgeFill;
elseif(strcmp(listboxstring,'Oligo Tracing'))
    tempMatrices = oligoHandler.NeuronPositionsSkeletonization;
end


[filename, pathname] = uigetfile('*.mat', 'Select MAT File where the Matrices should be copied.');
load([pathname filename]);

%csvHandler = handles.CSVCoordinates;
%neuronHandler = handles.NeuronCoordinates;
%imageHandler = handles.ImageHandler;
%optionHandler = handles.OptionHandler;
if(strcmp(trg, 'Manual Neurons 1'))
    csvHandler.ManualPositions1 = tempMatrices;
elseif(strcmp(trg,'Manual Neurons 2'))
    csvHandler.ManualPositions2 = tempMatrices;
elseif(strcmp(trg,'Manual Neurons 3'))
    csvHandler.ManualPositions3 = tempMatrices;
elseif(strcmp(trg,'Manual Neurons 4'))
    csvHandler.ManualPositions4 = tempMatrices;
elseif(strcmp(trg,'Neuronal Tracing'))    
    neuronHandler.NeuronPositionsSkeletonization = tempMatrices;
elseif(strcmp(trg,'Manual Neurons Main'))
    neuronHandler.ManualNeuronPositionsSparse = tempMatrices;
elseif(strcmp(trg,'Oligo Tracing'))
    oligoHandler.ManualNeuronPositionsSparse = tempMatrices;
end
handles.CSVCoordinates = csvHandler;
handles.NeuronCoordinates = neuronHandler;
handles.ImageHandler = imageHandler;
handles.OptionHandler = optionHandler;
if(isfield(handles,'OligoCoordinates') && isprop(handles.OligoCoordinates,'NeuronPositionsSkeletonization'))
    handles.OligoCoordinates = oligoHandler;
end

guidata(handles.figure1, handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end




% --------------------------------------------------------------------
function TestCurrentView_Callback(hObject, eventdata, handles)
% hObject    handle to TestCurrentView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function SetDensDistParam(controlWells,handles)
%Get histogram out of control wells
%Enter name of control wells at Batch
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
C=strread(controlWells,'%s','delimiter',';');
medHop = zeros(numel(C));
medHopArea = zeros(numel(C));
for(i=1:numel(C))
    %Check for each control Well histogram and calculate regarding
    selectedWell = C(i);
    selectedWell = selectedWell{1};
    selectedWellLong=selectedWell;
    if(length(selectedWell) == 2)
      selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    elseif(length(selectedWell) == 4)
        selectedWellLong = [selectedWell(1) '0' selectedWell(2:length(selectedWell))];
    end
    %Load nucleus image
    foldername = [imageHandler.Foldername '/ConvertedCellomics'];
    imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
    currentNucleus = imread(imagePathNucleusBig);
    [counts x] = imhist(currentNucleus);   
    firstDer = diff(counts);
    [maxval hop] = max(firstDer(10:60));  
    hop = hop + 13;
    medHop(i) = hop;
    
    [maxval hop] = max(firstDer(100:240));  
    hop = hop + 100 - 15;
    medHopArea(i) = hop;
end
meanHop = mean(medHop);
meanHopArea=mean(medHopArea);
optionHandler.MigDistLowerNucleusThreshold = meanHop(1);
%optionHandler.MigDistLowerFloodFillThreshold = meanHopArea(1);
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles);

function MultiParamAnalysis(handles)
%1. Get example picture List
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
%if(size(imageHandler.ExamplePicturesYIndices,1) <= 1)
    %SetRandomImageExcerptsForMultiParamAnalysis(handles);
    %handles = guidata(handles.figure1);
%end

h = waitbar(0,'Multiparametric Analysis is running. This could take a weekend. Please be patient.');
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles);
%Start multi Param analysis:
%1. Set all parameters to low value
%Relevant parameters: CompositeFill Overlap, CompositeFill Overlap Min, Min
%Neurite length, Skeleton: Distance Nuclei, Skeleton: Angle Nuclei,
%ManualLowThreshold
%Start 1D Optimization with 1. Param (EdgeFillNucleusAreaWithinNeurite
 multiParamStepContainer = MultiParamStepContainer();
 t=0;
%Threshold
exportText = sprintf('Changed Parameter;Parameter Value;DP;FP;Quality Index\r\n');
[optimum, t,exportText, multiParamStepContainer] = NewtonParamOptimization(exportText,0.7, 0.9, 0, 0, 0, 0, 0.001 ,t+1, multiParamStepContainer, 9, handles);
optionHandler = optimum;
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles);


lowNucThresh = optionHandler.NucleusThreshold - 5;
highNucThresh = optionHandler.NucleusThreshold + 5;

if(lowNucThresh <= 0)
    lowNucThresh = optionHandler.NucleusThreshold;
end


waitbar(1/7);

%Threshold Nuclei
[optimum, t,exportText, multiParamStepContainer] = NewtonParamOptimization(exportText,lowNucThresh, highNucThresh, 0, 0, 0, 0, 1 ,t+1, multiParamStepContainer, 13, handles);
optionHandler = optimum;
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles);
waitbar(2/7);

%Threshold Distance
%EF1
[optimum, t,exportText, multiParamStepContainer] = NewtonParamOptimization(exportText,0.4, 0.6, 0, 0, 0, 0, 0.02 ,t+1, multiParamStepContainer, 1, handles);
optionHandler = optimum;
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles);
waitbar(3/7);
[optimum, t,exportText, multiParamStepContainer] = NewtonParamOptimization(exportText, 0.2, optionHandler.EdgeFillNucleusAreaWithinNeurite, 0, 0, 0, 0, 0.02 ,t+1, multiParamStepContainer, 2, handles);
optionHandler = optimum;
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles);
waitbar(4/7);

%Composite Fill Look Around
[optimum, t,exportText, multiParamStepContainer] = NewtonParamOptimization(exportText, 0, 4, 0, 0, 0, 0, 1 ,t+1, multiParamStepContainer, 14, handles);
optionHandler = optimum;
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles);
waitbar(5/7);
%FuzzFilter Second
[optimum, t,exportText, multiParamStepContainer] = NewtonParamOptimization(exportText, 0, 400, 0, 0, 0, 0, 25 ,t+1, multiParamStepContainer, 12, handles);
optionHandler = optimum;
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles);
waitbar(6/7);
%FuzzFilter
[optimum, t, exportText, multiParamStepContainer] = NewtonParamOptimization(exportText, 1, 13, 0, 0, 0, 0, 0.1 ,t+1, multiParamStepContainer, 10, handles);
optionHandler = optimum;
handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles);
waitbar(7/7);
% [optimum, t,exportText, multiParamStepContainer] = NewtonParamOptimization(exportText,0.7, 0.9, 0, 0, 0, 0, 0.001 ,t+1, multiParamStepContainer, 9, handles);
%  optionHandler = optimum;
%  handles.OptionHandler = optionHandler;
%  guidata(handles.figure1, handles);


% Skeleton Min Neurite length
% [optimum, t, multiParamStepContainer] = NewtonParamOptimization(35, 50, 0, 0, 0, 0, 3 ,t+1, multiParamStepContainer, 3, handles);
% optionHandler = optimum;
% handles.OptionHandler = optionHandler;
% guidata(handles.figure1, handles);
% 
% Distance from Endpoint
% [optimum, t, multiParamStepContainer] = NewtonParamOptimization(15, 40, 0, 0, 0, 0, 5 ,t+1, multiParamStepContainer, 4, handles);
% imageHandler.MultiParamStepContainer = multiParamStepContainer;
% handles.ImageHandler = imageHandler;
% optionHandler = optimum;
% %Tolerance Angle
% [optimum, t, multiParamStepContainer] = NewtonParamOptimization(35, 60, 0, 0, 0, 0, 5 ,t+1, multiParamStepContainer, 5, handles);
% imageHandler.MultiParamStepContainer = multiParamStepContainer;
% handles.ImageHandler = imageHandler;
% optionHandler = optimum;
% 
% New parameters need to be analyzed:
% MinNumberStrongThresholdPixels = 18;
% [optimum, t, multiParamStepContainer] = NewtonParamOptimization(1, 20, 0, 0, 0, 0, 3 ,t+1, multiParamStepContainer, 6, handles);
% optionHandler = optimum;
% handles.OptionHandler = optionHandler;
% guidata(handles.figure1, handles);
% 
% 
% 
% MinDistanceBetweenNeurons = 23;
% [optimum, t, multiParamStepContainer] = NewtonParamOptimization(1, 20, 0, 0, 0, 0, 3 ,t+1, multiParamStepContainer, 7, handles);
% optionHandler = optimum;
% handles.OptionHandler = optionHandler;
% guidata(handles.figure1, handles);
% 
% MinSizeNeuriteArea = 2000;
% [optimum, t, multiParamStepContainer] = NewtonParamOptimization(500, 6000, 0, 0, 0, 0, 500 ,t+1, multiParamStepContainer, 8, handles);
% optionHandler = optimum;
% handles.OptionHandler = optionHandler;
% guidata(handles.figure1, handles);





%Distance and Angle again, as they should be 2d analyzed!
% [optimum, t, multiParamStepContainer] = NewtonParamOptimization(optionHandler.MaxDistanceFromEndpoint-15, optionHandler.MaxDistanceFromEndpoint+15, 0, 0, 0, 0, 3 ,t+1, multiParamStepContainer, 4, handles);
% imageHandler.MultiParamStepContainer = multiParamStepContainer;
% handles.ImageHandler = imageHandler;
% optionHandler = optimum;
% %Tolerance Angle
% [optimum, t, multiParamStepContainer] = NewtonParamOptimization(optionHandler.ToleranceAngleFromEndpoint-15, optionHandler.ToleranceAngleFromEndpoint+15, 0, 0, 0, 0, 3 ,t+1, multiParamStepContainer, 5, handles);
 imageHandler.MultiParamStepContainer = multiParamStepContainer;
 
[FileName,PathName] = uiputfile('MulPaResult.csv');
fileID = fopen(strcat(PathName,'/',FileName),'w');
fprintf(fileID,'%s',exportText);
fclose(fileID);

 handles.ImageHandler = imageHandler;
 optionHandler = optimum;
 
 handles.OptionHandler = optionHandler;
 guidata(handles.figure1, handles);
close(h);
%The previous parameters were 1D param and optimized separately. The
%remaining two are in relation to each other and must be optimized 2D!
%MaximizeDPRatio(handles,multiParamStepContainer);

%MinimizeFPRatio(handles,multiParamStepContainer);
%MaximizeDPRatio(handles,multiParamStepContainer);
%3. Evaluate and maximize DP ratio until DP is maximum and FP is lower than
%max FP

%optional: Export/Save multiParamStepContainer

%Recursive Optimization by Divide and Conquer between start and stop
function [optimum, t, exportText, multiParamStepContainer] = NewtonParamOptimization(exportText,start, stop, startRate, stopRate,startOptionHandler,stopOptionHandler,stopDiff,t, multiParamStepContainer, parameter, handles)
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
currentStep = MultiParamStep();
currentStep.OptionHandler = optionHandler;
if(startRate == 0 && stopRate == 0)
    bothZero=1;
else
    bothZero=0;
end
if(startRate == 0)
    optionHandler.SetValue(parameter,start);
    [TP FP TN]=CalculateCurrentOptionhandlerExperiment(parameter, optionHandler, handles);
    
    %Add it in MultiParamStepContainer
    startOptionHandler = Option(optionHandler);
    currentStep.OptionHandler = startOptionHandler;
    %FPRatioPercent = FP/(TP+TN);
    %TPRatioPercent = TP/(TP+TN);
    TPRatioPercent=TP;
    FPRatioPercent=FP;    
%     if(FPRatioPercent < 0.1)
%         startRate = TPRatioPercent;
%     elseif(FPRatioPercent < 0.15)
%         %Newton linear interpolation
%         FPFactor = (1/0.05)*(FPRatioPercent-0.1);
%         startRate = TPRatioPercent - (FPRatioPercent*FPFactor);
%     else
%         startRate = TPRatioPercent - FPRatioPercent;
%     end
    
    
    if(FPRatioPercent < (optionHandler.MaxAllowedFPFirst/100))
        startRate = TPRatioPercent;
    elseif(FPRatioPercent < (optionHandler.MaxAllowedFPSecond/100))
        %Newton linear interpolation
        FPFactor = (1/((optionHandler.MaxAllowedFPSecond/100) - (optionHandler.MaxAllowedFPFirst/100)))*(FPRatioPercent-(optionHandler.MaxAllowedFPFirst/100));
        startRate = TPRatioPercent - (FPRatioPercent*FPFactor);
    else
        startRate = TPRatioPercent - FPRatioPercent;
    end
    
    currentStep.T = t;
    currentStep.ChangedDim = parameter;
    currentStep.TPWells = TP;
    currentStep.FPWells = FP;  
    currentStep.TNWells = TN;
    multiParamStepContainer.MultiParamStepsDict(t) = currentStep;
    exportText = [exportText num2str(parameter) ';' num2str(start) ';' num2str(TPRatioPercent) ';' num2str(FPRatioPercent)  ';' num2str(startRate) ';' sprintf('\r\n')];
end
if(stopRate == 0)
    optionHandler.SetValue(parameter,stop);
    [TP FP TN]=CalculateCurrentOptionhandlerExperiment(parameter, optionHandler, handles);
    TPRatioPercent=TP;
    FPRatioPercent=FP;
    %FPRatioPercent = FP/(TP+TN);
    %TPRatioPercent = TP/(TP+TN);
%     if(FPRatioPercent < 0.1)
%         stopRate = TPRatioPercent;
%     elseif(FPRatioPercent < 0.15)
%         %Newton linear interpolation
%         FPFactor = (1/0.05)*(FPRatioPercent-0.1);
%         stopRate = TPRatioPercent - (FPRatioPercent*FPFactor);
%     else
%         stopRate = TPRatioPercent - FPRatioPercent;
%     end
    
    
    
    if(FPRatioPercent < (optionHandler.MaxAllowedFPFirst/100))
        stopRate = TPRatioPercent;
    elseif(FPRatioPercent < (optionHandler.MaxAllowedFPSecond/100))
        %Newton linear interpolation
        FPFactor = (1/((optionHandler.MaxAllowedFPSecond/100) - (optionHandler.MaxAllowedFPFirst/100)))*(FPRatioPercent-(optionHandler.MaxAllowedFPFirst/100));
        stopRate = TPRatioPercent - (FPRatioPercent*FPFactor);
    else
        stopRate = TPRatioPercent - FPRatioPercent;
    end
    stopOptionHandler = Option(optionHandler);
    currentStep.OptionHandler = stopOptionHandler;
    currentStep.T = t+bothZero;
    currentStep.ChangedDim = parameter;
    currentStep.TPWells = TP;
    currentStep.FPWells = FP;  
    currentStep.TNWells = TN;
    multiParamStepContainer.MultiParamStepsDict(t) = currentStep;
    exportText = [exportText num2str(parameter) ';' num2str(stop) ';' num2str(TPRatioPercent) ';' num2str(FPRatioPercent)  ';' num2str(stopRate) ';' sprintf('\r\n')];
end
if(stop-start > stopDiff && stopRate ~= startRate)
    %Three Cases:
    %1. Stop is bigger than Start
    %Look for HOP betwenn Mid and Stop
    if(stopRate > startRate)
        [optimum, t, exportText] = NewtonParamOptimization(exportText,start+((stop-start)/4),stop,0,stopRate,startOptionHandler,stopOptionHandler,stopDiff,t+1+bothZero, multiParamStepContainer, parameter, handles);
    else
        %2. Start is bigger than Stop
        %Look for HOP between Start and Mid
        [optimum, t, exportText] = NewtonParamOptimization(exportText,start,stop-((stop-start)/4),startRate,0,startOptionHandler,stopOptionHandler,stopDiff,t+1+bothZero, multiParamStepContainer, parameter, handles);
    end
%Terminate if difference between Start and Stop is <40% of initial
%difference
else
    if(startRate > stopRate)
        optimum = startOptionHandler;
    else
        optimum = stopOptionHandler;
    end
end


function [TP, FP, TN] = CalculateCurrentOptionhandlerExperiment(parameter, optionHandler, handles)
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
neuronHandler = handles.NeuronCoordinates;
csvHandler = handles.CSVCoordinates;
foldername = imageHandler.Foldername;
wellList = get(handles.lbWell, 'string');
%keyList = imageHandler.ExamplePictureMap.keys;
TP = 0;
FP = 0;
TN = 0;
[sizeY sizeX]=size(imageHandler.NucleusImage);
%for(j=numel(wellList):numel(wellList))

%Create Matlab job
parallel.defaultClusterProfile('local');
c = parcluster();
c.NumWorkers = optionHandler.NumberOfCores;
job = createJob(c,'JobData',optionHandler);
excludeList = optionHandler.ExcludedWells;
%One task per Well
for(j=optionHandler.StartWell:numel(wellList))
    selectedWell=wellList(j);
    selectedWell=selectedWell{1};
    if(numel(strfind(excludeList,selectedWell)) == 0)
        NucleusM = csvHandler.CellPosMatrix(selectedWell);
        if(size(NucleusM,2) > sizeX || size(NucleusM,1) > sizeY)
            NucleusM = NucleusM(1:sizeY,1:sizeX);
            csvHandler.CellPosMatrix(selectedWell) = NucleusM;
            handles.CSVCoordinates = csvHandler;
            guidata(handles.figure1, handles);
        end
    %     currentKey = keyList(j);
    %     currentKey = currentKey{1};
    %     wellLetter = imageHandler.ExamplePictureMap(currentKey);
    %     selectedWell = [wellLetter currentKey];
    %     if(mod(j,2) == 0)
    %             A = optionHandler.RectangleSizeA;
    %             B = optionHandler.RectangleSizeB;
    %     else
    %             B = optionHandler.RectangleSizeA;
    %             A = optionHandler.RectangleSizeB;
    %     end
        %if(parameter==1)
        %    [TPcur FPcur TNcur] = EdgeFillNeuronsAlgorithm(selectedWell, imageHandler.ExamplePicturesYIndices(selectedWell), imageHandler.ExamplePicturesYIndices(selectedWell)+A, imageHandler.ExamplePicturesXIndices(selectedWell), imageHandler.ExamplePicturesXIndices(selectedWell)+B, 0, handles);                   
            %Everything, not only image excerpts
            createTask(job, @ExecuteCompositeFillAndSkeleton, 4, {selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,0,'Neurite',0});
    end
end
submit(job);
wait(job);

    out=fetchOutputs(job);
    [TP FP] = CalculateTPPerWell(out);


    % TP = sum(cell2mat(out(:,2)));
    % FP = sum(cell2mat(out(:,3)));
    % TN = sum(cell2mat(out(:,4)));
    %ToDo: Merge NeuronHandler
    f=0;
    if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
        neuronHandler.NeuronPositionsSkeletonization = containers.Map();
    end
    if(numel(neuronHandler.AreaDictionary) ==0)
        neuronHandler.AreaDictionary = containers.Map();
    end
    if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
        neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
    end
    if(numel(neuronHandler.NeuronPositionsEdgeFill) ==0)
        neuronHandler.NeuronPositionsEdgeFill = containers.Map();
    end
    if(numel(neuronHandler.NeuronPositionsEdgeFillSecond) ==0)
        neuronHandler.NeuronPositionsEdgeFillSecond = containers.Map();
    end
    for(j=optionHandler.StartWell:numel(wellList))
    %Get current NeuronHandler result and write it in main Neuronhandler
        selectedWell=wellList(j);
        selectedWell=selectedWell{1}; 
        if(numel(strfind(excludeList,selectedWell)) == 0)
            f=f+1;   
            currentNeuronHandler = out(f,1);
            currentNeuronHandler=currentNeuronHandler{1};
            %Merge Neuron positions as well as neurite length matrix.
            
            neuronHandler.NeuronPositionsSkeletonization(selectedWell) = currentNeuronHandler.NeuronPositionsSkeletonization(selectedWell);
            neuronHandler.AreaDictionary(selectedWell) = currentNeuronHandler.AreaDictionary(selectedWell);
            %neuronHandler.NeuronStatMatrix(selectedWell) = currentNeuronHandler.NeuronStatMatrix(selectedWell);
            neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = currentNeuronHandler.NeuronPositionsSkelDeleted(selectedWell);
            neuronHandler.NeuronPositionsEdgeFill(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFill(selectedWell);
            neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
        end
    end

handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles);


function [TP FP] = CalculateTPPerWell(out)
    TPPerWell = cell2mat(out(:,2));
    FPPerWell = cell2mat(out(:,3));
    TNPerWell = cell2mat(out(:,4));
    %TPPerWell = TPPerWell ./ (TPPerWell+FPPerWell);
    %FPPerWell = FPPerWell ./ (TPPerWell+FPPerWell);
    for(j=1:numel(TPPerWell))
        if(TPPerWell(j)+FPPerWell(j)==0)
            FPPerWell(j)=FPPerWell(j);
            TPPerWell(j)=1;
        else
            FPPerWell(j)=(FPPerWell(j))/(TPPerWell(j)+TNPerWell(j));
            TPPerWell(j) = (TPPerWell(j))/(TPPerWell(j)+TNPerWell(j));
        end
    end

    TP = mean(TPPerWell);
    FP = mean(FPPerWell);

 function [neuronHandler, TP, FP, TN] = ExecuteCompositeFillAndSkeleton(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler)    
        %Get first filter from selected Well
        firstY=1;
        firstX=1;
        lastY=sizeY;
        lastX=sizeX;      
         if(optionHandler.FilterMulPa==1)
             path = [foldername '/ConvertedCellomics/' selectedWell '.mat'];
             TP=0;
             FP=0;
             TN=0;
             if(exist(path,'file')) 
                 load(path);
                 %For every positive Filter
                 filterCount = numel(FMask.PositiveFilters);
                 for(i=1:filterCount);
                    f1 = FMask.PositiveFilters{i};
                    firstInd=find(f1,1, 'first');
                    lastInd=find(f1,1, 'last');
                    [firstY firstX] = ind2sub([sizeY sizeX], firstInd);
                    [lastY lastX] = ind2sub([sizeY sizeX], lastInd);
                    [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithm(selectedWell,firstY,lastY, firstX, lastX,0,optionHandler,neuronHandler,NucleusM,foldername,channelName);                     
                    if(onlyCompFill == 0)
                        [neuronHandler TPSkel FPSkel TNSkel] = SkeletonizationNeuronsAlgorithm(selectedWell,firstY,lastY, firstX, lastX, 1, optionHandler,sizeY,sizeX,foldername,neuronHandler,NucleusM,1,channelName,competingHandler);
                        TP=TP+TPSkel;
                        FP=FP+FPSkel;
                        TN=TN+TNSkel;
                    end
                 end
             end
         else
             if(lastY - firstY > 10000 || lastX - firstX > 10000)        
                sizeXsmall = (lastX-firstX)/2;
                sizeYsmall = (lastY-firstY)/2;
                [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithm(selectedWell,firstY,sizeYsmall,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithm(selectedWell,firstY,sizeYsmall,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithm(selectedWell,sizeYsmall+1,lastY,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithm(selectedWell,sizeYsmall+1,lastY,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
              else
                [neuronHandler TP FP TN] = EdgeFillNeuronsAlgorithm(selectedWell,firstY,lastY, firstX, lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);  
              end
              if(onlyCompFill == 0)
                [neuronHandler TP FP TN] = SkeletonizationNeuronsAlgorithm(selectedWell,firstY,lastY, firstX, lastX, 1, optionHandler,sizeY,sizeX,foldername,neuronHandler,NucleusM,1,channelName,competingHandler);
              end
         end    
         
        % system('taskkill /fi "WINDOWTITLE eq FBird"');
        %Calculate only for filter 
 
        



function MinimizeFPRatio(handles,multiParamStepContainer)
%1. Analyze all sub pictures and count number of true positive, false
%positives and false negatives.
handles = guidata(handles.figure1);
wellList = get(handles.lbWell, 'string');
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
optionHandler.EdgeFillNucleusAreaWithinNeurite = 0.35;
optionHandler.EdgeFillNucleusAreaWithinNeuriteSecond = 0.2;
%optionHandler.SkeletonMinNeuriteLength = 25;
%optionHandler.MaxDistanceFromEndpoint = 20;
%optionHandler.ToleranceAngleFromEndpoint = 20;
t=0;
oldFPperDP = 999999;
currentFPperDP = 999998;
%while(currentFPperDP < oldFPperDP && currentFPperDP > 0.1)  %More than 10% FP)%FP/DP is decreasing
    %Try changing all 6 parameters and check which one reduces FP mostly.
    oldFPperDP = currentFPperDP;
    minFPIndex = 1;
    minFP = 999999;
    minFPperDP = 999999;
    multiParamStepContainerMethod = MultiParamStepContainer();
    handles.OptionHandler = optionHandler;
    guidata(handles.figure1, handles);
    multiParamStepContainer.MultiParamStepsDict(t) = multiParamStepContainerMethod.MultiParamStepsDict(minFPIndex);
    optionHandler = multiParamStepContainerMethod.MultiParamStepsDict(minFPIndex).OptionHandler;
    t=t+1;
    currentFP = multiParamStepContainerMethod.MultiParamStepsDict(minFPIndex).FPWells;
    currentTP = multiParamStepContainerMethod.MultiParamStepsDict(minFPIndex).TPWells;
    currentFPperDP = currentFP / currentTP;    
    
    %Do Newton Optimization for current Parameter with current
    %configuration
    [optimum t] = NewtonParamOptimization(handles, 1, t);
%end

handles.OptionHandler = optionHandler;
guidata(handles.figure1, handles);


function PlotMultiParamAnalysisResults(plotParam,handles)
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
multiParamStepContainer = imageHandler.MultiParamStepContainer;
keyList = multiParamStepContainer.MultiParamStepsDict.keys;
xAxis = zeros(0);
yTP = zeros(0);
yFP = zeros(0);
for(i=2:numel(keyList))
    currentKey = keyList(i);
    currentKey = currentKey{1};
    multiParamStep = multiParamStepContainer.MultiParamStepsDict(currentKey);
    if(multiParamStep.ChangedDim == plotParam)
        paramValue = multiParamStep.OptionHandler.GetValue(plotParam);        
        TPRatio = multiParamStep.TPWells / (multiParamStep.TPWells + multiParamStep.TNWells);
        FPRatio = multiParamStep.FPWells / (multiParamStep.TPWells + multiParamStep.TNWells);
        %Plot: X-Axis: Parameter Values
        xAxis = [xAxis paramValue];
        %Y-Axis: TP Ratio and FP Ratio
        yTP = [yTP TPRatio];
        yFP = [yFP FPRatio];
    end    
end
xAxis = xAxis.';
yTP = yTP.';
yFP = yFP.';
dataTP = [xAxis yTP];
dataFP = [xAxis yFP];
dataTP = sortrows(dataTP,1);
dataFP = sortrows(dataFP,1);
figure(1);
plot(dataTP(:,1),dataTP(:,2),'b');
hold on;
plot(dataFP(:,1),dataFP(:,2),'r');


function SetRandomImageExcerptsForMultiParamAnalysis(handles)
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
csvHandler = handles.CSVCoordinates;
wellList = get(handles.lbWell, 'string');
[sizeY sizeX] = size(imageHandler.NucleusImage);

sizeYScaled = uint16(sizeY/(1/ScalingFactor));
sizeXScaled = uint16(sizeX/(1/ScalingFactor));
imageHandler.ExamplePicturesYIndices = containers.Map();
imageHandler.ExamplePicturesXIndices = containers.Map();
densityWidth = optionHandler.MigrationDistanceDensityImageXSize;
densityHeight = optionHandler.MigrationDistanceDensityImageYSize;
%One excerpt for every wellnumber
% wellLettersForNumberMap = containers.Map();
     withinMarkerpoints=0;
%Iterate over all selected Wells
%for(i=optionHandler.StartWell:numel(wellList))
%for(i=1:numel(wellLettersForNumberMap.keys))
for(i=optionHandler.StartWell:numel(wellList))
    %Get random position within Markerpoints
    %Low right of rectangle has to be also within markerpoints!
    selectedWell=wellList(i);
    selectedWell=selectedWell{1};
        %Get density distribution to exclude too dense areas
        %if(saveCircle)    
        NucleusM = csvHandler.CellPosMatrix(selectedWell);
        [filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, densityWidth, densityHeight, optionHandler, NucleusM, imageHandler.Foldername, sizeY, sizeX);         
        hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./(1/ScalingFactor)));
        innerMask = logical(createMask(hInner));
        delete(hInner);
    
    while(~withinMarkerpoints)
        yPos = randi([1,sizeY]);
        xPos = randi([1,sizeX]);
        xPosOrig=xPos;
        yPosOrig=yPos;
        yPos = uint16(yPos/(1/ScalingFactor));
        xPos = uint16(xPos/(1/ScalingFactor));
        if(yPos < 1)
            yPos=1;
        end
        if(xPos <1)
            xPos=1;
        end         
        
        if(innerMask(yPos,xPos) == 1)
            %Random point is within inner circle
            continue;
        end          
                
        currentRing =logical(ones(int32(sizeY/(1/ScalingFactor)),int32(sizeX/(1/ScalingFactor))) - innerMask);

        hOuter = impoly(handles.axes2,double(markerPointCoordinates(num2str(ringNumber * 10))./(1/ScalingFactor)));
        outerMask = logical(createMask(hOuter));
        delete(hOuter);
        if(~outerMask(yPos,xPos)==1)            
            continue;
        end
        
        if(mod(i,2) == 0)
            A = optionHandler.RectangleSizeA;
            B = optionHandler.RectangleSizeB;
        else
            B = optionHandler.RectangleSizeA;
            A = optionHandler.RectangleSizeB;
        end
        %Same for Bottom Right Point
        hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./(1/ScalingFactor)));
        innerMask = logical(createMask(hInner));
        delete(hInner);
        if(uint16(yPos+(A/(1/ScalingFactor))) > sizeYScaled || uint16(xPos+(B/(1/ScalingFactor))) > sizeXScaled || (size(innerMask,1) > uint16(yPos+(A/(1/ScalingFactor))) && size(innerMask,2) > uint16(xPos+(B/(1/ScalingFactor))) && innerMask(uint16(yPos+(A/(1/ScalingFactor))),uint16(xPos+(B/(1/ScalingFactor)))) == 1))
            %Random point is within inner circle
            continue;
        end          
                
        currentRing =logical(ones(int32(sizeY/(1/ScalingFactor)),int32(sizeX/(1/ScalingFactor))) - innerMask);

        hOuter = impoly(handles.axes2,double(markerPointCoordinates(num2str(ringNumber * 10))./(1/ScalingFactor)));
        outerMask = logical(createMask(hOuter));
        delete(hOuter);
        if(~outerMask(uint16(yPos+(A/(1/ScalingFactor))),uint16(xPos+(B/(1/ScalingFactor))))==1)            
            continue;
        else
            withinMarkerpoints=1;
        end
    end
    imageHandler.ExamplePicturesYIndices(selectedWell) = yPosOrig;
    imageHandler.ExamplePicturesXIndices(selectedWell) = xPosOrig;
end
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function LoadBatch_Callback(hObject, eventdata, handles)
% hObject    handle to LoadBatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 1. Get global folder
try
foldername = uigetdir;
if(foldername ~=0)
% 2. Get name of experiment
expname = inputdlg('Experiment name');
if(numel(expname)>0)
expname = expname{1};
controlWells = inputdlg('Control Wells separated by ;');
if(numel(controlWells)>0)
    controlWells = controlWells{1};
    EightBitConversion = questdlg('Should the raw images automatically converted to 8 Bit?')
    if(~strcmp(EightBitConversion,'Cancel'))
        %FlipNuclei = questdlg('Should the position matrices flipped around the x Axis?')
        wellType = GetWellType(handles);
        if(wellType ~= WellType.Undefined)
            foldernameConverted=strcat(foldername,'/',expname);
            subfoldername = strcat(foldername, '/ConvertedCellomics');
            channel2=0;
            channel3=0;
            channel4=0;
            v2=0;
            v3=0;
            v4=0;
            str={'Neurons','Oligodendrocytes','Astrocytes','Channel Not Available'}
            [channel2,v2] = listdlg('PromptString','Which cell Type is Channel 2?',...
                    'SelectionMode','single',...
                    'ListString',str);

            [channel3,v3] = listdlg('PromptString','Which cell Type is Channel 3?',...
                    'SelectionMode','single',...
                    'ListString',str);    

            [channel4,v4] = listdlg('PromptString','Which cell Type is Channel 4?',...
                    'SelectionMode','single',...
                    'ListString',str); 
            if(channel2 == 4)
                v2=0;
            end
            if(channel3==4)
                v3=0;
            end
            if(channel4==4)
                v4=0;
            end       

            %Check if CSV files for channel exist
            pathNuclei = strcat(foldername,'/', expname, '-Kernkoordinaten.csv');
            pathNeurons = strcat(foldername,'/', expname, '.csv');

            pathNucleiMDB = strcat(foldername,'/', expname, '-Kernkoordinaten.mdb');
            pathNeuronsMDB = strcat(foldername,'/', expname, '.mdb');


            if(channel2~=0 && channel3~=0 && channel4~=0)
                choice = exist(pathNuclei,'file');
                choice2 = exist(pathNucleiMDB,'file');
                if(~strcmp(choice, 'Cancel'))
                    %Call 'Load Image Folder'
                    CreateOptionHandler(handles,controlWells);
                    handles = guidata(handles.figure1);
                    if(strcmp(EightBitConversion,'Yes'))
                        ConvertPicturesTo8Bit(handles,foldernameConverted);
                    end
                    LoadImageFolder(foldernameConverted,wellType,channel2,channel3,channel4,v2,v3,v4, handles); 

                    %Show question if Nuclei should be loaded by CSV or automatically
                    %calculated

                    % Handle response
                    if(choice)    
                           LoadCSV(['/' expname '-Kernkoordinaten.csv'],foldername,wellType,0,handles);
                    elseif(choice2)
                        LoadCSV(['/' expname '-Kernkoordinaten.csv'],foldername,wellType,1,handles);
                    else
                            %Create Nuclei positions by watershed        
                            ReloadWellList(handles,1);
                            watershedNucleusImage(wellType,controlWells,handles)
                            handles = guidata(handles.figure1);
                    end


                    %Call Load Neuron CSV
                    if(exist(pathNeurons,'file'))
                        LoadNeuronCSV(['/' expname '.csv'],foldername,wellType,0,handles);
                    elseif(exist(pathNeuronsMDB,'file'))
                        LoadNeuronCSV(['/' expname '.csv'],foldername,wellType,1,handles);
                    end
                    %Set correct thresholds for analyzing Migration Distance and Density
                    %Distribution
                    SetDensDistParam(controlWells,handles);
                    ReloadWellList(handles,0);
                    CalcContrastStretchIndices(handles);
                    handles = guidata(handles.figure1);
                    if(exist(pathNeurons,'file') || exist(pathNeuronsMDB,'file'))
                        %ToDo: Fix also Oligo and Astropositions
                        FixNeuronPositions(hObject, eventdata, handles);
                    end
                    guidata(handles.figure1, handles);
                    %Create Savegame
                    SaveFixedNeuronPositions(strcat(expname,'.mat'),foldername,handles);
                    %guidata(handles.figure1, handles);
                end
            end
        end
    end
end
end
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


function closeThreshFigure(source,callbackdata,handles)
    close(figure(source.Parent.Number));

function watershedNucleusImage(wellType,controlWells,handles)
    %Load Nucleus image
    %Return CSVHandler
    handles = guidata(handles.figure1);
    optionHandler = handles.OptionHandler;
    
    if(~isfield(handles,'CSVCoordinates'))
        csvHandler = CSVCoordinates();        
        csvHandler.ManualNeuronPositionsSparse = containers.Map();
    else
        csvHandler = handles.CSVCoordinates;
    end
    csvHandler.CellPosMatrix = containers.Map();
    imageHandler = handles.ImageHandler;
    wellList = get(handles.lbWell, 'string');
    foldername = [imageHandler.Foldername '/ConvertedCellomics'];
    wellList = get(handles.lbWell, 'string'); 
    waitbarHandle = waitbar(0,'Searching for Nuclei', 'WindowStyle', 'modal');
        %Set user defined threshold
    f=figure('units','normalized','outerposition',[0 0.05 1 0.95]);
    %Create slider in figure
    ax = axes('Units','pixels');
    %slider = uicontrol('Style','slider','Min',1,'Max',255,'Position',[10 10 500 30],'Value',15,'Callback',{@refreshNucThresh,controlWells,handles});
    
    slider = uicontrol('Style','slider','Units','normalized','Min',1,'Max',255,'Position',[0.125 0.05 0.78 0.02],'Value',15,'Callback',{@refreshNucThresh,controlWells,handles});
    label = uicontrol('Style','text','Units','normalized','Position',[0 0.015 0.07 0.05],'String','Threshold');
    editBox = uicontrol('Style','edit','Units','normalized','String','15','Position',[0.065 0.05 0.05 0.02],'Value',15,'Callback',{@refreshNucThresh,controlWells,handles});
    applyBtn = uicontrol('Units','normalized','Position',[0.125 0.01 0.05 0.04],'String','Apply','Callback',{@closeThreshFigure,handles})
    refreshNucThresh(slider,'empty',controlWells,handles);
    uiwait(f);
    handles = guidata(handles.figure1);
    optionHandler = handles.OptionHandler;
    threshold = optionHandler.NucleusThreshold;
    if(numel(csvHandler.ManualNeuronPositionsSparse) ==0)
        csvHandler.ManualNeuronPositionsSparse = containers.Map();
    end
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
        if(wellType ~= 0)
            if(wellType == WellType.Well96 || wellType == WellType.Well9616  || wellType == WellType.Well9616255)
                csvHandler.CellPosMatrix(selectedWell) = sparse(7168, 7168);
                csvHandler.ManualNeuronPositionsSparse(selectedWell) = sparse(7168, 7168);
            elseif(wellType == WellType.Well96484)
                csvHandler.CellPosMatrix(selectedWell) = sparse(11264,11264);
                csvHandler.ManualNeuronPositionsSparse(selectedWell) = sparse(11264,11264);
            elseif(wellType == WellType.Well48 ||wellType == WellType.Well4816)
                csvHandler.CellPosMatrix(selectedWell) = sparse(18432, 18432);
                csvHandler.ManualNeuronPositionsSparse(selectedWell) = sparse(18432, 18432);
            elseif(wellType == WellType.OTSingle || wellType == WellType.OTSingle16)
                csvHandler.CellPosMatrix(selectedWell) = sparse(13824, 10752);
                csvHandler.ManualNeuronPositionsSparse(selectedWell) = sparse(13824, 10752);
            elseif(wellType == WellType.OTWhole || wellType == WellType.OTWhole16)
                csvHandler.CellPosMatrix(selectedWell) = sparse(22016, 10240);
                csvHandler.ManualNeuronPositionsSparse(selectedWell) = sparse(22016, 10240);
            end %if
        else
            [sizeY sizeX] = size(imageHandler.NucleusImage);
            csvHandler.CellPosMatrix(selectedWell) = sparse(sizeY, sizeX);
        end
    
    
        imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
        imagePathNucleusSmall = [foldername '/' selectedWellLong 'NucleusSmall' optionHandler.FileEnding];
        imagePathNeuriteBig = [foldername '/' selectedWellLong 'NeuriteBig' optionHandler.FileEnding];
        imagePathNeuriteSmall = [foldername '/' selectedWellLong 'NeuriteSmall' optionHandler.FileEnding];
        NeuriteImage = imread(imagePathNeuriteBig);
        [sizeYorig sizeXorig] = size(NeuriteImage);
        [sizeY sizeX] = size(NeuriteImage);
        %if(sav)
            BinaryImage = logical(zeros(sizeYorig, sizeXorig));
            SkeletonImage = logical(zeros(sizeYorig, sizeXorig));
        %end
        %else
        %    BinaryImage(yStartPos:yEndPos,xStartPos:xEndPos) = logical(zeros(sizeY, sizeX));
        %    SkeletonImage(yStartPos:yEndPos,xStartPos:xEndPos) = logical(zeros(sizeY, sizeX));
        %end
        %Get density distribution to exclude too dense areas
        %filepath = [imageHandler.Foldername '/ConvertedCellomics/' selectedWell 'BinaryCut.tif'];
        %if(~exist(filepath,'file'))        
        nucleusPic = imread(imagePathNucleusBig);
        
%         %Another (hopefully more intelligent) way to get initial Nucleus
%         %Threshold:
%         [sizeY sizeX] = size(nucleusPic);
%         GW = sizeY * sizeX;
%         PW = GW * 0.022;
%         PW = GW - PW;
%         threshold = 1;
%         pixelCount = nnz(find(nucleusPic <= threshold));
%         while(pixelCount < PW)
%             threshold = threshold + 1;
%             pixelCount = nnz(find(nucleusPic <= threshold));
%         end
%         if(threshold<5)
%             threshold=5;
%         end
        %nucleusPic = CutOutCircles(NucleusImage,selectedWell,1,0,1,optionHandler, foldername, sizeY, sizeX, NucleusM, NeuronM,0);
        nucleusPic(nucleusPic<threshold) = 0;%double(double(imageHandler.NucleusImage) ./ 255);
        nucleusPic(nucleusPic>=threshold) = 1;
        nucleusPic = logical(nucleusPic);        
        nucleusPic = (xor(bwareaopen(nucleusPic,1),  bwareaopen(nucleusPic,15000)));
        nucleusPicBinBigNuclei = logical(xor(bwareaopen(nucleusPic,250),  bwareaopen(nucleusPic,15000)));        
        nucleusPic=uint8(nucleusPic);
        nucleusPic(nucleusPic==0)=255;
        nucleusPic(nucleusPic==1)=0;
        
        %Invoke ImageJ for watershed as ImageJ's watershed is far better than
        %Matlabs watershed!    
    %Do Watershed
    javaaddpath('lib/ij.jar');
        javaaddpath('lib/mij.jar');
        MIJ.start(false);
        MIJ.createImage(nucleusPic);
        MIJ.run('Options...', 'iterations=1 count=1 edm=Overwrite do=Nothing');
        MIJ.run('Watershed');
        %nucleusPicBin = logical(MIJ.getCurrentImage);
        MIJ.run('Save', 'save=temp.tif');
        MIJ.closeAllWindows;
        MIJ.exit;
        nucleusPic = imread('temp.tif');
        delete temp.tif;
        nucleusPic = logical(1-nucleusPic);
        imagePathWatershed = [foldername '/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
        imwrite(nucleusPic,imagePathWatershed);
        
        
        
%         D = bwdist(~nucleusPicBinBigNuclei);
%         D = -D;
%         D(~nucleusPic) = -Inf;
%         L = logical(watershed(D));
%         D=0;
%         nucleusPic(find(~logical(L))) = 0;
        nucleusPic = xor(bwareaopen(nucleusPic,35),  bwareaopen(nucleusPic,15000));
        
        %Treat every area of nucleus image as nucleus and save it to CSV
        %Handler
        
        balancePointStructs = regionprops(nucleusPic, 'Centroid');
        
%         L = bwlabel(nucleusPic);
%         labelCount = numel(max(L(:)));
%         for(j=1:labelCount);
%             %Get balance point of current nucleus
%         end

        for(j=1:numel(balancePointStructs))
            currentCentroid = balancePointStructs(j);
            y = round(currentCentroid.Centroid(2));
            x = round(currentCentroid.Centroid(1));
            %Save data from currentCentroid to CSVHandler
            currentMatrix = csvHandler.CellPosMatrix(selectedWell);
            currentMatrix(y,x) = 1;
            csvHandler.CellPosMatrix(selectedWell) = currentMatrix;
        end        
    end
    close(waitbarHandle);
    handles.CSVCoordinates = csvHandler;
guidata(handles.figure1, handles);



    
    
function refreshNucThresh(source,callbackdata,controlWells,handles)
    
    %Reload Nucleus Image
    
    optionHandler = handles.OptionHandler;
    controlWells=strread(controlWells,'%s','delimiter',';');
    %controlWells = strsplit(controlWells,';');
    nucPrev=nucleusPreview(handles,controlWells);
    [sizeY sizeX] = size(nucPrev);
    bwImage = uint8(zeros(sizeY,sizeX));
    %Make image black and white    
    if(strcmp(source.Style, 'edit'))
        threshold = str2num(get(source,'String'));
    else
        threshold = get(source,'value');
    end
    source.Parent.Children(4).Value = threshold;
    source.Parent.Children(2).Value = threshold;
    source.Parent.Children(2).String = num2str(threshold);
    bwImage(nucPrev>=threshold)=0;
    bwImage(nucPrev<threshold)=255;
    
    %Do Watershed
    javaaddpath('lib/ij.jar');
        javaaddpath('lib/mij.jar');
        MIJ.start(false);
        MIJ.createImage(bwImage);
        MIJ.run('Options...', 'iterations=1 count=1 edm=Overwrite do=Nothing');
        MIJ.run('Watershed');
        %nucleusPicBin = logical(MIJ.getCurrentImage);
        MIJ.run('Save', 'save=temp.tif');
        MIJ.closeAllWindows;
        MIJ.exit;
        bwImage = imread('temp.tif');
        delete temp.tif;
        bwImage = logical(1-bwImage);
    
    

%     D = bwdist(~nucleusPicBinBigNuclei);
%     D = -D;
%     D(~bwImage) = -Inf;
%     L = watershed(D);
%     bwImage(find(~logical(L))) = 0;
%     bwImage = xor(bwareaopen(bwImage,35),  bwareaopen(bwImage,5000));
    
    
    
    figure(source.Parent.Number);
    imshow(bwImage);
    optionHandler.NucleusThreshold = threshold;
    handles.optionHandler=optionHandler;
    guidata(handles.figure1, handles);
    
function totalImage=nucleusPreview(handles,controlWells)
    %Load 6 Nuclei Images into one well
    optionHandler = handles.OptionHandler;
    imageHandler = handles.ImageHandler;
    i = optionHandler.StartWell;
    wellList = get(handles.lbWell, 'string');
    selectedWell = wellList{i};
    foldername = [imageHandler.Foldername '/ConvertedCellomics'];
    
    if(length(selectedWell) == 2)
      selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    elseif(length(selectedWell) == 4)
        selectedWellLong = [selectedWell(1) '0' selectedWell(2:length(selectedWell))];
    else
        selectedWellLong = selectedWell;
    end
    imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
    nucImage = imread(imagePathNucleusBig);
    [sizeY sizeX] = size(nucImage);
    
    %Take only 1st control well for speed purposes
    controlWellCount = numel(controlWells);
    maxCount = controlWellCount;
%     if(controlWellCount > 9)
%         %Drei Reihen
%         sizeXTotal = sizeX * 3;
%         sizeYTotal = sizeY * 3;
%         maxCount=9;
%     elseif(controlWellCount <= 9 && controlWellCount > 6)
%         %Ebenfalls drei Reihen
%         sizeXTotal = sizeX * 3;
%         sizeYTotal = sizeY * 3;
%     elseif(controlWellCount <= 6 && controlWellCount >3)
%         %Zwei Reihen
%         sizeXTotal = sizeX * 3;
%         sizeYTotal = sizeY * 2;
%     else
%         %Eine Reihe
%         sizeXTotal = sizeX * controlWellCount;
%         sizeYTotal = sizeY;
%     end  
    sizeYTotal = sizeY;
    sizeXTotal = sizeX;
    
    totalImage = uint8(zeros(sizeYTotal,sizeXTotal));
    
    %Iteriere �ber KontrolWells
    startY = 1;
    startX = 1;
    
    %Take currently only 1st control well
    for(i=1:1)
        selectedWell=controlWells(i);
        selectedWell=selectedWell{1};
        if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
        elseif(length(selectedWell) == 4)
            selectedWellLong = [selectedWell(1) '0' selectedWell(2:length(selectedWell))];
        else
            selectedWellLong = selectedWell;
        end
        %selectedWellLong=selectedWellLong{1};
        imagePathNucleusBig = strcat(foldername, '/', selectedWellLong, 'NucleusBig', optionHandler.FileEnding);
        nucImage = imread(imagePathNucleusBig);
        totalImage(startY:startY+sizeY-1,startX:startX+sizeX-1) = nucImage;
        if(mod(i,3)==0)
            %New line
            startY = startY + sizeY;
            startX = 1;
        else
            startX = startX + sizeX;
        end
    end    

% --------------------------------------------------------------------
function MultiParamAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to MultiParamAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
MultiParamAnalysis(handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --------------------------------------------------------------------
function PlotMultiParamAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to PlotMultiParamAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
PlotMultiParamAnalysisResults(9,handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
    


% --------------------------------------------------------------------
function CalcDensDistParam_Callback(hObject, eventdata, handles)
% hObject    handle to CalcDensDistParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
controlWells = inputdlg('Control Wells separated by ;');
controlWells = controlWells{1};
SetDensDistParam(controlWells,handles);
ReloadWellList(handles,0);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --------------------------------------------------------------------
function ExportManualStatistics_Callback(hObject, eventdata, handles)
% hObject    handle to ExportManualStatistics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
%Iterate over manual counted Neurons
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
neuronHandler = handles.NeuronCoordinates;
csvHandler = handles.CSVCoordinates;
imageHandler = handles.ImageHandler;
wellList = get(handles.lbWell, 'string');
foldername = imageHandler.Foldername;
exportText = sprintf('Well;X;Y;Neurite Detected(1=Neurite on Position, 2=Neurite nearby 0=No Neurite found);Detection Method Skeleton(Manual Neurite = 1, Skeleton Neurite = 2, Both = 3);Composite Fill Detected? (1=Yes, 0=No);NeuriteLength;# Branching Points;Neurite Area\r\n');
for i=optionHandler.StartWell:numel(wellList)
    disp(num2str(i));
    selectedWell = wellList{i};
    NucleusM = csvHandler.CellPosMatrix(selectedWell);
    NeuronM = neuronHandler.CellPosMatrix(selectedWell);
    CompositeFillM = logical(neuronHandler.NeuronPositionsEdgeFill(selectedWell));
    if(length(selectedWell) == 2)
      selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    elseif(length(selectedWell) == 4)
        selectedWellLong = [selectedWell(1) '0' selectedWell(2:length(selectedWell))];
    else
        selectedWellLong = selectedWell;
    end
    %Iterate over all manual counted Neurons for current Well
    NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
    SkelM = neuronHandler.NeuronPositionsSkeletonization(selectedWell);
    [ManualIndicesY ManualIndicesX] = find(NeuronManualM);
    [SkeletonIndicesY SkeletonIndicesX] = find(SkelM);
    
    %Prepare nucleus and neurite pic like within Skeletonization
    imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
    imagePathNucleusSmall = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusSmall' optionHandler.FileEnding];
    imagePathNeuriteBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NeuriteBig' optionHandler.FileEnding];
    imagePathNeuriteSmall = [foldername '/ConvertedCellomics/' selectedWellLong 'NeuriteSmall' optionHandler.FileEnding];
    NeuriteImage = imread(imagePathNeuriteBig);
    [sizeYorig sizeXorig] = size(NeuriteImage);
    [sizeY sizeX] = size(NeuriteImage);
    %if(sav)
        BinaryImage = logical(zeros(sizeYorig, sizeXorig));
        SkeletonImage = logical(zeros(sizeYorig, sizeXorig));
    %end
    %else
    %    BinaryImage(yStartPos:yEndPos,xStartPos:xEndPos) = logical(zeros(sizeY, sizeX));
    %    SkeletonImage(yStartPos:yEndPos,xStartPos:xEndPos) = logical(zeros(sizeY, sizeX));
    %end
    %Get density distribution to exclude too dense areas
    %filepath = [imageHandler.Foldername '/ConvertedCellomics/' selectedWell 'BinaryCut.tif'];
    %if(~exist(filepath,'file'))
    if(optionHandler.CutOutSphereCore == 1)
        brightNeurites = CutOutCircles(NeuriteImage,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM, NeuronM,0);
    else
        brightNeurites=NeuriteImage;
    end
    NucleusImage = imread(imagePathNucleusBig);
    if(optionHandler.CutOutSphereCore == 1)
        nucleusPic = CutOutCircles(NucleusImage,selectedWell,1,0,1,optionHandler, foldername, sizeY, sizeX, NucleusM, NeuronM,0);
    else
        nucleusPic=NucleusImage;
    end
    nucleusPic(nucleusPic<optionHandler.NucleusThreshold) = 0;%double(double(imageHandler.NucleusImage) ./ 255);
    nucleusPic(nucleusPic>=optionHandler.NucleusThreshold) = 1;
    nucleusPic = logical(nucleusPic);
    nucleusPic = logical(xor(bwareaopen(nucleusPic,1),  bwareaopen(nucleusPic,15000)));
    nucleusPicBinBigNuclei = xor(bwareaopen(nucleusPic,250),  bwareaopen(nucleusPic,15000));
    D = bwdist(~nucleusPicBinBigNuclei);
    D = -D;
    D(~nucleusPic) = -Inf;
    L = logical(watershed(D));
    nucleusPic(find(~logical(L))) = 0;
    nucleusPic = xor(bwareaopen(nucleusPic,35),  bwareaopen(nucleusPic,15000));
    stretchlimLow = neuronHandler.StretchlimResult;
    [BW BWStrong lMin] = ThresholdPic(brightNeurites,optionHandler,stretchlimLow,sizeYorig, sizeXorig,0, 'Neurite');
    %BW = imfill(BW);
    %ToDo: Fill holes and imdilate
    SE = strel('disk', 1);
    BW = imdilate(BW,SE);
    %BW=BW(yStartPos:yEndPos,xStartPos:xEndPos);
    %BWStrong=BWStrong(yStartPos:yEndPos,xStartPos:xEndPos);
    %nucleusPic = nucleusPic(yStartPos:yEndPos,xStartPos:xEndPos);
    %ToDo: Filtere Hintergrundrauschen aus Neuritenbild heraus.
    %Idee: Errechne hell zu dunkel-Verh�ltnis f�r verschiedene Bereiche des
    %Neuritenbildes (Schwarze Fl�chen, graue Fl�chen und wei�e Fl�chen
    %1. hist3 mit sehr hellen Fl�chen
    %2. hist3 mit sehr dunlen Fl�chen
    % -> Wenn beide hist3 hinreichend klein -> Setze Neuritenbild auf schwarz
    % f�r den jeweiligen Bereich
    % -> Wenn beide hist3 gro� genug sind, handelt es sich um valide Neuriten  
    [thresholdedRowsBig thresholdedColsBig] = find(NeuriteImage > 25);
    [thresholdedRowsSmall thresholdedColsSmall] = find(NeuriteImage > 11);
    thresholdedRowsBig = [thresholdedRowsBig;sizeYorig];
    thresholdedColsBig = [thresholdedColsBig;sizeXorig];
    thresholdedRowsBig = [thresholdedRowsBig;0];
    thresholdedColsBig = [thresholdedColsBig;0];
    thresholdedRowsSmall = [thresholdedRowsSmall;sizeYorig];
    thresholdedColsSmall = [thresholdedColsSmall;sizeXorig];
    thresholdedRowsSmall = [thresholdedRowsSmall;0];
    thresholdedColsSmall = [thresholdedColsSmall;0];
    BrightDensity = (hist3([double(thresholdedRowsBig),double(thresholdedColsBig)],[128,128]));
    LowDensity = (hist3([double(thresholdedRowsSmall),double(thresholdedColsSmall)],[128,128]));
    %LowDensity = imcomplement(LowDensity);
    %ToDo: Save Neurite objects with single pictures to trace
    %Check if already a Neuron by EdgeFillNeurons in this area available
    %If not, trace Neurite and check if Neuron possible on one or
    %another side  

    %The same maybe faster:
    L = uint16(bwlabel(BW,8));
    NotList = [0];
    while 1
        %Get label with highest occurance in L:
        mostOccuredLabel = mode(L(~ismember(L,NotList)));
        NotList = [NotList mostOccuredLabel];
        currentNeuriteArea = logical(sparse(zeros(size(BW,1), size(BW,2))));
        currentNeuriteArea(find(L == mostOccuredLabel))=1;
        areaSize = bwarea(full(currentNeuriteArea));
        if(areaSize > optionHandler.MinSizeNeuriteArea)
            %Apply harder threshold in relevant part of BW
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
            %Try to split up neurite
            NeuriteSubRegion = NeuriteImage(newNeurite.cutYPosStart-1:newNeurite.cutYPosEnd-1,newNeurite.cutXPosStart-1:newNeurite.cutXPosEnd-1);
            NeuriteSubRegion = imcomplement(NeuriteSubRegion);
            subHist = imhist(NeuriteSubRegion);
            thresholdSub = IJIsoData(subHist);
            currentNeuriteArea = im2bw(NeuriteSubRegion,thresholdSub/255);
            BW(newNeurite.cutYPosStart:newNeurite.cutYPosEnd,newNeurite.cutXPosStart:newNeurite.cutXPosEnd) = 1-currentNeuriteArea;
        else
            break;
        end
    end
    %Again: Remove to small areas of BW and Cut out Circles again
    LB = 75;
    %UB = 4000;
    %CuttedIndices = BW;
    %BW = xor(bwareaopen(BW,LB),  bwareaopen(BW,UB));

    BW = bwareaopen(BW,LB);
    if(optionHandler.CutOutSphereCore == 1)
        BW = logical(CutOutCircles(BW,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM, NeuronM,0));
    else
        BW=logical(BW);
    end
    %CC = bwconncomp(BW);
    Label = uint16(bwlabel(BW));    
    [D, IDX]=bwdist(BW);
    IDX = uint32(IDX);
    D=0;
    for(j=1:numel(ManualIndicesY))
        %Check for each Neuron:
        %If it has Neurite below
        xKern=ManualIndicesX(j);
        yKern=ManualIndicesY(j);
        x=xKern;
        y=yKern;
        
        if(BW(y,x) > 0)  
            exportText=SkeletonizeForBranchingPointAnalysis(BW,BWStrong,Label,sizeY,sizeX,sizeYorig,sizeXorig,y,x,exportText,optionHandler,selectedWell,1,BrightDensity,LowDensity,SkelM,1,CompositeFillM,yKern,xKern);
        else
            %If no Neurite below, check for Neurite nearby.                        
            [y,x] = ind2sub([sizeY,sizeX],IDX(sub2ind([sizeY,sizeX],y,x)));
            exportText=SkeletonizeForBranchingPointAnalysis(BW,BWStrong,Label,sizeY,sizeX,sizeYorig,sizeXorig,y,x,exportText,optionHandler,selectedWell,2,BrightDensity,LowDensity,SkelM,1,CompositeFillM,yKern,xKern);
            %exportText = [exportText selectedWell ';' num2str(x) ';' num2str(y) ';' num2str(2) sprintf(';;;;\r\n')];
        end
    end
    
    
    
    for(j=1:numel(SkeletonIndicesY))
        %Check for each Neuron:
        %If it has Neurite below
        xKern=SkeletonIndicesX(j);
        yKern=SkeletonIndicesY(j);
        x=xKern;
        y=yKern;
        
        if(BW(y,x) > 0)  
            exportText=SkeletonizeForBranchingPointAnalysis(BW,BWStrong,Label,sizeY,sizeX,sizeYorig,sizeXorig,y,x,exportText,optionHandler,selectedWell,1,BrightDensity,LowDensity,NeuronManualM,2,CompositeFillM,yKern,xKern);
        else
            %If no Neurite below, check for Neurite nearby.                        
            [y,x] = ind2sub([sizeY,sizeX],IDX(sub2ind([sizeY,sizeX],y,x)));
            exportText=SkeletonizeForBranchingPointAnalysis(BW,BWStrong,Label,sizeY,sizeX,sizeYorig,sizeXorig,y,x,exportText,optionHandler,selectedWell,2,BrightDensity,LowDensity,NeuronManualM,2,CompositeFillM,yKern,xKern);
            %exportText = [exportText selectedWell ';' num2str(x) ';' num2str(y) ';' num2str(2) sprintf(';;;;\r\n')];
        end
    end
end
[FileName,PathName] = uiputfile('NeuronStatistics.csv');
fileID = fopen(strcat(PathName,'/',FileName),'w');
fprintf(fileID,'%s',exportText);
fclose(fileID);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


function [NeuriteLengthMatrix NeuriteObjectMatrix] = SkeletonizeForBranchingPointAnalysis(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,NeuriteLengthMatrix,foldername,AnalyzeM,RawString,NeuriteObjectMatrix)
NucleusM = logical(NucleusM);
mkdir([foldername '/ConvertedCellomics/' selectedWell '_' RawString]);
SkelDeletedMat = logical(sparse(double(sizeY), double(sizeX)));
SkelNeurons = logical(sparse(double(sizeY), double(sizeX)));
NeuronPositionsEdgeFillNeurite = logical(sparse(double(sizeY), double(sizeX)));
ringNumber = optionHandler.DensityDistributionRingNumber;
selectedWellLong=selectedWell;
neuriteAreaCount=0;
neuriteCounter=0;
TP = 0;
FP = 0;
TN=0;
yStartPos=1;
xStartPos=1;
yEndPos=sizeY;
xEndPos=sizeX;
%NucleusM = NucleusM(yStartPos:yEndPos,xStartPos:xEndPos);
TruePositives = zeros(0,0);
%Just for testing purposes.
mkdir([foldername '/ConvertedCellomics/'],selectedWell);
%Get density distribution to exclude too dense areas
%if(saveCircle)     
densityWidth = optionHandler.MigrationDistanceDensityImageXSize;
densityHeight = optionHandler.MigrationDistanceDensityImageYSize;
[filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, densityWidth, densityHeight, optionHandler, NucleusM, foldername, sizeY, sizeX); 

innerMask = ones(uint16(sizeY/(1/ScalingFactor)),uint16(sizeX/(1/ScalingFactor)));

if(markerPointCoordinates ~= 0)
    mk = markerPointCoordinates('10')./(1/ScalingFactor);
    innerMask = roipoly(innerMask,mk(:,1),mk(:,2));
    currentRing =logical(ones(int32(sizeY/(1/ScalingFactor)),int32(sizeX/(1/ScalingFactor))) - innerMask);
    currentRing =imresize(currentRing,[sizeY,sizeX]);

    if(isKey(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
        ManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
    else
        ManualM = sparse(logical(zeros(sizeY,sizeX))); 
        neuronHandler.ManualNeuronPositionsSparse(selectedWell) =ManualM;
    end
    
    
    
    if(yStartPos<=1)
        ManualM = ManualM.*currentRing;
        AnalyzeM= AnalyzeM.*currentRing;
    end
    ManualM = ManualM(yStartPos:yEndPos,xStartPos:xEndPos);
    TN = nnz(ManualM);
    [D Ind]= bwdist(full(logical(NucleusM)));
    
    if(length(selectedWell) == 2)
              selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    end
    imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
    imagePathNucleusSmall = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusSmall' optionHandler.FileEnding];
    imagePathNeuriteBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NeuriteBig' optionHandler.FileEnding];
    imagePathNeuriteSmall = [foldername '/ConvertedCellomics/' selectedWellLong 'NeuriteSmall' optionHandler.FileEnding];
    NeuriteImage = imread(imagePathNeuriteBig);
    BinaryImage = logical(zeros(sizeY, sizeX));
    SkeletonImage = logical(zeros(sizeY, sizeX));

    if(optionHandler.CutOutSphereCore == 1)
        brightNeurites = CutOutCircles(NeuriteImage,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM, 0);
    else
        brightNeurites=NeuriteImage;
    end

    %Implement Top hat filter with rolling ball
    %SE = strel('ball',8,7);
    %brightNeurites = imtophat(brightNeurites,SE);

    stretchlimLow = neuronHandler.StretchlimResult;
    [BW BWStrong lMin] = ThresholdPic(brightNeurites,optionHandler,stretchlimLow,sizeY, sizeX,0, 'Neurite');
    %BW = imfill(BW);
    %ToDo: Fill holes and imdilate
    SE = strel('disk', 1);
    BW = imdilate(BW,SE);
    BW=BW(yStartPos:yEndPos,xStartPos:xEndPos);
    BWStrong=BWStrong(yStartPos:yEndPos,xStartPos:xEndPos);
    NucleusImage = imread(imagePathNucleusBig);
    if(optionHandler.CutOutSphereCore == 1)
        nucleusPic = CutOutCircles(NucleusImage,selectedWell,1,0,1,optionHandler, foldername, sizeY, sizeX, NucleusM, 0);
    else
        nucleusPic=NucleusImage;
    end
    nucleusPic(nucleusPic<optionHandler.NucleusThreshold) = 0;%double(double(imageHandler.NucleusImage) ./ 255);
    nucleusPic(nucleusPic>=optionHandler.NucleusThreshold) = 1;
    nucleusPic = logical(nucleusPic);
    nucleusPic2 = nucleusPic(yStartPos:yEndPos,xStartPos:xEndPos);
    overlayAreaIndices = find(nucleusPic2 > 0 & BW > 0);
    overlayArea=logical(zeros(sizeY,sizeX));
    overlayArea(overlayAreaIndices)=1;

    nucleusPic = logical(xor(bwareaopen(nucleusPic,1),  bwareaopen(nucleusPic,15000)));
    nucleusPicBinBigNuclei = xor(bwareaopen(nucleusPic,250),  bwareaopen(nucleusPic,15000));
    D = bwdist(~nucleusPicBinBigNuclei);
    D = -D;
    D(~nucleusPic) = -Inf;
    L = watershed(D);
    nucleusPic(find(~logical(L))) = 0;
    nucleusPic = xor(bwareaopen(nucleusPic,35),  bwareaopen(nucleusPic,15000));
    nucleusPic = nucleusPic(yStartPos:yEndPos,xStartPos:xEndPos);

    %ToDo: Save Neurite objects with single pictures to trace
    %Check if already a Neuron by EdgeFillNeurons in this area available
    %If not, trace Neurite and check if Neuron possible on one or
    %another side  

    %The same maybe faster:
    L = uint16(bwlabel(BW,8));
    NotList = [0];

    %Again: Remove to small areas of BW and Cut out Circles again
    LB = 75;
    %UB = 4000;
    %CuttedIndices = BW;
    %BW = xor(bwareaopen(BW,LB),  bwareaopen(BW,UB));

    BW = bwareaopen(BW,LB);
    if(optionHandler.CutOutSphereCore == 1)
        BW = logical(CutOutCircles(BW,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM, 0));
    else
        BW=logical(BW);
    end
    [DistanceNeurites IndNeurites] = bwdist(BW);
    labeledNeurites = uint16(bwlabel(BW));

    %[whitePointRows, whitePointCols] = find(BW);
    alreadyVisited = logical(ones(size(BW,1), size(BW,2)));
    neuriteAreasDict = containers.Map();
    i=0;
    %BW2= logical(zeros(size(BW,1), size(BW,2)));
    secondThreshold = 20;
    exportText='abc';    
    
    %Now take Cellomics Neurons and look for next Neurite starting there.
    CellomicsNeuronIndices = find(AnalyzeM);
    
    posYList=zeros(numel(CellomicsNeuronIndices),1);
    posXList=zeros(numel(CellomicsNeuronIndices),1);
    neuronIndicesList = zeros(numel(CellomicsNeuronIndices),2);
    %DEBUG START
    for(i=1:numel(CellomicsNeuronIndices))
       [posY posX] = ind2sub([7168 7168],CellomicsNeuronIndices(i));
       posYList(i)=posY;
       posXList(i)=posX;
       neuronIndicesList(i,1) = CellomicsNeuronIndices(i); 
    end    
    %DEBUG END
    
    
    for(i=1:numel(CellomicsNeuronIndices))        
        %DEBUG START
        disp(num2str(i));        
        [currentY currentX] = ind2sub([sizeY sizeX],CellomicsNeuronIndices(i));
        nucY=currentY;
        nucX = currentX;
        neuriteInd = IndNeurites(CellomicsNeuronIndices(i));
        currentNeuriteArea = logical(zeros(sizeY,sizeX));
        %Get whole area (FloodFill)
        labeledInd = labeledNeurites(neuriteInd);
        currentImage=logical(zeros(sizeY,sizeX));
        currentNeuriteArea(labeledNeurites==labeledInd)=1;
        
        newNeuriteBig = Neurite();
        [rowIndices colIndices] = find(currentNeuriteArea);
        yIndex = min(rowIndices);
        xIndex = min(colIndices);
        cutYPosStartRef = min(rowIndices);
        cutXPosStartRef = min(colIndices);
        newNeuriteBig.cutYPosStart = min([rowIndices;currentY]) - 2;
        newNeuriteBig.cutXPosStart = min([colIndices;currentX]) - 2;    

        newNeuriteBig.cutYPosEnd = max([rowIndices;currentY]) + 2;
        newNeuriteBig.cutXPosEnd = max([colIndices;currentX]) + 2;

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
        currentNeuriteAreaSmall=currentNeuriteArea(newNeuriteBig.cutYPosStart:newNeuriteBig.cutYPosEnd,newNeuriteBig.cutXPosStart:newNeuriteBig.cutXPosEnd);
        currentNeuriteAreaSmall=1-currentNeuriteAreaSmall;
        yMid = mean([newNeuriteBig.cutYPosStart newNeuriteBig.cutYPosEnd]);
        xMid = mean([newNeuriteBig.cutXPosStart newNeuriteBig.cutXPosEnd]);
        neuriteIndex = sub2ind([sizeY sizeX], yIndex, xIndex);
        NucAreaSmall = nucleusPic(newNeuriteBig.cutYPosStart:newNeuriteBig.cutYPosEnd,newNeuriteBig.cutXPosStart:newNeuriteBig.cutXPosEnd);  
        SubNucM = AnalyzeM(newNeuriteBig.cutYPosStart:newNeuriteBig.cutYPosEnd,newNeuriteBig.cutXPosStart:newNeuriteBig.cutXPosEnd);
        SE = strel('disk', 3);
        %NucAreaSmallDil = imdilate(NucAreaSmall,SE);
        %AreaResultSmallLabel = bwlabel(NucAreaSmallDil);
        AreaResultSmallLabel = bwlabel(NucAreaSmall);
        [sizeYskel sizeXskel] = size(NucAreaSmall);
        AreaResultSmall = logical(zeros(sizeYskel,sizeXskel));
        currentY=currentY - newNeuriteBig.cutYPosStart-1;
        currentX=currentX - newNeuriteBig.cutXPosStart-1;
        if(currentY==0)
            currentY=1;
        end
        if(currentX==0)
            currentX=1;
        end
        
        
        %ToDo: All other Cellomics Nuclei have to be 1 also in Nucleus
        %Image. Therefore: Check on which NucleiAreaLabels in current area,
        %there is a regarding Cellomics Neuron
        [CellomicsNeuronsInAreaY CellomicsNeuronsInAreaX] = find(AnalyzeM(newNeuriteBig.cutYPosStart:newNeuriteBig.cutYPosEnd,newNeuriteBig.cutXPosStart:newNeuriteBig.cutXPosEnd));
        %Get number of TPs, FPs and TNs in Area
        AnalyzeMSmall=AnalyzeM(newNeuriteBig.cutYPosStart:newNeuriteBig.cutYPosEnd,newNeuriteBig.cutXPosStart:newNeuriteBig.cutXPosEnd);                
        ManualMSmall=ManualM(newNeuriteBig.cutYPosStart:newNeuriteBig.cutYPosEnd,newNeuriteBig.cutXPosStart:newNeuriteBig.cutXPosEnd);                
        TPs = numel(find(AnalyzeMSmall.*ManualMSmall));
        FPs = numel(find(AnalyzeMSmall.*(~ManualMSmall)));
        TNs = numel(find(ManualMSmall.*(~AnalyzeMSmall)));
        
        [DD IIDX] = bwdist(AreaResultSmallLabel);
        [currentY currentX] = ind2sub([sizeYskel sizeXskel],IIDX(sub2ind([sizeYskel sizeXskel],currentY,currentX)));        
        nucLabel=AreaResultSmallLabel(currentY,currentX); %yStart und xStart: take positions of next given Cellomics Nucleus        
        
        if(nucLabel~=0)
            AreaResultSmall(AreaResultSmallLabel~=nucLabel)=0;
            AreaResultSmall(AreaResultSmallLabel==nucLabel)=1;
        end
        for(zz=1:numel(CellomicsNeuronsInAreaY))            
            %Same as above: Look for index!
            [curY curX] = ind2sub([sizeYskel sizeXskel],IIDX(sub2ind([sizeYskel sizeXskel],CellomicsNeuronsInAreaY(zz),CellomicsNeuronsInAreaX(zz))));        
            nucLabel=AreaResultSmallLabel(curY,curX);
            if(nucLabel~=0)
                AreaResultSmall(AreaResultSmallLabel==nucLabel)=1;
            end
        end
        newNeuriteBig.nucleusImage = AreaResultSmall;
        %Skeletonize
        newNeuriteBig.image = currentNeuriteAreaSmall;
        ignoreNeurite=0;
        if(sizeYskel > 550 || sizeXskel > 550)
            ignoreNeurite=1;
            newNeuriteBig.skeletonImage=zeros(sizeYskel,sizeXskel);
        else
            try
                newNeuriteBig.CalculateBaiSkeleton(16);
                newNeuriteBig.FindBranchingCrossingPoints();
                newNeuriteBig.KillShortBranches();
                newNeuriteBig.CalculateNeuriteLength();
            catch
                ignoreNeurite=1;
                newNeuriteBig.skeletonImage=zeros(sizeYskel,sizeXskel);
            end
        end
        
        newNeuriteBig.nucleusImage = AreaResultSmall;        
        %Same Postprocessing as in Skeletonization Algorithm
        oldSkel = newNeuriteBig.skeletonImage;            
      
        bothIndices = find(logical(newNeuriteBig.skeletonImage) & logical(newNeuriteBig.nucleusImage));

        %Doppelt gemoppelt h�lt besser :-D
        if(ignoreNeurite==0)
        newNeuriteBig.FindBranchingCrossingPoints2();           
        newNeuriteBig.FindBranchingCrossingPoints2(); 
        %Finde �berdeckungsl�nge zwischen Nucleus und Neuritenbild

        %minNeuriteLength = optionHandler.SkeletonMinNeuriteLength - numel(bothIndices);
        %if(minNeuriteLength < 11)
            minNeuriteLength=11;
        %end
        
        subtractImage = newNeuriteBig.skeletonImage - newNeuriteBig.nucleusImage;
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
        newNeuriteBig.skeletonImage(labelSubtract==maxIndex)=1;
        subtractImage(labelSubtract==maxIndex) = 1;
        %ToDo: Ensure that at least one subtractImageLabel remains and 
        %ensure that bothIndices are counted before finding new
        %BranchingCrossingPoints.


        %Do this only if there is at least no Nucleus Point, touching a
        %Neurite Point

        [DSubtract,IDX] = bwdist(subtractImage);
        [DNuc,IDXNuc] = bwdist(newNeuriteBig.nucleusImage);
        neuriteCounter=neuriteCounter+1;

        %Get minimum of DSubtract on Indices where
        %newNeuriteBig.nucleusImage is 1
        
        minDist = min(DSubtract(newNeuriteBig.nucleusImage==1));
        if(minDist >= 2)
            %1. Get mid of Nucleus and find shortest path to next Neurite
            [Ay Ax]= find(newNeuriteBig.nucleusImage>=1);
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
        newNeuriteBig.skeletonImage = subtractImage;
        %newNeuriteBig.FindBranchingCrossingPoints2();           
        end


        %Create figure with overlayed Nucleus, Skeleton and original
        %Neurite in /SkeletonExport/selectedWell/ID.png
        
        currentObject = zeros(15,20);
        if(numel(find(newNeuriteBig.skeletonImage==0))==0)
             newNeuriteBig.skeletonImage=zeros(sizeYskel,sizeXskel);
             ignoreNeurite=1;
        end
        %if(ignoreNeurite==0)
            newNeuriteBig.skeletonImage = parsiSkel(newNeuriteBig.skeletonImage);
            branchImage = logical(bwmorph(newNeuriteBig.skeletonImage, 'branchpoints'));
            endpointImage = logical(bwmorph(newNeuriteBig.skeletonImage, 'endpoints'));
            %2 Endpoints: No branches, maybe just circles
            newNeuriteBig.xBranchingPoints = zeros(0);
            newNeuriteBig.yBranchingPoints = zeros(0);
            [yBranchingPoints xBranchingPoints] = find(branchImage==1);
            [yBranchingPoints xBranchingPoints] = KillCloseBranches(yBranchingPoints,xBranchingPoints);
            newNeuriteBig.yBranchingPoints = yBranchingPoints;
            newNeuriteBig.xBranchingPoints = xBranchingPoints;
            
            %Reset skeleton image if it is empty.
            
               
            
            newNeuriteBig.skeletonImage = parsiSkel(newNeuriteBig.skeletonImage);


            newNeuriteBig.nucleusImage = AreaResultSmall;
            %figure(3);
            %imshow(imfuse(newNeuriteBig.skeletonImage,newNeuriteBig.nucleusImage));
            exportRGB = zeros(sizeYskel,sizeXskel,3);
            exportRGB(:,:,3) = 1-newNeuriteBig.image;
            exportRGB(:,:,1) = newNeuriteBig.nucleusImage;
            exportRGB(:,:,2) = newNeuriteBig.skeletonImage;
            %figure(1);
            %imshow(exportRGB);
            
            %positionID has to be mid of Neurite in big picture
            
            positionID = sub2ind([sizeY sizeX],cutYPosStartRef, cutXPosStartRef);
            %Check if Neurite is already available in Dictionary            
            connections=SplitUpNeuronsAndNuclei(SubNucM,newNeuriteBig.nucleusImage,newNeuriteBig.skeletonImage,newNeuriteBig,sizeY,sizeX);
            
            %Calculate and insert distance for every neuron from
            %sphere core
            mkFirst = markerPointCoordinates('10');
            %Get comma separated Nucleus IDs
            for(t=1:size(connections,1))
                nucIDXList = connections(t,4:end);
                %Split list up by ","
                u=1;                
                globalID = nucIDXList(u);
                while((u<numel(nucIDXList)) && (globalID ~= 0))
                    [nucY nuX]=ind2sub([sizeY sizeX],globalID);
                    if(~isKey(NeuriteObjectMatrix,[selectedWell ';' num2str(globalID)]))
                        %Calculate distance from sphere core
                        %Iterate over all markerpointcoordinates
                        minDist = 99999999;
                        minIndex = 0;
                        for(v=1:size(mkFirst,1))
                            currentYSphereCore = mkFirst(v,1);
                            currentXSphereCore = mkFirst(v,2);
                            distance = pdist([nucY nuX;currentYSphereCore currentXSphereCore]);
                            if (distance < minDist)
                                minDist=distance;
                                minIndex = v;
                            end                            
                        end
                        NeuriteObjectMatrix([selectedWell ';' num2str(globalID)]) = minDist;                            
                    end
                    u=u+1;
                    globalID = nucIDXList(u);
                end
            end
            
            
            labelSubtract = bwlabel(newNeuriteBig.skeletonImage);   
            alreadyThere=0;
            for(t=1:size(connections,1))
                %Save to Excel:
                %Values for Sub Neurite
                %0. Neurite ID
                currentObject(1,1)=connections(t,1);
                %1. Subneurite Length
                currentObject(2,1)=connections(t,2);
                %2. Subeurite Branching Points
                currentObject(3,1)=connections(t,3);
                %3. Comma separated Nucleus IDs
                currentObject(4,1:end)=connections(t,1:end);
                %3. Values for whole Neurite
                %1. Number Nuclei
                currentObject(5,1) = TPs+FPs;
                %2. Number total Branching Points
                currentObject(6,1) = numel(newNeuriteBig.yBranchingPoints);
                %3. Neurite length
                currentObject(7,1) = nnz(newNeuriteBig.skeletonImage);
                %4. Number of distinct Neurons per Skeleton
                currentObject(8,1) = max(labelSubtract(:));
                %ID shown to the user
                currentObject(9,1) = neuriteCounter;
                %5. Export skeleton picture with ID
                %ToDo: Calculate number of TPs, FPs and TNs on current
                %Neurite.
                %1. Get all manual positions on Neurite Area
                currentObject(10,1) = TPs;
                currentObject(11,1) = FPs;
                currentObject(12,1) = TNs;
                currentObject(13,1) = bwarea(newNeuriteBig.nucleusImage);
                currentObject(14,1) = ignoreNeurite;
                currentObject(15,1) = CellomicsNeuronIndices(i);

                %DEBUG START
                neuronIndicesList(i,2)=1;
                %DEBUG END

                %2. Get all positions of current method on Neurite Area
                %3. Check which are identical

                %Furthermore: Get total size of all Nuclei, detected by
                %current method.
                if(~isKey(NeuriteLengthMatrix,[selectedWell ';' num2str(positionID) ';' num2str(t)]))
                    NeuriteLengthMatrix([selectedWell ';' num2str(positionID) ';' num2str(t)]) = currentObject;
                else
                    obj = NeuriteLengthMatrix([selectedWell ';' num2str(positionID) ';' num2str(t)]);
                    position = nnz(obj(15,:))+1;
                    obj(15,position) = CellomicsNeuronIndices(i);
                    NeuriteLengthMatrix([selectedWell ';' num2str(positionID) ';' num2str(t)]) = obj;
                 %Add Nucleus to list
                    alreadyThere=1;
                end
            end
            if(~alreadyThere)
                imwrite(exportRGB,[foldername '/ConvertedCellomics/' selectedWell '_' RawString '/' num2str(positionID) '.png']);
            end        
    end
end


% --- Executes on button press in cbPreprocessedImages.
function cbPreprocessedImages_Callback(hObject, eventdata, handles)
% hObject    handle to cbPreprocessedImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
% Hint: get(hObject,'Value') returns toggle state of cbPreprocessedImages
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on button press in cbOligoPicture.
function cbOligoPicture_Callback(hObject, eventdata, handles)
% hObject    handle to cbOligoPicture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
% Hint: get(hObject,'Value') returns toggle state of cbOligoPicture
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --- Executes on button press in cbAstroPicture.
function cbAstroPicture_Callback(hObject, eventdata, handles)
% hObject    handle to cbAstroPicture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
if(imageHandler.ZoomState)
    RefreshZoomedImage(handles);
else
    RefreshUnzoomedImage(handles);
end
% Hint: get(hObject,'Value') returns toggle state of cbAstroPicture
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --------------------------------------------------------------------
function CheckSphereCoreAvailability_Callback(hObject, eventdata, handles)
% hObject    handle to CheckSphereCoreAvailability (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
%Check for each well if Sphere core is available
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
csvHandler = handles.CSVCoordinates;
handles = guidata(handles.figure1);
wellList = get(handles.lbWell, 'string');
for i=optionHandler.StartWell:numel(wellList)
    selectedWell=wellList(i);
    selectedWell=selectedWell{1};
    markerPointCoordinates=-1;
    filterDistance = -1;
    nonFilterDistance=-1;
    imageHandler = handles.ImageHandler;
    csvHandler = handles.CSVCoordinates;
    neuronHandler = handles.NeuronCoordinates;
    [sizeY sizeX] = size(imageHandler.NucleusImage);
    SphereArea=-1;
    optionHandler = handles.OptionHandler;
    path = [imageHandler.Foldername '/ConvertedCellomics/MarkerPointCoordinates-' selectedWell '.mat'];
    ringNumber = optionHandler.DensityDistributionRingNumber;
    NucleusM = logical(csvHandler.CellPosMatrix(selectedWell));
    SphereAreaSizeX = optionHandler.MigrationDistanceDensityImageXSize;
    SphereAreaSizeY = optionHandler.MigrationDistanceDensityImageYSize;
    [filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, optionHandler, NucleusM, imageHandler.Foldername, sizeY, sizeX); 

    if(markerPointCoordinates~=0)
        hInner = impoly(handles.axes2,double(markerPointCoordinates('0')./(1/ScalingFactor)));
        innerMask = createMask(hInner);
        disp([selectedWell ': ' num2str(nnz(innerMask))]);    
    else
        disp([selectedWell ': ' num2str(0)]);    
    end
    delete(hInner);
end
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --------------------------------------------------------------------
function ExportNeuriteLengthStat_Callback(hObject, eventdata, handles)
% hObject    handle to ExportNeuriteLengthStat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Iterate over Neurite length Matrix
try
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
neuronHandler = handles.NeuronCoordinates;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
foldername = [imageHandler.Foldername '/ConvertedCellomics'];

counter=0;
for(j=optionHandler.StartWell:numel(wellList))
    exportText = sprintf('Well;Neurite ID;SubNeurite ID;SubNeurite Length;# Branching Points Subneurite;Connections To Nuclei IDs;#Nuclei;#Branching Points Total;Total Neurite Length;#Distinct Neurites per Skeleton;\r\n');
    counter=counter+1;
    selectedWell=wellList(j);
    selectedWell=selectedWell{1};    
    mapObj = neuronHandler.NeuriteLengthMatrix(selectedWell);
    keyset = mapObj.keys;
    keyOne = keyset(1);
    keyOne=keyOne{1};
    mapObj=mapObj(keyOne);
    keyset = mapObj.keys;
    for(i=1:numel(keyset))
        currentKey = keyset(i);
        currentKey = currentKey{1};
        currentValue = mapObj(currentKey);
        %Column 3:
        %Get last NonZero element.
        [lastNonZeroNumber] = max(find(currentValue(4,1:end)));
        exportText = [exportText currentKey ';' num2str(currentValue(2)) ';' num2str(currentValue(3)) ';' num2str(currentValue(4,4:lastNonZeroNumber)) ';' num2str(currentValue(5)) ';' num2str(currentValue(6)) ';' num2str(currentValue(7)) ';' num2str(currentValue(8)) ';' sprintf('\r\n')];
    end
    fileID = fopen(strcat(foldername,'/',selectedWell,'DataOmnisphero.csv'),'w');
    fprintf(fileID,'%s',exportText);
    fclose(fileID);
end
    catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


function ExportNeuriteLengthStat(AnalyzeDict,exportname,figtitle,answer,handles)
handles = guidata(handles.figure1);
QuickFixManualM(handles);
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
neuronHandler = handles.NeuronCoordinates;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};



% %1. Calculate
NucleusM = csvHandler.CellPosMatrix(selectedWell);
foldername = imageHandler.Foldername;
[sizeY sizeX] = size(imageHandler.NeuriteImage);
NeuriteLengthMatrix = containers.Map();
NeuriteObjectMatrix = containers.Map();

 parallel.defaultClusterProfile('local');
maxCores =optionHandler.NumberOfCores;
% if(maxCores>1)
    c = parcluster();
    c.NumWorkers = maxCores;
    job = createJob(c,'JobData',optionHandler);
% else
%     out = zeros(0);
%     outNew = zeros(0);
% end
 excludeList = optionHandler.ExcludedWells;
% 
% %DEBUG CHANGE
  %for(j=20:numel(wellList))
  waitbarHandle = waitbar(0,'Please wait. I am currently skeletonizing for the Morphology. You can view my progress in the image folder.', 'WindowStyle', 'modal');
  waitbar(0.5,waitbarHandle);
  for(j=optionHandler.StartWell:numel(wellList))
      selectedWell=wellList(j);
      selectedWell=selectedWell{1};
      %selectedWell='E7';
      isKey(AnalyzeDict, selectedWell)
      if(numel(strfind(excludeList,selectedWell)) == 0 && isKey(AnalyzeDict, selectedWell))
            NucleusM = csvHandler.CellPosMatrix(selectedWell);
          %First check if markerpoints are created correctly   
          path = [imageHandler.Foldername '/ConvertedCellomics/MarkerPointCoordinates-' selectedWell '.mat']; 
          AnalyzeM = AnalyzeDict(selectedWell);
          %DEBUG!!!
          %out = SkeletonizeForBranchingPointAnalysis(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,NeuriteLengthMatrix,foldername,AnalyzeM,exportname,NeuriteObjectMatrix);
%           if(maxCores>1)
            createTask(job, @SkeletonizeForBranchingPointAnalysis, 2, {selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,NeuriteLengthMatrix,foldername,AnalyzeM,exportname,NeuriteObjectMatrix});
%           else
%               [NeuriteLengthMatrix NeuriteObjectMatrix] = SkeletonizeForBranchingPointAnalysis(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,NeuriteLengthMatrix,foldername,AnalyzeM,exportname,NeuriteObjectMatrix);
%               out = [out;NeuriteLengthMatrix];
%               outNew = [outNew;NeuriteObjectMatrix];
%           end
      end
  end
  
%   if(maxCores>1)
      submit(job);
      wait(job);
      close(waitbarHandle);
      outComplete=fetchOutputs(job);  
      out = outComplete(:,1);
      outNew = outComplete(:,2);
%   end
%job1013=c.Jobs(3);
%out=fetchOutputs(job1013);

%Add experiment to database
%Fix
%javaclasspath(postgresql-9.4-1203.jdbc42);
props=java.util.Properties;
props.setProperty('user',answer{4});
props.setProperty('password', answer{5});
driver=org.postgresql.Driver;
url=strcat('jdbc:postgresql://',answer{1},':',answer{2},'/',answer{3});
conn=driver.connect(url,props);
if(strcmp(exportname,'DataManual'))
    typeID=1;
elseif(strcmp(exportname,'DataCellomics'))
    typeID=2;
elseif(strcmp(exportname,'DataOmnisphero'))
    typeID=3;
end




csvHandler.NeuriteLengthMatrix=out;
guidata(handles.figure1, handles);
counter=0;
for(j=optionHandler.StartWell:numel(out))
    %Export separately for every well
    counter=counter+1;
     selectedWell=wellList(j);
     selectedWell=selectedWell{1};
    insertExpStr = strcat('INSERT INTO OMNISPHERO_EXPERIMENT (NAME, WELL, TYPE) VALUES(''',figtitle,''',''',selectedWell,''',',num2str(typeID),')');    
    ps=conn.prepareStatement(insertExpStr);
    ps.execute();
    ps.close();
    selectExperimentIDStr='SELECT CURRVAL(''OMNISPHERO_EXPERIMENT_ID_SEQ'')';
    ps=conn.prepareStatement(selectExperimentIDStr);
    rs=ps.executeQuery();
    rs.next();
    experimentID=str2num(char(rs.getString(1)));
    ps.close();
    rs.close();
    exportText = sprintf('Well;Neurite ID;SubNeurite ID;Position ID;SubNeurite Length;# Branching Points Subneurite;Connections To Nuclei IDs;#Nuclei;#Branching Points Total;Total Neurite Length;#Distinct Neurites per Skeleton;#TP;#FP;#TN;Nucleus Area Total\r\n');
    mapObj = out{counter};    
    keyset = mapObj.keys;
    
    
    %DEBUG START
    %Check for every Cellomics Neuron Index, if it is in DB
%     ManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
%     CellomicsNeuronIndices=find(ManualM);
%     for(i=1:numel(CellomicsNeuronIndices))
%         [posY posX] = ind2sub([7168 7168], CellomicsNeuronIndices(i));
%         select = strcat({'SELECT ID FROM OMNISPHERO_NUCLEUS WHERE positiony BETWEEN '},num2str(posY-1),{' AND '},{num2str(posY+1)}, {' AND positionx BETWEEN '},num2str(posX-1),{' AND '},{num2str(posX+1)},{' AND EXPERIMENT_ID=478'});
%         ps=conn.prepareStatement(select);
%         rs=ps.executeQuery();
%         if(rs.next())
%         else
%             test='ERROR';
%         end
%         ps.close();
%         rs.close();
%     end
    
    %DEBUG END
    
    for(i=1:numel(keyset))
        currentKey = keyset(i);
        currentKey = currentKey{1};
        currentValue = mapObj(currentKey);    
        
        %Column 3:
        %Get last NonZero element.
        splitKey = strsplit(currentKey,';');
        well = splitKey{1};
        positionID = splitKey{2};
        positionID = str2num(splitKey{2});
        [positionY positionX] = ind2sub([7168 7168],positionID);
       
        neuriteID = num2str(currentValue(9));
        subID = num2str(splitKey{3});
        [lastNonZeroNumber] = max(find(currentValue(4,1:end)));
        
        %Check if Neurite already exists.
        selectNeuStr = strcat('SELECT ID AS NEURITE_ID FROM OMNISPHERO_NEURITE WHERE POSITIONY = ',num2str(positionY),' AND POSITIONX = ',num2str(positionX), ' AND EXPERIMENT_ID =', num2str(experimentID));
        ps=conn.prepareStatement(selectNeuStr);
        rs=ps.executeQuery();
        if(rs.next())
            curNeuID=str2num(char(rs.getString(1)));
            ps.close();
            rs.close();
        else
            ps.close();
            rs.close();
            insertNeuStr = strcat('INSERT INTO OMNISPHERO_NEURITE (POSITIONY,POSITIONX,TOTALLENGTH,TOTALBRANCHINGPOINTS,FPCOUNT,TNCOUNT,EXPERIMENT_ID,NUCLEUS_AREA_TOTAL, IMAGE_ID, REJECTED) VALUES (',num2str(positionY),',',num2str(positionX),',',num2str(currentValue(7)),',',num2str(currentValue(6)),',',num2str(currentValue(11)),',',num2str(currentValue(12)),',',num2str(experimentID),',', num2str(currentValue(13)),',',num2str(positionID),',',num2str(currentValue(14)),')');
            ps=conn.prepareStatement(insertNeuStr);
            ps.execute();
            ps.close();
            curNeuIDSTR = 'SELECT CURRVAL(''OMNISPHERO_NEURITE_ID_SEQ'')';
            ps=conn.prepareStatement(curNeuIDSTR);
            rs=ps.executeQuery();
            rs.next();
            curNeuID = str2num(char(rs.getString(1)));
            ps.close();
            rs.close();
        end
        
        %Insert Nuclei separately
        for(zz=1:nnz(currentValue(15,:)))
            nucPosInd = currentValue(15,zz);
            [nucPosY nucPosX] = ind2sub([7168 7168],nucPosInd);
            selectNucStr = strcat({'SELECT ID AS NUCLEUS_ID FROM OMNISPHERO_NUCLEUS WHERE POSITIONY BETWEEN '}, {num2str(nucPosY-1)}, {' AND '}, {num2str(nucPosY+1)}, {' AND POSITIONX BETWEEN '}, num2str(nucPosX-1),  {' AND '}, {num2str(nucPosX+1)}, {' AND EXPERIMENT_ID ='}, num2str(experimentID));
            ps=conn.prepareStatement(selectNucStr);
            rs=ps.executeQuery();
            if(rs.next())
                curNucID=str2num(char(rs.getString(1)));
                ps.close();
                rs.close();
            else            
                ps.close();
                rs.close();
                insertNucStr = strcat('INSERT INTO OMNISPHERO_NUCLEUS (POSITIONY,POSITIONX,EXPERIMENT_ID) VALUES (',num2str(nucPosY),',',num2str(nucPosX),',',num2str(experimentID),')');
                ps=conn.prepareStatement(insertNucStr);
                ps.execute();
                ps.close();
                curNucIDSTR = 'SELECT CURRVAL(''OMNISPHERO_NUCLEUS_ID_SEQ'')';
            end
        end
        insertSubNeuStr = strcat('INSERT INTO OMNISPHERO_SUBNEURITE (LENGTH, BRANCHINGPOINTS,NEURITEID) VALUES (',num2str(currentValue(2)),',',num2str(currentValue(3)),',',num2str(curNeuID),')');
        ps=conn.prepareStatement(insertSubNeuStr);
        ps.execute();
        ps.close();
        curSubNeuIDSTR = 'SELECT CURRVAL(''OMNISPHERO_SUBNEURITE_ID_SEQ'')';
        ps=conn.prepareStatement(curSubNeuIDSTR);
        rs=ps.executeQuery();
        rs.next();
        curSubNeuID = str2num(char(rs.getString(1)));
        ps.close();
        rs.close();
        %Nucleus To Do!
        nucleiIDsList = currentValue(4,4:lastNonZeroNumber);
        for(j=1:numel(nucleiIDsList))
            currentNucleusPosID = nucleiIDsList(j);
            [posY posX] = ind2sub([7168 7168], currentNucleusPosID);
            checkIfNucExist = strcat({'SELECT ID AS NUCLEUS_ID FROM OMNISPHERO_NUCLEUS WHERE POSITIONY BETWEEN '},{num2str(posY-1)},{' AND '},{num2str(posY+1)},{' AND POSITIONX BETWEEN '},{num2str(posX-1)},{' AND '},{num2str(posX+1)},{' AND EXPERIMENT_ID = '},num2str(experimentID));
            %checkIfNucExist = strcat('SELECT nuc.id as NUCLEUS_ID, neu.id AS NEURITE_ID, sub.id AS SUBNEURITE_ID FROM OMNISPHERO_NEURITE neu INNER JOIN OMNISPHERO_EXPERIMENT ex ON neu.experiment_id = ex.id INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.neuriteid = neu.id INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON map.subneuriteid = sub.id INNER JOIN OMNISPHERO_NUCLEUS nuc ON nuc.id = map.nucleusid WHERE ex.name =''', figtitle, ''' AND ex.well =''', selectedWell, {''' AND ex.type='}, {num2str(typeID)}, {' AND nuc.positiony BETWEEN '}, {num2str(posY-1)}, {' AND '}, {num2str(posY+1)}, {' AND nuc.positionx BETWEEN '}, {num2str(posX-1)}, {' AND '}, {num2str(posX+1)});
            ps=conn.prepareStatement(checkIfNucExist);
            rs=ps.executeQuery();
            if(~rs.next)
                ps.close();
                rs.close();
                %Nucleus doesn't exist. Add it.
                insertCommand = strcat('INSERT INTO OMNISPHERO_NUCLEUS (POSITIONX, POSITIONY, EXPERIMENT_ID) VALUES (',num2str(posX),',',num2str(posY),',',num2str(experimentID),')');
                ps=conn.prepareStatement(insertCommand);
                ps.execute();
                ps.close();
                %Get Nucleus ID
                curSubNeuIDSTR = 'SELECT CURRVAL(''OMNISPHERO_NUCLEUS_ID_SEQ'')';
                ps=conn.prepareStatement(curSubNeuIDSTR);
                rs=ps.executeQuery();
                rs.next();
                nucID = str2num(char(rs.getString(1)));
                ps.close();
                rs.close();
            else
                %Nucleus already exit. Just get id and add it to mapping.
                nucID = str2num(char(rs.getString(1)));
                ps.close();
                rs.close();
            end
            %Add mapping
            insertMap = strcat('INSERT INTO OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING (NUCLEUSID, SUBNEURITEID) VALUES (',num2str(nucID),',',num2str(curSubNeuID),')');
            ps=conn.prepareStatement(insertMap);
            ps.execute();
            ps.close();        
        end        
        exportText = [exportText well ';' neuriteID ';' subID ';' positionID ';' num2str(currentValue(2)) ';' num2str(currentValue(3)) ';' num2str(currentValue(4,4:lastNonZeroNumber)) ';' num2str(currentValue(5)) ';' num2str(currentValue(6)) ';' num2str(currentValue(7)) ';' num2str(currentValue(8)) ';' num2str(currentValue(10)) ';' num2str(currentValue(11)) ';' num2str(currentValue(12)) ';' num2str(currentValue(13)) ';' sprintf('\r\n')];
    end
    %[FileName,PathName] = uiputfile('NeuronStatistics.csv');
    %fileID = fopen(strcat(foldername,'/',selectedWell,exportname,'.csv'),'w');
    %fprintf(fileID,'%s',exportText);
    %fclose(fileID);
end

counter=0;
for(j=optionHandler.StartWell:numel(outNew))
    counter=counter+1;
    selectedWell=wellList(j);
    selectedWell=selectedWell{1};
    selectExpStr = strcat('SELECT ID FROM OMNISPHERO_EXPERIMENT WHERE NAME = ''',figtitle,''' AND WELL = ''',selectedWell,''' AND TYPE = ',num2str(typeID));
    ps=conn.prepareStatement(selectExpStr);
    rs=ps.executeQuery();
    rs.next();
    experimentID=str2num(char(rs.getString(1)));
    ps.close();
    rs.close();
    mapObj = outNew{counter};    
    keyset = mapObj.keys;
    
    
    for(i=1:numel(keyset))
        currentKey = keyset(i);
        currentKey = currentKey{1};
        currentValue = mapObj(currentKey);

        curSubNeuIDSTR = 'SELECT nextval(''OMNISPHERO_Neuron_ID_SEQ'')';
        ps=conn.prepareStatement(curSubNeuIDSTR);
        rs=ps.executeQuery();
        rs.next();
        neuID = str2num(char(rs.getString(1)));
        ps.close();
        rs.close();        
        insertNucStr = strcat('INSERT INTO OMNISPHERO_NEURON (experiment_FK,ID,DISTANCE_TO_SPHERECORE) VALUES (',num2str(experimentID),',',num2str(neuID),',',num2str(currentValue),')');
        ps=conn.prepareStatement(insertNucStr);
        ps.execute();
        ps.close();        
        
    end
end
conn.close();



% --------------------------------------------------------------------
function ExportNeuriteLengthStat_Cellomics_Callback(hObject, eventdata, handles)
% hObject    handle to ExportNeuriteLengthStat_Cellomics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
neuronHandler = handles.NeuronCoordinates;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
h=get(gca,'Title');
figtitle=get(h,'String');
% %1. Calculate
NucleusM = csvHandler.CellPosMatrix(selectedWell);
foldername = [imageHandler.Foldername '/ConvertedCellomics'];
[sizeY sizeX] = size(imageHandler.NeuriteImage);
NeuriteLengthMatrix = containers.Map();
figtitle = inputdlg('Please enter the experiment name!');
%Ask user for DB connection details
list=dir('*jdbc*');  
jdbcFolder = pwd;
if(numel(list)<=0)
        %errordlg('No Postgre JDBC driver found within working directory. Please download and save the Postgre JDBC driver to the working directory!' ,'Error');
        jdbcFolder = uigetdir('','Select JDBC driver folder');
        list=dir(strcat(jdbcFolder,'/','*jdbc*'));     
end
   % [FileName,PathName,FilterIndex] = uigetfile('*.jdbc*','Select Postgre JDBC Driver');
    FileName = list(1).name;
    PathName = jdbcFolder;
    unzip(strcat(PathName,'/',FileName),PathName);
    %javaclasspath(strcat(FileName(1:length(FileName)-4)))
    javaclasspath(jdbcFolder);
    %Ask user for password
    prompt = {'Postgre Hostname', 'Postgre Port', 'DB Name', 'Postgre Username', 'Password'};
    dlg_title='Postgre Connection';
    num_lines=1;
    defaultans={'localhost','5432','postgres','postgres','postgres'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    waitbarHandle=waitbar(0.1,'Calculating Neuronal Morphology. Please be patient.', 'WindowStyle', 'modal');
    props=java.util.Properties;
    props.setProperty('user',answer{4});
    props.setProperty('password', answer{5});
    driver=org.postgresql.Driver;
    url=strcat('jdbc:postgresql://',answer{1},':',answer{2},'/',answer{3});
    conn=driver.connect(url,props);
    %Check if necessary tables exist. If not, create them.   
    CreateNonExistingTable(conn);
    conn.close();
    %Check if empty
    
    if(numel(neuronHandler.ManualNeuronPositionsSparse) ~=0)
        ExportNeuriteLengthStat(neuronHandler.ManualNeuronPositionsSparse,'DataManual',figtitle,answer,handles);
    end
    waitbar(0.3)
    if(numel(neuronHandler.CellPosMatrix) ~=0)
        ExportNeuriteLengthStat(neuronHandler.CellPosMatrix,'DataCellomics',figtitle,answer,handles); 
    end
    waitbar(0.6);
    if(numel(neuronHandler.NeuronPositionsSkeletonization) ~=0)
        ExportNeuriteLengthStat(neuronHandler.NeuronPositionsSkeletonization,'DataOmnisphero',figtitle,answer,handles); 
    end
    close(waitbarHandle);
catch errorObj
    % If there is a problem, we display the error message
    %close(waitbarHandle);
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function CreateNonExistingTable(conn)
    sqlString = strcat('SELECT EXISTS ( SELECT 1 FROM information_schema.tables WHERE  table_name = ''omnisphero_datasource'')');
    sqlStatement = conn.prepareStatement(sqlString);
    sqlReader = sqlStatement.executeQuery();
    sqlReader.next();
    sqlResult = char(sqlReader.getString(1));
    if(sqlResult=='f')
        %Create Table
        
       insertSql = ['CREATE TABLE omnisphero_datasource(id integer NOT NULL, name character varying, CONSTRAINT omnisphero_datasource_pkey PRIMARY KEY (id)) WITH (OIDS=FALSE);'...
'CREATE TABLE omnisphero_experiment(id serial NOT NULL, name character varying, well character varying, type integer, CONSTRAINT experiment_pkey PRIMARY KEY (id))WITH (  OIDS=FALSE);'...
'CREATE TABLE omnisphero_neurite(id serial NOT NULL, positiony integer, positionx integer, totallength integer, totalbranchingpoints integer, fpcount integer, tncount integer, experiment_id integer, nucleus_area_total bigint, image_id integer, rejected integer, CONSTRAINT omnisphero_neurite_pkey PRIMARY KEY (id), CONSTRAINT omnisphero_neurite_experiment_id_fkey FOREIGN KEY (experiment_id) REFERENCES omnisphero_experiment (id) MATCH SIMPLE ON UPDATE NO ACTION ON DELETE NO ACTION) WITH (OIDS=FALSE);'...
'CREATE TABLE omnisphero_nucleus(id serial NOT NULL,positionx integer,positiony integer,area integer,experiment_id integer, CONSTRAINT omnisphero_nucleus_pkey PRIMARY KEY (id))WITH (  OIDS=FALSE);'...
'CREATE TABLE omnisphero_subneurite(id serial NOT NULL, length integer, branchingpoints integer, neuriteid integer, CONSTRAINT omnisphero_subneurite_pkey PRIMARY KEY (id), CONSTRAINT omnisphero_subneurite_neuriteid_fkey FOREIGN KEY (neuriteid) REFERENCES omnisphero_neurite (id) MATCH SIMPLE  ON UPDATE NO ACTION ON DELETE NO ACTION) WITH (OIDS=FALSE);'...
'CREATE TABLE omnisphero_subneurite_nucleus_mapping(nucleusid integer, subneuriteid integer, CONSTRAINT omnisphero_subneurite_nucleus_mapping_nucleusid_fkey FOREIGN KEY (nucleusid) REFERENCES omnisphero_nucleus (id) MATCH SIMPLE ON UPDATE NO ACTION ON DELETE NO ACTION, CONSTRAINT omnisphero_subneurite_nucleus_mapping_subneuriteid_fkey FOREIGN KEY (subneuriteid) REFERENCES omnisphero_subneurite (id) MATCH SIMPLE ON UPDATE NO ACTION ON DELETE NO ACTION) WITH (OIDS=FALSE);'
'CREATE TABLE omnisphero_neuron(experiment_FK integer, id serial NOT NULL, distance_to_spherecore integer, CONSTRAINT omnisphero_neuron_pkey PRIMARY KEY (id), CONSTRAINT omnisphero_neurite_experiment_FK_fkey FOREIGN KEY (experiment_FK) REFERENCES omnisphero_experiment (id) MATCH SIMPLE ON UPDATE NO ACTION ON DELETE NO ACTION) WITH (OIDS=FALSE);'];

        insertSql = strsplit(insertSql,';');
        for(i=1:numel(insertSql))
            sql=insertSql{i};
            ps=conn.prepareStatement(sql);
            ps.execute();           
            ps.close();
        end    
    end



% --------------------------------------------------------------------
function ExportNucleiSizes_Callback(hObject, eventdata, handles)
% hObject    handle to ExportNucleiSizes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
%Iterate over binary Nuclei Pictures of all Wells. Export Total Pixel Size of
%Neurospheres and Total Size of Nuclei
foldername = imageHandler.Foldername;
SphereAreaSizeX = optionHandler.MigrationDistanceDensityImageXSize;
SphereAreaSizeY = optionHandler.MigrationDistanceDensityImageYSize;
exportText = sprintf('Well;Migration Area Total;Nucleus Area Total;\r\n');
[sizeY sizeX] = size(imageHandler.NucleusImage);
wellList = get(handles.lbWell, 'string');
for i=optionHandler.StartWell:numel(wellList)
    selectedWell=wellList(i);
    selectedWell=selectedWell{1};
    selectedWellLong=selectedWell;AnalyzeM
     if(length(selectedWell) == 2)
         selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    end
    path = [foldername '/ConvertedCellomics/MarkerPointCoordinates-' selectedWell '.mat'];
    NucleusM = csvHandler.CellPosMatrix(selectedWell);
    ringNumber = optionHandler.DensityDistributionRingNumber;
    %Get density distribution to exclude too dense areas
    %if(saveCircle)    
    [filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, optionHandler, NucleusM, foldername, sizeY, sizeX); 
    
    %MarkerPointCoordinates are available. Get second ring:
    innerMask = ones(uint16(sizeY/(1/ScalingFactor)),uint16(sizeX/(1/ScalingFactor)));
    if(markerPointCoordinates ~=0 && numel(markerPointCoordinates('10') > 0))
        mk = markerPointCoordinates('10')./(1/ScalingFactor);
        innerMask = roipoly(innerMask,mk(:,1),mk(:,2));
        outerMask = ones(uint16(sizeY/(1/ScalingFactor)),uint16(sizeX/(1/ScalingFactor)));
        mk = markerPointCoordinates(num2str(ringNumber * 10))./(1/ScalingFactor);
        outerMask = roipoly(outerMask,mk(:,1),mk(:,2));        
        currentRing =logical(outerMask - innerMask);    
        currentRing = imresize(currentRing, [sizeY, sizeX]);
    else
        currentRing = logical(ones(sizeY,sizeX));
    end
    outerMask = logical(imresize(outerMask, [sizeY, sizeX]));
    innerMask = logical(imresize(innerMask, [sizeY, sizeX]));
    %Get density distribution to exclude too dense areas
    %if(saveCircle)    1 
   
    imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
    NucleusImage = imread(imagePathNucleusBig);
    if(optionHandler.CutOutSphereCore == 1)
        nucleusPic = CutOutCircles(NucleusImage,selectedWell,1,0,1,optionHandler, foldername, sizeY, sizeX, NucleusM, 0);
    else
        nucleusPic=NucleusImage;
    end
    nucleusPic(nucleusPic<optionHandler.NucleusThreshold) = 0;%double(double(imageHandler.NucleusImage) ./ 255);
    nucleusPic(nucleusPic>=optionHandler.NucleusThreshold) = 1;
    nucleusPic = logical(nucleusPic);
    nucleusPic = logical(xor(bwareaopen(nucleusPic,1),  bwareaopen(nucleusPic,15000)));
    nucleusPicBinBigNuclei = xor(bwareaopen(nucleusPic,250),  bwareaopen(nucleusPic,15000));
    D = bwdist(~nucleusPicBinBigNuclei);
    D = -D;
    D(~nucleusPic) = -Inf;
    L = watershed(D);
    nucleusPic(find(~logical(L))) = 0;
    D=0;
    L=0;
    nucleusPic = xor(bwareaopen(nucleusPic,35),  bwareaopen(nucleusPic,15000));
    totalSize=nnz(outerMask)-nnz(innerMask);
    nucleusWithoutOuter=logical(nucleusPic.*outerMask);
    totalSizeNuclei=nnz(nucleusWithoutOuter);
    exportText = [exportText selectedWell ';' num2str(totalSize) ';' num2str(totalSizeNuclei) ';' sprintf('\r\n')];
    disp(selectedWell);
end
[FileName,PathName] = uiputfile('NucleiSizeStatistic.csv');
fileID = fopen(strcat(PathName,'/',FileName),'w');
fprintf(fileID,'%s',exportText);
fclose(fileID);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --------------------------------------------------------------------
function ExportNucleiDensity_Callback(hObject, eventdata, handles)
% hObject    handle to ExportNucleiDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
%Calculate Density Image of Nuclei and export it. As image and as CSV
%Also calculate Derivation of it.
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
wellList = get(handles.lbWell, 'string');
exportText = sprintf('Well;Mean Gradient Nuclei Positions;Mean Gradient Nucleus Image;\r\n');
for i=optionHandler.StartWell:numel(wellList)
        selectedWell=wellList(i);
        selectedWell=selectedWell{1};
        selectedWellLong=selectedWell;
         if(length(selectedWell) == 2)
             selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
        end

    %Iterate over binary Nuclei Pictures of all Wells. Export Total Pixel Size of
    %Neurospheres and Total Size of Nuclei
    foldername = [imageHandler.Foldername '/ConvertedCellomics'];
    SphereAreaSizeX = optionHandler.MigrationDistanceDensityImageXSize;
    SphereAreaSizeY = optionHandler.MigrationDistanceDensityImageYSize;
    NucleusM = csvHandler.CellPosMatrix(selectedWell);
    [sizeY sizeX] = size(NucleusM);
    [DensityMSized] = getDensityDistributionsFromCSV(NucleusM,SphereAreaSizeX, SphereAreaSizeY);

    if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    elseif(length(selectedWell) == 4)
            selectedWellLong = [selectedWell(1) '0' selectedWell(2:length(selectedWell))];
    end
    %Read Nucleus Image
    imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
    NucleusImage = imread(imagePathNucleusBig);
    nucleusPic=NucleusImage;
    %nucleusPic = CutOutCircles(NucleusImage,selectedWell,1,0,1,optionHandler, foldername, sizeY, sizeX, NucleusM, NeuronM,0);
    nucleusPic(nucleusPic<optionHandler.NucleusThreshold) = 0;%double(double(imageHandler.NucleusImage) ./ 255);
    nucleusPic(nucleusPic>=optionHandler.NucleusThreshold) = 1;
    nucleusPic = logical(nucleusPic);
    nucleusPic = logical(xor(bwareaopen(nucleusPic,1),  bwareaopen(nucleusPic,15000)));
    nucleusPicBinBigNuclei = xor(bwareaopen(nucleusPic,250),  bwareaopen(nucleusPic,15000));
    D = bwdist(~nucleusPicBinBigNuclei);
    D = -D;
    D(~nucleusPic) = -Inf;
    L = watershed(D);
    nucleusPic(find(~logical(L))) = 0;
    nucleusPic = xor(bwareaopen(nucleusPic,35),  bwareaopen(nucleusPic,15000));
    [row col] = find(nucleusPic);
    DensityNucImage = uint8(hist3([row col],[SphereAreaSizeY SphereAreaSizeX]));
    gradientImage=imgradient(DensityMSized);
    gradientNucImage = imgradient(DensityNucImage);
    %Calculate mean gradient within Migration area on nonzero-pixels.
    path = [foldername '/MarkerPointCoordinates-' selectedWell '.mat'];
    NucleusM = csvHandler.CellPosMatrix(selectedWell);   
    ringNumber = optionHandler.DensityDistributionRingNumber;
    %Get density distribution to exclude too dense areas
    %if(saveCircle)        
    [filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, optionHandler, NucleusM, imageHandler.FolderName, sizeY, sizeX); 
    
   
    innerMask = ones(uint16(sizeY/(1/ScalingFactor)),uint16(sizeX/(1/ScalingFactor)));
    if(markerPointCoordinates ~=0 && numel(markerPointCoordinates('10') > 0))
        mk = markerPointCoordinates('10')./(1/ScalingFactor);
        innerMask = roipoly(innerMask,mk(:,1),mk(:,2));
        outerMask = ones(uint16(sizeY/(1/ScalingFactor)),uint16(sizeX/(1/ScalingFactor)));
        mk = markerPointCoordinates(num2str(ringNumber * 10))./(1/ScalingFactor);
        outerMask = roipoly(outerMask,mk(:,1),mk(:,2));        
        currentRing =logical(outerMask - innerMask);    
        currentRing = imresize(currentRing, [sizeY, sizeX]);
    else
        currentRing = logical(ones(sizeY,sizeX));
    end
    currentRingSmall = imresize(currentRing,[SphereAreaSizeY SphereAreaSizeX]);
    gradientImage = gradientImage.*currentRingSmall;
    gradientNucImage = gradientNucImage.*currentRingSmall;
    meanGrad = mean(gradientImage(find(gradientImage>0)));
    meanGradNucImage = mean(gradientNucImage(find(gradientNucImage>0)));
    disp(selectedWell);
    disp(strcat('Mean Gradient Positions: ',num2str(meanGrad)));
    disp(strcat('Mean Gradient Image: ',num2str(meanGradNucImage)));
    exportText = [exportText selectedWell ';' num2str(meanGrad) ';' num2str(meanGradNucImage) ';' sprintf('\r\n')];
end
[FileName,PathName] = uiputfile('StarDetectionResult.csv');
fileID = fopen(strcat(PathName,'/',FileName),'w');
fprintf(fileID,'%s',exportText);
fclose(fileID);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% --------------------------------------------------------------------
function MorphologyDBExport_Callback(hObject, eventdata, handles)
% hObject    handle to MorphologyDBExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try    
    list=dir('*jdbc*');  
    jdbcFolder = pwd;
    if(numel(list)<=0)
        %errordlg('No Postgre JDBC driver found within working directory. Please download and save the Postgre JDBC driver to the working directory!' ,'Error');
        jdbcFolder = uigetdir('','Select JDBC driver folder');
        list=dir(strcat(jdbcFolder,'/','*jdbc*')); 
    end    
   % [FileName,PathName,FilterIndex] = uigetfile('*.jdbc*','Select Postgre JDBC Driver');
    FileName = list(1).name;
    PathName = jdbcFolder;
    unzip(strcat(PathName,'/',FileName),PathName);
    %javaclasspath(strcat(FileName(1:length(FileName)-4)))
    javaclasspath(jdbcFolder);
        props=java.util.Properties;
        %Ask user for password
        prompt = {'Postgre Hostname', 'Postgre Port', 'DB Name', 'Postgre Username', 'Password'};
        dlg_title='Postgre Connection';
        num_lines=1;
        defaultans={'localhost','5432','postgres','postgres','postgres'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

        props.setProperty('user',answer{4});
        props.setProperty('password', answer{5});
        driver=org.postgresql.Driver;
        url=strcat('jdbc:postgresql://',answer{1},':',answer{2},'/',answer{3});
        conn=driver.connect(url,props);


        csvText=sprintf('Experiment ID;Experiment Name;Well;Type (1=Manual, 2=Cellomics, 3=Omnisphero);Total Subneuritelength per Nucleus;Average Subneuritelength per Nucleus;Mean Number Branchingpoints;Mean Subneurite Count\r\n');


        selectExps = 'SELECT ID, name, well, type FROM OMNISPHERO_EXPERIMENT';
        ps = conn.prepareStatement(selectExps);
        rs=ps.executeQuery();
        while(rs.next())
            expID = char(rs.getString(1));
            expName = char(rs.getString(2));
            currentWell = char(rs.getString(3));
            typeID = char(rs.getString(4));
            
            selectOne = ['select expname,typeid,wellname, avg(sublengthtotal) as Total_Subneuritelength_PerNucleus FROM(SELECT nuc.id,ex.name as expname,ex.well as wellname,ex.type as typeid,'...
' sum(sub.length) as sublengthtotal FROM OMNISPHERO_NUCLEUS nuc INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON map.nucleusid=nuc.id	INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid'...
	' INNER JOIN OMNISPHERO_NEURITE neu on neu.id = sub.neuriteid INNER JOIN OMNISPHERO_EXPERIMENT ex ON neu.experiment_id = ex.id'...
   ' WHERE ((neu.totallength > 10 AND neu.totallength <= 100 AND neu.totalbranchingpoints <3 ) OR (neu.totallength > 100 AND neu.totallength <= 200 AND neu.totalbranchingpoints <4 ) OR (neu.totallength > 200 AND neu.totallength <= 300 AND neu.totalbranchingpoints <4)OR (neu.totallength > 300 AND neu.totallength <= 360 AND neu.totalbranchingpoints <5))'...
	 'AND sub.id not in	(SELECT subid FROM (SELECT sub.id as subid, ex.name,ex.type,ex.well,count(distinct nuc.id) as distinctNuclei FROM OMNISPHERO_NUCLEUS nuc INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON'...
    ' map.nucleusid=nuc.id INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid	INNER JOIN OMNISPHERO_NEURITE neu ON sub.neuriteid=neu.id INNER JOIN OMNISPHERO_EXPERIMENT ex on neu.experiment_id = ex.id'...    
			' GROUP BY ex.name,ex.type, ex.well, sub.id order by sub.id asc) as innerResult WHERE innerResult.distinctNuclei >= 3) AND ex.name=''' expName ''' AND ex.well = ''' currentWell '''  AND ex.type=' typeID ' GROUP BY nuc.id, ex.name, ex.well,ex.type) as mainneurite group by expname, wellname, typeid order by expname,typeid, wellname'];
        
        selectTwo = ['select expname,typeid,wellname, avg(sublength) as Average_Subneuritelength_PerNucleus FROM(SELECT nuc.id,ex.name as expname,ex.well as wellname,ex.type as typeid, avg(sub.length) as sublength FROM OMNISPHERO_NUCLEUS nuc'...
	' INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON map.nucleusid=nuc.id INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid	INNER JOIN OMNISPHERO_NEURITE neu on neu.id = sub.neuriteid'...
	' INNER JOIN OMNISPHERO_EXPERIMENT ex ON neu.experiment_id = ex.id WHERE ((neu.totallength > 10 AND neu.totallength <= 100 AND neu.totalbranchingpoints <3 ) OR (neu.totallength > 100 AND neu.totallength <= 200 AND neu.totalbranchingpoints <4 ) OR (neu.totallength > 200 AND neu.totallength <= 300 AND neu.totalbranchingpoints <4)OR (neu.totallength > 300 AND neu.totallength <= 360 AND neu.totalbranchingpoints <5))'...
	' AND sub.id not in (SELECT subid FROM (SELECT sub.id as subid, ex.name,ex.type,ex.well,count(distinct nuc.id) as distinctNuclei FROM OMNISPHERO_NUCLEUS nuc'...
			' INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON map.nucleusid=nuc.id INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid	INNER JOIN OMNISPHERO_NEURITE neu ON sub.neuriteid=neu.id'...
			' INNER JOIN OMNISPHERO_EXPERIMENT ex on neu.experiment_id = ex.id GROUP BY ex.name,ex.type, ex.well, sub.id order by sub.id asc) as innerResult WHERE innerResult.distinctNuclei >= 3)'...
            ' AND ex.name=''' expName ''' AND ex.well = ''' currentWell '''  AND ex.type=' typeID...
	' GROUP BY nuc.id, ex.name, ex.well,ex.type) as mainneurite group by expname, wellname, typeid order by expname,typeid, wellname'];

    selectThree = ['select expname,typeid,wellname, avg(sublengthtotal) as Total_Subneuritelength_PerNucleus, avg(numberbranchingpoints) as numberbranchingpoints,avg(distinctsubneurites) as distinctsubneurites FROM'...
    ' (SELECT nuc.id,ex.name as expname,ex.well as wellname,ex.type as typeid,sum(sub.length) as sublengthtotal, sum(sub.branchingpoints) as numberbranchingpoints, count(distinct map.subneuriteid) as distinctsubneurites FROM OMNISPHERO_NUCLEUS nuc'...
	' INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON map.nucleusid=nuc.id INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid	INNER JOIN OMNISPHERO_NEURITE neu on neu.id = sub.neuriteid'...
	' INNER JOIN OMNISPHERO_EXPERIMENT ex ON neu.experiment_id = ex.id'...
	' WHERE ((neu.totallength > 10 AND neu.totallength <= 100 AND neu.totalbranchingpoints <3 ) OR (neu.totallength > 100 AND neu.totallength <= 200 AND neu.totalbranchingpoints <4 ) OR (neu.totallength > 200 AND neu.totallength <= 300 AND neu.totalbranchingpoints <4)OR (neu.totallength > 300 AND neu.totallength <= 360 AND neu.totalbranchingpoints <5))'...
    ' AND sub.id not in (SELECT subid FROM (SELECT sub.id as subid, ex.name,ex.type,ex.well,count(distinct nuc.id) as distinctNuclei FROM OMNISPHERO_NUCLEUS nuc	INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON map.nucleusid=nuc.id'...
	' INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid INNER JOIN OMNISPHERO_NEURITE neu ON sub.neuriteid=neu.id INNER JOIN OMNISPHERO_EXPERIMENT ex on neu.experiment_id = ex.id'...	
    ' GROUP BY ex.name,ex.type, ex.well, sub.id order by sub.id asc) as innerResult WHERE innerResult.distinctNuclei >= 3) AND ex.name=''' expName ''' AND ex.well = ''' currentWell '''  AND ex.type=' typeID ' GROUP BY nuc.id, ex.name, ex.well,ex.type) as mainneurite group by expname, wellname, typeid order by expname,typeid, wellname'];

            
%             selectAVG = ['select avg(totallength) as meanneuritelength, avg(totallength/mainneurite.subneurite_count::float) AS MEANSINGLENEURITELENGTH, avg(totalbranchingpoints) as meannumberbranchingpoints, avg(mainneurite.SUBNEURITE_COUNT) as meanSubneuriteCount, avg(mainneurite.subneurite_count/mainneurite.NUCLEUS_COUNT::float) AS Subneurites_Per_Nucleus from('...
%             ' SELECT count(distinct(map.nucleusid)) AS NUCLEUS_COUNT, count(distinct(map.subneuriteid)) AS SUBNEURITE_COUNT, neu.image_id as image_id, neu.id as neuriteid, neu.nucleus_area_total as nucleus_area_total, neu.tncount as tncount, neu.fpcount as fpcount, neu.totalbranchingpoints as totalbranchingpoints, neu.totallength as totallength, neu.positionx as positionx, neu.positiony as positiony '...
%              ' FROM OMNISPHERO_NEURITE neu'...
%               ' INNER JOIN OMNISPHERO_EXPERIMENT ex ON neu.experiment_id = ex.id'...
%               ' INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.neuriteid = neu.id'...
%               ' INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON map.subneuriteid = sub.id'...
%               ' WHERE ex.name=''' expName ''' AND ex.well = ''' currentWell '''  AND ex.type=' typeID...
%               ' GROUP BY neu.image_id, neu.nucleus_area_total, neu.tncount, neu.fpcount, neu.totalbranchingpoints, neu.totallength, neu.positionx, neu.positiony, neu.id'...
%         ') as mainneurite'];
%             avgStatement = conn.prepareStatement(selectAVG);
%             avgReader = avgStatement.executeQuery();
%             avgReader.next();
            
            oneStatement = conn.prepareStatement(selectOne); 
            oneReader = oneStatement.executeQuery();
            oneReader.next();
            meanNeuriteLength = char(oneReader.getString(4));
            oneStatement.close();
            oneReader.close();
            
            twoStatement = conn.prepareStatement(selectTwo);
            twoReader = twoStatement.executeQuery();
            twoReader.next();
            meanSingleNeuriteLength = char(twoReader.getString(4));
            twoStatement.close();
            twoReader.close();
            
            threeStatement = conn.prepareStatement(selectThree);
            threeReader = threeStatement.executeQuery();
            threeReader.next();
            meanNumberBranchingPoints = char(threeReader.getString(5));
            meanSubneuriteCount = char(threeReader.getString(6));
            threeStatement.close();
            threeReader.close();            
            csvText = [csvText expID ';' expName ';' currentWell ';' typeID ';' meanNeuriteLength ';' meanSingleNeuriteLength ';' meanNumberBranchingPoints ';' meanSubneuriteCount ';' sprintf('\r\n')];
        end
        ps.close();
        rs.close();
        conn.close();
        [FileName,PathName] = uiputfile('NeuronalMorphologyTotal.csv');
        fileID = fopen(strcat(PathName,'/',FileName),'w');
        fprintf(fileID,'%s',csvText);
        fclose(fileID);    
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end



% --------------------------------------------------------------------
function MulPa_Small_Callback(hObject, eventdata, handles)
% hObject    handle to MulPa_Small (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function QuantifyNuclei_Callback(hObject, eventdata, handles)
% hObject    handle to QuantifyNuclei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
csvHandler = handles.CSVCoordinates;
%Print warning if nuclei are already quantified
wellList = get(handles.lbWell, 'string');
selectedWellNumber = get(handles.lbWell,'Value');
selectedWell = wellList{selectedWellNumber};
WellNucleusDict = csvHandler.CellPosMatrix;
buttonResult = 'Yes';
if(isKey(WellNucleusDict, selectedWell))
    NucleusM = WellNucleusDict(selectedWell);
    if(nnz(NucleusM) > 0)
        %Warning
        buttonResult = questdlg('Nuclei are already quantified. Are you sure you want to overwrite the current Nuclei quantification?');
    end
end
if(strcmp(buttonResult,'Yes'))
    %Quantify nuclei
    watershedNucleusImage(0,selectedWell,handles);
end


% --------------------------------------------------------------------
function CompFillAllWells_Callback(hObject, eventdata, handles)
% hObject    handle to CompFillAllWells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
%foldername = uigetdir;
%Call EdgeFill Neurons and Skeleton Neurons and play with parameters in
%Option Handler
ExecuteOptionExperiment(handles,1,'Neurite');
guidata(handles.figure1, handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --------------------------------------------------------------------
function PreprocessWatershedNuclei_Callback(hObject, eventdata, handles)
% hObject    handle to PreprocessWatershedNuclei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
WatershedNuclei(handles);


function WatershedNuclei(handles)
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
wellList = get(handles.lbWell, 'string');
waitbarHandle = waitbar(0,'Watersheding Nuclei. Please wait.', 'WindowStyle', 'modal');
foldername = [imageHandler.Foldername '/ConvertedCellomics'];
for i=optionHandler.StartWell:numel(wellList)
    waitbar((i/numel(wellList)),waitbarHandle);
        selectedWell=wellList(i);
        selectedWell=selectedWell{1};
        selectedWellLong=selectedWell;
         if(length(selectedWell) == 2)
             selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
        end
    imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
    nucleusPic = imread(imagePathNucleusBig);
    nucleusPic(nucleusPic<optionHandler.NucleusThreshold) = 0;
    nucleusPic(nucleusPic>=optionHandler.NucleusThreshold) = 255;
    nucleusPic=imcomplement(nucleusPic);
    nucleusPic = uint8(nucleusPic);

    %Do Watershed
    %Check if there is already a watershed image. If not, create it.
    imagePathWatershed = [foldername '/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
    javaaddpath('Omnisphero src/lib/ij.jar');
    javaaddpath('Omnisphero src/lib/mij.jar');
    javaaddpath('lib/ij.jar');
    javaaddpath('lib/mij.jar');
    MIJ.start(false);
    MIJ.createImage(nucleusPic);
    MIJ.run('Options...', 'iterations=1 count=1 edm=Overwrite do=Nothing');
    MIJ.run('Watershed');
    %nucleusPicBin = logical(MIJ.getCurrentImage);
    MIJ.run('Save', 'save=temp.tif');
    MIJ.closeAllWindows;
    MIJ.exit;
    nucleusPic = imread('temp.tif');
    delete temp.tif;
    nucleusPic = logical(1-nucleusPic);
    imwrite(nucleusPic,imagePathWatershed);
end
close(waitbarHandle);

% --------------------------------------------------------------------
function CalculateHighThroughput_Callback(hObject, eventdata, handles)
% hObject    handle to CalculateHighThroughput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
set(handles.cbNucleusPicture ,'Value',1);
set(handles.cbNeuritePicture ,'Value',1);

needSphereCutter = questdlg('Do you have an experiment with multiple spheres in one OT?');
additionalOligo = questdlg('Do you have Oligodendrocytes in your experiment?');
additionalAstro = questdlg('Do you have Astrocytes in your experiment?');
additionalOligo = strcmp(additionalOligo,'Yes');
additionalAstro = strcmp(additionalAstro,'Yes');
%2. Select folder with Save files
%3. Iterate over all save files
htFolder = uigetdir('','Select Batch folder');
htFolder = strcat(htFolder,'/');
list=dir(strcat(htFolder,'*.zip'));

imFolder = uigetdir('','Select image parent folder');
imFolder = strcat(imFolder,'/');

%Get Config xml file
configfilename=strcat(htFolder,'Config.xml');
optionFields = zeros(0);
if(exist(configfilename,'file'))
    struct = xml2struct(configfilename);
    %Get child objects
    configStruct = struct.Configuration;
    optionFields = fieldnames(configStruct);
end

%Create subdirectory for Save files
mkdir(strcat(htFolder,'Results'));
resultFolder = strcat(htFolder,'Results');

for(i=1:numel(list))
    selectedExp = list(i).name(1:length(list(i).name)-8);    
    %4. Load current file
    Load(htFolder,list(i).name,handles);  
    %Get handlers after loading
    selectedExp = list(i).name(1:length(list(i).name)-8);  
    LoadImageFolder(strcat(imFolder,selectedExp),WellType.Undefined,0,0,0,0,0,0,handles); 
    handles = guidata(handles.figure1);
    imageHandler = handles.ImageHandler;
    optionHandler = handles.OptionHandler;
    csvHandler = handles.CSVCoordinates;
    if(strcmp(needSphereCutter,'Yes'))
        %Do Watershed
        %Check if there is already a watershed image. If not, create it.
        %Iterate over all wells
        wellList = get(handles.lbWell, 'string');
        %waitbarHandle = waitbar(0,'Watersheding Nuclei. Please wait.', 'WindowStyle', 'modal');
        foldername = [imageHandler.Foldername '/ConvertedCellomics'];
        for i=optionHandler.StartWell:numel(wellList)
            %waitbar((i/numel(wellList)),waitbarHandle);
            selectedWell=wellList(i);
            selectedWell=selectedWell{1};
            selectedWellLong=selectedWell;
            if(length(selectedWell) == 2)
                 selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
            end
            imagePathWatershed = [foldername '/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
            if(exist(imagePathWatershed,'file'))
                nucleusPic = imread(imagePathWatershed);
            else
                imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
                nucleusPic = imread(imagePathNucleusBig);
                nucleusPic(nucleusPic<optionHandler.NucleusThreshold) = 0;
                nucleusPic(nucleusPic>=optionHandler.NucleusThreshold) = 255;
                nucleusPic=imcomplement(nucleusPic);
                nucleusPic = uint8(nucleusPic);
                javaaddpath('lib/ij.jar');
                javaaddpath('lib/mij.jar');
                MIJ.start(false);
                MIJ.createImage(nucleusPic);
                MIJ.run('Options...', 'iterations=1 count=1 edm=Overwrite do=Nothing');
                MIJ.run('Watershed');
                %nucleusPicBin = logical(MIJ.getCurrentImage);
                MIJ.run('Save', 'save=temp.tif');
                MIJ.closeAllWindows;
                MIJ.exit;
                nucleusPic = imread('temp.tif');
                delete temp.tif;
                nucleusPic = logical(1-nucleusPic);
                imwrite(nucleusPic,imagePathWatershed);
            end
        end
        [newCSVHandler wellList] = Imagecutter_V5_9(1, additionalOligo, additionalAstro, imageHandler.Foldername, csvHandler, optionHandler);
        PostProcessImageCutter(wellList,newCSVHandler,handles);
        
        %handles.CSVCoordinates = csvHandler;
        %guidata(handles.figure1, handles);
    end
    QuickFixManualPositions(handles);
    %Set Image folder
    %LoadImageFolder(strcat(imFolder,selectedExp),WellType.Undefined,0,0,0,0,0,0,handles);    
    %Load first well that imageHandler is set    
    %5. Set changed parameters    
    handles = guidata(handles.figure1);
    optionHandler = handles.OptionHandler;
    for(j=1:numel(optionFields))
        fieldName = optionFields{j};
        fieldValue = configStruct.(optionFields{j}).Text;
        optionHandler.SetStringValue(fieldName,fieldValue);
    end
    %6. Run "Neuronal Quantification on all Wells"
    handles.OptionHandler = optionHandler;
    guidata(handles.figure1, handles);
    handles = guidata(handles.figure1);
    
    
%Call EdgeFill Neurons and Skeleton Neurons and play with parameters in
%Option Handler
    ChangeWell(1,handles);
    ExecuteOptionExperiment(handles,0,'Neurite');
    if(strcmp(additionalOligo,'Yes'))
        handles = guidata(handles.figure1);
        ExecuteOptionExperiment(handles,0,'Oligo');
    end
    handles = guidata(handles.figure1);
    %Save file automatically in SubFolder
    SaveFixedNeuronPositions(list(i).name(1:length(list(i).name)-4),strcat(resultFolder,'/'),handles);    
end
catch errorObj
    % If there is a problem, we display the error message
    disp(getReport(errorObj,'extended','hyperlinks','off'));
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


function PostProcessImageCutter(wellList,newCSVHandler,handles)
    csvHandler = handles.CSVCoordinates;
        
    %Delete and rewrite maps
    newCellPosMatrix = containers.Map();
    csvHandler.ManualNeuronPositionsSparse = containers.Map();
    for(i=1:numel(wellList))
        selectedWell =wellList{i};
        newCellPosMatrix(selectedWell) = newCSVHandler.CellPosMatrix(selectedWell);
        csvHandler.ManualNeuronPositionsSparse(selectedWell) = sparse(logical(zeros(7168,7168)));
    end
    csvHandler.CellPosMatrix = newCellPosMatrix;
    handles.CSVCoordinates = csvHandler;
    guidata(handles.figure1, handles);
    ReloadWellList(handles,0);
    
    

function QuickFixManualM(handles)
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
[sizeY sizeX] = size(imageHandler.NucleusImage);
wellList = get(handles.lbWell, 'string');

for currentWellIndex=optionHandler.StartWell:numel(wellList)
    selectedWell = wellList(currentWellIndex);
    selectedWell = selectedWell{1};
    if(isKey(neuronHandler.ManualNeuronPositionsSparse, selectedWell))
        ManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
        if(size(ManualM,2)>sizeX)
             ManualM=ManualM(:,1:sizeX);
             neuronHandler.ManualNeuronPositionsSparse(selectedWell)=ManualM;
        end

        if(size(ManualM,1)>sizeY)
             ManualM=ManualM(1:sizeY,:);
             neuronHandler.ManualNeuronPositionsSparse(selectedWell)=ManualM;
        end
    end
end
handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles);


function QuickFixManualPositions(handles)
%Sometimes the ManualPositionsMatrix is accidently not a sparse matrix
%This will lead to Saving problems
%This function changes every ManualPositionMatrix into a sparse matrix
handles = guidata(handles.figure1);
neuronHandler = handles.NeuronCoordinates;
if(numel(neuronHandler.ManualNeuronPositionsSparse) ==0) 
    neuronHandler.ManualNeuronPositionsSparse = containers.Map()
end
keyList = keys(neuronHandler.ManualNeuronPositionsSparse);
for(i=1:numel(keyList))
    currentKey=keyList{i};
    neuronHandler.ManualNeuronPositionsSparse(currentKey) = sparse(neuronHandler.ManualNeuronPositionsSparse(currentKey));
end
handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles);

function TotalExport_Callback(hObject, eventdata, handles)
try
    neuronHandler = handles.NeuronCoordinates;
    figtitle = inputdlg('Please enter the experiment name!');
    list=dir('*jdbc*');  
    jdbcFolder = pwd;
    if(numel(list)<=0)
        %errordlg('No Postgre JDBC driver found within working directory. Please download and save the Postgre JDBC driver to the working directory!' ,'Error');
        jdbcFolder = uigetdir('','Select JDBC driver folder');
        list=dir(strcat(jdbcFolder,'/','*jdbc*')); 
    end    
    % [FileName,PathName,FilterIndex] = uigetfile('*.jdbc*','Select Postgre JDBC Driver');
    FileName = list(1).name;
    PathName = jdbcFolder;
    percNeurons = inputdlg('Enter percentage of Neurons for Morphological Analysis','Percentage Neurons',1,{'100'});
    percNeurons = percNeurons{1};
    percNeurons = str2num(percNeurons);
    percNeurons = percNeurons/100;
    [FileName,PathName] = uiputfile('migdist.csv');
    %unzip(strcat(PathName,'/',FileName),PathName);
    %javaclasspath(strcat(FileName(1:length(FileName)-4)))
    javaclasspath(jdbcFolder);
    props=java.util.Properties;
    %Ask user for password
    prompt = {'Postgre Hostname', 'Postgre Port', 'DB Name', 'Postgre Username', 'Password'};
    dlg_title='Postgre Connection';
    num_lines=1;
    defaultans={'localhost','5432','postgres','postgres','postgres'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

    props.setProperty('user',answer{4});
    props.setProperty('password', answer{5});
    driver=org.postgresql.Driver;
    url=strcat('jdbc:postgresql://',answer{1},':',answer{2},'/',answer{3});
    
    
    ExportNeuriteLengthStat(neuronHandler.NeuronPositionsSkeletonization,'DataOmnisphero',figtitle,answer,handles);
    conn=driver.connect(url,props);
    TotalExport(FileName,PathName,percNeurons,conn,handles);
    conn.close();
catch errorObj
    % If there is a problem, we display the error message    
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
    conn.close();
end



function TotalExport(FileName,PathName,PercentageNeurons,conn,handles)
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
if(isfield(handles,'OligoCoordinates') && isprop(handles.OligoCoordinates,'NeuronPositionsSkeletonization'))
    oligoHandler = handles.OligoCoordinates;
end
expName=get(gcf,'name');
expName=expName(1:length(expName)-4);
fileID = fopen(strcat(PathName,'/',FileName),'w');
foldername = [imageHandler.Foldername '/ConvertedCellomics'];
handles = guidata(handles.figure1);
ringCount = optionHandler.DensityDistributionRingNumber;
densityWidth = optionHandler.MigrationDistanceDensityImageXSize;
densityHeight = optionHandler.MigrationDistanceDensityImageYSize
wellList = get(handles.lbWell, 'string');
waitbarHandle = waitbar(0,'Calculating and exporting everything for all wells', 'WindowStyle', 'modal');
selectedWellNumber = get(handles.lbWell,'Value');
exportText = sprintf('Well;Ring Number;Migration Distance(px);Migration Radius (px);Manual measured Migration Distance (px);Migration Distance without filter (px);Number of Pixels;Number Nuclei;Cellomics Neurons;Manual Neurons;Skeleton Neurons;Cellomics Neurons per Nuclei;Manual Neurons per Nuclei;Skeleton Neurons per Nuclei;Nuclei per Pixels;;Composite Fill Oligos;Composite Fill Oligos per Nuclei;Skeleton Oligos;Skeleton Oligos per Nuclei;;Mean Gradient Nuclei Positions;Mean Gradient Nucleus Image;;Total Subneuritelength per Nucleus unlimited;Average Subneuritelength per Nucleus unlimited;Mean Number Branchingpoints unlimited;Mean Subneurite Count unlimited;;Total Subneuritelength per Nucleus limited;Average Subneuritelength per Nucleus limited;Mean Number Branchingpoints limited;Mean Subneurite Count limited;;Neurons Omnisphero & Manual;Additional Manual;Additional Omnisphero;;Oligos CompFill Omnisphero & Manual;Additional Manual;Additional CompFill Oligo;;Oligos Omnisphero & Manual;Additional Manual;Additional Omnisphero\r\n');

%1. Create Neuron distribution Sheet

ExportDistanceDistribution2(PathName,expName,neuronHandler.NeuronPositionsSkeletonization,'Neuro',handles);
if(isfield(handles,'OligoCoordinates') && isprop(handles.OligoCoordinates,'NeuronPositionsSkeletonization'))
    ExportDistanceDistribution2(PathName,expName,oligoHandler.NeuronPositionsSkeletonization,'OligoSkel',handles);
    everyThingAvailable = ExportDistanceDistribution2(PathName,expName,oligoHandler.NeuronPositionsEdgeFill,'OligoComp',handles);
end

[sizeY sizeX] = size(imageHandler.NucleusImage);
WellNucleusDict = csvHandler.CellPosMatrix;

%if(0)
for currentWellIndex=optionHandler.StartWell:numel(wellList)
    selectedWell = wellList(currentWellIndex);
    selectedWell = selectedWell{1};
    disp(selectedWell);
    NeuronM = neuronHandler.NeuronPositionsSkeletonization(selectedWell);
    NucleusM = WellNucleusDict(selectedWell);
    %Get percentage of Neurons
    [filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, densityWidth, densityHeight, optionHandler, NucleusM, imageHandler.Foldername, sizeY, sizeX); 
    currentDistancesList = zeros(0);
    innerMask = ones(uint16(sizeY/(1/ScalingFactor)),uint16(sizeX/(1/ScalingFactor)));       
    mk = markerPointCoordinates('10')./(1/ScalingFactor);
    innerMask = roipoly(innerMask,mk(:,1),mk(:,2));
    currentRing =logical(ones(int32(sizeY/(1/ScalingFactor)),int32(sizeX/(1/ScalingFactor))) - innerMask);
    currentRing =imresize(currentRing,[sizeY,sizeX]);
    NeuronM = NeuronM .* currentRing;
    numberNeurons = nnz(NeuronM);
    limitFactor = PercentageNeurons * numberNeurons;
    limitFactor = round(limitFactor);
    %1. Migration Distance and Density Distribution
    [distFilter, distWithoutFilter, sphereRadius, SphereArea, markerPointCoordinates densDist densDistCompFill neuronHandler]=calculateDensityDistribution(selectedWell,densityWidth,densityHeight,neuronHandler,handles);
    if(numel(densDist) == 1)
        %Fill up densDist with dummy values
        densDist = zeros(ringCount,9);         
    end
    
    %Oligos
    if(isfield(handles,'OligoCoordinates') && isprop(handles.OligoCoordinates,'NeuronPositionsSkeletonization') && everyThingAvailable)
        oligoHandler = handles.OligoCoordinates;
        [distFilterOligos, distWithoutFilterOligos, sphereRadiusOligos, SphereAreaOligos, markerPointCoordinatesOligos densDistOligos densDistCompFillOligos oligoHandler]=calculateDensityDistribution(selectedWell,densityWidth,densityHeight,oligoHandler,handles);
        if(numel(densDistOligos) == 1)
            %Fill up densDist with dummy values
            densDistOligos = zeros(ringCount,9);        
        end
        if(numel(densDistCompFillOligos) == 1)
            %Fill up densDist with dummy values
            densDistCompFillOligos = zeros(ringCount,9);        
        end
    end
    
    if(distFilter >-1)
        
        
        
        
        
        
    %2. Add gradient to output
    SphereAreaSizeX = optionHandler.MigrationDistanceDensityImageXSize;
    SphereAreaSizeY = optionHandler.MigrationDistanceDensityImageYSize;
    NucleusM = csvHandler.CellPosMatrix(selectedWell);
    [sizeY sizeX] = size(NucleusM);
    [DensityMSized] = getDensityDistributionsFromCSV(NucleusM,SphereAreaSizeX, SphereAreaSizeY);
    selectedWellLong = selectedWell;
    if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    elseif(length(selectedWell) == 4)
            selectedWellLong = [selectedWell(1) '0' selectedWell(2:length(selectedWell))];
    end
    %Read Nucleus Image
    imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
    NucleusImage = imread(imagePathNucleusBig);
    nucleusPic=NucleusImage;
    %nucleusPic = CutOutCircles(NucleusImage,selectedWell,1,0,1,optionHandler, foldername, sizeY, sizeX, NucleusM, NeuronM,0);
    nucleusPic(nucleusPic<optionHandler.NucleusThreshold) = 0;%double(double(imageHandler.NucleusImage) ./ 255);
    nucleusPic(nucleusPic>=optionHandler.NucleusThreshold) = 1;
    nucleusPic = logical(nucleusPic);
    nucleusPic = logical(xor(bwareaopen(nucleusPic,1),  bwareaopen(nucleusPic,15000)));
    nucleusPicBinBigNuclei = xor(bwareaopen(nucleusPic,250),  bwareaopen(nucleusPic,15000));
    D = bwdist(~nucleusPicBinBigNuclei);
    D = -D;
    D(~nucleusPic) = -Inf;
    L = watershed(D);
    nucleusPic(find(~logical(L))) = 0;
    nucleusPic = xor(bwareaopen(nucleusPic,35),  bwareaopen(nucleusPic,15000));
    [row col] = find(nucleusPic);
    DensityNucImage = uint8(hist3([row col],[SphereAreaSizeY SphereAreaSizeX]));
    gradientImage=imgradient(DensityMSized);
    gradientNucImage = imgradient(DensityNucImage);
    %Calculate mean gradient within Migration area on nonzero-pixels.
    
    
    
    
    
    NucleusM = csvHandler.CellPosMatrix(selectedWell);   
    ringNumber = optionHandler.DensityDistributionRingNumber;
    [filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, optionHandler, NucleusM, imageHandler.Foldername, sizeY, sizeX)
    %Get density distribution to exclude too dense areas
    %if(saveCircle)   
   
    innerMask = ones(uint16(sizeY/(1/ScalingFactor)),uint16(sizeX/(1/ScalingFactor)));
    if(markerPointCoordinates ~=0 && numel(markerPointCoordinates('10') > 0))
        mk = markerPointCoordinates('10')./(1/ScalingFactor);
        innerMask = roipoly(innerMask,mk(:,1),mk(:,2));
        outerMask = ones(uint16(sizeY/(1/ScalingFactor)),uint16(sizeX/(1/ScalingFactor)));
        mk = markerPointCoordinates(num2str(ringNumber * 10))./(1/ScalingFactor);
        outerMask = roipoly(outerMask,mk(:,1),mk(:,2));        
        currentRing =logical(outerMask - innerMask);    
        currentRing = imresize(currentRing, [sizeY, sizeX]);
    else
        currentRing = logical(ones(sizeY,sizeX));
    end
    currentRingSmall = imresize(currentRing,[SphereAreaSizeY SphereAreaSizeX]);
    gradientImage = gradientImage.*currentRingSmall;
    gradientNucImage = gradientNucImage.*currentRingSmall;
    meanGrad = mean(gradientImage(find(gradientImage>0)));
    meanGradNucImage = mean(gradientNucImage(find(gradientNucImage>0)));
    disp(selectedWell);
    disp(strcat('Mean Gradient Positions: ',num2str(meanGrad)));
    disp(strcat('Mean Gradient Image: ',num2str(meanGradNucImage)));
        
    
    %3. Extract morphology from DB
    selectExps = ['SELECT ID, name, well, type FROM OMNISPHERO_EXPERIMENT ex WHERE ex.well = ''' selectedWell ''' AND ex.name=''' expName ''' AND ex.type=3'];
    ps = conn.prepareStatement(selectExps);
    rs=ps.executeQuery();
    meanSubneuriteCount=0;
    meanNumberBranchingPoints=0;
    meanSingleNeuriteLength=0;
    meanNeuriteLength=0;
    %limitFactor = 99999;
    while(rs.next())
        expID = char(rs.getString(1));
        %expName = char(rs.getString(2));
        currentWell = char(rs.getString(3));
        typeID = char(rs.getString(4));

        selectOne = ['select expname,typeid,wellname, avg(sublengthtotal) as Total_Subneuritelength_PerNucleus FROM(SELECT nuc.id,ex.name as expname,ex.well as wellname,ex.type as typeid,'...
' sum(sub.length) as sublengthtotal FROM OMNISPHERO_NUCLEUS nuc INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON map.nucleusid=nuc.id	INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid'...
' INNER JOIN OMNISPHERO_NEURITE neu on neu.id = sub.neuriteid INNER JOIN OMNISPHERO_EXPERIMENT ex ON neu.experiment_id = ex.id'...
' WHERE ((neu.totallength > 10 AND neu.totallength <= 100 AND neu.totalbranchingpoints <3 ) OR (neu.totallength > 100 AND neu.totallength <= 200 AND neu.totalbranchingpoints <4 ) OR (neu.totallength > 200 AND neu.totallength <= 300 AND neu.totalbranchingpoints <4)OR (neu.totallength > 300 AND neu.totallength <= 360 AND neu.totalbranchingpoints <5))'...
 'AND sub.id not in	(SELECT subid FROM (SELECT sub.id as subid, ex.name,ex.type,ex.well,count(distinct nuc.id) as distinctNuclei FROM OMNISPHERO_NUCLEUS nuc INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON'...
' map.nucleusid=nuc.id INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid	INNER JOIN OMNISPHERO_NEURITE neu ON sub.neuriteid=neu.id INNER JOIN OMNISPHERO_EXPERIMENT ex on neu.experiment_id = ex.id'...    
        ' GROUP BY ex.name,ex.type, ex.well, sub.id order by sub.id asc) as innerResult WHERE innerResult.distinctNuclei >= 3) AND ex.name=''' expName ''' AND ex.well = ''' currentWell '''  AND ex.type=' typeID ' GROUP BY nuc.id, ex.name, ex.well,ex.type ORDER BY sublengthtotal DESC LIMIT ' num2str(limitFactor) ' ) as mainneurite group by expname, wellname, typeid order by expname,typeid, wellname'];

    selectTwo = ['select expname,typeid,wellname, avg(sublength) as Average_Subneuritelength_PerNucleus FROM(SELECT nuc.id,ex.name as expname,ex.well as wellname,ex.type as typeid, avg(sub.length) as sublength, sum(sub.length) as sublengthtotal FROM OMNISPHERO_NUCLEUS nuc'...
' INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON map.nucleusid=nuc.id INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid	INNER JOIN OMNISPHERO_NEURITE neu on neu.id = sub.neuriteid'...
' INNER JOIN OMNISPHERO_EXPERIMENT ex ON neu.experiment_id = ex.id WHERE ((neu.totallength > 10 AND neu.totallength <= 100 AND neu.totalbranchingpoints <3 ) OR (neu.totallength > 100 AND neu.totallength <= 200 AND neu.totalbranchingpoints <4 ) OR (neu.totallength > 200 AND neu.totallength <= 300 AND neu.totalbranchingpoints <4)OR (neu.totallength > 300 AND neu.totallength <= 360 AND neu.totalbranchingpoints <5))'...
' AND sub.id not in (SELECT subid FROM (SELECT sub.id as subid, ex.name,ex.type,ex.well,count(distinct nuc.id) as distinctNuclei FROM OMNISPHERO_NUCLEUS nuc'...
        ' INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON map.nucleusid=nuc.id INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid	INNER JOIN OMNISPHERO_NEURITE neu ON sub.neuriteid=neu.id'...
        ' INNER JOIN OMNISPHERO_EXPERIMENT ex on neu.experiment_id = ex.id GROUP BY ex.name,ex.type, ex.well, sub.id order by sub.id asc) as innerResult WHERE innerResult.distinctNuclei >= 3)'...
        ' AND ex.name=''' expName ''' AND ex.well = ''' currentWell '''  AND ex.type=' typeID...
' GROUP BY nuc.id, ex.name, ex.well,ex.type ORDER BY sublengthtotal DESC LIMIT ' num2str(limitFactor) ' ) as mainneurite group by expname, wellname, typeid order by expname,typeid, wellname'];

selectThree = ['select expname,typeid,wellname, avg(sublengthtotal) as Total_Subneuritelength_PerNucleus, avg(numberbranchingpoints) as numberbranchingpoints,avg(distinctsubneurites) as distinctsubneurites FROM'...
' (SELECT nuc.id,ex.name as expname,ex.well as wellname,ex.type as typeid,sum(sub.length) as sublengthtotal, sum(sub.branchingpoints) as numberbranchingpoints, count(distinct map.subneuriteid) as distinctsubneurites FROM OMNISPHERO_NUCLEUS nuc'...
' INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON map.nucleusid=nuc.id INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid	INNER JOIN OMNISPHERO_NEURITE neu on neu.id = sub.neuriteid'...
' INNER JOIN OMNISPHERO_EXPERIMENT ex ON neu.experiment_id = ex.id'...
' WHERE ((neu.totallength > 10 AND neu.totallength <= 100 AND neu.totalbranchingpoints <3 ) OR (neu.totallength > 100 AND neu.totallength <= 200 AND neu.totalbranchingpoints <4 ) OR (neu.totallength > 200 AND neu.totallength <= 300 AND neu.totalbranchingpoints <4)OR (neu.totallength > 300 AND neu.totallength <= 360 AND neu.totalbranchingpoints <5))'...
' AND sub.id not in (SELECT subid FROM (SELECT sub.id as subid, ex.name,ex.type,ex.well,count(distinct nuc.id) as distinctNuclei FROM OMNISPHERO_NUCLEUS nuc	INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON map.nucleusid=nuc.id'...
' INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid INNER JOIN OMNISPHERO_NEURITE neu ON sub.neuriteid=neu.id INNER JOIN OMNISPHERO_EXPERIMENT ex on neu.experiment_id = ex.id'...	
' GROUP BY ex.name,ex.type, ex.well, sub.id order by sub.id asc) as innerResult WHERE innerResult.distinctNuclei >= 3) AND ex.name=''' expName ''' AND ex.well = ''' currentWell '''  AND ex.type=' typeID ' GROUP BY nuc.id, ex.name, ex.well,ex.type ORDER BY sublengthtotal DESC LIMIT ' num2str(limitFactor) ' ) as mainneurite group by expname, wellname, typeid order by expname,typeid, wellname'];

 selectOneUnlimited = ['select expname,typeid,wellname, avg(sublengthtotal) as Total_Subneuritelength_PerNucleus FROM(SELECT nuc.id,ex.name as expname,ex.well as wellname,ex.type as typeid,'...
' sum(sub.length) as sublengthtotal FROM OMNISPHERO_NUCLEUS nuc INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON map.nucleusid=nuc.id	INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid'...
' INNER JOIN OMNISPHERO_NEURITE neu on neu.id = sub.neuriteid INNER JOIN OMNISPHERO_EXPERIMENT ex ON neu.experiment_id = ex.id'...
' WHERE ((neu.totallength > 10 AND neu.totallength <= 100 AND neu.totalbranchingpoints <3 ) OR (neu.totallength > 100 AND neu.totallength <= 200 AND neu.totalbranchingpoints <4 ) OR (neu.totallength > 200 AND neu.totallength <= 300 AND neu.totalbranchingpoints <4)OR (neu.totallength > 300 AND neu.totallength <= 360 AND neu.totalbranchingpoints <5))'...
 'AND sub.id not in	(SELECT subid FROM (SELECT sub.id as subid, ex.name,ex.type,ex.well,count(distinct nuc.id) as distinctNuclei FROM OMNISPHERO_NUCLEUS nuc INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON'...
' map.nucleusid=nuc.id INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid	INNER JOIN OMNISPHERO_NEURITE neu ON sub.neuriteid=neu.id INNER JOIN OMNISPHERO_EXPERIMENT ex on neu.experiment_id = ex.id'...    
        ' GROUP BY ex.name,ex.type, ex.well, sub.id order by sub.id asc) as innerResult WHERE innerResult.distinctNuclei >= 3) AND ex.name=''' expName ''' AND ex.well = ''' currentWell '''  AND ex.type=' typeID ' GROUP BY nuc.id, ex.name, ex.well,ex.type ORDER BY sublengthtotal DESC LIMIT ' num2str(99999) ' ) as mainneurite group by expname, wellname, typeid order by expname,typeid, wellname'];

    selectTwoUnlimited = ['select expname,typeid,wellname, avg(sublength) as Average_Subneuritelength_PerNucleus FROM(SELECT nuc.id,ex.name as expname,ex.well as wellname,ex.type as typeid, avg(sub.length) as sublength, sum(sub.length) as sublengthtotal FROM OMNISPHERO_NUCLEUS nuc'...
' INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON map.nucleusid=nuc.id INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid	INNER JOIN OMNISPHERO_NEURITE neu on neu.id = sub.neuriteid'...
' INNER JOIN OMNISPHERO_EXPERIMENT ex ON neu.experiment_id = ex.id WHERE ((neu.totallength > 10 AND neu.totallength <= 100 AND neu.totalbranchingpoints <3 ) OR (neu.totallength > 100 AND neu.totallength <= 200 AND neu.totalbranchingpoints <4 ) OR (neu.totallength > 200 AND neu.totallength <= 300 AND neu.totalbranchingpoints <4)OR (neu.totallength > 300 AND neu.totallength <= 360 AND neu.totalbranchingpoints <5))'...
' AND sub.id not in (SELECT subid FROM (SELECT sub.id as subid, ex.name,ex.type,ex.well,count(distinct nuc.id) as distinctNuclei FROM OMNISPHERO_NUCLEUS nuc'...
        ' INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON map.nucleusid=nuc.id INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid	INNER JOIN OMNISPHERO_NEURITE neu ON sub.neuriteid=neu.id'...
        ' INNER JOIN OMNISPHERO_EXPERIMENT ex on neu.experiment_id = ex.id GROUP BY ex.name,ex.type, ex.well, sub.id order by sub.id asc) as innerResult WHERE innerResult.distinctNuclei >= 3)'...
        ' AND ex.name=''' expName ''' AND ex.well = ''' currentWell '''  AND ex.type=' typeID...
' GROUP BY nuc.id, ex.name, ex.well,ex.type ORDER BY sublengthtotal DESC LIMIT ' num2str(99999) ' ) as mainneurite group by expname, wellname, typeid order by expname,typeid, wellname'];

selectThreeUnlimited = ['select expname,typeid,wellname, avg(sublengthtotal) as Total_Subneuritelength_PerNucleus, avg(numberbranchingpoints) as numberbranchingpoints,avg(distinctsubneurites) as distinctsubneurites FROM'...
' (SELECT nuc.id,ex.name as expname,ex.well as wellname,ex.type as typeid,sum(sub.length) as sublengthtotal, sum(sub.branchingpoints) as numberbranchingpoints, count(distinct map.subneuriteid) as distinctsubneurites FROM OMNISPHERO_NUCLEUS nuc'...
' INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON map.nucleusid=nuc.id INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid	INNER JOIN OMNISPHERO_NEURITE neu on neu.id = sub.neuriteid'...
' INNER JOIN OMNISPHERO_EXPERIMENT ex ON neu.experiment_id = ex.id'...
' WHERE ((neu.totallength > 10 AND neu.totallength <= 100 AND neu.totalbranchingpoints <3 ) OR (neu.totallength > 100 AND neu.totallength <= 200 AND neu.totalbranchingpoints <4 ) OR (neu.totallength > 200 AND neu.totallength <= 300 AND neu.totalbranchingpoints <4)OR (neu.totallength > 300 AND neu.totallength <= 360 AND neu.totalbranchingpoints <5))'...
' AND sub.id not in (SELECT subid FROM (SELECT sub.id as subid, ex.name,ex.type,ex.well,count(distinct nuc.id) as distinctNuclei FROM OMNISPHERO_NUCLEUS nuc	INNER JOIN OMNISPHERO_SUBNEURITE_NUCLEUS_MAPPING map ON map.nucleusid=nuc.id'...
' INNER JOIN OMNISPHERO_SUBNEURITE sub on sub.id = map.subneuriteid INNER JOIN OMNISPHERO_NEURITE neu ON sub.neuriteid=neu.id INNER JOIN OMNISPHERO_EXPERIMENT ex on neu.experiment_id = ex.id'...	
' GROUP BY ex.name,ex.type, ex.well, sub.id order by sub.id asc) as innerResult WHERE innerResult.distinctNuclei >= 3) AND ex.name=''' expName ''' AND ex.well = ''' currentWell '''  AND ex.type=' typeID ' GROUP BY nuc.id, ex.name, ex.well,ex.type ORDER BY sublengthtotal DESC LIMIT ' num2str(99999) ' ) as mainneurite group by expname, wellname, typeid order by expname,typeid, wellname'];

        oneStatement = conn.prepareStatement(selectOne); 
        oneReader = oneStatement.executeQuery();
        if(oneReader.next())        
            meanNeuriteLength = char(oneReader.getString(4));
        else
            meanNeuriteLength = 'NaN';
        end
        oneStatement.close();
        oneReader.close();

        twoStatement = conn.prepareStatement(selectTwo);
        twoReader = twoStatement.executeQuery();
        
        if(twoReader.next())
            meanSingleNeuriteLength = char(twoReader.getString(4));
        else
            meanSingleNeuriteLength = 'NaN';
        end
        twoStatement.close();
        twoReader.close();

        threeStatement = conn.prepareStatement(selectThree);
        threeReader = threeStatement.executeQuery();
        
        if(threeReader.next())
            meanNumberBranchingPoints = char(threeReader.getString(5));
            meanSubneuriteCount = char(threeReader.getString(6));
        else
            meanNumberBranchingPoints='NaN';
            meanSubneuriteCount='NaN';
        end
        threeStatement.close();
        threeReader.close();   
        
        
        oneStatementUnlimited = conn.prepareStatement(selectOneUnlimited); 
        oneReaderUnlimited = oneStatementUnlimited.executeQuery();
        if(oneReaderUnlimited.next())        
            meanNeuriteLengthUnlimited = char(oneReaderUnlimited.getString(4));
        else
            meanNeuriteLengthUnlimited = 'NaN';
        end
        oneStatementUnlimited.close();
        oneReaderUnlimited.close();

        twoStatementUnlimited = conn.prepareStatement(selectTwoUnlimited);
        twoReaderUnlimited = twoStatementUnlimited.executeQuery();
        
        if(twoReaderUnlimited.next())
            meanSingleNeuriteLengthUnlimited = char(twoReaderUnlimited.getString(4));
        else
            meanSingleNeuriteLengthUnlimited = 'NaN';
        end
        twoStatementUnlimited.close();
        twoReaderUnlimited.close();

        threeStatementUnlimited = conn.prepareStatement(selectThreeUnlimited);
        threeReaderUnlimited = threeStatementUnlimited.executeQuery();
        
        if(threeReaderUnlimited.next())
            meanNumberBranchingPointsUnlimited = char(threeReaderUnlimited.getString(5));
            meanSubneuriteCountUnlimited = char(threeReaderUnlimited.getString(6));
        else
            meanNumberBranchingPointsUnlimited='NaN';
            meanSubneuriteCountUnlimited='NaN';
        end
        threeStatement.close();
        threeReader.close(); 
        
    end
    ps.close();
    rs.close();
    
    
    
    
    
    
    
    
    
    %Normal Quality Check
    if(isKey(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
        ManualNeurons = logical(neuronHandler.ManualNeuronPositionsSparse(selectedWell) .* currentRing);
    else
        ManualNeurons = sparse(logical(zeros(sizeY,sizeX)));
    end


filterDistance=0;
nonFilterDistance=0;
SphereArea=0;
NucleusArea64=0;
markerPointCoordinates=0;
ManualNeuronsInvert = logical(full(logical(~ManualNeurons)));

if(~isa(neuronHandler.NeuronPositionsEdgeFill,'containers.Map'))
    neuronHandler.NeuronPositionsEdgeFill = containers.Map();
end
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
CellomicsNeuronsInvertWithoutCut = logical(full(~CellomicsNeuronsWithoutCut));
ManualNeuronsInvertWithoutCut = logical(full(~ManualNeuronsWithoutCut));
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

if(isfield(handles,'OligoCoordinates') && isprop(handles.OligoCoordinates,'NeuronPositionsSkeletonization') && everyThingAvailable)
        distFilterOligos = num2str(distFilterOligos,'%f');
        distWithoutFilterOligos = num2str(distWithoutFilterOligos,'%f');
        %densDistOligos = zeros(100);
        
        
        
       SkeletonOligos = logical(oligoHandler.NeuronPositionsSkeletonization(selectedWell) .* currentRing);
       CompFillOligos = logical(oligoHandler.NeuronPositionsEdgeFill(selectedWell) .* currentRing);
       OligosDeleted = oligoHandler.NeuronPositionsSkelDeleted(selectedWell);

        SkeletonOligos = SkeletonOligos - (OligosDeleted);
        ind = SkeletonOligos<0;
        SkeletonOligos(ind) = 0;   
        bothManualSkeletonizedOligos = nnz(SkeletonOligos .* ManualNeurons);        
        AdditionalSkeletonManualOligos = nnz(SkeletonOligos .* ManualNeuronsInvert);

        bothManualCompositeFillOligos = nnz(CompFillOligos .* ManualNeurons);
        AdditionalCompositeFillManualOligos = nnz(CompFillOligos .* ManualNeuronsInvert);
    
        SkeletonOligosInvert = ~SkeletonOligos;
        CompFillOligosInvert = ~CompFillOligos;
        
        AdditionalManualSkeletonOligos = nnz(ManualNeurons .* SkeletonOligosInvert);
        AdditionalManualCompFillOligos = nnz(ManualNeurons .* CompFillOligosInvert);
else
    densDistOligos = zeros(100);
    densDistCompFillOligos = zeros(100);
    distFilterOligos = 0;
    distWithoutFilterOligos = 0;
    AdditionalManualSkeletonOligos=0;
    AdditionalSkeletonManualOligos=0;
    bothManualSkeletonizedOligos=0;
    bothManualCompositeFillOligos = 0;
    AdditionalCompositeFillManualOligos=0;
    AdditionalManualCompFillOligos=0;
    
end
currentRing=0;

        distFilter = num2str(distFilter,'%f');
        distWithoutFilter = num2str(distWithoutFilter,'%f');
        sphereRadius = num2str(sphereRadius,'%f');
        for i=1:optionHandler.DensityDistributionRingNumber        
          if(isKey(imageHandler.MigDistWellMapping,selectedWell))
            manualDist = num2str(imageHandler.MigDistWellMapping(selectedWell), '%f')
            manualDist = strrep(manualDist,'.',',');
            exportText = [exportText selectedWell ';' num2str(i) ';' distFilter ';' sphereRadius ';' manualDist ';' distWithoutFilter ';' num2str(densDist(i,1)),'.',',' ';' num2str(densDist(i,2)) ';' num2str(densDist(i,3)) ';' num2str(densDist(i,4)) ';' num2str(densDist(i,5)) ';' num2str(densDist(i,6)) ';' num2str(densDist(i,7)) ';' num2str(densDist(i,8)) ';' num2str(densDist(i,9)) ';;' num2str(densDistCompFillOligos(i,5)) ';' num2str(densDistCompFillOligos(i,8)) ';' num2str(densDistOligos(i,5)) ';' num2str(densDistOligos(i,8)) ';;' num2str(meanGrad) ';' num2str(meanGradNucImage) ';;' meanNeuriteLengthUnlimited ';' meanSingleNeuriteLengthUnlimited ';' meanNumberBranchingPointsUnlimited ';' meanSubneuriteCountUnlimited ';;' meanNeuriteLength ';' meanSingleNeuriteLength ';' meanNumberBranchingPoints ';' meanSubneuriteCount ';;' num2str(bothManualSkeletonized) ';' num2str(AdditionalManualSkeleton) ';' num2str(AdditionalSkeletonManual) ';;' num2str(bothManualCompositeFillOligos) ';' num2str(AdditionalManualCompFillOligos) ';' num2str(AdditionalCompositeFillManualOligos) ';;' num2str(bothManualSkeletonizedOligos) ';' num2str(AdditionalManualSkeletonOligos) ';' num2str(AdditionalSkeletonManualOligos) ';' sprintf('\r\n')];
          else
            exportText = [exportText selectedWell ';' num2str(i) ';' distFilter ';' sphereRadius ';;' distWithoutFilter ';' num2str(densDist(i,1)) ';' num2str(densDist(i,2)) ';' num2str(densDist(i,3)) ';' num2str(densDist(i,4)) ';' num2str(densDist(i,5)) ';' num2str(densDist(i,6)) ';' num2str(densDist(i,7)) ';' num2str(densDist(i,8)) ';' num2str(densDist(i,9)) ';;' num2str(densDistCompFillOligos(i,5)) ';' num2str(densDistCompFillOligos(i,8)) ';' num2str(densDistOligos(i,5)) ';' num2str(densDistOligos(i,8)) ';;' num2str(meanGrad) ';' num2str(meanGradNucImage) ';;' meanNeuriteLengthUnlimited ';' meanSingleNeuriteLengthUnlimited ';' meanNumberBranchingPointsUnlimited ';' meanSubneuriteCountUnlimited ';;' meanNeuriteLength ';' meanSingleNeuriteLength ';' meanNumberBranchingPoints ';' meanSubneuriteCount ';;' num2str(bothManualSkeletonized) ';' num2str(AdditionalManualSkeleton) ';' num2str(AdditionalSkeletonManual) ';;' num2str(bothManualCompositeFillOligos) ';' num2str(AdditionalManualCompFillOligos) ';' num2str(AdditionalCompositeFillManualOligos) ';;' num2str(bothManualSkeletonizedOligos) ';' num2str(AdditionalManualSkeletonOligos) ';' num2str(AdditionalSkeletonManualOligos) ';' sprintf('\r\n')];
          end
        end
    else
        for i=1:optionHandler.DensityDistributionRingNumber
            exportText = [exportText selectedWell ';' num2str(i) ';' num2str(0) ';' num2str(0) ';;' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';;' num2str(0) ';' num2str(0) ';' num2str(0) ';' num2str(0) ';;' num2str(meanGrad) ';' num2str(meanGradNucImage) ';;' meanNeuriteLengthUnlimited ';' meanSingleNeuriteLengthUnlimited ';' meanNumberBranchingPointsUnlimited ';' meanSubneuriteCountUnlimited ';;' meanNeuriteLength ';' meanSingleNeuriteLength ';' meanNumberBranchingPoints ';' meanSubneuriteCount ';;' num2str(bothManualSkeletonized) ';' num2str(AdditionalManualSkeleton) ';' num2str(AdditionalSkeletonManual) ';;' num2str(bothManualCompositeFillOligos) ';' num2str(AdditionalManualCompFillOligos) ';' num2str(AdditionalCompositeFillManualOligos) ';;' num2str(bothManualSkeletonizedOligos) ';' num2str(AdditionalManualSkeletonOligos) ';' num2str(AdditionalSkeletonManualOligos) ';' sprintf('\r\n')];
        end
    end
    waitbar((currentWellIndex/numel(wellList)),waitbarHandle);
end
conn.close();
%Write exportText to File
fprintf(fileID,'%s',exportText);
close(waitbarHandle);
 


function conn=OpenPostgreConnection(answer)
list=dir('*jdbc*');  
jdbcFolder = pwd;
if(numel(list)<=0)
        %errordlg('No Postgre JDBC driver found within working directory. Please download and save the Postgre JDBC driver to the working directory!' ,'Error');
        jdbcFolder = uigetdir('','Select JDBC driver folder');
        list=dir(strcat(jdbcFolder,'/','*jdbc*'));     
end
% [FileName,PathName,FilterIndex] = uigetfile('*.jdbc*','Select Postgre JDBC Driver');
FileName = list(1).name;
PathName = jdbcFolder;
unzip(strcat(PathName,'/',FileName),PathName);

%javaclasspath(strcat(FileName(1:length(FileName)-4)))
javaclasspath(jdbcFolder);

props=java.util.Properties;
props.setProperty('user',answer{4});
props.setProperty('password', answer{5});
driver=org.postgresql.Driver;
url=strcat('jdbc:postgresql://',answer{1},':',answer{2},'/',answer{3});
conn=driver.connect(url,props);

% --------------------------------------------------------------------
function High_Throughput_Export_Callback(hObject, eventdata, handles)
% hObject    handle to High_Throughput_Export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    %Create JDBC connection

%2. Select folder with Save files
htFolder = uigetdir('','Select Batch folder');
htFolder = strcat(htFolder,'/');
list=dir(strcat(htFolder,'*.zip'));

imFolder = uigetdir('','Select image parent folder');
imFolder = strcat(imFolder,'/');
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
% %1. Calculate
NeuriteLengthMatrix = containers.Map();
%Ask user for DB connection details
%Ask user for password
prompt = {'Postgre Hostname', 'Postgre Port', 'DB Name', 'Postgre Username', 'Password'};
dlg_title='Postgre Connection';
num_lines=1;
defaultans={'localhost','5432','postgres','postgres','postgres'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
percNeurons = inputdlg('Enter percentage of Neurons for Morphological Analysis','Percentage Neurons',1,{'100'});
percNeurons = percNeurons{1};
percNeurons = str2num(percNeurons);
percNeurons = percNeurons/100;
conn = OpenPostgreConnection(answer);
%Check if necessary tables exist. If not, create them.   
CreateNonExistingTable(conn);
conn.close();


set(handles.cbNucleusPicture ,'Value',1);
set(handles.cbNeuritePicture ,'Value',1);

PathName= strcat(htFolder,'/Results');
%Iterate over all save files in folder
mkdir(strcat(htFolder,'Results'));
resultFolder = strcat(htFolder,'Results');


%Get Config xml file
configfilename=strcat(htFolder,'Config.xml');
optionFields = zeros(0);
if(exist(configfilename,'file'))
    struct = xml2struct(configfilename);
    %Get child objects
    configStruct = struct.Configuration;
    optionFields = fieldnames(configStruct);
end


for(i=1:numel(list))
    selectedWell = list(i).name(1:length(list(i).name)-8);    
    %4. Load current file
    Load(htFolder,list(i).name,handles);  
    handles = guidata(handles.figure1);
    neuronHandler = handles.NeuronCoordinates;
    QuickFixManualPositions(handles);
    %Set Image folder    
    LoadImageFolder(strcat(imFolder,selectedWell),WellType.Undefined,0,0,0,0,0,0,handles); 
    
    %5. Set changed parameters    
    handles = guidata(handles.figure1);
    optionHandler = handles.OptionHandler;
    for(j=1:numel(optionFields))
        fieldName = optionFields{j};
        fieldValue = configStruct.(optionFields{j}).Text;
        optionHandler.SetStringValue(fieldName,fieldValue);
    end
    %6. Run "Neuronal Quantification on all Wells"
    handles.OptionHandler = optionHandler;
    guidata(handles.figure1, handles);
    handles = guidata(handles.figure1);
    
    
    ChangeWell(1,handles);
    %Calculate Neuronal Morphology
  %DEBUG COMMENT OUT
    ExportNeuriteLengthStat(neuronHandler.NeuronPositionsSkeletonization,'DataOmnisphero',selectedWell,answer,handles);
    FileName = strcat(selectedWell,'.csv');
    conn = OpenPostgreConnection(answer);
    TotalExport(FileName,PathName,percNeurons,conn,handles); 
    conn.close();
    %Make Total Export for current Save
end

    
function [filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, optionHandler, NucleusM, foldername, sizeY, sizeX)
%Get density distribution to exclude too dense areas
%if(saveCircle)
path = [foldername '/ConvertedCellomics/MarkerPointCoordinates-' selectedWell '.mat'];
if(exist(path,'file'))        
    load(path);
else
     [filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = calculateMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, optionHandler, NucleusM, [foldername '/ConvertedCellomics'], sizeY, sizeX); 
     %Save hInner to file         
     save('-v7.3',strcat(foldername,'/ConvertedCellomics/MarkerPointCoordinates-',selectedWell),'markerPointCoordinates', 'filterDistance', 'nonFilterDistance', 'sphereRadius', 'SphereArea');
end


% --------------------------------------------------------------------
function NeuronalOligoQuantification_Callback(hObject, eventdata, handles)
% hObject    handle to NeuronalOligoQuantification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------
waitbarHandle = waitbar(0.5,'Executing Skeletonization. This could take up to 45 minutes.', 'WindowStyle', 'modal');
try
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
csvHandler = handles.CSVCoordinates;
if(isfield(handles,'OligoCoordinates') && isprop(handles.OligoCoordinates,'NeuronPositionsSkeletonization'))
   oligoHandler = handles.OligoCoordinates;
else
   %oligoHandler = CSVCoordinates();
end
foldername = imageHandler.Foldername;
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
selectedWell = wellList{selectedWellNumber};
NucleusM = csvHandler.CellPosMatrix(selectedWell);
[sizeY sizeX] = size(imageHandler.NeuriteImage);
if(imageHandler.ZoomState)
    xMin = uint32(imageHandler.ZoomState(1));
    yMin = uint32(imageHandler.ZoomState(2));
    xMax = uint32(imageHandler.ZoomState(3) + xMin);
    yMax = uint32(imageHandler.ZoomState(4) + yMin);
    [neuronHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell, yMin, yMax, xMin, xMax, 0, optionHandler, neuronHandler, NucleusM,foldername,'Neurite');
    if(isfield(handles,'OligoCoordinates') && isprop(handles.OligoCoordinates,'NeuronPositionsSkeletonization'))
        [oligoHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell, yMin, yMax, xMin, xMax, 0, optionHandler, oligoHandler, NucleusM,foldername,'Oligo');
    end
    [neuronHandler a b c] = SkeletonizationNeuronsAlgorithm(selectedWell,yMin, yMax, xMin, xMax, 0,optionHandler, sizeY,sizeX,foldername,neuronHandler,NucleusM,0,'Neurite',0);
    [neuronHandler a b c] = SkeletonizationNeuronsAlgorithm(selectedWell,yMin, yMax, xMin, xMax, 0,optionHandler, sizeY,sizeX,foldername,neuronHandler,NucleusM,0,'Oligo',neuronHandler);
    RefreshZoomedImage(handles);
else
    %Check size of image handler. If image is bigger than 10.000 pixels,
    %split it up into four pics
    if(size(imageHandler.NucleusImage,1) > 10000 || size(imageHandler.NucleusImage,2) > 10000)        
        sizeXsmall = sizeX/2;
        sizeYsmall = sizeY/2;
        [neuronHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell,1,sizeYsmall,1,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,'Neurite');
        [neuronHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell,1,sizeYsmall,sizeXsmall+1,sizeX,1,optionHandler,neuronHandler,NucleusM,foldername,'Neurite');
        [neuronHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell,sizeYsmall+1,sizeY,sizeXsmall+1,sizeX,1,optionHandler,neuronHandler,NucleusM,foldername,'Neurite');
        [neuronHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell,sizeYsmall+1,sizeY,1,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,'Neurite');
        if(isfield(handles,'OligoCoordinates') && isprop(handles.OligoCoordinates,'NeuronPositionsSkeletonization'))
            [oligoHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell,1,sizeYsmall,1,sizeXsmall,1,optionHandler,oligoHandler,NucleusM,foldername,'Oligo');
            [oligoHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell,1,sizeYsmall,sizeXsmall+1,sizeX,1,optionHandler,oligoHandler,NucleusM,foldername,'Oligo');
            [oligoHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell,sizeYsmall+1,sizeY,sizeXsmall+1,sizeX,1,optionHandler,oligoHandler,NucleusM,foldername,'Oligo');
            [oligoHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell,sizeYsmall+1,sizeY,1,sizeXsmall,1,optionHandler,oligoHandler,NucleusM,foldername,'Oligo');
        end
    else
        [neuronHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell,1,sizeY,1,sizeX,1,optionHandler,neuronHandler,NucleusM,foldername,'Neurite');
        if(isfield(handles,'OligoCoordinates') && isprop(handles.OligoCoordinates,'NeuronPositionsSkeletonization'))
            [oligoHandler a b c] = EdgeFillNeuronsAlgorithm(selectedWell,1,sizeY,1,sizeX,1,optionHandler,oligoHandler,NucleusM,foldername,'Oligo');
        end
    end
    [neuronHandler a b c] = SkeletonizationNeuronsAlgorithm(selectedWell,1, sizeY, 1, sizeX, 1,optionHandler, sizeY,sizeX,foldername,neuronHandler,NucleusM,1,'Neurite',0);    
    if(isfield(handles,'OligoCoordinates') && isprop(handles.OligoCoordinates,'NeuronPositionsSkeletonization'))
        [oligoHandler a b c] = SkeletonizationNeuronsAlgorithm(selectedWell,1, sizeY, 1, sizeX, 1,optionHandler, sizeY,sizeX,foldername,oligoHandler,NucleusM,1,'Oligo',neuronHandler);    
    end
    RefreshUnzoomedImage(handles);
end
handles.NeuronCoordinates = neuronHandler;
if(isfield(handles,'OligoCoordinates') && isprop(handles.OligoCoordinates,'NeuronPositionsSkeletonization'))
    handles.OligoCoordinates = oligoHandler;
end
guidata(handles.figure1, handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end
close(waitbarHandle);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over popupPlus1.
function popupPlus1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to popupPlus1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SplitMultiSphereWells_Callback(hObject, eventdata, handles)
% hObject    handle to SplitMultiSphereWells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
neuronHandler = handles.NeuronCoordinates;
optionHandler = handles.OptionHandler;
csvHandler = handles.CSVCoordinates;

neuronSplitter = questdlg('Would you like to split up the Neuron Channel?');
oligoSplitter = questdlg('Would you like to split up the Oligo Channel?');
astroSplitter = questdlg('Would you like to split up the Astro Channel?');

if(strcmp(neuronSplitter,'Yes'))
    neuronSplitter = 1;
else
    neuronSplitter = 0;
end

if(strcmp(oligoSplitter,'Yes'))
    oligoSplitter = 1;
    %oligoHandler = handles.OligoCoordinates;
else
    oligoSplitter = 0;
end

if(strcmp(astroSplitter,'Yes'))
    astroSplitter = 1;
else
    astroSplitter = 0;
end





%Do Watershed
%Check if there is already a watershed image. If not, create it.
%Iterate over all wells
wellList = get(handles.lbWell, 'string');
%waitbarHandle = waitbar(0,'Watersheding Nuclei. Please wait.', 'WindowStyle', 'modal');
foldern = imageHandler.Foldername;
foldername = [imageHandler.Foldername '/ConvertedCellomics'];
imageHandler.NucleusImage = 0;
imageHandler.NeuriteImage = 0;
for i=optionHandler.StartWell:numel(wellList)
    %waitbar((i/numel(wellList)),waitbarHandle);
    selectedWell=wellList(i);
    selectedWell=selectedWell{1};
    selectedWellLong=selectedWell;
    if(length(selectedWell) == 2)
         selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    end
    imagePathWatershed = [foldername '/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
    if(exist(imagePathWatershed,'file'))
        nucleusPic = imread(imagePathWatershed);
    else
        imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
        nucleusPic = imread(imagePathNucleusBig);
        nucleusPic(nucleusPic<optionHandler.NucleusThreshold) = 0;
        nucleusPic(nucleusPic>=optionHandler.NucleusThreshold) = 255;
        nucleusPic=imcomplement(nucleusPic);
        nucleusPic = uint8(nucleusPic);
        javaaddpath('lib/ij.jar');
        javaaddpath('lib/mij.jar');
        MIJ.start(false);
        MIJ.createImage(nucleusPic);
        MIJ.run('Options...', 'iterations=1 count=1 edm=Overwrite do=Nothing');
        MIJ.run('Watershed');
        %nucleusPicBin = logical(MIJ.getCurrentImage);
        MIJ.run('Save', 'save=temp.tif');
        MIJ.closeAllWindows;
        MIJ.exit;
        nucleusPic = imread('temp.tif');
        delete temp.tif;
        nucleusPic = logical(1-nucleusPic);
        imwrite(nucleusPic,imagePathWatershed);
    end
end
nucleusPic = 0;
[newCSVHandler wellList] = Imagecutter_V5_9(neuronSplitter, oligoSplitter, astroSplitter, foldern, csvHandler, optionHandler);
PostProcessImageCutter(wellList,newCSVHandler,handles);
ChangeWell(1,handles);
%handles.CSVCoordinates = csvHandler;
%guidata(handles.figure1, handles);


function everyThingAvailable=ExportDistanceDistribution2(PathName,expName,currentHandler,modeStr,handles)
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
csvHandler = handles.CSVCoordinates;

SphereAreaSizeX = optionHandler.MigrationDistanceDensityImageXSize;
SphereAreaSizeY = optionHandler.MigrationDistanceDensityImageYSize;

[sizeY sizeX] = size(imageHandler.NucleusImage);
%Iterate over all Well numbers
wellNumberList = zeros(0);
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
checkMap = containers.Map();
everyThingAvailable=1;
%exportText = sprintf('Well;Ring Number;Migration Distance(px);Migration Radius (px);Manual measured Migration Distance (px);Migration Distance without filter (px);Number of Pixels;Number Nuclei;Cellomics Neurons;Manual Neurons;Skeleton Neurons;Cellomics Neurons per Nuclei;Manual Neurons per Nuclei;Skeleton Neurons per Nuclei;Nuclei per Pixels;;Mean Gradient Nuclei Positions;Mean Gradient Nucleus Image;;Total Subneuritelength per Nucleus;Average Subneuritelength per Nucleus;Mean Number Branchingpoints;Mean Subneurite Count;;Neurons Omnisphero & Manual;Additional Manual;Additional Omnisphero;\r\n');
for currentWellIndex=optionHandler.StartWell:numel(wellList)
    selectedWell = wellList(currentWellIndex);
    selectedWell = selectedWell{1};
    if(isKey(currentHandler,selectedWell))
        NucleusM = csvHandler.CellPosMatrix(selectedWell);    
        NeuronM = currentHandler(selectedWell);


        currentWellNumber = str2num(GetWellNumberFromString(selectedWell));
        if(~any(wellNumberList == currentWellNumber))
            wellNumberList = [wellNumberList currentWellNumber];
        end
        [filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, SphereAreaSizeX, SphereAreaSizeY, optionHandler, NucleusM, imageHandler.Foldername, sizeY, sizeX); 
        currentDistancesList = zeros(0);


        innerMask = ones(uint16(sizeY/(1/ScalingFactor)),uint16(sizeX/(1/ScalingFactor)));   

        mk = markerPointCoordinates('10')./(1/ScalingFactor);
        innerMask = roipoly(innerMask,mk(:,1),mk(:,2));
        currentRing =logical(ones(int32(sizeY/(1/ScalingFactor)),int32(sizeX/(1/ScalingFactor))) - innerMask);
        currentRing =imresize(currentRing,[sizeY,sizeX]);
        NeuronM = NeuronM .* currentRing;
        %Original Code:
        mkFirst = markerPointCoordinates('10');
        mkFirst=double(mkFirst);
        %Iterate over all Nuclei for current well
        [nucleiY nucleiX] = find(NeuronM);
        for(i=1:numel(nucleiY))
            currentY = nucleiY(i);
            currentX = nucleiX(i);
            minDist = 99999999;
            minIndex = 0;
            for(v=1:size(mkFirst,1))
                currentYSphereCore = mkFirst(v,1);
                currentXSphereCore = mkFirst(v,2);
                distance = pdist([currentY currentX;currentYSphereCore currentXSphereCore]);
                if (distance < minDist)
                    minDist=distance;
                    minIndex = v;
                end                            
            end
            currentDistancesList = [currentDistancesList minDist];
        end
        checkMap(selectedWell) = currentDistancesList;    
    else
        everyThingAvailable=0;
    end
end 

if(everyThingAvailable)
    distancesList = uint16(zeros(numel(wellNumberList),5000));    
    maxJ = 0;
    for (i=1:numel(wellNumberList))
        %Create list of all elements in current WellNumber    
        currentWellNumberList = zeros(0);
        if(isKey(checkMap,strcat('A',num2str(wellNumberList(i)))))
            smallList = checkMap(strcat('A',num2str(wellNumberList(i))));
            currentWellNumberList = [currentWellNumberList smallList];
        end
        if(isKey(checkMap,strcat('B',num2str(wellNumberList(i)))))
            smallList = checkMap(strcat('B',num2str(wellNumberList(i))));
            currentWellNumberList = [currentWellNumberList smallList];
        end
        if(isKey(checkMap,strcat('C',num2str(wellNumberList(i)))))
            smallList = checkMap(strcat('C',num2str(wellNumberList(i))));
            currentWellNumberList = [currentWellNumberList smallList];
        end
        if(isKey(checkMap,strcat('D',num2str(wellNumberList(i)))))
            smallList = checkMap(strcat('D',num2str(wellNumberList(i))));
            currentWellNumberList = [currentWellNumberList smallList];
        end
        if(isKey(checkMap,strcat('E',num2str(wellNumberList(i)))))
            smallList = checkMap(strcat('E',num2str(wellNumberList(i))));
            currentWellNumberList = [currentWellNumberList smallList];
        end
        if(isKey(checkMap,strcat('F',num2str(wellNumberList(i)))))
            smallList = checkMap(strcat('F',num2str(wellNumberList(i))));
            currentWellNumberList = [currentWellNumberList smallList];
        end
        if(isKey(checkMap,strcat('G',num2str(wellNumberList(i)))))
            smallList = checkMap(strcat('G',num2str(wellNumberList(i))));
            currentWellNumberList = [currentWellNumberList smallList];
        end
        if(isKey(checkMap,strcat('H',num2str(wellNumberList(i)))))
            smallList = checkMap(strcat('H',num2str(wellNumberList(i))));
            currentWellNumberList = [currentWellNumberList smallList];
        end
        currentWellNumberList = sort(currentWellNumberList);
        %for(j=1:numel(currentWellNumberList))        
        distancesList(i,1:numel(currentWellNumberList)) = currentWellNumberList;
        %end    
        if(numel(currentWellNumberList) > maxJ)
            maxJ=numel(currentWellNumberList);
        end
    end

    exportText = '';
    for (i=1:numel(wellNumberList))
        exportText = [exportText num2str(wellNumberList(i)) ';'];
    end
    exportText = [exportText sprintf('\r\n')];


    fileID = fopen(strcat(PathName,'/',expName,'_',modeStr,'_DistanceDistribution.csv'),'w');
    foldername = [imageHandler.Foldername '/ConvertedCellomics'];
    for(line=1:maxJ)
        column=1;
        firstWellNumber=wellNumberList(1);
        while(column<=numel(wellNumberList))        
            if(distancesList(column,line) == 0)
                exportText = [exportText ';'];
            else
                exportText = [exportText num2str(distancesList(column,line)) ';'];
            end
            column = column +1;
        end
        exportText = [exportText ';' sprintf('\r\n')];
    end
    %Write exportText to File
    fprintf(fileID,'%s',exportText);
    fclose(fileID);
end

function ExportDistanceDistribution(PathName,conn,expName,handles)
imageHandler = handles.ImageHandler;
optionHandler = handles.OptionHandler;
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
%expName=get(gcf,'name');
%expName=expName(1:length(expName)-4);
fileID = fopen(strcat(PathName,'/',expName,'DistanceDistribution.csv'),'w');
foldername = [imageHandler.Foldername '/ConvertedCellomics'];
%Iterate over all Well numbers
wellNumberList = zeros(0);
selectedWellNumber = get(handles.lbWell,'Value');
wellList = get(handles.lbWell, 'string');
%exportText = sprintf('Well;Ring Number;Migration Distance(px);Migration Radius (px);Manual measured Migration Distance (px);Migration Distance without filter (px);Number of Pixels;Number Nuclei;Cellomics Neurons;Manual Neurons;Skeleton Neurons;Cellomics Neurons per Nuclei;Manual Neurons per Nuclei;Skeleton Neurons per Nuclei;Nuclei per Pixels;;Mean Gradient Nuclei Positions;Mean Gradient Nucleus Image;;Total Subneuritelength per Nucleus;Average Subneuritelength per Nucleus;Mean Number Branchingpoints;Mean Subneurite Count;;Neurons Omnisphero & Manual;Additional Manual;Additional Omnisphero;\r\n');
for currentWellIndex=optionHandler.StartWell:numel(wellList)
    selectedWell = wellList(currentWellIndex);
    selectedWell = selectedWell{1};
    currentWellNumber = str2num(GetWellNumberFromString(selectedWell));
    if(~any(wellNumberList == currentWellNumber))
        wellNumberList = [wellNumberList currentWellNumber];
    end    
end
exportText = '';
for (i=1:numel(wellNumberList))
    exportText = [exportText num2str(wellNumberList(i)) ';'];
end
exportText = [exportText sprintf('\r\n')];
distancesList = uint16(zeros(numel(wellNumberList),5000));
maxJ = 0;
for (i=1:numel(wellNumberList))
    selectString = ['select DISTANCE_TO_SPHERECORE from omnisphero_experiment exp INNER JOIN OMNISPHERO_NEURON neu ON exp.id = neu.experiment_fk WHERE exp.name LIKE ''' expName ''' AND exp.well Like ''%' num2str(wellNumberList(i)) ''' order by 1'];
    ps = conn.prepareStatement(selectString);
    rs=ps.executeQuery();    
    j=1;
    while(rs.next())
        currentDistance = char(rs.getString(1));
        distancesList(i,j) = str2num(currentDistance);
        j=j+1;
    end
    if(j>maxJ-1)
        maxJ = j-1;
    end
    ps.close();
    rs.close();
end
for(line=1:maxJ)
    column=1;
    firstWellNumber=wellNumberList(1);
    while(column<=numel(wellNumberList))        
        if(distancesList(column,line) == 0)
            exportText = [exportText ';'];
        else
            exportText = [exportText num2str(distancesList(column,line)) ';'];
        end
        column = column +1;
    end
    exportText = [exportText ';' sprintf('\r\n')];
end
%Write exportText to File
fprintf(fileID,'%s',exportText);
fclose(fileID);


function wellNumber=GetWellNumberFromString(wellName)
     if(length(wellName)==2)
         wellNumber = wellName(2);
     elseif(length(wellName) == 3)
        if(strcmp(wellName(2),'0'))
            wellNumber = wellName(3);
        else
            wellNumber = wellName(2:3);
        end
     end
         

% --------------------------------------------------------------------
function LoadZebrafish_Callback(hObject, eventdata, handles)
% hObject    handle to LoadZebrafish (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
neuronHandler = CSVCoordinates();
csvHandler = CSVCoordinates();
 imageHandler = ImageHandler();
 optionHandler = Option();
 handles.OptionHandler = optionHandler;
optionHandler.FileEnding = '.png';
[wellList folder ScalingFactor csvHandler] = ZebraFish(csvHandler);
optionHandler.ScalingFactor = ScalingFactor;
 imageHandler.Foldername = folder;
set(handles.lbWell, 'String', wellList);
handles.CSVCoordinates = csvHandler;
handles.NeuronCoordinates = neuronHandler;
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles);


% --------------------------------------------------------------------
function LoadGFAPS100bbraintissue_Callback(hObject, eventdata, handles)
% hObject    handle to LoadGFAPS100bbraintissue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
neuronHandler = CSVCoordinates();
csvHandler = CSVCoordinates();
 imageHandler = ImageHandler();
 optionHandler = Option();
 handles.OptionHandler = optionHandler;
optionHandler.FileEnding = '.png';
[wellList folder ScalingFactor csvHandler neuronHandler] = GFAPS100B(csvHandler,neuronHandler);
optionHandler.ScalingFactor = ScalingFactor;
 imageHandler.Foldername = folder;
set(handles.lbWell, 'String', wellList);
handles.CSVCoordinates = csvHandler;
handles.NeuronCoordinates = neuronHandler;
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles);
% --------------------------------------------------------------------
function LoadGFAPbraintissue_Callback(hObject, eventdata, handles)
% hObject    handle to LoadGFAPbraintissue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
neuronHandler = CSVCoordinates();
neuronHandler.ManualNeuronPositions = containers.Map();
csvHandler = CSVCoordinates();
 imageHandler = ImageHandler();
 optionHandler = Option();
 handles.OptionHandler = optionHandler;
optionHandler.FileEnding = '.png';
[wellList folder ScalingFactor csvHandler] = GFAPonly(csvHandler);
optionHandler.ScalingFactor = ScalingFactor;
 imageHandler.Foldername = folder;
set(handles.lbWell, 'String', wellList);
handles.CSVCoordinates = csvHandler;
handles.NeuronCoordinates = neuronHandler;
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function LoadImageCorticalNeurons_Callback(hObject, eventdata, handles)
% hObject    handle to LoadImageCorticalNeurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
neuronHandler = CSVCoordinates();
neuronHandler.ManualNeuronPositions = containers.Map();
csvHandler = CSVCoordinates();
 imageHandler = ImageHandler();
 optionHandler = Option();
 optionHandler.MinSizeNeuriteArea = 25000;
 handles.OptionHandler = optionHandler;
optionHandler.FileEnding = '.png';
[wellList folder ScalingFactor csvHandler ] = CorticalNeurons(csvHandler);
optionHandler.ScalingFactor = ScalingFactor;
 imageHandler.Foldername = folder;
set(handles.lbWell, 'String', wellList);
handles.CSVCoordinates = csvHandler;
handles.NeuronCoordinates = neuronHandler;
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles);


% --------------------------------------------------------------------




% --------------------------------------------------------------------
function LoadCorticalNeuronsNucleiSelf_Callback(hObject, eventdata, handles)
% hObject    handle to LoadCorticalNeuronsNucleiSelf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
neuronHandler = CSVCoordinates();
csvHandler = CSVCoordinates();
 imageHandler = ImageHandler();
 optionHandler = Option();
 optionHandler.MinSizeNeuriteArea = 25000;
 handles.OptionHandler = optionHandler;
optionHandler.FileEnding = '.png';
[wellList folder ScalingFactor csvHandler neuronHandler] = CorticalNeuronsNucleusSkeletonownNeuronDetect(csvHandler, neuronHandler);
optionHandler.ScalingFactor = ScalingFactor;
 imageHandler.Foldername = folder;
set(handles.lbWell, 'String', wellList);
handles.CSVCoordinates = csvHandler;
handles.NeuronCoordinates = neuronHandler;
%NeuronPositionsSkeletonization = neuronHandler.NeuronPositionsSkeletonization(wellname);
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles);





function EdgeCompositAlgorithmOligos(selectedWell,handles)
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
imagePathNucleusBig = [foldername '/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
imagePathNucleusSmall = [foldername '/' selectedWellLong 'NucleusSmall' optionHandler.FileEnding];
imagePathOligoBig = [foldername '/' selectedWellLong 'OligoBig' optionHandler.FileEnding];
imagePathOligoSmall = [foldername '/' selectedWellLong 'OligoSmall' optionHandler.FileEnding];

%Check if NeuronStatMatrix exists
if(~isa(neuronHandler.NeuronStatMatrix,'containers.Map'))
    neuronHandler.NeuronStatMatrix = containers.Map();
end
%if(isKey(neuronHandler.NeuronStatMatrix,selectedWell))
%    currentNeuronStat = neuronHandler.NeuronStatMatrix(selectedWell);
%else
%    currentNeuronStat = containers.Map();
%end

imageHandler.ResizedNeuriteImage = imread(imagePathOligoSmall);
imageHandler.ResizedNucleusImage = imread(imagePathNucleusSmall);    
imageHandler.NeuriteImage = imread(imagePathOligoBig);
imageHandler.NucleusImage = imread(imagePathNucleusBig);


[filterDistance nonFilterDistance sphereRadius SphereArea markerPointCoordinates] = loadMigrationDistance(selectedWell, densityWidth, densityHeight, optionHandler, NucleusM, imageHandler.Foldername, sizeY, sizeX); 
%MarkerPointCoordinates are available. Get second ring:
if(markerPointCoordinates ~= 0)
    hInner = impoly(handles.axes2,double(markerPointCoordinates('10')./(1/ScalingFactor)));
    innerMask = logical(createMask(hInner));
    hOuter = impoly(handles.axes2,double(markerPointCoordinates(num2str(ringNumber * 10))./(1/ScalingFactor)));
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
        if(startY-50>0 && startX-50>0 && startX+50 < sizeXsmall && startY + 50 < sizeYsmall)
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
            [xNuc,yNuc,distance] = csvHandler.FindNucleusForNeuron(xFocus,yFocus,NucleusM,75,0,0);
            EdgeCompositPositionList(yNuc,xNuc) = 1;
        end
        
    end
end
%Check if Neuronhandler already has Dictionary of NeuronPositionsEdgeComposite
if(numel(neuronHandler.OligoPositionsEdgeComposite) ==0)
    neuronHandler.OligoPositionsEdgeComposite = containers.Map();
end
neuronHandler.OligoPositionsEdgeComposite(selectedWell) = EdgeCompositPositionList;
handles.NeuronCoordinates = neuronHandler;
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function CompositeOligos_Callback(hObject, eventdata, handles)
% hObject    handle to CompositeOligos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
oligoHandler = CSVCoordinates();
%foldername = uigetdir;
%Call EdgeFill Neurons and Skeleton Neurons and play with parameters in
%Option Handler
ExecuteOptionExperimentOligo(handles,1,'Neurite');
handles.OligoCoordinates = oligoHandler;
guidata(handles.figure1, handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function ExecuteOptionExperimentOligo(handles,onlyCompFill,channelName)
h = waitbar(0.5,'Neuronal quantification in progress. This could take a few hours. Please be patient.', 'WindowStyle', 'modal');
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
optionHandler = handles.OptionHandler;
oligoHandler = CSVCoordinates();
if(strcmp(channelName,'Oligo'))
    competingHandler = handles.NeuronCoordinates;
    if(isfield(handles,'OligoCoordinates'))
       neuronHandler = handles.OligoCoordinates;
    else
       neuronHandler = CSVCoordinates();
    end
else
    neuronHandler = handles.NeuronCoordinates;
    competingHandler = 0;
end
wellList = get(handles.lbWell, 'string');
foldername = imageHandler.Foldername;
[sizeY sizeX] = size(imageHandler.NucleusImage);
%Select STARTWELL

if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
if(numel(neuronHandler.AreaDictionary) ==0)
    neuronHandler.AreaDictionary = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
    neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFill) ==0)
    neuronHandler.NeuronPositionsEdgeFill = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFillSecond) ==0)
    neuronHandler.NeuronPositionsEdgeFillSecond = containers.Map();
end
if(numel(neuronHandler.NeuronStatMatrix) ==0)
    neuronHandler.NeuronStatMatrix = containers.Map();
end

%Create Matlab job
%Leave one Core out to ensure that GUI isn't freezing
maxCores =optionHandler.NumberOfCores;
if(maxCores > 1)
    c = parcluster();
    c.NumWorkers = maxCores;
    job = createJob(c,'JobData',optionHandler);
end

for(j=optionHandler.StartWell:numel(wellList))
    selectedWell=wellList(j);
    selectedWell=selectedWell{1};
    selectedWellLong=selectedWell;
    if(length(selectedWell) == 2)
              selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    end
    NucleusM = csvHandler.CellPosMatrix(selectedWell);    
    if(size(NucleusM,2) > sizeX || size(NucleusM,1) > sizeY)
        NucleusM = NucleusM(1:sizeY,1:sizeX);
        csvHandler.CellPosMatrix(selectedWell) = NucleusM;
        handles.CSVCoordinates = csvHandler;
        guidata(handles.figure1, handles);
    end
    if(maxCores > 1)
        %Check if watershed images are available. If not create them before
        %processing as ImageJ operation can't be done in parallel.
        %Check if there is already a watershed image. If not, create it.
        imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
        if(~exist(imagePathWatershed,'file'))
            WatershedNuclei(handles);
        end
        createTask(job, @ExecuteCompositeFillAndSkeletonOligo, 4, {selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler});
    else
       [oligoHandler, TP, FP, TN] = ExecuteCompositeFillAndSkeletonOligo(selectedWell,sizeY,sizeX,optionHandler,oligoHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler);
    end
end
if(maxCores > 1)
    submit(job);
    wait(job);
    %job = out(2);
    out=fetchOutputs(job);
    %ToDo: Merge NeuronHandler
    TP = sum(cell2mat(out(:,2)));
    FP = sum(cell2mat(out(:,3)));
    TN = sum(cell2mat(out(:,4)));
end
% TPPerWell = cell2mat(out(:,2));
% FPPerWell = cell2mat(out(:,3));
% TNPerWell = cell2mat(out(:,4));
% TPPerWell = TPPerWell ./ (TPPerWell+FPPerWell);
% FPPerWell = FPPerWell ./ (TPPerWell+FPPerWell);
i=0;

if(maxCores > 1)
    for(j=optionHandler.StartWell:numel(wellList))
     %Get current NeuronHandler result and write it in main Neuronhandler
        i=i+1;
        selectedWell=wellList(j);
        selectedWell=selectedWell{1};    
    
        currentNeuronHandler = out(i,1);
        currentNeuronHandler=currentNeuronHandler{1};
        %Merge Neuron positions as well as neurite length matrix.
        neuronHandler.NeuronPositionsSkeletonization(selectedWell) = currentNeuronHandler.NeuronPositionsSkeletonization(selectedWell);
        neuronHandler.AreaDictionary(selectedWell) = currentNeuronHandler.AreaDictionary(selectedWell);
        neuronHandler.NeuronPositionsEdgeFill(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFill(selectedWell);
        neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
        neuronHandler.NeuronStatMatrix(selectedWell) = currentNeuronHandler.NeuronStatMatrix(selectedWell);
        neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = currentNeuronHandler.NeuronPositionsSkelDeleted(selectedWell);
        %neuronHandler.NeuriteLengthMatrix(selectedWell) = currentNeuronHandler.NeuriteLengthMatrix;
    end
end
close(h);
if(strcmp(channelName,'Oligo'))
     handles.OligoCoordinates = neuronHandler;
else
    handles.NeuronCoordinates = neuronHandler;
end
guidata(handles.figure1, handles);


 function [oligoHandler, TP, FP, TN] = ExecuteCompositeFillAndSkeletonOligo(selectedWell,sizeY,sizeX,optionHandler,oligoHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler)    
        %Get first filter from selected Well
        firstY=1;
        firstX=1;
        lastY=sizeY;
        lastX=sizeX;      
         if(optionHandler.FilterMulPa==1)
             path = [foldername '/ConvertedCellomics/' selectedWell '.mat'];
             TP=0;
             FP=0;
             TN=0;
             if(exist(path,'file')) 
                 load(path);
                 %For every positive Filter
                 filterCount = numel(FMask.PositiveFilters);
                 for(i=1:filterCount);
                    f1 = FMask.PositiveFilters{i};
                    firstInd=find(f1,1, 'first');
                    lastInd=find(f1,1, 'last');
                    [firstY firstX] = ind2sub([sizeY sizeX], firstInd);
                    [lastY lastX] = ind2sub([sizeY sizeX], lastInd);
                    [oligoHandler TP FP TN] = EdgeFillNeuronsAlgorithmOligo(selectedWell,firstY,lastY, firstX, lastX,0,optionHandler,oligoHandler,NucleusM,foldername,channelName);                     
                end
             end
         else
             if(lastY - firstY > 10000 || lastX - firstX > 10000)        
                sizeXsmall = (lastX-firstX)/2;
                sizeYsmall = (lastY-firstY)/2;
                [oligoHandler TP FP TN] = EdgeFillNeuronsAlgorithmOligo(selectedWell,firstY,sizeYsmall,firstX,sizeXsmall,1,optionHandler,oligoHandler,NucleusM,foldername,channelName);
                [oligoHandler TP FP TN] = EdgeFillNeuronsAlgorithmOligo(selectedWell,firstY,sizeYsmall,sizeXsmall+1,lastX,1,optionHandler,oligoHandler,NucleusM,foldername,channelName);
                [oligoHandler TP FP TN] = EdgeFillNeuronsAlgorithmOligo(selectedWell,sizeYsmall+1,lastY,sizeXsmall+1,lastX,1,optionHandler,oligoHandler,NucleusM,foldername,channelName);
                [oligoHandler TP FP TN] = EdgeFillNeuronsAlgorithmOligo(selectedWell,sizeYsmall+1,lastY,firstX,sizeXsmall,1,optionHandler,oligoHandler,NucleusM,foldername,channelName);
              else
                [oligoHandler TP FP TN] = EdgeFillNeuronsAlgorithmOligo(selectedWell,firstY,lastY, firstX, lastX,1,optionHandler,oligoHandler,NucleusM,foldername,channelName);  
              end
              
         end    
         
        % system('taskkill /fi "WINDOWTITLE eq FBird"');
        %Calculate only for filter 
        
function [oligoHandler TP FP TN] = EdgeFillNeuronsAlgorithmOligo(selectedWell, yStartPos, yEndPos, xStartPos, xEndPos, save, optionHandler, oligoHandler, NucleusM, foldername, channelName)
[sizeYorig sizeXorig] = size(NucleusM);
NucleusM = NucleusM(yStartPos:yEndPos,xStartPos:xEndPos);

%[D IDX] = bwdist(full(NucleusM));
TN=0;
if(isfield(oligoHandler.ManualNeuronPositionsSparse,selectedWell))
    NeuronManualM = oligoHandler.ManualNeuronPositionsSparse(selectedWell);
    NeuronManualM = NeuronManualM(yStartPos:yEndPos,xStartPos:xEndPos);
    TN = nnz(NeuronManualM);
else
    NeuronManualM = sparse(logical(zeros(sizeYorig,sizeXorig)));
    if(isa(oligoHandler.ManualNeuronPositionsSparse,'containers.Map') == 0)
        oligoHandler.ManualNeuronPositionsSparse = containers.Map();
    end
    oligoHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
end
%Load neurite and nucleus image from file 
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end
imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
imagePathNeuriteBig = [foldername '/ConvertedCellomics/' selectedWellLong 'OligoBig' optionHandler.FileEnding];

%Check if NeuronStatMatrix exists
if(~isa(oligoHandler.NeuronStatMatrix,'containers.Map'))
    oligoHandler.NeuronStatMatrix = containers.Map();
end
if(~isa(oligoHandler.OmniNucleusPositions,'containers.Map'))
    oligoHandler.OmniNucleusPositions = containers.Map();
end
currentNeuronStat = containers.Map();
neuritePic = imread(imagePathNeuriteBig);
nucleusImage = imread(imagePathNucleusBig);
sizeY = yEndPos - yStartPos+1;
sizeX = xEndPos - xStartPos+1;
[sizeYOrig sizeXOrig]= size(neuritePic);

%Threshold Neurite Picture the known way
    
[sizeY sizeX] = size(neuritePic);
if(optionHandler.FilterMulPa~=1 && optionHandler.CutOutSphereCore == 1)
    neuritePic = CutOutCircles(neuritePic,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM, 0);
end
[neuritePic neuritePicStrong lMin] = ThresholdPic(neuritePic,optionHandler, oligoHandler.StretchlimResult, sizeY, sizeX,0,channelName);

SE = strel('disk', 3);
neuritePic = imdilate(neuritePic,SE);
neuritePic = neuritePic(yStartPos:yEndPos,xStartPos:xEndPos);
neuritePicStrong = neuritePicStrong(yStartPos:yEndPos,xStartPos:xEndPos);
nucleusPic=nucleusImage;
if(optionHandler.FilterMulPa~=1 && optionHandler.CutOutSphereCore == 1)
    nucleusPic = CutOutCircles(nucleusPic,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM,0, 0);
end
nucleusPic(nucleusPic<optionHandler.NucleusThreshold) = 0;
nucleusPic(nucleusPic>=optionHandler.NucleusThreshold) = 255;
nucleusPic=imcomplement(nucleusPic);
nucleusPic = uint8(nucleusPic);

%Do Watershed
%Check if there is already a watershed image. If not, create it.
imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
if(exist(imagePathWatershed,'file'))
    nucleusPic = imread(imagePathWatershed);
else
    javaaddpath('lib/ij.jar');
    javaaddpath('lib/mij.jar');
    MIJ.start(false);
    MIJ.createImage(nucleusPic);
    MIJ.run('Options...', 'iterations=1 count=1 edm=Overwrite do=Nothing');
    MIJ.run('Watershed');
    %nucleusPicBin = logical(MIJ.getCurrentImage);
    MIJ.run('Save', 'save=temp.tif');
    MIJ.closeAllWindows;
    MIJ.exit;
    nucleusPic = imread('temp.tif');
    delete temp.tif;
    nucleusPic = logical(1-nucleusPic);
    imwrite(nucleusPic,imagePathWatershed);
end
nucleusPic = nucleusPic(yStartPos:yEndPos,xStartPos:xEndPos);
[sizeYsmall sizeXsmall] = size(nucleusPic);
fusedPic = imfuse(nucleusPic, neuritePic);
%figure(1);
%imshow(fusedPic);
neuronPointIndices = fusedPic(:,:,1) > 250 & fusedPic(:,:,2) > 250 & fusedPic(:,:,3) > 250;
[sizeYneuronPoints sizeXneuronPoints] = size(neuronPointIndices);
%For every Nucleus:
%Check if Neurite with edges around.
%If yes: Count overlap of Nucleus and Neurites
%Count how many points of Nucleus are within the Neurite edges
EdgeFillNeuronPositionList = sparse(double(sizeY),double(sizeX));
EdgeFillNeuronPositionListSecond = sparse(double(sizeY),double(sizeX));
EdgeFillNeuronPositionList = EdgeFillNeuronPositionList(yStartPos:yEndPos,xStartPos:xEndPos);
EdgeFillNeuronPositionListSecond = EdgeFillNeuronPositionListSecond(yStartPos:yEndPos,xStartPos:xEndPos);
%stringCSVStatisticExport = sprintf('Position X;Position Y;Overlap (percent)\r\n');
[NonZeroY, NonZeroX] = find(NucleusM);
NonZeroY = uint16(NonZeroY);
NonZeroX = uint16(NonZeroX);

dbgYesCount = 0;
dbgNoCount = 0;
dbgNoStartFound = 0;
TP = 0;
FP = 0;
meanSize=zeros(0);
meanBrightness=zeros(0);
count=0;
nucleusPicBin = logical(nucleusPic);
% Get all connected components and label them
LabelComps=bwlabel(neuronPointIndices,8);
LabelCompsNeu = bwlabel(neuritePic);
%Check for each connected component, if also a place on NucleusPic is given
LabelCompsNuc = bwlabel(nucleusPicBin);
for i=1:numel(NonZeroY)    
    if(mod(i,10) == 0)
       disp([num2str(i) ' of ' num2str(numel(NonZeroY))]);
    end
    startX=0;
    startY=0;
    maxDistWhiteNuc = optionHandler.EdgeFillLookAround;
    maxDistWhiteNucIndex=-9;
    maxDistWhiteNucIndexSaved=-9;
    for j=-maxDistWhiteNuc:maxDistWhiteNuc
        for k=-maxDistWhiteNuc:maxDistWhiteNuc
            maxDistWhiteNucIndex = abs(abs(j)+abs(k));
            currentX = NonZeroX(i) + k;
            currentY = NonZeroY(i) + j;
            %If point or one neighbour is white in Overlaypic.            
            if(currentY < sizeY && currentX < sizeX && currentY>0 && currentX>0 && currentX<sizeXneuronPoints && currentY < sizeYneuronPoints && neuronPointIndices(currentY,currentX) > 0 && (startX ==0 || maxDistWhiteNucIndex < maxDistWhiteNucIndexSaved))
                startX= currentX;
                startY = currentY;                
            end
        end
    end
    if(startX > 0)
        %FloodFill white Area from current point. Check also nucleus points
        %stack = java.util.Stack();
        %Idea to be more efficient:
        %Cut out NucleusArea from startX - 50 until startX + 50
        %Do FloodFill Operation only on this subset
        if(startY-150>0 && startX-150>0 && startX+150 < sizeXsmall && startY + 150 < sizeYsmall)
            SubArea = neuronPointIndices(startY-150:startY+150,startX-150:startX+150);
            SubNucleusM = NucleusM(startY-150:startY+150,startX-150:startX+150);
            SubNucleusPic = nucleusPic(startY-150:startY+150,startX-150:startX+150);
            SubNucleusImage = nucleusImage(startY-150:startY+150,startX-150:startX+150);
            
            SubLabelComps = LabelComps(startY-150:startY+150,startX-150:startX+150);
            SubLabelCompsNuc = LabelCompsNuc(startY-150:startY+150,startX-150:startX+150);
            SubLabelCompsNeu = LabelCompsNeu(startY-150:startY+150,startX-150:startX+150);
            conCompIndex = SubLabelComps(150, 150);
            conCompIndexNuc = SubLabelCompsNuc(150, 150);
            conCompIndexNeu = SubLabelCompsNeu(150, 150);
            
            if(conCompIndex == 0)
               conCompIndex = SubLabelComps(151, 151); 
               if(conCompIndex == 0)
                   conCompIndex = SubLabelComps(149, 149); 
                   if(conCompIndex == 0)
                     conCompIndex = SubLabelComps(151, 149); 
                        if(conCompIndex == 0)
                            conCompIndex = SubLabelComps(149, 151); 
                        end
                   end
               end
            end
            if(conCompIndexNuc == 0)
               conCompIndexNuc = SubLabelCompsNuc(151, 151); 
               if(conCompIndexNuc == 0)
                   conCompIndexNuc = SubLabelCompsNuc(149, 149); 
                   if(conCompIndexNuc == 0)
                     conCompIndexNuc = SubLabelCompsNuc(151, 149); 
                        if(conCompIndexNuc == 0)
                            conCompIndexNuc = SubLabelCompsNuc(149, 151); 
                        end
                   end
               end
            end
            if(conCompIndexNeu == 0)
               conCompIndexNeu = SubLabelCompsNeu(151, 151); 
               if(conCompIndexNeu == 0)
                   conCompIndexNeu = SubLabelCompsNeu(149, 149); 
                   if(conCompIndexNeu == 0)
                     conCompIndexNeu = SubLabelCompsNeu(151, 149); 
                        if(conCompIndexNeu == 0)
                            conCompIndexNeu = SubLabelCompsNeu(149, 151); 
                        end
                   end
               end
            end
            
            
            %Extract area
            numberNeuritePixels = bwarea(SubArea == 1 & SubLabelCompsNuc == conCompIndexNuc);
            numberNucleusPixels = bwarea(SubLabelCompsNuc == conCompIndexNuc);            
            numberNeuritePixelsTotal = bwarea(conCompIndexNeu == SubLabelCompsNeu);
            
            count=count+1;
            meanSize = [meanSize numberNucleusPixels];
            meanBrightness = [meanBrightness CalculateNucleusBrightness(SubNucleusPic, SubNucleusImage)];
            [maxY, maxX] = size(SubNucleusPic);
            %Recalculate startX and startY            
         %   stack.push([50 50]);
            %FloodFill SubArea as well as Binary NucleusPic
            %Count number of pixels of both during FloodFill operation
      %      [numberNeuritePixels, numberNucleusPixels] = FloodFillForEdgeFill(SubNucleusPic,SubArea,maxX,maxY,stack);
            %WhiteArea2 = imfill(SubArea,[50 50]);         
            %[yIndices xIndices] = find(WhiteArea2);
            %yIndices = yIndices + (startY-50);
            %xIndices = xIndices + (startX-50);
        else
            %Extract connected component at given point
            conCompIndex = LabelComps(startY, startX);
            conCompIndexNuc = LabelCompsNuc(startY, startX);
            conCompIndexNeu = LabelCompsNeu(startY, startX);
            %Extract area
            numberNeuritePixels = bwarea(LabelCompsNuc == conCompIndexNuc & neuronPointIndices == 1);
            numberNucleusPixels = bwarea(LabelCompsNuc == conCompIndexNuc);
            numberNeuritePixelsTotal = bwarea(LabelCompsNeu == conCompIndexNeu);
%            stack.push([startX startY]);
            WhiteArea = zeros(sizeY,sizeX);
    %        [numberNeuritePixels, numberNucleusPixels] = FloodFillForEdgeFill(nucleusPic,neuronPointIndices,sizeX,sizeY,stack);
             %Schwerpunkt aller Punkte berechnen.
            %[yIndices xIndices] = find(WhiteArea);
        end     
        %Check if Area has not too less and not too much pixels

        %Check relation between numberNeuritePixels and numberNucleusPixels
        
        %Add neurite independent of it's neuritePixels per
        %NucleusPixels. Save that value.
        %Data structure: Matrix with positions AND relation
        
        %Mean calculation of brightness and NucleusArea
        
        area = optionHandler.EdgeFillNucleusAreaWithinNeurite;
        areaSecond = optionHandler.EdgeFillNucleusAreaWithinNeuriteSecond;
        areaAround = optionHandler.CompositeFillNeuriteAreaWithinNucleus;
        areaAroundMax = optionHandler.CompositeFillNeuriteAreaWithinNucleusMax;
        edgeFillMarked = 0;
        edgeFillMarkedSecond = 0;
        
        %Find next position of NucleusM from current NucleusOM
        %[curY curX] = ind2sub([sizeY,sizeX],IDX(NonZeroY(i),NonZeroX(i)));
        
        if(numberNeuritePixels/numberNucleusPixels >= area && numberNucleusPixels/numberNeuritePixelsTotal >= areaAround && numberNucleusPixels/numberNeuritePixelsTotal <= areaAroundMax && numberNeuritePixelsTotal > optionHandler.MinSizeNeuriteAreaCompFill) 
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
        
        if(NeuronManualM(uint16(yStartPos)-1+uint16(NonZeroY(i)),uint16(xStartPos)-1+uint16(NonZeroX(i))))
            manualMarked=1;            
            if(edgeFillMarked)
                TP = TP + 1;
                TN = TN - 1;
            end
        else
            manualMarked=0;
            if(edgeFillMarked)
                FP = FP + 1;
            end
        end
        if(save)
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
            eFillFirst = 0;
            eFillSecond=0;
        else
            eFillFirst = EdgeFillNeuronPositionList;
            eFillSecond=EdgeFillNeuronPositionListSecond;
        end
    else
        dbgNoStartFound = dbgNoStartFound +1;
    end    
end

meanSize = median(meanSize);
meanBrightness = median(meanBrightness);

if(numel(oligoHandler.NucleusBrightnessDict) ==0)% || ~isKey(oligoHandler.NucleusBrightnessDict, selectedWell)
    oligoHandler.NucleusBrightnessDict = containers.Map();
    oligoHandler.NucleusAreaDict = containers.Map();
    oligoHandler.NucleusBrightnessDict(selectedWell) = meanBrightness;
    oligoHandler.NucleusAreaDict(selectedWell) = meanSize;
else
    oligoHandler.NucleusAreaDict(selectedWell) = meanBrightness;
    oligoHandler.NucleusAreaDict(selectedWell) = meanSize;
end
if(numel(oligoHandler.NeuronPositionsEdgeFill) ==0)
    oligoHandler.NeuronPositionsEdgeFill = containers.Map();
end
if(~isKey(oligoHandler.NeuronPositionsEdgeFill, selectedWell)) 
    oligoHandler.NeuronPositionsEdgeFill(selectedWell) = sparse(logical(zeros(sizeYOrig,sizeXOrig)));
    fullMatrix = oligoHandler.NeuronPositionsEdgeFill(selectedWell);
    fullMatrix(yStartPos:yEndPos,xStartPos:xEndPos) = EdgeFillNeuronPositionList;
    oligoHandler.NeuronPositionsEdgeFill(selectedWell) = fullMatrix;
    %oligoHandler.NeuronPositionsEdgeFill(selectedWell) = EdgeFillNeuronPositionList;
else
    fullMatrix = oligoHandler.NeuronPositionsEdgeFill(selectedWell);
    fullMatrix(yStartPos:yEndPos,xStartPos:xEndPos) = EdgeFillNeuronPositionList;
    oligoHandler.NeuronPositionsEdgeFill(selectedWell) = fullMatrix;
end
if(numel(oligoHandler.NeuronPositionsEdgeFillSecond) ==0)
    oligoHandler.NeuronPositionsEdgeFillSecond = containers.Map();
end
if(~isKey(oligoHandler.NeuronPositionsEdgeFillSecond, selectedWell))
    oligoHandler.NeuronPositionsEdgeFillSecond(selectedWell) = sparse(logical(zeros(sizeYOrig,sizeXOrig)));
    fullMatrix = oligoHandler.NeuronPositionsEdgeFillSecond(selectedWell);
    fullMatrix(yStartPos:yEndPos,xStartPos:xEndPos) = EdgeFillNeuronPositionListSecond;
    oligoHandler.NeuronPositionsEdgeFillSecond(selectedWell) = fullMatrix;
    %oligoHandler.NeuronPositionsEdgeFillSecond(selectedWell) = EdgeFillNeuronPositionListSecond;
else
    fullMatrix = oligoHandler.NeuronPositionsEdgeFillSecond(selectedWell);
    fullMatrix(yStartPos:yEndPos,xStartPos:xEndPos) = EdgeFillNeuronPositionListSecond;
    oligoHandler.NeuronPositionsEdgeFillSecond(selectedWell) = fullMatrix;
end
oligoHandler.NeuronStatMatrix(selectedWell) = currentNeuronStat;
%handles.NeuronCoordinates = oligoHandler;
%handles.ImageHandler = imageHandler;
%guidata(handles.figure1, handles);    
%Check if oligoHandler already has Dictionary of NeuronPositionsEdgeComposite        
        


% --------------------------------------------------------------------
function CorticalNucMorphology_Callback(hObject, eventdata, handles)
% hObject    handle to CorticalNucMorphology (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
% Get variables from opened GUI    
handles = guidata(handles.figure1);
% Get variable of Options from set seetings in GUI
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
foldername = imageHandler.Foldername;
ExecuteOptionExperimentCorticalNuc(handles,1,'Neurite');
guidata(handles.figure1, handles);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
[summary] = Measurement (foldername, csvHandler, neuronHandler)
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function ExecuteOptionExperimentCorticalNuc(handles,onlyCompFill,channelName)
h = waitbar(0.5,'Neuronal morpholgy analyzed. This could take a few hours. Please be patient.', 'WindowStyle', 'modal');
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
optionHandler = handles.OptionHandler;

if(strcmp(channelName,'Oligo'))
    competingHandler = handles.NeuronCoordinates;
    if(isfield(handles,'OligoCoordinates'))
       neuronHandler = handles.OligoCoordinates;
    else
       neuronHandler = CSVCoordinates();
    end
else
    neuronHandler = handles.NeuronCoordinates;
    competingHandler = 0;
end
wellList = get(handles.lbWell, 'string');
foldername = imageHandler.Foldername;
[sizeY sizeX] = size(imageHandler.NucleusImage);
%Select STARTWELL

if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
if(numel(neuronHandler.AreaDictionary) ==0)
    neuronHandler.AreaDictionary = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
    neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFill) ==0)
    neuronHandler.NeuronPositionsEdgeFill = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFillSecond) ==0)
    neuronHandler.NeuronPositionsEdgeFillSecond = containers.Map();
end
if(numel(neuronHandler.NeuronStatMatrix) ==0)
    neuronHandler.NeuronStatMatrix = containers.Map();
end

%Create Matlab job
%Leave one Core out to ensure that GUI isn't freezing
maxCores =optionHandler.NumberOfCores;
if(maxCores > 1)
    c = parcluster();
    c.NumWorkers = maxCores;
    job = createJob(c,'JobData',optionHandler);
end

for(j=optionHandler.StartWell:numel(wellList))
    selectedWell=wellList(j);
    selectedWell=selectedWell{1};
    selectedWellLong=selectedWell;
    if(length(selectedWell) == 2)
              selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    end
    
    NucleusM = csvHandler.CellPosMatrix(selectedWell);    
    if(size(NucleusM,2) > sizeX || size(NucleusM,1) > sizeY)
        NucleusM = NucleusM(1:sizeY,1:sizeX);
        csvHandler.CellPosMatrix(selectedWell) = NucleusM;
        handles.CSVCoordinates = csvHandler;
        guidata(handles.figure1, handles);
    end
    if(maxCores > 1)
        %Check if watershed images are available. If not create them before
        %processing as ImageJ operation can't be done in parallel.
        %Check if there is already a watershed image. If not, create it.
        imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
        if(~exist(imagePathWatershed,'file'))
            WatershedNuclei(handles);
        end
        createTask(job, @ExecuteCorticalNuc, 4, {selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler});
    else
       [neuronHandler,  summary] = ExecuteCorticalNuc(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler);
    end
end
if(maxCores > 1)
    submit(job);
    wait(job);
    %job = out(2);
    out=fetchOutputs(job);
    %ToDo: Merge NeuronHandler
    TP = sum(cell2mat(out(:,2)));
    FP = sum(cell2mat(out(:,3)));
    TN = sum(cell2mat(out(:,4)));
end
% TPPerWell = cell2mat(out(:,2));
% FPPerWell = cell2mat(out(:,3));
% TNPerWell = cell2mat(out(:,4));
% TPPerWell = TPPerWell ./ (TPPerWell+FPPerWell);
% FPPerWell = FPPerWell ./ (TPPerWell+FPPerWell);
i=0;

if(maxCores > 1)
    for(j=optionHandler.StartWell:numel(wellList))
     %Get current NeuronHandler result and write it in main Neuronhandler
        i=i+1;
        selectedWell=wellList(j);
        selectedWell=selectedWell{1};    
    
        currentNeuronHandler = out(i,1);
        currentNeuronHandler=currentNeuronHandler{1};
        %Merge Neuron positions as well as neurite length matrix.
        neuronHandler.NeuronPositionsSkeletonization(selectedWell) = currentNeuronHandler.NeuronPositionsSkeletonization(selectedWell);
        neuronHandler.AreaDictionary(selectedWell) = currentNeuronHandler.AreaDictionary(selectedWell);
        neuronHandler.NeuronPositionsEdgeFill(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFill(selectedWell);
        neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
        neuronHandler.NeuronStatMatrix(selectedWell) = currentNeuronHandler.NeuronStatMatrix(selectedWell);
        neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = currentNeuronHandler.NeuronPositionsSkelDeleted(selectedWell);
        %neuronHandler.NeuriteLengthMatrix(selectedWell) = currentNeuronHandler.NeuriteLengthMatrix;
    end
end
close(h);
if(strcmp(channelName,'Oligo'))
     handles.OligoCoordinates = neuronHandler;
else
    handles.NeuronCoordinates = neuronHandler;
end
guidata(handles.figure1, handles);


 function [neuronHandler, summary] = ExecuteCorticalNuc(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler)    
        %Get first filter from selected Well
        selectedWellLong=selectedWell;
        if(length(selectedWell) == 2)
                  selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
        end
        imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
        nucleusImage = imread(imagePathNucleusBig);
        [sizeYorig sizeXorig] = size(nucleusImage);
        firstY=1;
        firstX=1;
        lastY=sizeYorig;
        lastX=sizeXorig;
        sizeY = sizeYorig;
        sizeX = sizeXorig;
         if(optionHandler.FilterMulPa==1)
             path = [foldername '/ConvertedCellomics/' selectedWell '.mat'];
             TP=0;
             FP=0;
             TN=0;
             if(exist(path,'file')) 
                 load(path);
                 %For every positive Filter
                 filterCount = numel(FMask.PositiveFilters);
                 for(i=1:filterCount);
                    f1 = FMask.PositiveFilters{i};
                    firstInd=find(f1,1, 'first');
                    lastInd=find(f1,1, 'last');
                    [firstY firstX] = ind2sub([sizeY sizeX], firstInd);
                    [lastY lastX] = ind2sub([sizeY sizeX], lastInd);
                    [neuronHandler summary] = MorphologyCorticalNuc(selectedWell,firstY,lastY, firstX, lastX,0,optionHandler,neuronHandler,NucleusM,foldername,channelName);                     
                end
             end
         else
             if(lastY - firstY > 10000 || lastX - firstX > 10000)        
                sizeXsmall = (lastX-firstX)/2;
                sizeYsmall = (lastY-firstY)/2;
                [neuronHandler summary] = MorphologyCorticalNuc(selectedWell,firstY,sizeYsmall,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyCorticalNuc(selectedWell,firstY,sizeYsmall,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyCorticalNuc(selectedWell,sizeYsmall+1,lastY,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyCorticalNuc(selectedWell,sizeYsmall+1,lastY,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
              else
                [neuronHandler summary] = MorphologyCorticalNuc(selectedWell,firstY,lastY, firstX, lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);  
              end
              
         end    
         
        % system('taskkill /fi "WINDOWTITLE eq FBird"');
        %Calculate only for filter 
        
function [neuronHandler summary] = MorphologyCorticalNuc(selectedWell, yStartPos, yEndPos, xStartPos, xEndPos, save, optionHandler, neuronHandler, NucleusM, foldername, channelName)
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end
imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
nucleusImage = imread(imagePathNucleusBig);

[sizeYorig sizeXorig] = size(nucleusImage);
%yStartPos = 1;
%yEndPos= sizeYorig;
%xStartPos = 1;
%xEndPos = sizeXorig;
NucleusM = NucleusM(yStartPos:yEndPos,xStartPos:xEndPos);

%[D IDX] = bwdist(full(NucleusM));
TN=0;
if(isfield(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
    NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
    NeuronManualM = NeuronManualM(yStartPos:yEndPos,xStartPos:xEndPos);
    TN = nnz(NeuronManualM);
else
    NeuronManualM = sparse(logical(zeros(sizeYorig,sizeXorig)));
    if(isa(neuronHandler.ManualNeuronPositionsSparse,'containers.Map') == 0)
        neuronHandler.ManualNeuronPositionsSparse = containers.Map();
    end
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
end
%Load neurite and nucleus image from file 
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end

% From here onwards starts the actual algorithm!
%These are for now set fixed parameters, but can later be included into
%Option menue
minOverlapp = 0.3;
A = 15000;
B = 4;
N = 50;
C = 250;
D=10;
E = 15000;
NeuronBody=200;
NeuriteThreshold = 0.7;

NeuronNucleiArea = 500;
SE= strel('disk',8);
RatioFill = 0.6;
Maxlengthsquare = 100;

AxonMeasurement=0;

% Image locations:
imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
imagePathBinary = [foldername '/ConvertedCellomics/' selectedWellLong 'Binary' optionHandler.FileEnding];
if(AxonMeasurement==1);
imagePathBinaryHard = [foldername '/ConvertedCellomics/' selectedWellLong 'BinaryHard' optionHandler.FileEnding];
end
%Check if NeuronStatMatrix exists
if(~isa(neuronHandler.NeuronStatMatrix,'containers.Map'))
    neuronHandler.NeuronStatMatrix = containers.Map();
end
if(~isa(neuronHandler.OmniNucleusPositions,'containers.Map'))
    neuronHandler.OmniNucleusPositions = containers.Map();
end
currentNeuronStat = containers.Map();
binaryimgNeurite = imread(imagePathBinary);
if(AxonMeasurement==1);
binaryimgNeuriteHard = imread(imagePathBinaryHard);
end
nucleusPicBinary = imread(imagePathWatershed);

sizeY = yEndPos - yStartPos+1;
sizeX = xEndPos - xStartPos+1;
[sizeYOrig sizeXOrig]= size(binaryimgNeurite);

%Threshold Neurite Picture the known way
    
[sizeY sizeX] = size(binaryimgNeurite);
if(optionHandler.FilterMulPa~=1 && optionHandler.CutOutSphereCore == 1)
    neuritePic = CutOutCircles(binaryimgNeurite,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM, 0);
end

if(optionHandler.FilterMulPa~=1 && optionHandler.CutOutSphereCore == 1)
    nucleusPic = CutOutCircles(nucleusPicBinary,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM,0, 0);
end

nucleusPicBinary = nucleusPicBinary(yStartPos:yEndPos,xStartPos:xEndPos);
[sizeYsmall sizeXsmall] = size(nucleusPicBinary);

binaryimgNeurite = binaryimgNeurite(yStartPos:yEndPos,xStartPos:xEndPos);
if(AxonMeasurement==1);
binaryimgNeuriteHard = binaryimgNeuriteHard(yStartPos:yEndPos,xStartPos:xEndPos);
end
binaryimgNeurite=logical(binaryimgNeurite);
%binaryimgNeurite=uint8(binaryimgNeurite);  
%Algorithm to smooth edges:
%1) Find edges:
%E = edge(binaryimgNeurite,'canny');
%E = logical(E);
%2) Add edge image to normal image
%binaryimgNeurite = binaryimgNeurite + E;
%3) This image will have holes which have to be filled. Since crossing
%Neurite will create 'big' holes, the small 'holes' have to be extracted!
% Fill holes in image
%Holes = imfill(binaryimgNeurite,'holes');
% Get only filled areas
%Holes = Holes - binaryimgNeurite;
% Extract 'big' holes
%HolesBig = bwareaopen(Holes, 50);
% Get small holes, by supstraction the big holes from all holes
%Holes = Holes - HolesBig;
% Fill small holes:
%binaryimgNeurite= binaryimgNeurite + Holes;
%binaryimgNeurite = bwmorph(binaryimgNeurite,'majority');
%binaryimgNeurite = bwmorph(binaryimgNeurite,'majority');
SkeletonImage = bwmorph(binaryimgNeurite,'thin',Inf);
SkeletonImage = bwmorph(SkeletonImage,'spur',4);
        %binaryimgNeurite(binaryimgNeurite==0)=255;
        %binaryimgNeurite(binaryimgNeurite==1)=0;
        %javaaddpath('lib/ij.jar');
        %javaaddpath('lib/mij.jar');
        %MIJ.start(false);
        %MIJ.createImage(binaryimgNeurite);
        %MIJ.run('Options...', 'iterations=1 count=1 edm=Overwrite do=Nothing');
        %MIJ.run('Skeletonize');
        %nucleusPicBin = logical(MIJ.getCurrentImage);
        %MIJ.run('Save', 'save=temp.tif');
        %MIJ.closeAllWindows;
        %MIJ.exit;
        %SkeletonImage = imread('temp.tif');
        %delete temp.tif;
        %SkeletonImage = logical(1-SkeletonImage);
        	
        %SkeletonImage = bwareaopen(SkeletonImage,50);
        if(AxonMeasurement==1);
        binaryimgNeuriteHard=uint8(binaryimgNeuriteHard);
        
        binaryimgNeuriteHard(binaryimgNeuriteHard==0)=255;
        binaryimgNeuriteHard(binaryimgNeuriteHard==1)=0;
        javaaddpath('lib/ij.jar');
        javaaddpath('lib/mij.jar');
        MIJ.start(false);
        MIJ.createImage(binaryimgNeuriteHard);
        MIJ.run('Options...', 'iterations=1 count=1 edm=Overwrite do=Nothing');
        MIJ.run('Skeletonize');
        %nucleusPicBin = logical(MIJ.getCurrentImage);
        MIJ.run('Save', 'save=temp.tif');
        MIJ.closeAllWindows;
        MIJ.exit;
        SkeletonImageHard = imread('temp.tif');
        delete temp.tif;
        SkeletonImageHard = logical(1-SkeletonImageHard);
        SkeletonImageHard = bwmorph(SkeletonImageHard,'spur',10);
        SkeletonImageHard = bwareaopen(SkeletonImageHard,50);
        end
        
        
        
        
        
        
        
        
        
        
        
        
        %binaryimgNeurite = logical(1-binaryimgNeurite);
        %labelBinaryNucleus = bwlabel(nucleusPicBinary);
        %Count number of nuclei present
        %k = max(max(labelBinaryNucleus));
        % Function to extract only neuron nuclei. 
        %Create empty image with size equal to original image
        %[r,c] = size(labelBinaryNucleus);
        %neuronNuclei = zeros(r,c);
        %neuronNuclei = uint8(neuronNuclei);
        %labelSingleNucleus = zeros(r,c);
        %count = 0;
        % Now extract each single lable (one separate Nuclei) and check
        % whether the overlapp with the neurite channel is sufficient to be
        % a neuron nuclei. If so keep it, else delete it.
        %for(w=1:k)
        % This loop iterates over the image for a certain lable and extracts it. E.g. Nuclei with lable one is extracted in the first step.        
        %   for(u=1:r)
        %   for(v=1:c)
        %       if(labelBinaryNucleus(u,v)==w)
        %         labelSingleNucleus(u,v)=1;
        %      else
        %         labelSingleNucleus(u,v)=0; 
        %      end   
        %   end
        %   end
        %   % Add Neurite Image to the single nuclei image.
        %   labelSingleNucleus = uint8(labelSingleNucleus);
        %   binaryimgNeurite = uint8(binaryimgNeurite);
        %   CloseConnectedComponent = labelSingleNucleus +  binaryimgNeurite;
           % Count number of Pixels with value 2 (equals the close
           % connected component of nuclei and neurons staining!)
        %   numberPixelOverlay = length(CloseConnectedComponent(CloseConnectedComponent==2));
           % Count number of pixel of nuclei
        %   numberPixelNuclei = length(labelSingleNucleus(labelSingleNucleus==1));
           % Check if overlapp is bigger than minOverlapp:
        %   if(numberPixelOverlay/numberPixelNuclei>minOverlapp)
        %      neuronNuclei = neuronNuclei + labelSingleNucleus;
              %count = count +1;
        %   end
        %   tell = w;
        %end
        %neuronNuclei = logical(neuronNuclei);
        %sn = regionprops(neuronNuclei,'centroid');
        %sndimensions = size(sn);
        % Generating the Neuron matrix:
        %NeuronM = sparse(zeros(size(S)));
        %if(sndimensions > 0)
        %    for(b=1:sndimensions)
        %    xn_centroid = sn(b).Centroid(1);
        %    yn_centroid = sn(b).Centroid(2);
        %    NeuronM(uint16(yn_centroid),uint16(xn_centroid)) = 1;
        %    end
        %end
        %neuronNuclei = uint8(neuronNuclei);
        %newfile = [selectedWellLong 'NeuronNuclei' optionHandler.FileEnding];
        %imwrite(neuronNuclei,[foldername '/ConvertedCellomics/' newfile] );
        %NumberofNeurons = sndimensions(1);
        %Substract Cell soma:
        SkeletonImage = uint8(SkeletonImage);
        if(AxonMeasurement==1);
        SkeletonImageHard = uint8(SkeletonImageHard);
        end
        SkeletonImage = SkeletonImage;
        if(AxonMeasurement==1);
        SkeletonImageHard = SkeletonImageHard;
        end
        %Save SkeletonImage without Nuclei
        Dummy = zeros(sizeYOrig,sizeXOrig);
        Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = SkeletonImage;
        SkeletonImage = Dummy;
        newfile = [selectedWellLong 'Skeleton' optionHandler.FileEnding];
        imwrite(SkeletonImage,[foldername '/ConvertedCellomics/' newfile] );
        if(AxonMeasurement==1);
        Dummy = zeros(sizeYOrig,sizeXOrig);
        Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = SkeletonImageHard;
        SkeletonImageHard = Dummy;
        newfile = [selectedWellLong 'AstroBig' optionHandler.FileEnding];
        imwrite(SkeletonImageHard,[foldername '/ConvertedCellomics/' newfile] );
        newfile = [selectedWellLong 'AstroSmall' optionHandler.FileEnding];
        imwrite(SkeletonImageHard,[foldername '/ConvertedCellomics/' newfile] );
        end
        Dummy = zeros(sizeYOrig,sizeXOrig);
        Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = binaryimgNeurite;
        binaryimgNeurite = Dummy;
        newfile = [selectedWellLong 'Binary' optionHandler.FileEnding];
        imwrite(binaryimgNeurite,[foldername '/ConvertedCellomics/' newfile] );
        %Create Matrix with dimension of maximum images
        NeuronM = NucleusM;
        %Dummy = zeros(sizeYOrig,sizeXOrig);
        %Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = NeuronM;
        %NeuronM = Dummy;
        summary=0;





neuronHandler.NeuronPositionsEdgeFill(selectedWell) = NeuronM;
%handles.NeuronCoordinates = neuronHandler;
%handles.ImageHandler = imageHandler;
%guidata(handles.figure1, handles);    
%Check if Neuronhandler already has Dictionary of NeuronPositionsEdgeComposite  
    %function [Summary] = SummaryOut(n,summary{n,1},summary{n,2},summary{n,3},summary{n,4},summary{n,5} ,summary{n,6})
       
    
    

%This function is used to postprocess binary images from the load cortical
%neuron function
% --------------------------------------------------------------------
function CorticalMorphology_Callback(hObject, eventdata, handles)
% hObject    handle to CorticalMorphology (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
% Defines handles for GUI    
handles = guidata(handles.figure1);
% Sets optionHandler to variable currently set in GUI under option menu
optionHandler = handles.OptionHandler;
% ImageHandler contains recent loaded images
imageHandler = handles.ImageHandler;
%image path is obtained from GUI
foldername = imageHandler.Foldername;
% Executes function: ExecuteOptionExperimentCorticalTrans
ExecuteOptionExperimentCorticalTrans(handles,1,'Neurite');
% Analysis results are returned to the GUI
guidata(handles.figure1, handles);
% csvHandler is obtained from GUI
csvHandler = handles.CSVCoordinates;
% csvHandler is obtained from GUI
neuronHandler = handles.NeuronCoordinates;
% Execution of MeasurementTran
[summary Summary2 Summary3 Summary4 Summary5 ExpName csvHandler] = MeasurementTrans (foldername, csvHandler, neuronHandler)
% This function is used to write summary - Summary 4 in Excel files
Replace = 0;
[r c] = size(Summary2)
for(z=1:(c-1))
    if(Summary2{2,(z+1)} > 0)
        Replace = 1;
    end
    for(p=1:(r-1))
        if(Replace == 1)
        if(isempty(Summary2{(p+1),(z+1)}) == 1)
            Summary2{(p+1),(z+1)} = 0;
        end
        end
    end
    Replace = 0;
end

[r c] = size(Summary3)
for(z=1:(c-1))
    if(Summary3{2,(z+1)} > 0)
        Replace = 1;
    end
    for(p=1:(r-1))
        if(Replace == 1)
        if(isempty(Summary3{(p+1),(z+1)}) == 1)
            Summary3{(p+1),(z+1)} = 0;
        end
        end
    end
    Replace = 0;
end

[r c] = size(Summary4)

for(z=1:(c-1))
    if(Summary4{2,(z+1)} > 0)
        Replace = 1;
    end
    for(p=1:(r-1))
        if(Replace == 1)
        if(isempty(Summary4{(p+1),(z+1)}) == 1)
            Summary4{(p+1),(z+1)} = 0;
        end
        end
    end
    Replace = 0;
end
% Bounces summary - Summary4 back to command window
summary = summary
Summary2 = Summary2
Summary3 = Summary3
Summary4 = Summary4
Summary5 = Summary5
% Filename from MeasurementTrans is used as name for Excell files
filename = regexprep(ExpName,'.png','.xlsx','ignorecase');  
%All summaries are written in the same excel file under different sheeds
xlswrite(filename,summary);
sheet = 2;
xlswrite(filename,Summary2,sheet);
sheet = 3;
xlswrite(filename,Summary3,sheet);
sheet = 4;
xlswrite(filename,Summary4,sheet);
sheet = 5;
xlswrite(filename,Summary5,sheet);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% This function is mainly used to identify whether only neurons or also
% other cell types are evaluated and whether the job is performed on a
% single CPU core or on mutiple cores
function ExecuteOptionExperimentCorticalTrans(handles,onlyCompFill,channelName)
h = waitbar(0.5,'Neuronal morpholgy analyzed. This could take a few hours. Please be patient.', 'WindowStyle', 'modal');
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
optionHandler = handles.OptionHandler;
%Here is checked whether Oligo data should be analyzed
if(strcmp(channelName,'Oligo'))
    competingHandler = handles.NeuronCoordinates;
    if(isfield(handles,'OligoCoordinates'))
       neuronHandler = handles.OligoCoordinates;
    else
       neuronHandler = CSVCoordinates();
    end
else
    neuronHandler = handles.NeuronCoordinates;
    competingHandler = 0;
end
wellList = get(handles.lbWell, 'string');
foldername = imageHandler.Foldername;
[sizeY sizeX] = size(imageHandler.NucleusImage);
%Select STARTWELL

if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
if(numel(neuronHandler.AreaDictionary) ==0)
    neuronHandler.AreaDictionary = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
    neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFill) ==0)
    neuronHandler.NeuronPositionsEdgeFill = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFillSecond) ==0)
    neuronHandler.NeuronPositionsEdgeFillSecond = containers.Map();
end
if(numel(neuronHandler.NeuronStatMatrix) ==0)
    neuronHandler.NeuronStatMatrix = containers.Map();
end

%Create Matlab job
%Leave one Core out to ensure that GUI isn't freezing
maxCores =optionHandler.NumberOfCores;
if(maxCores > 1)
    c = parcluster();
    c.NumWorkers = maxCores;
    job = createJob(c,'JobData',optionHandler);
end

for(j=optionHandler.StartWell:numel(wellList))
    selectedWell=wellList(j);
    selectedWell=selectedWell{1};
    selectedWellLong=selectedWell;
    if(length(selectedWell) == 2)
              selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    end
    NucleusM = csvHandler.CellPosMatrix(selectedWell); 
    if(size(csvHandler.ManualPositions1)>0)
        try
    NeuronM = csvHandler.ManualPositions1(selectedWell);
        catch
    NeuronM = 0;
        end
    else
    NeuronM = 0;  
    end
    if(size(NucleusM,2) > sizeX || size(NucleusM,1) > sizeY)
        NucleusM = NucleusM(1:sizeY,1:sizeX);
        csvHandler.CellPosMatrix(selectedWell) = NucleusM;
        handles.CSVCoordinates = csvHandler;
        guidata(handles.figure1, handles);
    end
    if(maxCores > 1)
        %Check if watershed images are available. If not create them before
        %processing as ImageJ operation can't be done in parallel.
        %Check if there is already a watershed image. If not, create it.
        imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
        if(~exist(imagePathWatershed,'file'))
            WatershedNuclei(handles);
        end
        createTask(job, @ExecuteCorticalTrans, 4, {selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler});
    else
        % This is usally started on a normal PC and obtains all relevant
        % data. For filters the dimensions are given to the next function
        % as sizeY,sizeX
       [NeuronM, neuronHandler,  summary] = ExecuteCorticalTrans(NeuronM, selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler);
    end
end
if(maxCores > 1)
    submit(job);
    wait(job);
    %job = out(2);
    out=fetchOutputs(job);
    %ToDo: Merge NeuronHandler
    TP = sum(cell2mat(out(:,2)));
    FP = sum(cell2mat(out(:,3)));
    TN = sum(cell2mat(out(:,4)));
end
% TPPerWell = cell2mat(out(:,2));
% FPPerWell = cell2mat(out(:,3));
% TNPerWell = cell2mat(out(:,4));
% TPPerWell = TPPerWell ./ (TPPerWell+FPPerWell);
% FPPerWell = FPPerWell ./ (TPPerWell+FPPerWell);
i=0;

if(maxCores > 1)
    for(j=optionHandler.StartWell:numel(wellList))
     %Get current NeuronHandler result and write it in main Neuronhandler
        i=i+1;
        selectedWell=wellList(j);
        selectedWell=selectedWell{1};    
    
        currentNeuronHandler = out(i,1);
        currentNeuronHandler=currentNeuronHandler{1};
        %Merge Neuron positions as well as neurite length matrix.
        neuronHandler.NeuronPositionsSkeletonization(selectedWell) = currentNeuronHandler.NeuronPositionsSkeletonization(selectedWell);
        neuronHandler.AreaDictionary(selectedWell) = currentNeuronHandler.AreaDictionary(selectedWell);
        neuronHandler.NeuronPositionsEdgeFill(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFill(selectedWell);
        neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
        neuronHandler.NeuronStatMatrix(selectedWell) = currentNeuronHandler.NeuronStatMatrix(selectedWell);
        neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = currentNeuronHandler.NeuronPositionsSkelDeleted(selectedWell);
        %neuronHandler.NeuriteLengthMatrix(selectedWell) = currentNeuronHandler.NeuriteLengthMatrix;
    end
end
close(h);
if(strcmp(channelName,'Oligo'))
     handles.OligoCoordinates = neuronHandler;
else
    handles.NeuronCoordinates = neuronHandler;
end
guidata(handles.figure1, handles);

% This function is used to determine whether the whole umage is evaluated or only a filter 
 function [NeuronM, neuronHandler, summary] = ExecuteCorticalTrans(NeuronM, selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler)    
        %Get first filter from selected Well
        selectedWellLong=selectedWell;
        if(length(selectedWell) == 2)
                  selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
        end
        imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
        nucleusImage = imread(imagePathNucleusBig);
        % Gives the original dimension of the image and defines x,y
        % coordinates of image edges
        [sizeYorig sizeXorig] = size(nucleusImage);
        firstY=1;
        firstX=1;
        lastY=sizeYorig;
        lastX=sizeXorig;
        sizeY = sizeYorig;
        sizeX = sizeXorig;
        % Here it is checked whether filter analysis is chosen in the GUI
        % or not
         if(optionHandler.FilterMulPa==1)
             path = [foldername '/ConvertedCellomics/' selectedWell '.mat'];
             TP=0;
             FP=0;
             TN=0;
             if(exist(path,'file')) 
                 load(path);
                 %For every positive Filter
                 filterCount = numel(FMask.PositiveFilters);
                 for(i=1:filterCount);
                    % Create respective x,y coordinates for filter
                    f1 = FMask.PositiveFilters{i};
                    firstInd=find(f1,1, 'first');
                    lastInd=find(f1,1, 'last');
                    [firstY firstX] = ind2sub([sizeY sizeX], firstInd);
                    [lastY lastX] = ind2sub([sizeY sizeX], lastInd);
                    % Forwards previous variables and new x,y coordinates
                    % of filter to MorphologyCorticalTrans
                    [neuronHandler summary] = MorphologyCorticalTrans(selectedWell,firstY,lastY, firstX, lastX,0,optionHandler,neuronHandler,NucleusM,foldername,channelName);                     
                end
             end
         else
             % If no filter is chosen, the complete image is evluated in
             % MorphologyCorticalTrans. If the image is bigger than 10000,
             % the image is divided to get no out of memory. This should
             % not be the case for any cortical neuron image
             if(lastY - firstY > 10000 || lastX - firstX > 10000)        
                sizeXsmall = (lastX-firstX)/2;
                sizeYsmall = (lastY-firstY)/2;
                [neuronHandler summary] = MorphologyCorticalTrans(selectedWell,firstY,sizeYsmall,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyCorticalTrans(selectedWell,firstY,sizeYsmall,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyCorticalTrans(selectedWell,sizeYsmall+1,lastY,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyCorticalTrans(selectedWell,sizeYsmall+1,lastY,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
             else
                 % if the image is smaller than 10000, original x,y
                 % coordinates are used
                [NeuronM neuronHandler summary] = MorphologyCorticalTrans(NeuronM, selectedWell,firstY,lastY, firstX, lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);  
              end
              
         end    
         
        % system('taskkill /fi "WINDOWTITLE eq FBird"');
        %Calculate only for filter 
  
% This function smothes the binary image obtained from load cortical neurons and generates a skeleton from the smoothed binary image. This function can be run on the complete image or only on filters depending on the settings in the GUI options        
function [NeuronM neuronHandler summary] = MorphologyCorticalTrans(NeuronM, selectedWell, yStartPos, yEndPos, xStartPos, xEndPos, save, optionHandler, neuronHandler, NucleusM, foldername, channelName)
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end
imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
nucleusImage = imread(imagePathNucleusBig);

[sizeYorig sizeXorig] = size(nucleusImage);
%yStartPos = 1;
%yEndPos= sizeYorig;
%xStartPos = 1;
%xEndPos = sizeXorig;
NucleusM = NucleusM(yStartPos:yEndPos,xStartPos:xEndPos);

%[D IDX] = bwdist(full(NucleusM));
TN=0;
if(isfield(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
    NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
    NeuronManualM = NeuronManualM(yStartPos:yEndPos,xStartPos:xEndPos);
    TN = nnz(NeuronManualM);
else
    NeuronManualM = sparse(logical(zeros(sizeYorig,sizeXorig)));
    if(isa(neuronHandler.ManualNeuronPositionsSparse,'containers.Map') == 0)
        neuronHandler.ManualNeuronPositionsSparse = containers.Map();
    end
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
end
%Load neurite and nucleus image from file 
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end

% From here onwards starts the actual algorithm!
%These are for now set fixed parameters, but can later be included into
%Option menue


% This variable is used in case we want to analyse primary structures. Till
% now this feature is not implemented, but might get implemented in a later
% version of OMnisphero
AxonMeasurement=0;

% Image locations:
imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
imagePathBinary = [foldername '/ConvertedCellomics/' selectedWellLong 'Binary' optionHandler.FileEnding];
if(AxonMeasurement==1);
imagePathBinaryHard = [foldername '/ConvertedCellomics/' selectedWellLong 'BinaryHard' optionHandler.FileEnding];
end
%Check if NeuronStatMatrix exists
if(~isa(neuronHandler.NeuronStatMatrix,'containers.Map'))
    neuronHandler.NeuronStatMatrix = containers.Map();
end
if(~isa(neuronHandler.OmniNucleusPositions,'containers.Map'))
    neuronHandler.OmniNucleusPositions = containers.Map();
end
currentNeuronStat = containers.Map();
%Read images from imagepath
binaryimgNeurite = imread(imagePathBinary);




if(AxonMeasurement==1);
binaryimgNeuriteHard = imread(imagePathBinaryHard);
end
nucleusPicBinary = imread(imagePathWatershed);

%On this instance we will delte user delted binary particles, if
%applicable (both neuronal cell bodies and binary nuclei):
%1) We will generate a matrix containing all manual counted Neurons!
NeuronDummy = zeros(size(binaryimgNeurite));
NeuronM = NeuronM + NeuronDummy;
%Now we label all particles in binaryimgNeurite and nucleusPicBinary
binaryimgNeuriteLabel = bwlabel(binaryimgNeurite);
nucleusPicBinaryLabel = bwlabel(nucleusPicBinary);
ParticleCount = max(max(binaryimgNeuriteLabel));
DummyNeuron = zeros(size(binaryimgNeurite));
DummyNuclei = zeros(size(binaryimgNeurite));
if(sum(sum(NeuronM)) >0)
for(u=1:ParticleCount)
    binaryimgNeuriteOne = binaryimgNeuriteLabel;
    ii = binaryimgNeuriteOne == u;
    binaryimgNeuriteOne(ii) = 255;
    binaryimgNeuriteOne = binaryimgNeuriteOne - 250;
    binaryimgNeuriteOne = uint8(binaryimgNeuriteOne);
    binaryimgNeuriteOne = logical(binaryimgNeuriteOne);
    nucleusPicBinaryOne = nucleusPicBinaryLabel;
    ii = nucleusPicBinaryOne == u;
    nucleusPicBinaryOne(ii) = 255;
    nucleusPicBinaryOne = nucleusPicBinaryOne - 250;
    nucleusPicBinaryOne = uint8(nucleusPicBinaryOne);
    nucleusPicBinaryOne = logical(nucleusPicBinaryOne);
    TestNeuron = binaryimgNeuriteOne + NeuronM;
    TestNuclei = nucleusPicBinaryOne + NeuronM;
    if(max(max(TestNeuron))==2)
        binaryimgNeuriteOne = zeros(size(binaryimgNeurite));
    end
    if(max(max(TestNuclei))==2)
        nucleusPicBinaryOne = zeros(size(binaryimgNeurite));
    end
DummyNeuron= DummyNeuron+binaryimgNeuriteOne;
DummyNuclei = DummyNuclei +nucleusPicBinaryOne;

end
binaryimgNeurite = DummyNeuron;
nucleusPicBinary = DummyNuclei;
newfile = [selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
imwrite(nucleusPicBinary,[foldername '/ConvertedCellomics/' newfile] );
end
sizeY = yEndPos - yStartPos+1;
sizeX = xEndPos - xStartPos+1;
[sizeYOrig sizeXOrig]= size(binaryimgNeurite);

%Threshold Neurite Picture the known way
    
[sizeY sizeX] = size(binaryimgNeurite);
% Cut of filter neuronal image
if(optionHandler.FilterMulPa~=1 && optionHandler.CutOutSphereCore == 1)
    neuritePic = CutOutCircles(binaryimgNeurite,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM, 0);
end
% Cut of filter neuronal image
if(optionHandler.FilterMulPa~=1 && optionHandler.CutOutSphereCore == 1)
    nucleusPic = CutOutCircles(nucleusPicBinary,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM,0, 0);
end
% Cut of filter neuronal image
nucleusPicBinary = nucleusPicBinary(yStartPos:yEndPos,xStartPos:xEndPos);
[sizeYsmall sizeXsmall] = size(nucleusPicBinary);
% Cut of filter neuronal image
binaryimgNeurite = binaryimgNeurite(yStartPos:yEndPos,xStartPos:xEndPos);
if(AxonMeasurement==1);
binaryimgNeuriteHard = binaryimgNeuriteHard(yStartPos:yEndPos,xStartPos:xEndPos);
end
binaryimgNeurite=logical(binaryimgNeurite);
binaryimgNeuriteOrig = binaryimgNeurite;
  
%Algorithm to smooth edges of primary binary image from load cortical neurons:
%1) Find edges using canny method:
E = edge(binaryimgNeurite,'canny');
% Convert image to logical
E = logical(E);
%2) Add edge image to normal image
binaryimgNeurite = binaryimgNeurite + E;
%3) This image will have holes which have to be filled. Since crossing
%Neurite will create 'big' holes, the small 'holes' have to be extracted!
% Fill holes in image
Holes = imfill(binaryimgNeurite,'holes');
% Get only filled areas
Holes = Holes - binaryimgNeurite;
% Extract 'big' holes
HolesBig = bwareaopen(Holes, 10);
% Get small holes, by supstraction the big holes from all holes
Holes = Holes - HolesBig;
% Fill small holes:
binaryimgNeurite= binaryimgNeurite + Holes;
%Try to close small gaps with bridge function:
binaryimgNeuriteTry = bwmorph(binaryimgNeurite,'majority');
% In some instances the thickening of the binary image will combine two
% binary images resulting in two nuclei but only one binary image.
% Therefore we need to check here, wheter we reduce the number of binary
% particles and if yes to seperate them again:
%1) Get number of particles:
binaryimgNeuriteOrigLabel = bwlabel(binaryimgNeuriteOrig);
binaryimgNeuriteTryLabel = bwlabel(binaryimgNeuriteTry);
nParticlebinaryOrig = max(max(binaryimgNeuriteOrigLabel));
nParticlebinary = max(max(binaryimgNeuriteTryLabel));
if(nParticlebinary==nParticlebinaryOrig)
    %Create image difference
    binaryimgNeurite= bwmorph(binaryimgNeurite,'majority');
    
end    
%Creae skeleton by infinite thinning of the binary image
SkeletonImage = bwmorph(binaryimgNeurite,'thin',inf);
% Spur removes the last 4 pixel from endpoints. Thereby very small
% processes are eliminated
SkeletonImage = bwmorph(SkeletonImage,'spur',4);
        % only applies for AxonMeasurement
        if(AxonMeasurement==1);
        binaryimgNeuriteHard=uint8(binaryimgNeuriteHard);
        binaryimgNeuriteHard(binaryimgNeuriteHard==0)=255;
        binaryimgNeuriteHard(binaryimgNeuriteHard==1)=0;
        javaaddpath('lib/ij.jar');
        javaaddpath('lib/mij.jar');
        MIJ.start(false);
        MIJ.createImage(binaryimgNeuriteHard);
        MIJ.run('Options...', 'iterations=1 count=1 edm=Overwrite do=Nothing');
        MIJ.run('Skeletonize');
        %nucleusPicBin = logical(MIJ.getCurrentImage);
        MIJ.run('Save', 'save=temp.tif');
        MIJ.closeAllWindows;
        MIJ.exit;
        SkeletonImageHard = imread('temp.tif');
        delete temp.tif;
        SkeletonImageHard = logical(1-SkeletonImageHard);
        SkeletonImageHard = bwmorph(SkeletonImageHard,'spur',10);
        SkeletonImageHard = bwareaopen(SkeletonImageHard,50);
        end
        
        %Convert image to 8-bit
        SkeletonImage = uint8(SkeletonImage);
        if(AxonMeasurement==1);
        SkeletonImageHard = uint8(SkeletonImageHard);
        end
        %Images are saved now. In case a filter was chosen, only the area
        %of the filter is saved
        Dummy = zeros(sizeYOrig,sizeXOrig);
        Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = SkeletonImage;
        SkeletonImage = Dummy;
        newfile = [selectedWellLong 'Skeleton' optionHandler.FileEnding];
        imwrite(SkeletonImage,[foldername '/ConvertedCellomics/' newfile] );
        if(AxonMeasurement==1);
        Dummy = zeros(sizeYOrig,sizeXOrig);
        Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = SkeletonImageHard;
        SkeletonImageHard = Dummy;
        newfile = [selectedWellLong 'AstroBig' optionHandler.FileEnding];
        imwrite(SkeletonImageHard,[foldername '/ConvertedCellomics/' newfile] );
        newfile = [selectedWellLong 'AstroSmall' optionHandler.FileEnding];
        imwrite(SkeletonImageHard,[foldername '/ConvertedCellomics/' newfile] );
        end
        Dummy = zeros(sizeYOrig,sizeXOrig);
        Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = binaryimgNeurite;
        binaryimgNeurite = Dummy;
        newfile = [selectedWellLong 'Binary' optionHandler.FileEnding];
        imwrite(binaryimgNeurite,[foldername '/ConvertedCellomics/' newfile] );
        Dummy = zeros(sizeYOrig,sizeXOrig);
        Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = NeuronM;
        NeuronM = Dummy;
        summary=0;








% --------------------------------------------------------------------
function GFAPQuantification_Callback(hObject, eventdata, handles)
% hObject    handle to GFAPQuantification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
%foldername = uigetdir;
%Call EdgeFill Neurons and Skeleton Neurons and play with parameters in
%Option Handler
imageHandler = handles.ImageHandler;
foldername = imageHandler.Foldername;
ExecuteOptionExperimentGFAPOnly(handles,1,'Neurite');
guidata(handles.figure1, handles);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
[summary] = MeasurementGFAPOnly (foldername, csvHandler)
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function ExecuteOptionExperimentGFAPOnly(handles,onlyCompFill,channelName)
h = waitbar(0.5,'Neuronal morpholgy analyzed. This could take a few hours. Please be patient.', 'WindowStyle', 'modal');
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
optionHandler = handles.OptionHandler;

if(strcmp(channelName,'Oligo'))
    competingHandler = handles.NeuronCoordinates;
    if(isfield(handles,'OligoCoordinates'))
       neuronHandler = handles.OligoCoordinates;
    else
       neuronHandler = CSVCoordinates();
    end
else
    neuronHandler = handles.NeuronCoordinates;
    competingHandler = 0;
end
wellList = get(handles.lbWell, 'string');
foldername = imageHandler.Foldername;
[sizeY sizeX] = size(imageHandler.NucleusImage);
%Select STARTWELL

if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
if(numel(neuronHandler.AreaDictionary) ==0)
    neuronHandler.AreaDictionary = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
    neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFill) ==0)
    neuronHandler.NeuronPositionsEdgeFill = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFillSecond) ==0)
    neuronHandler.NeuronPositionsEdgeFillSecond = containers.Map();
end
if(numel(neuronHandler.NeuronStatMatrix) ==0)
    neuronHandler.NeuronStatMatrix = containers.Map();
end

%Create Matlab job
%Leave one Core out to ensure that GUI isn't freezing
maxCores =optionHandler.NumberOfCores;
if(maxCores > 1)
    c = parcluster();
    c.NumWorkers = maxCores;
    job = createJob(c,'JobData',optionHandler);
end

for(j=optionHandler.StartWell:numel(wellList))
    selectedWell=wellList(j);
    selectedWell=selectedWell{1};
    selectedWellLong=selectedWell;
    if(length(selectedWell) == 2)
              selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    end
    NucleusM = csvHandler.CellPosMatrix(selectedWell);    
    if(size(NucleusM,2) > sizeX || size(NucleusM,1) > sizeY)
        NucleusM = NucleusM(1:sizeY,1:sizeX);
        csvHandler.CellPosMatrix(selectedWell) = NucleusM;
        handles.CSVCoordinates = csvHandler;
        guidata(handles.figure1, handles);
    end
    if(maxCores > 1)
        %Check if watershed images are available. If not create them before
        %processing as ImageJ operation can't be done in parallel.
        %Check if there is already a watershed image. If not, create it.
        imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
        if(~exist(imagePathWatershed,'file'))
            WatershedNuclei(handles);
        end
        createTask(job, @ExecuteGFAPOnly, 4, {selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler});
    else
       [neuronHandler,  summary] = ExecuteGFAPOnly(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler);
    end
end
if(maxCores > 1)
    submit(job);
    wait(job);
    %job = out(2);
    out=fetchOutputs(job);
    %ToDo: Merge NeuronHandler
    TP = sum(cell2mat(out(:,2)));
    FP = sum(cell2mat(out(:,3)));
    TN = sum(cell2mat(out(:,4)));
end
% TPPerWell = cell2mat(out(:,2));
% FPPerWell = cell2mat(out(:,3));
% TNPerWell = cell2mat(out(:,4));
% TPPerWell = TPPerWell ./ (TPPerWell+FPPerWell);
% FPPerWell = FPPerWell ./ (TPPerWell+FPPerWell);
i=0;

if(maxCores > 1)
    for(j=optionHandler.StartWell:numel(wellList))
     %Get current NeuronHandler result and write it in main Neuronhandler
        i=i+1;
        selectedWell=wellList(j);
        selectedWell=selectedWell{1};    
    
        currentNeuronHandler = out(i,1);
        currentNeuronHandler=currentNeuronHandler{1};
        %Merge Neuron positions as well as neurite length matrix.
        neuronHandler.NeuronPositionsSkeletonization(selectedWell) = currentNeuronHandler.NeuronPositionsSkeletonization(selectedWell);
        neuronHandler.AreaDictionary(selectedWell) = currentNeuronHandler.AreaDictionary(selectedWell);
        neuronHandler.NeuronPositionsEdgeFill(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFill(selectedWell);
        neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
        neuronHandler.NeuronStatMatrix(selectedWell) = currentNeuronHandler.NeuronStatMatrix(selectedWell);
        neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = currentNeuronHandler.NeuronPositionsSkelDeleted(selectedWell);
        %neuronHandler.NeuriteLengthMatrix(selectedWell) = currentNeuronHandler.NeuriteLengthMatrix;
    end
end
close(h);
if(strcmp(channelName,'Oligo'))
     handles.OligoCoordinates = neuronHandler;
else
    handles.NeuronCoordinates = neuronHandler;
end
guidata(handles.figure1, handles);


 function [neuronHandler, summary] = ExecuteGFAPOnly(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler)    
        %Get first filter from selected Well
        selectedWellLong=selectedWell;
        if(length(selectedWell) == 2)
                  selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
        end
        imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
        nucleusImage = imread(imagePathNucleusBig);
        [sizeYorig sizeXorig] = size(nucleusImage);
        firstY=1;
        firstX=1;
        lastY=sizeYorig;
        lastX=sizeXorig;
        sizeY = sizeYorig;
        sizeX = sizeXorig;
         if(optionHandler.FilterMulPa==1)
             path = [foldername '/ConvertedCellomics/' selectedWell '.mat'];
             TP=0;
             FP=0;
             TN=0;
             if(exist(path,'file')) 
                 load(path);
                 %For every positive Filter
                 filterCount = numel(FMask.PositiveFilters);
                 for(i=1:filterCount);
                    f1 = FMask.PositiveFilters{i};
                    firstInd=find(f1,1, 'first');
                    lastInd=find(f1,1, 'last');
                    [firstY firstX] = ind2sub([sizeY sizeX], firstInd);
                    [lastY lastX] = ind2sub([sizeY sizeX], lastInd);
                    [neuronHandler summary] = MorphologyGFAPOnly(selectedWell,firstY,lastY, firstX, lastX,0,optionHandler,neuronHandler,NucleusM,foldername,channelName);                     
                end
             end
         else
             if(lastY - firstY > 10000 || lastX - firstX > 10000)        
                sizeXsmall = (lastX-firstX)/2;
                sizeYsmall = (lastY-firstY)/2;
                [neuronHandler summary] = MorphologyGFAPOnly(selectedWell,firstY,sizeYsmall,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyGFAPOnly(selectedWell,firstY,sizeYsmall,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyGFAPOnly(selectedWell,sizeYsmall+1,lastY,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyGFAPOnly(selectedWell,sizeYsmall+1,lastY,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
              else
                [neuronHandler summary] = MorphologyGFAPOnly(selectedWell,firstY,lastY, firstX, lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);  
              end
              
         end    
         
        % system('taskkill /fi "WINDOWTITLE eq FBird"');
        %Calculate only for filter 
        
function [neuronHandler summary] = MorphologyGFAPOnly(selectedWell, yStartPos, yEndPos, xStartPos, xEndPos, save, optionHandler, neuronHandler, NucleusM, foldername, channelName)
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end
imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
nucleusImage = imread(imagePathNucleusBig);


% Load filter if xStart is unequal 1
if(xStartPos > 1)
filterNumber = selectedWellLong(1,3);
filterNumber = str2num(filterNumber);
path = [foldername '/ConvertedCellomics/' selectedWellLong '.mat'];
if(exist(path,'file'))
        load(path);
        imageHandler.PositiveFilters = FMask.PositiveFilters;
        imageHandler.NegativeFilters = FMask.NegativeFilters;
end
FMask=0;

    %Get Filter
    filterList = imageHandler.PositiveFilters;
    filterNumber = size(filterList);
    filter = filterNumber(2);
    filterMasks = zeros(size(nucleusImage));
    for(u=1:filter)
    filterNumber = u;    
    filterHandle = filterList(filterNumber);
    filterHandle = filterHandle{1};
    filterMask = (filterHandle);
    %Remove excluded parts of picture
    if(iscell(imageHandler.NegativeFilters) > 0)
        negativeMask = imageHandler.NegativeFilters;
        for i=1:numel(negativeMask)
          workMask = negativeMask{i};
          workMask = int8(workMask -1);
          workMask = logical(workMask .* (-1));
          workMask = sparse(logical(workMask));
          filterMask = (logical((filterMask .* workMask)));
        end
    end
    filterMasks = filterMasks + filterMask;
    end
    workMask=0;
    filterMask = filterMasks;
    %Multiply mask with image
    nucleusImage = uint8(nucleusImage).*uint8(full(filterMask));
    % With the help of the filterMask we can now determine the closest
    % sourround box coordinates
    Coordinates = regionprops(filterMask,'BoundingBox');
    xStartPos = Coordinates(1).BoundingBox(1);
    xEndPos = xStartPos + Coordinates(1).BoundingBox(3);
    yStartPos = Coordinates(1).BoundingBox(2);
    yEndPos = yStartPos + Coordinates(1).BoundingBox(4);
end






[sizeYorig sizeXorig] = size(nucleusImage);
%yStartPos = 1;
%yEndPos= sizeYorig;
%xStartPos = 1;
%xEndPos = sizeXorig;
NucleusM = NucleusM(yStartPos:yEndPos,xStartPos:xEndPos);

%[D IDX] = bwdist(full(NucleusM));
TN=0;
if(isfield(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
    NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
    NeuronManualM = NeuronManualM(yStartPos:yEndPos,xStartPos:xEndPos);
    TN = nnz(NeuronManualM);
else
    NeuronManualM = sparse(logical(zeros(sizeYorig,sizeXorig)));
    if(isa(neuronHandler.ManualNeuronPositionsSparse,'containers.Map') == 0)
        neuronHandler.ManualNeuronPositionsSparse = containers.Map();
    end
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
end
%Load neurite and nucleus image from file 
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end

% From here onwards starts the actual algorithm!
%These are for now set fixed parameters, but can later be included into
%Option menue


% Image locations:
imagePathBinary = [foldername '/ConvertedCellomics/' selectedWellLong 'Binary' optionHandler.FileEnding];
binaryimgNeurite = imread(imagePathBinary);
if(xStartPos > 1)
binaryimgNeurite = uint8(binaryimgNeurite).*uint8(full(filterMask));
end


sizeY = yEndPos - yStartPos+1;
sizeX = xEndPos - xStartPos+1;
[sizeYOrig sizeXOrig]= size(binaryimgNeurite);

binaryimgNeurite = binaryimgNeurite(yStartPos:yEndPos,xStartPos:xEndPos);
Dummy = zeros(sizeYOrig,sizeXOrig);
Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = binaryimgNeurite;
binaryimgNeurite = Dummy;
newfile = [selectedWellLong 'Binary' optionHandler.FileEnding];
imwrite(binaryimgNeurite,[foldername '/ConvertedCellomics/' newfile] );

%Check if NeuronStatMatrix exists
if(~isa(neuronHandler.NeuronStatMatrix,'containers.Map'))
    neuronHandler.NeuronStatMatrix = containers.Map();
end
if(~isa(neuronHandler.OmniNucleusPositions,'containers.Map'))
    neuronHandler.OmniNucleusPositions = containers.Map();
end
currentNeuronStat = containers.Map();


        
        NeuronM = NucleusM;
       
        Dummy = zeros(sizeYorig ,sizeXorig);
        Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = NeuronM;
        NeuronM = Dummy;
        summary=0;





neuronHandler.NeuronPositionsEdgeFill(selectedWell) = NeuronM;


% --------------------------------------------------------------------
function LoadSphereCorona_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSphereCorona (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
neuronHandler = CSVCoordinates();
neuronHandler.ManualNeuronPositions = containers.Map();
csvHandler = CSVCoordinates();
 imageHandler = ImageHandler();
 optionHandler = Option();
 optionHandler.MinSizeNeuriteArea = 25000;
 handles.OptionHandler = optionHandler;
optionHandler.FileEnding = '.png';
[wellList folder ScalingFactor csvHandler neuronHandler] = SphereOutgrowth(csvHandler, neuronHandler);
optionHandler.ScalingFactor = ScalingFactor;
 imageHandler.Foldername = folder;
set(handles.lbWell, 'String', wellList);
handles.CSVCoordinates = csvHandler;
handles.NeuronCoordinates = neuronHandler;
%NeuronPositionsSkeletonization = neuronHandler.NeuronPositionsSkeletonization(wellname);
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles);




% --------------------------------------------------------------------
function SphereOutgrowth_Callback(hObject, eventdata, handles)
% hObject    handle to SphereOutgrowth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
%foldername = uigetdir;
%Call EdgeFill Neurons and Skeleton Neurons and play with parameters in
%Option Handler
imageHandler = handles.ImageHandler;
foldername = imageHandler.Foldername;
csvHandler = handles.CSVCoordinates;
[summary] = MeasurementSphereOutgrowth(foldername, csvHandler);
guidata(handles.figure1, handles);
neuronHandler = handles.NeuronCoordinates;
[summary] = MeasurementGFAPOnly (foldername, csvHandler)
fileName = 'Results';
xlswrite(fileName,summary)
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end






% --------------------------------------------------------------------
function SphereRadius_Callback(hObject, eventdata, handles)
% hObject    handle to SphereRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to SphereOutgrowth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
%foldername = uigetdir;
%Call EdgeFill Neurons and Skeleton Neurons and play with parameters in
%Option Handler
imageHandler = handles.ImageHandler;
foldername = imageHandler.Foldername;
csvHandler = handles.CSVCoordinates;
wellList = get(handles.lbWell, 'string');
WatershedNuclei(handles);
[summary] = SphereRadius(wellList, foldername, csvHandler);

guidata(handles.figure1, handles);
neuronHandler = handles.NeuronCoordinates;
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end







% --------------------------------------------------------------------
function LoadIPSSynapse_Callback(hObject, eventdata, handles)
% hObject    handle to LoadIPSSynapse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
neuronHandler = CSVCoordinates();
csvHandler = CSVCoordinates();
 imageHandler = ImageHandler();
 optionHandler = Option();
 optionHandler.MinSizeNeuriteArea = 25000;
 handles.OptionHandler = optionHandler;
optionHandler.FileEnding = '.png';
[wellList folder ScalingFactor csvHandler neuronHandler] = IPSSynapseCount(csvHandler, neuronHandler);
optionHandler.ScalingFactor = ScalingFactor;
 imageHandler.Foldername = folder;
set(handles.lbWell, 'String', wellList);
handles.CSVCoordinates = csvHandler;
handles.NeuronCoordinates = neuronHandler;

%NeuronPositionsSkeletonization = neuronHandler.NeuronPositionsSkeletonization(wellname);
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles);


% --------------------------------------------------------------------
function IPSSynapeAnalyze_Callback(hObject, eventdata, handles)
% hObject    handle to IPSSynapeAnalyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
% Get variables from opened GUI    
handles = guidata(handles.figure1);
% Get variable of Options from set seetings in GUI
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
foldername = imageHandler.Foldername;
ExecuteOptionExperimentIPSSynapse(handles,1,'Neurite');
guidata(handles.figure1, handles);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
[summary] = MeasurementIPS (foldername, csvHandler, neuronHandler)
filename = [foldername '/Results.xlsx'];
xlswrite(filename,summary);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function ExecuteOptionExperimentIPSSynapse(handles,onlyCompFill,channelName)
h = waitbar(0.5,'Neuronal morpholgy analyzed. This could take a few hours. Please be patient.', 'WindowStyle', 'modal');
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
optionHandler = handles.OptionHandler;

if(strcmp(channelName,'Oligo'))
    competingHandler = handles.NeuronCoordinates;
    if(isfield(handles,'OligoCoordinates'))
       neuronHandler = handles.OligoCoordinates;
    else
       neuronHandler = CSVCoordinates();
    end
else
    neuronHandler = handles.NeuronCoordinates;
    competingHandler = 0;
end
wellList = get(handles.lbWell, 'string');
foldername = imageHandler.Foldername;
[sizeY sizeX] = size(imageHandler.NucleusImage);
%Select STARTWELL

if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
if(numel(neuronHandler.AreaDictionary) ==0)
    neuronHandler.AreaDictionary = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
    neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFill) ==0)
    neuronHandler.NeuronPositionsEdgeFill = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFillSecond) ==0)
    neuronHandler.NeuronPositionsEdgeFillSecond = containers.Map();
end
if(numel(neuronHandler.NeuronStatMatrix) ==0)
    neuronHandler.NeuronStatMatrix = containers.Map();
end

%Create Matlab job
%Leave one Core out to ensure that GUI isn't freezing
maxCores =optionHandler.NumberOfCores;
if(maxCores > 1)
    c = parcluster();
    c.NumWorkers = maxCores;
    job = createJob(c,'JobData',optionHandler);
end

for(j=optionHandler.StartWell:numel(wellList))
    selectedWell=wellList(j);
    selectedWell=selectedWell{1};
    selectedWellLong=selectedWell;
    if(length(selectedWell) == 2)
              selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    end
    
    NucleusM = csvHandler.CellPosMatrix(selectedWell);    
    if(size(NucleusM,2) > sizeX || size(NucleusM,1) > sizeY)
        NucleusM = NucleusM(1:sizeY,1:sizeX);
        csvHandler.CellPosMatrix(selectedWell) = NucleusM;
        handles.CSVCoordinates = csvHandler;
        guidata(handles.figure1, handles);
    end
    if(maxCores > 1)
        %Check if watershed images are available. If not create them before
        %processing as ImageJ operation can't be done in parallel.
        %Check if there is already a watershed image. If not, create it.
        imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
        if(~exist(imagePathWatershed,'file'))
            WatershedNuclei(handles);
        end
        createTask(job, @ExecuteIPSSynapes, 4, {selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler});
    else
       [neuronHandler,  summary] = ExecuteIPSSynapes(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler);
    end
end
if(maxCores > 1)
    submit(job);
    wait(job);
    %job = out(2);
    out=fetchOutputs(job);
    %ToDo: Merge NeuronHandler
    TP = sum(cell2mat(out(:,2)));
    FP = sum(cell2mat(out(:,3)));
    TN = sum(cell2mat(out(:,4)));
end
% TPPerWell = cell2mat(out(:,2));
% FPPerWell = cell2mat(out(:,3));
% TNPerWell = cell2mat(out(:,4));
% TPPerWell = TPPerWell ./ (TPPerWell+FPPerWell);
% FPPerWell = FPPerWell ./ (TPPerWell+FPPerWell);
i=0;

if(maxCores > 1)
    for(j=optionHandler.StartWell:numel(wellList))
     %Get current NeuronHandler result and write it in main Neuronhandler
        i=i+1;
        selectedWell=wellList(j);
        selectedWell=selectedWell{1};    
    
        currentNeuronHandler = out(i,1);
        currentNeuronHandler=currentNeuronHandler{1};
        %Merge Neuron positions as well as neurite length matrix.
        neuronHandler.NeuronPositionsSkeletonization(selectedWell) = currentNeuronHandler.NeuronPositionsSkeletonization(selectedWell);
        neuronHandler.AreaDictionary(selectedWell) = currentNeuronHandler.AreaDictionary(selectedWell);
        neuronHandler.NeuronPositionsEdgeFill(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFill(selectedWell);
        neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
        neuronHandler.NeuronStatMatrix(selectedWell) = currentNeuronHandler.NeuronStatMatrix(selectedWell);
        neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = currentNeuronHandler.NeuronPositionsSkelDeleted(selectedWell);
        %neuronHandler.NeuriteLengthMatrix(selectedWell) = currentNeuronHandler.NeuriteLengthMatrix;
    end
end
close(h);
if(strcmp(channelName,'Oligo'))
     handles.OligoCoordinates = neuronHandler;
else
    handles.NeuronCoordinates = neuronHandler;
end
guidata(handles.figure1, handles);


 function [neuronHandler, summary] = ExecuteIPSSynapes(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler)    
        %Get first filter from selected Well
        selectedWellLong=selectedWell;
        if(length(selectedWell) == 2)
                  selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
        end
        imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
        nucleusImage = imread(imagePathNucleusBig);
        [sizeYorig sizeXorig] = size(nucleusImage);
        firstY=1;
        firstX=1;
        lastY=sizeYorig;
        lastX=sizeXorig;
        sizeY = sizeYorig;
        sizeX = sizeXorig;
         if(optionHandler.FilterMulPa==1)
             path = [foldername '/ConvertedCellomics/' selectedWell '.mat'];
             TP=0;
             FP=0;
             TN=0;
             if(exist(path,'file')) 
                 load(path);
                 %For every positive Filter
                 filterCount = numel(FMask.PositiveFilters);
                 for(i=1:filterCount);
                    f1 = FMask.PositiveFilters{i};
                    firstInd=find(f1,1, 'first');
                    lastInd=find(f1,1, 'last');
                    [firstY firstX] = ind2sub([sizeY sizeX], firstInd);
                    [lastY lastX] = ind2sub([sizeY sizeX], lastInd);
                    [neuronHandler summary] = MorphologyCorticalNuc(selectedWell,firstY,lastY, firstX, lastX,0,optionHandler,neuronHandler,NucleusM,foldername,channelName);                     
                end
             end
         else
             if(lastY - firstY > 10000 || lastX - firstX > 10000)        
                sizeXsmall = (lastX-firstX)/2;
                sizeYsmall = (lastY-firstY)/2;
                [neuronHandler summary] = MorphologyIPSSynapses(selectedWell,firstY,sizeYsmall,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyIPSSynapses(selectedWell,firstY,sizeYsmall,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyIPSSynapses(selectedWell,sizeYsmall+1,lastY,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyIPSSynapses(selectedWell,sizeYsmall+1,lastY,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
              else
                [neuronHandler summary] = MorphologyIPSSynapses(selectedWell,firstY,lastY, firstX, lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);  
              end
              
         end    
         
        % system('taskkill /fi "WINDOWTITLE eq FBird"');
        %Calculate only for filter 
        
function [neuronHandler summary] = MorphologyIPSSynapses(selectedWell, yStartPos, yEndPos, xStartPos, xEndPos, save, optionHandler, neuronHandler, NucleusM, foldername, channelName)
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end
imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
nucleusImage = imread(imagePathNucleusBig);

[sizeYorig sizeXorig] = size(nucleusImage);
%yStartPos = 1;
%yEndPos= sizeYorig;
%xStartPos = 1;
%xEndPos = sizeXorig;
NucleusM = NucleusM(yStartPos:yEndPos,xStartPos:xEndPos);

%[D IDX] = bwdist(full(NucleusM));
TN=0;
if(isfield(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
    NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
    NeuronManualM = NeuronManualM(yStartPos:yEndPos,xStartPos:xEndPos);
    TN = nnz(NeuronManualM);
else
    NeuronManualM = sparse(logical(zeros(sizeYorig,sizeXorig)));
    if(isa(neuronHandler.ManualNeuronPositionsSparse,'containers.Map') == 0)
        neuronHandler.ManualNeuronPositionsSparse = containers.Map();
    end
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
end
%Load neurite and nucleus image from file 
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end

% From here onwards starts the actual algorithm!
%These are for now set fixed parameters, but can later be included into
%Option menue

        summary=0;








% --------------------------------------------------------------------
function LoadImageCorticalNeuronsRat_Callback(hObject, eventdata, handles)
% hObject    handle to LoadImageCorticalNeuronsRat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to LoadImageCorticalNeurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
neuronHandler = CSVCoordinates();
neuronHandler.ManualNeuronPositions = containers.Map();
csvHandler = CSVCoordinates();
csvHandler.ManualPositions1 = containers.Map();
csvHandler.ManualPositions2 = containers.Map();
 imageHandler = ImageHandler();
 optionHandler = Option();
 optionHandler.MinSizeNeuriteArea = 25000;
 handles.OptionHandler = optionHandler;
optionHandler.FileEnding = '.png';
Question = questdlg('Do you want to use an automated preselection?')
if(strcmp(Question,'No'))
    Auto = 0;
else
    Auto = 1;
end

[wellList folder ScalingFactor csvHandler ] = CorticalNeuronsRat(csvHandler, Auto);
optionHandler.ScalingFactor = ScalingFactor;
imageHandler.Foldername = folder;
set(handles.lbWell, 'String', wellList);
handles.CSVCoordinates = csvHandler;
handles.NeuronCoordinates = neuronHandler;
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles);
% Now we will save the file:
%Filename will be the foldername:
[token,remain] = strtok(folder, '\');
lengthRemain = length(remain);
while(length(remain)>0)
[token,remain] = strtok(remain, '\');
end
a{1} = token;
a{2} = '_Raw';
a{3}= '.mat';
FileName = [a{1} a{2} a{3}];
%Pathname is folder and a '\'
a{1} = folder;
a{2} = '\';
PathName = [a{1} a{2}];
try
SaveFixedNeuronPositions(FileName,PathName,handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --------------------------------------------------------------------
function CorticalMorphologyRat_Callback(hObject, eventdata, handles)
% hObject    handle to CorticalMorphologyRat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
% Defines handles for GUI    
handles = guidata(handles.figure1);
% Sets optionHandler to variable currently set in GUI under option menu
optionHandler = handles.OptionHandler;
% ImageHandler contains recent loaded images
imageHandler = handles.ImageHandler;
%image path is obtained from GUI
foldername = imageHandler.Foldername;
Question = questdlg('Do you want Omnisphero to identify Neurons for you?')
if(strcmp(Question,'No'))
    Stats = 0;
else
    Stats = 1;
end
Question = questdlg('Is this a reanalyzation?')
if(strcmp(Question,'No'))
    ReAnalysis = 0;
else
    ReAnalysis = 1;
end
% Executes function: ExecuteOptionExperimentCorticalTrans
ExecuteOptionExperimentCorticalTransRat(handles,1,'Neurite',Stats,ReAnalysis);
% Analysis results are returned to the GUI
guidata(handles.figure1, handles);
% csvHandler is obtained from GUI
csvHandler = handles.CSVCoordinates;
% csvHandler is obtained from GUI
neuronHandler = handles.NeuronCoordinates;



% Execution of MeasurementTran
[summary Summary2 Summary3 Summary4 Summary5 ExpName csvHandler neuronHandler PercentInd ] = MeasurementTransRat (foldername, csvHandler, neuronHandler, Stats,ReAnalysis)
% This function is used to write summary - Summary 4 in Excel files
% Now we will save the file:
%Pathname is folder and a '\'
a{1} = foldername;
a{2} = '\';
PathName = [a{1} a{2}];
%Filename will be the foldername:
[token,remain] = strtok(foldername, '\');
lengthRemain = length(remain);
while(length(remain)>0)
[token,remain] = strtok(remain, '\');
end
a{1} = token;
a{2} = '_Analyzed';
a{3} = '_0';
a{4}= '.mat';
a{5}='.zip';
FileName = [a{1} a{2} a{3} a{4}];
FileNameExist = [a{1} a{2} a{3} a{4} a{5}];
filename = regexprep(ExpName,'.png','.xlsx','ignorecase'); 
TempDir = dir(PathName);
TempDir=struct2cell(TempDir);
if(any(ismember(TempDir(1,:),FileNameExist)))
    a{1} = token;
    a{2} = '_Analyzed';
    a{3} = '_1';
    a{4}= '.mat';
    FileName = [a{1} a{2} a{3} a{4}];
    filename = regexprep(ExpName,'.png','_1.xlsx','ignorecase'); 
end

try
SaveFixedNeuronPositions(FileName,PathName,handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

Replace = 0;
[r c] = size(Summary2)
for(z=1:(c-1))
    if(Summary2{2,(z+1)} > 0)
        Replace = 1;
    end
    for(p=1:(r-1))
        if(Replace == 1)
        if(isempty(Summary2{(p+1),(z+1)}) == 1)
            Summary2{(p+1),(z+1)} = 0;
        end
        end
    end
    Replace = 0;
end

[r c] = size(Summary3)
for(z=1:(c-1))
    if(Summary3{2,(z+1)} > 0)
        Replace = 1;
    end
    for(p=1:(r-1))
        if(Replace == 1)
        if(isempty(Summary3{(p+1),(z+1)}) == 1)
            Summary3{(p+1),(z+1)} = 0;
        end
        end
    end
    Replace = 0;
end

[r c] = size(Summary4)

for(z=1:(c-1))
    if(Summary4{2,(z+1)} > 0)
        Replace = 1;
    end
    for(p=1:(r-1))
        if(Replace == 1)
        if(isempty(Summary4{(p+1),(z+1)}) == 1)
            Summary4{(p+1),(z+1)} = 0;
        end
        end
    end
    Replace = 0;
end
% Bounces summary - Summary4 back to command window
summary = summary
Summary2 = Summary2
Summary3 = Summary3
Summary4 = Summary4
Summary5 = Summary5
% Filename from MeasurementTrans is used as name for Excell files
 
%All summaries are written in the same excel file under different sheeds
xlswrite(filename,summary);
sheet = 2;
%xlswrite(filename,SummarySSE,sheet);
xlswrite(filename,Summary2,sheet);
sheet = 3;
xlswrite(filename,Summary3,sheet);
sheet = 4;
xlswrite(filename,Summary4,sheet);
sheet = 5;
xlswrite(filename,Summary5,sheet);
sheet = 6;
xlswrite(filename,PercentInd,sheet);

catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

% This function is mainly used to identify whether only neurons or also
% other cell types are evaluated and whether the job is performed on a
% single CPU core or on mutiple cores
function ExecuteOptionExperimentCorticalTransRat(handles,onlyCompFill,channelName, Stats, ReAnalysis)
h = waitbar(0.5,'Neuronal morpholgy analyzed. This could take a few hours. Please be patient.', 'WindowStyle', 'modal');
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
optionHandler = handles.OptionHandler;
%Here is checked whether Oligo data should be analyzed
if(strcmp(channelName,'Oligo'))
    competingHandler = handles.NeuronCoordinates;
    if(isfield(handles,'OligoCoordinates'))
       neuronHandler = handles.OligoCoordinates;
    else
       neuronHandler = CSVCoordinates();
    end
else
    neuronHandler = handles.NeuronCoordinates;
    competingHandler = 0;
end
wellList = get(handles.lbWell, 'string');
foldername = imageHandler.Foldername;
[sizeY sizeX] = size(imageHandler.NucleusImage);
%Select STARTWELL

if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
if(numel(neuronHandler.AreaDictionary) ==0)
    neuronHandler.AreaDictionary = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
    neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFill) ==0)
    neuronHandler.NeuronPositionsEdgeFill = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFillSecond) ==0)
    neuronHandler.NeuronPositionsEdgeFillSecond = containers.Map();
end
if(numel(neuronHandler.NeuronStatMatrix) ==0)
    neuronHandler.NeuronStatMatrix = containers.Map();
end

%Create Matlab job
%Leave one Core out to ensure that GUI isn't freezing
maxCores =optionHandler.NumberOfCores;
if(maxCores > 1)
    c = parcluster();
    c.NumWorkers = maxCores;
    job = createJob(c,'JobData',optionHandler);
end

for(j=optionHandler.StartWell:numel(wellList))
    selectedWell=wellList(j);
    selectedWell=selectedWell{1};
    selectedWellLong=selectedWell;
    if(length(selectedWell) == 2)
              selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    end
    NucleusM = csvHandler.CellPosMatrix(selectedWell); 
    if(size(csvHandler.ManualPositions1)>1)
        try
    %NeuronM = csvHandler.ManualPositions1(selectedWell);
    NeuronM = 0;
        catch
    NeuronM = 0;
        end
    else
    NeuronM = 0;  
    if(Stats ==1)
    else    
    %NucleusM = zeros(size(csvHandler.CellPosMatrix(selectedWell)));
    end
    end
    if(size(NucleusM,2) > sizeX || size(NucleusM,1) > sizeY)
        NucleusM = NucleusM(1:sizeY,1:sizeX);
        csvHandler.CellPosMatrix(selectedWell) = NucleusM;
        handles.CSVCoordinates = csvHandler;
        guidata(handles.figure1, handles);
    end
    if(maxCores > 1)
        %Check if watershed images are available. If not create them before
        %processing as ImageJ operation can't be done in parallel.
        %Check if there is already a watershed image. If not, create it.
        imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
        if(~exist(imagePathWatershed,'file'))
            WatershedNuclei(handles);
        end
        createTask(job, @ExecuteCorticalTransRat, 4, {selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler});
    else
        % This is usally started on a normal PC and obtains all relevant
        % data. For filters the dimensions are given to the next function
        % as sizeY,sizeX
       [NeuronM, neuronHandler,  summary] = ExecuteCorticalTransRat(NeuronM, selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler, Stats, ReAnalysis);
    end
end
if(maxCores > 1)
    submit(job);
    wait(job);
    %job = out(2);
    out=fetchOutputs(job);
    %ToDo: Merge NeuronHandler
    TP = sum(cell2mat(out(:,2)));
    FP = sum(cell2mat(out(:,3)));
    TN = sum(cell2mat(out(:,4)));
end
% TPPerWell = cell2mat(out(:,2));
% FPPerWell = cell2mat(out(:,3));
% TNPerWell = cell2mat(out(:,4));
% TPPerWell = TPPerWell ./ (TPPerWell+FPPerWell);
% FPPerWell = FPPerWell ./ (TPPerWell+FPPerWell);
i=0;

if(maxCores > 1)
    for(j=optionHandler.StartWell:numel(wellList))
     %Get current NeuronHandler result and write it in main Neuronhandler
        i=i+1;
        selectedWell=wellList(j);
        selectedWell=selectedWell{1};    
    
        currentNeuronHandler = out(i,1);
        currentNeuronHandler=currentNeuronHandler{1};
        %Merge Neuron positions as well as neurite length matrix.
        neuronHandler.NeuronPositionsSkeletonization(selectedWell) = currentNeuronHandler.NeuronPositionsSkeletonization(selectedWell);
        neuronHandler.AreaDictionary(selectedWell) = currentNeuronHandler.AreaDictionary(selectedWell);
        neuronHandler.NeuronPositionsEdgeFill(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFill(selectedWell);
        neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
        neuronHandler.NeuronStatMatrix(selectedWell) = currentNeuronHandler.NeuronStatMatrix(selectedWell);
        neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = currentNeuronHandler.NeuronPositionsSkelDeleted(selectedWell);
        %neuronHandler.NeuriteLengthMatrix(selectedWell) = currentNeuronHandler.NeuriteLengthMatrix;
    end
end
close(h);
if(strcmp(channelName,'Oligo'))
     handles.OligoCoordinates = neuronHandler;
else
    handles.NeuronCoordinates = neuronHandler;
end
guidata(handles.figure1, handles);

% This function is used to determine whether the whole umage is evaluated or only a filter 
 function [NeuronM, neuronHandler, summary, Stats, ReAnalysis] = ExecuteCorticalTransRat(NeuronM, selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler, Stats, ReAnalysis)    
        %Get first filter from selected Well
        selectedWellLong=selectedWell;
        if(length(selectedWell) == 2)
                  selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
        end
        imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
        nucleusImage = imread(imagePathNucleusBig);
        % Gives the original dimension of the image and defines x,y
        % coordinates of image edges
        [sizeYorig sizeXorig] = size(nucleusImage);
        firstY=1;
        firstX=1;
        lastY=sizeYorig;
        lastX=sizeXorig;
        sizeY = sizeYorig;
        sizeX = sizeXorig;
        % Here it is checked whether filter analysis is chosen in the GUI
        % or not
         if(optionHandler.FilterMulPa==1)
             path = [foldername '/ConvertedCellomics/' selectedWell '.mat'];
             TP=0;
             FP=0;
             TN=0;
             if(exist(path,'file')) 
                 load(path);
                 %For every positive Filter
                 filterCount = numel(FMask.PositiveFilters);
                 for(i=1:filterCount);
                    % Create respective x,y coordinates for filter
                    f1 = FMask.PositiveFilters{i};
                    firstInd=find(f1,1, 'first');
                    lastInd=find(f1,1, 'last');
                    [firstY firstX] = ind2sub([sizeY sizeX], firstInd);
                    [lastY lastX] = ind2sub([sizeY sizeX], lastInd);
                    % Forwards previous variables and new x,y coordinates
                    % of filter to MorphologyCorticalTrans
                    [neuronHandler summary] = MorphologyCorticalTransRat(selectedWell,firstY,lastY, firstX, lastX,0,optionHandler,neuronHandler,NucleusM,foldername,channelName,Stats,ReAnalysis);                     
                end
             end
         else
             % If no filter is chosen, the complete image is evluated in
             % MorphologyCorticalTrans. If the image is bigger than 10000,
             % the image is divided to get no out of memory. This should
             % not be the case for any cortical neuron image
             if(lastY - firstY > 10000 || lastX - firstX > 10000)        
                sizeXsmall = (lastX-firstX)/2;
                sizeYsmall = (lastY-firstY)/2;
                [neuronHandler summary] = MorphologyCorticalTransRat(selectedWell,firstY,sizeYsmall,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName,Stats,ReAnalysis);
                [neuronHandler summary] = MorphologyCorticalTransRat(selectedWell,firstY,sizeYsmall,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName,Stats,ReAnalysis);
                [neuronHandler summary] = MorphologyCorticalTransRat(selectedWell,sizeYsmall+1,lastY,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName,Stats,ReAnalysis);
                [neuronHandler summary] = MorphologyCorticalTransRat(selectedWell,sizeYsmall+1,lastY,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName,Stats,ReAnalysis);
             else
                 % if the image is smaller than 10000, original x,y
                 % coordinates are used
                [NeuronM neuronHandler summary] = MorphologyCorticalTransRat(NeuronM, selectedWell,firstY,lastY, firstX, lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName,Stats,ReAnalysis);  
              end
              
         end    
         
        % system('taskkill /fi "WINDOWTITLE eq FBird"');
        %Calculate only for filter 
  
% This function smothes the binary image obtained from load cortical neurons and generates a skeleton from the smoothed binary image. This function can be run on the complete image or only on filters depending on the settings in the GUI options        
function [NeuronM neuronHandler summary] = MorphologyCorticalTransRat(NeuronM, selectedWell, yStartPos, yEndPos, xStartPos, xEndPos, save, optionHandler, neuronHandler, NucleusM, foldername, channelName, Stats,ReAnalysis)
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end
imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
nucleusImage = imread(imagePathNucleusBig);

[sizeYorig sizeXorig] = size(nucleusImage);
%yStartPos = 1;
%yEndPos= sizeYorig;
%xStartPos = 1;
%xEndPos = sizeXorig;
NucleusM = NucleusM(yStartPos:yEndPos,xStartPos:xEndPos);

%[D IDX] = bwdist(full(NucleusM));
TN=0;
if(isfield(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
    NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
    NeuronManualM = NeuronManualM(yStartPos:yEndPos,xStartPos:xEndPos);
    TN = nnz(NeuronManualM);
else
    NeuronManualM = sparse(logical(zeros(sizeYorig,sizeXorig)));
    if(isa(neuronHandler.ManualNeuronPositionsSparse,'containers.Map') == 0)
        neuronHandler.ManualNeuronPositionsSparse = containers.Map();
    end
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
end
%Load neurite and nucleus image from file 
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end

% From here onwards starts the actual algorithm!
%These are for now set fixed parameters, but can later be included into
%Option menue


% This variable is used in case we want to analyse primary structures. Till
% now this feature is not implemented, but might get implemented in a later
% version of OMnisphero
AxonMeasurement=0;

% Image locations:
imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
imagePathBinary = [foldername '/ConvertedCellomics/' selectedWellLong 'Binary' optionHandler.FileEnding];
imagePathSkel = [foldername '/ConvertedCellomics/' selectedWellLong 'Skeleton' optionHandler.FileEnding];
if(AxonMeasurement==1);
imagePathBinaryHard = [foldername '/ConvertedCellomics/' selectedWellLong 'BinaryHard' optionHandler.FileEnding];
end
%Check if NeuronStatMatrix exists
if(~isa(neuronHandler.NeuronStatMatrix,'containers.Map'))
    neuronHandler.NeuronStatMatrix = containers.Map();
end
if(~isa(neuronHandler.OmniNucleusPositions,'containers.Map'))
    neuronHandler.OmniNucleusPositions = containers.Map();
end
currentNeuronStat = containers.Map();
%Read images from imagepath
binaryimgNeurite = imread(imagePathBinary);
SkeletonImage = imread(imagePathSkel);



if(AxonMeasurement==1);
binaryimgNeuriteHard = imread(imagePathBinaryHard);
end
nucleusPicBinary = imread(imagePathWatershed);

% %On this instance we will delte user delted binary particles, if
% %applicable (both neuronal cell bodies and binary nuclei):
% %1) We will generate a matrix containing all manual counted Neurons!
% NeuronDummy = zeros(size(binaryimgNeurite));
% NeuronM = NeuronM + NeuronDummy;
% if(Stats == 1)
%     NeuronM = NucleusM + NeuronDummy;
% end    
% %Now we label all particles in binaryimgNeurite and nucleusPicBinary
% binaryimgNeuriteLabel = bwlabel(binaryimgNeurite);
% nucleusPicBinaryLabel = bwlabel(nucleusPicBinary);
% ParticleCount = max(max(binaryimgNeuriteLabel));
% DummyNeuron = zeros(size(binaryimgNeurite));
% DummyNuclei = zeros(size(binaryimgNeurite));
% if(sum(sum(NeuronM)) >0)
% for(u=1:ParticleCount)
%     binaryimgNeuriteOne = binaryimgNeuriteLabel;
%     ii = binaryimgNeuriteOne == u;
%     binaryimgNeuriteOne(ii) = 255;
%     binaryimgNeuriteOne = binaryimgNeuriteOne - 250;
%     binaryimgNeuriteOne = uint8(binaryimgNeuriteOne);
%     binaryimgNeuriteOne = logical(binaryimgNeuriteOne);
%     nucleusPicBinaryOne = nucleusPicBinaryLabel;
%     ii = nucleusPicBinaryOne == u;
%     nucleusPicBinaryOne(ii) = 255;
%     nucleusPicBinaryOne = nucleusPicBinaryOne - 250;
%     nucleusPicBinaryOne = uint8(nucleusPicBinaryOne);
%     nucleusPicBinaryOne = logical(nucleusPicBinaryOne);
%     TestNeuron = binaryimgNeuriteOne + NeuronM;
%     TestNuclei = nucleusPicBinaryOne + NeuronM;
%     if(max(max(TestNeuron))==2)
%         binaryimgNeuriteOne = binaryimgNeuriteOne;
%     else    
%         binaryimgNeuriteOne = zeros(size(binaryimgNeurite));
%     end
%     if(max(max(TestNuclei))==2)
%         nucleusPicBinaryOne = nucleusPicBinaryOne;
%     else    
%         nucleusPicBinaryOne = zeros(size(binaryimgNeurite));
%     end
% DummyNeuron= DummyNeuron+binaryimgNeuriteOne;
% DummyNuclei = DummyNuclei +nucleusPicBinaryOne;
% 
% end
% binaryimgNeurite = DummyNeuron;
% nucleusPicBinary = DummyNuclei;
% newfile = [selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
% imwrite(nucleusPicBinary,[foldername '/ConvertedCellomics/' newfile] );
% end
% Here we need to delete everything of images if NucleusM is 0:
if(sum(sum(NucleusM)) == 0);
    nucleusPicBinary = zeros(size(binaryimgNeurite));
    newfile = [selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
    imwrite(nucleusPicBinary,[foldername '/ConvertedCellomics/' newfile] );
    binaryimgNeurite = zeros(size(binaryimgNeurite));
end  
% In case we have no neurons we also need to delete the binary images!!!
% if(sum(sum(NeuronM)) == 0);
%     nucleusPicBinary = zeros(size(binaryimgNeurite));
%     newfile = [selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
%     imwrite(nucleusPicBinary,[foldername '/ConvertedCellomics/' newfile] );
%     binaryimgNeurite = zeros(size(binaryimgNeurite));
% end 


sizeY = yEndPos - yStartPos+1;
sizeX = xEndPos - xStartPos+1;
[sizeYOrig sizeXOrig]= size(binaryimgNeurite);

%Threshold Neurite Picture the known way
    
[sizeY sizeX] = size(binaryimgNeurite);
% Cut of filter neuronal image
if(optionHandler.FilterMulPa~=1 && optionHandler.CutOutSphereCore == 1)
    neuritePic = CutOutCircles(binaryimgNeurite,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM, 0);
end
% Cut of filter neuronal image
if(optionHandler.FilterMulPa~=1 && optionHandler.CutOutSphereCore == 1)
    nucleusPic = CutOutCircles(nucleusPicBinary,selectedWell,1,0,0,optionHandler, foldername, sizeY, sizeX, NucleusM,0, 0);
end
% Cut of filter neuronal image
nucleusPicBinary = nucleusPicBinary(yStartPos:yEndPos,xStartPos:xEndPos);
[sizeYsmall sizeXsmall] = size(nucleusPicBinary);
% Cut of filter neuronal image
binaryimgNeurite = binaryimgNeurite(yStartPos:yEndPos,xStartPos:xEndPos);
if(AxonMeasurement==1);
binaryimgNeuriteHard = binaryimgNeuriteHard(yStartPos:yEndPos,xStartPos:xEndPos);
end
binaryimgNeurite=logical(binaryimgNeurite);
binaryimgNeuriteOrig = binaryimgNeurite;
  
%Images are saved now. In case a filter was chosen, only the area
        %of the filter is saved
        Dummy = zeros(sizeYOrig,sizeXOrig);
        Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = SkeletonImage;
        SkeletonImage = Dummy;
        newfile = [selectedWellLong 'Skeleton' optionHandler.FileEnding];
        imwrite(SkeletonImage,[foldername '/ConvertedCellomics/' newfile] );
        Dummy = zeros(sizeYOrig,sizeXOrig);
        Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = binaryimgNeurite;
        binaryimgNeurite = Dummy;
        newfile = [selectedWellLong 'Binary' optionHandler.FileEnding];
        imwrite(binaryimgNeurite,[foldername '/ConvertedCellomics/' newfile] );
        Dummy = zeros(sizeYOrig,sizeXOrig);
        Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = NeuronM;
        NeuronM = Dummy;
        summary=0;


% --------------------------------------------------------------------
function LUHMES_Morphology_Callback(hObject, eventdata, handles)
% hObject    handle to LUHMES_Morphology (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
neuronHandler = CSVCoordinates();
neuronHandler.ManualNeuronPositions = containers.Map();
csvHandler = CSVCoordinates();
 imageHandler = ImageHandler();
 optionHandler = Option();
 optionHandler.MinSizeNeuriteArea = 25000;
 handles.OptionHandler = optionHandler;
optionHandler.FileEnding = '.png';
[wellList folder ScalingFactor csvHandler  Summary SummaryName Summary2 Summary2Name Summary3Name Summary3 Summary4 SummaryEnd SummaryEndTop Name Stats] = LUHMES_Morphology(csvHandler);
optionHandler.ScalingFactor = ScalingFactor;
 imageHandler.Foldername = folder;
set(handles.lbWell, 'String', wellList);
handles.CSVCoordinates = csvHandler;
handles.NeuronCoordinates = neuronHandler;
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles);
Name2 = '.xlsx';
File = [Name Name2];
if(Stats==0)
xlswrite(File,Summary); 
xlswrite(File,SummaryName);
     
sheet = 2;
     
xlswrite(File,Summary2,sheet); 
xlswrite(File,SummaryEnd, sheet);
xlswrite(File,SummaryEndTop, sheet);
sheet = 3;
     
     xlswrite(File,Summary3,sheet); 
     xlswrite(File,SummaryEnd, sheet);
     xlswrite(File,SummaryEndTop, sheet);
     
     sheet = 4;
     
     xlswrite(File,Summary4,sheet); 
     xlswrite(File,SummaryEnd, sheet);
     xlswrite(File,SummaryEndTop, sheet);
else
    xlswrite(File,Summary2);
    xlswrite(File,Summary);
    xlswrite(File,SummaryName);
end    


% --------------------------------------------------------------------
function Ryr_mTor_FKBP12_CorticalNeurons_Callback(hObject, eventdata, handles)
% hObject    handle to Ryr_mTor_FKBP12_CorticalNeurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
neuronHandler = CSVCoordinates();
csvHandler = CSVCoordinates();
 imageHandler = ImageHandler();
 optionHandler = Option();
 optionHandler.MinSizeNeuriteArea = 25000;
 handles.OptionHandler = optionHandler;
optionHandler.FileEnding = '.png';
[wellList folder ScalingFactor csvHandler neuronHandler] = RainbowCells(csvHandler, neuronHandler);
optionHandler.ScalingFactor = ScalingFactor;
 imageHandler.Foldername = folder;
set(handles.lbWell, 'String', wellList);
handles.CSVCoordinates = csvHandler;
handles.NeuronCoordinates = neuronHandler;
%NeuronPositionsSkeletonization = neuronHandler.NeuronPositionsSkeletonization(wellname);
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles);


% --------------------------------------------------------------------
function LoadCaParticlesBrain_Callback(hObject, eventdata, handles)
% hObject    handle to LoadCaParticlesBrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
neuronHandler = CSVCoordinates();
neuronHandler.ManualNeuronPositions = containers.Map();
csvHandler = CSVCoordinates();
 imageHandler = ImageHandler();
 optionHandler = Option();
 handles.OptionHandler = optionHandler;
optionHandler.FileEnding = '.png';
[wellList folder ScalingFactor csvHandler] = CaOnly(csvHandler);
optionHandler.ScalingFactor = ScalingFactor;
 imageHandler.Foldername = folder;
set(handles.lbWell, 'String', wellList);
handles.CSVCoordinates = csvHandler;
handles.NeuronCoordinates = neuronHandler;
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles);



% --------------------------------------------------------------------
function AnalyzeCaBrain_Callback(hObject, eventdata, handles)
% hObject    handle to AnalyzeCaBrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to AnalyzeCaBrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
%foldername = uigetdir;
%Call EdgeFill Neurons and Skeleton Neurons and play with parameters in
%Option Handler
imageHandler = handles.ImageHandler;
foldername = imageHandler.Foldername;
ExecuteOptionExperimentCaOnly(handles,1,'Neurite');
guidata(handles.figure1, handles);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
[summary] = MeasurementCaOnly (foldername, csvHandler)
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function ExecuteOptionExperimentCaOnly(handles,onlyCompFill,channelName)
h = waitbar(0.5,'Neuronal morpholgy analyzed. This could take a few hours. Please be patient.', 'WindowStyle', 'modal');
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
optionHandler = handles.OptionHandler;

if(strcmp(channelName,'Oligo'))
    competingHandler = handles.NeuronCoordinates;
    if(isfield(handles,'OligoCoordinates'))
       neuronHandler = handles.OligoCoordinates;
    else
       neuronHandler = CSVCoordinates();
    end
else
    neuronHandler = handles.NeuronCoordinates;
    competingHandler = 0;
end
wellList = get(handles.lbWell, 'string');
foldername = imageHandler.Foldername;
[sizeY sizeX] = size(imageHandler.NucleusImage);
%Select STARTWELL

if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
if(numel(neuronHandler.AreaDictionary) ==0)
    neuronHandler.AreaDictionary = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
    neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFill) ==0)
    neuronHandler.NeuronPositionsEdgeFill = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFillSecond) ==0)
    neuronHandler.NeuronPositionsEdgeFillSecond = containers.Map();
end
if(numel(neuronHandler.NeuronStatMatrix) ==0)
    neuronHandler.NeuronStatMatrix = containers.Map();
end

%Create Matlab job
%Leave one Core out to ensure that GUI isn't freezing
maxCores =optionHandler.NumberOfCores;
if(maxCores > 1)
    c = parcluster();
    c.NumWorkers = maxCores;
    job = createJob(c,'JobData',optionHandler);
end

for(j=optionHandler.StartWell:numel(wellList))
    selectedWell=wellList(j);
    selectedWell=selectedWell{1};
    selectedWellLong=selectedWell;
    if(length(selectedWell) == 2)
              selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    end
    NucleusM = csvHandler.CellPosMatrix(selectedWell);    
    if(size(NucleusM,2) > sizeX || size(NucleusM,1) > sizeY)
        NucleusM = NucleusM(1:sizeY,1:sizeX);
        csvHandler.CellPosMatrix(selectedWell) = NucleusM;
        handles.CSVCoordinates = csvHandler;
        guidata(handles.figure1, handles);
    end
    if(maxCores > 1)
        %Check if watershed images are available. If not create them before
        %processing as ImageJ operation can't be done in parallel.
        %Check if there is already a watershed image. If not, create it.
        imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
        if(~exist(imagePathWatershed,'file'))
            WatershedNuclei(handles);
        end
        createTask(job, @ExecuteCaOnly, 4, {selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler});
    else
       [neuronHandler,  summary] = ExecuteCaOnly(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler);
    end
end
if(maxCores > 1)
    submit(job);
    wait(job);
    %job = out(2);
    out=fetchOutputs(job);
    %ToDo: Merge NeuronHandler
    TP = sum(cell2mat(out(:,2)));
    FP = sum(cell2mat(out(:,3)));
    TN = sum(cell2mat(out(:,4)));
end
% TPPerWell = cell2mat(out(:,2));
% FPPerWell = cell2mat(out(:,3));
% TNPerWell = cell2mat(out(:,4));
% TPPerWell = TPPerWell ./ (TPPerWell+FPPerWell);
% FPPerWell = FPPerWell ./ (TPPerWell+FPPerWell);
i=0;

if(maxCores > 1)
    for(j=optionHandler.StartWell:numel(wellList))
     %Get current NeuronHandler result and write it in main Neuronhandler
        i=i+1;
        selectedWell=wellList(j);
        selectedWell=selectedWell{1};    
    
        currentNeuronHandler = out(i,1);
        currentNeuronHandler=currentNeuronHandler{1};
        %Merge Neuron positions as well as neurite length matrix.
        neuronHandler.NeuronPositionsSkeletonization(selectedWell) = currentNeuronHandler.NeuronPositionsSkeletonization(selectedWell);
        neuronHandler.AreaDictionary(selectedWell) = currentNeuronHandler.AreaDictionary(selectedWell);
        neuronHandler.NeuronPositionsEdgeFill(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFill(selectedWell);
        neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
        neuronHandler.NeuronStatMatrix(selectedWell) = currentNeuronHandler.NeuronStatMatrix(selectedWell);
        neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = currentNeuronHandler.NeuronPositionsSkelDeleted(selectedWell);
        %neuronHandler.NeuriteLengthMatrix(selectedWell) = currentNeuronHandler.NeuriteLengthMatrix;
    end
end
close(h);
if(strcmp(channelName,'Oligo'))
     handles.OligoCoordinates = neuronHandler;
else
    handles.NeuronCoordinates = neuronHandler;
end
guidata(handles.figure1, handles);


 function [neuronHandler, summary] = ExecuteCaOnly(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler)    
        %Get first filter from selected Well
        selectedWellLong=selectedWell;
        if(length(selectedWell) == 2)
                  selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
        end
        imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
        nucleusImage = imread(imagePathNucleusBig);
        [sizeYorig sizeXorig] = size(nucleusImage);
        firstY=1;
        firstX=1;
        lastY=sizeYorig;
        lastX=sizeXorig;
        sizeY = sizeYorig;
        sizeX = sizeXorig;
         if(optionHandler.FilterMulPa==1)
             path = [foldername '/ConvertedCellomics/' selectedWell '.mat'];
             TP=0;
             FP=0;
             TN=0;
             if(exist(path,'file')) 
                 load(path);
                 %For every positive Filter
                 filterCount = numel(FMask.PositiveFilters);
                 for(i=1:filterCount);
                    f1 = FMask.PositiveFilters{i};
                    firstInd=find(f1,1, 'first');
                    lastInd=find(f1,1, 'last');
                    [firstY firstX] = ind2sub([sizeY sizeX], firstInd);
                    [lastY lastX] = ind2sub([sizeY sizeX], lastInd);
                    [neuronHandler summary] = MorphologyCaOnly(selectedWell,firstY,lastY, firstX, lastX,0,optionHandler,neuronHandler,NucleusM,foldername,channelName);                     
                end
             end
         else
             if(lastY - firstY > 10000 || lastX - firstX > 10000)        
                sizeXsmall = (lastX-firstX)/2;
                sizeYsmall = (lastY-firstY)/2;
                [neuronHandler summary] = MorphologyCaOnly(selectedWell,firstY,sizeYsmall,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyCaOnly(selectedWell,firstY,sizeYsmall,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyCaOnly(selectedWell,sizeYsmall+1,lastY,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyCaOnly(selectedWell,sizeYsmall+1,lastY,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
              else
                [neuronHandler summary] = MorphologyCaOnly(selectedWell,firstY,lastY, firstX, lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);  
              end
              
         end    
         
        % system('taskkill /fi "WINDOWTITLE eq FBird"');
        %Calculate only for filter 
        
function [neuronHandler summary] = MorphologyCaOnly(selectedWell, yStartPos, yEndPos, xStartPos, xEndPos, save, optionHandler, neuronHandler, NucleusM, foldername, channelName)
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end
imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
nucleusImage = imread(imagePathNucleusBig);


% Load filter if xStart is unequal 1
if(xStartPos > 1)
filterNumber = selectedWellLong(1,3);
filterNumber = str2num(filterNumber);
path = [foldername '/ConvertedCellomics/' selectedWellLong '.mat'];
if(exist(path,'file'))
        load(path);
        imageHandler.PositiveFilters = FMask.PositiveFilters;
        imageHandler.NegativeFilters = FMask.NegativeFilters;
end
FMask=0;

    %Get Filter
    filterList = imageHandler.PositiveFilters;
    filterNumber = size(filterList);
    filter = filterNumber(2);
    filterMasks = zeros(size(nucleusImage));
    for(u=1:filter)
    filterNumber = u;    
    filterHandle = filterList(filterNumber);
    filterHandle = filterHandle{1};
    filterMask = (filterHandle);
    %Remove excluded parts of picture
    if(iscell(imageHandler.NegativeFilters) > 0)
        negativeMask = imageHandler.NegativeFilters;
        for i=1:numel(negativeMask)
          workMask = negativeMask{i};
          workMask = int8(workMask -1);
          workMask = logical(workMask .* (-1));
          workMask = sparse(logical(workMask));
          filterMask = (logical((filterMask .* workMask)));
        end
    end
    filterMasks = filterMasks + filterMask;
    end
    workMask=0;
    filterMask = filterMasks;
    %Save mask:
    f{1} = selectedWellLong;
    f{2} ='Mask.png';
    newFileM = [f{1} f{2}];
    imwrite(filterMask,[foldername '/ConvertedCellomics/' newFileM]);
    %Multiply mask with image
    nucleusImage = uint8(nucleusImage).*uint8(full(filterMask));
    % With the help of the filterMask we can now determine the closest
    % sourround box coordinates
    Coordinates = regionprops(filterMask,'BoundingBox');
    xStartPos = Coordinates(1).BoundingBox(1);
    xEndPos = xStartPos + Coordinates(1).BoundingBox(3);
    yStartPos = Coordinates(1).BoundingBox(2);
    yEndPos = yStartPos + Coordinates(1).BoundingBox(4);
end






[sizeYorig sizeXorig] = size(nucleusImage);
%yStartPos = 1;
%yEndPos= sizeYorig;
%xStartPos = 1;
%xEndPos = sizeXorig;
NucleusM = NucleusM(yStartPos:yEndPos,xStartPos:xEndPos);

%[D IDX] = bwdist(full(NucleusM));
TN=0;
if(isfield(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
    NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
    NeuronManualM = NeuronManualM(yStartPos:yEndPos,xStartPos:xEndPos);
    TN = nnz(NeuronManualM);
else
    NeuronManualM = sparse(logical(zeros(sizeYorig,sizeXorig)));
    if(isa(neuronHandler.ManualNeuronPositionsSparse,'containers.Map') == 0)
        neuronHandler.ManualNeuronPositionsSparse = containers.Map();
    end
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
end
%Load neurite and nucleus image from loadcaparticlesbrain 
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end

% From here onwards starts the actual algorithm!
%These are for now set fixed parameters, but can later be included into
%Option menue


% Image locations:
imagePathBinary = [foldername '/ConvertedCellomics/' selectedWellLong 'Binary' optionHandler.FileEnding];
binaryimgNeurite = imread(imagePathBinary);
if(xStartPos > 1)
binaryimgNeurite = uint8(binaryimgNeurite).*uint8(full(filterMask));
end


sizeY = yEndPos - yStartPos+1;
sizeX = xEndPos - xStartPos+1;
[sizeYOrig sizeXOrig]= size(binaryimgNeurite);

binaryimgNeurite = binaryimgNeurite(yStartPos:yEndPos,xStartPos:xEndPos);
Dummy = zeros(sizeYOrig,sizeXOrig);
Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = binaryimgNeurite;
binaryimgNeurite = Dummy;
newfile = [selectedWellLong 'Binary' optionHandler.FileEnding];
imwrite(binaryimgNeurite,[foldername '/ConvertedCellomics/' newfile] );

%Check if NeuronStatMatrix exists
if(~isa(neuronHandler.NeuronStatMatrix,'containers.Map'))
    neuronHandler.NeuronStatMatrix = containers.Map();
end
if(~isa(neuronHandler.OmniNucleusPositions,'containers.Map'))
    neuronHandler.OmniNucleusPositions = containers.Map();
end
currentNeuronStat = containers.Map();


        
        NeuronM = NucleusM;
       
        Dummy = zeros(sizeYorig ,sizeXorig);
        Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = NeuronM;
        NeuronM = Dummy;
        summary=0;





neuronHandler.NeuronPositionsEdgeFill(selectedWell) = NeuronM;


% --------------------------------------------------------------------
function LoadMicrogliaMorphology_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMicrogliaMorphology (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
neuronHandler = CSVCoordinates();
neuronHandler.ManualNeuronPositions = containers.Map();
csvHandler = CSVCoordinates();
 imageHandler = ImageHandler();
 optionHandler = Option();
 handles.OptionHandler = optionHandler;
optionHandler.FileEnding = '.png';
[wellList folder ScalingFactor csvHandler] = MicrogliaMorphology(csvHandler);
optionHandler.ScalingFactor = ScalingFactor;
 imageHandler.Foldername = folder;
set(handles.lbWell, 'String', wellList);
handles.CSVCoordinates = csvHandler;
handles.NeuronCoordinates = neuronHandler;
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles);


% --------------------------------------------------------------------
function AnalyzeMicrogliaActivation_Callback(hObject, eventdata, handles)
% hObject    handle to AnalyzeMicrogliaActivation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
expname = inputdlg('Experiment name');
Question = questdlg('Are you analyzing image stacks?')
if(strcmp(Question,'No'))
    CD68 = 0;
else
    CD68 = 1;
end
Question = questdlg('Is this a manual validation')
if(strcmp(Question,'No'))
    RAND = 0;
else
    RAND = 1;
end
token = expname{1};
a{1} = token;
a{2} = '.xlsx';
ExcelName = [a{1} a{2}];
%foldername = uigetdir;
%Call EdgeFill Neurons and Skeleton Neurons and play with parameters in
%Option Handler
imageHandler = handles.ImageHandler;
foldername = imageHandler.Foldername;
ExecuteAnalyzeMicrogliaActivation(handles,1,'Neurite');
guidata(handles.figure1, handles);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
[summary] = MeasureMicroGliaActivation (foldername, csvHandler, ExcelName, CD68, RAND)
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function ExecuteAnalyzeMicrogliaActivation(handles,onlyCompFill,channelName)
h = waitbar(0.5,'Neuronal morpholgy analyzed. This could take a few hours. Please be patient.', 'WindowStyle', 'modal');
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
optionHandler = handles.OptionHandler;

if(strcmp(channelName,'Oligo'))
    competingHandler = handles.NeuronCoordinates;
    if(isfield(handles,'OligoCoordinates'))
       neuronHandler = handles.OligoCoordinates;
    else
       neuronHandler = CSVCoordinates();
    end
else
    neuronHandler = handles.NeuronCoordinates;
    competingHandler = 0;
end
wellList = get(handles.lbWell, 'string');
foldername = imageHandler.Foldername;
[sizeY sizeX] = size(imageHandler.NucleusImage);
%Select STARTWELL

if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
if(numel(neuronHandler.AreaDictionary) ==0)
    neuronHandler.AreaDictionary = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
    neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFill) ==0)
    neuronHandler.NeuronPositionsEdgeFill = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFillSecond) ==0)
    neuronHandler.NeuronPositionsEdgeFillSecond = containers.Map();
end
if(numel(neuronHandler.NeuronStatMatrix) ==0)
    neuronHandler.NeuronStatMatrix = containers.Map();
end

%Create Matlab job
%Leave one Core out to ensure that GUI isn't freezing
maxCores =optionHandler.NumberOfCores;
if(maxCores > 1)
    c = parcluster();
    c.NumWorkers = maxCores;
    job = createJob(c,'JobData',optionHandler);
end

for(j=optionHandler.StartWell:numel(wellList))
    selectedWell=wellList(j);
    selectedWell=selectedWell{1};
    selectedWellLong=selectedWell;
    if(length(selectedWell) == 2)
              selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    end
    NucleusM = csvHandler.CellPosMatrix(selectedWell);    
    if(size(NucleusM,2) > sizeX || size(NucleusM,1) > sizeY)
        NucleusM = NucleusM(1:sizeY,1:sizeX);
        csvHandler.CellPosMatrix(selectedWell) = NucleusM;
        handles.CSVCoordinates = csvHandler;
        guidata(handles.figure1, handles);
    end
    if(maxCores > 1)
        %Check if watershed images are available. If not create them before
        %processing as ImageJ operation can't be done in parallel.
        %Check if there is already a watershed image. If not, create it.
        imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
        if(~exist(imagePathWatershed,'file'))
            WatershedNuclei(handles);
        end
        createTask(job, @ExecuteMeasureMicroGliaActivation, 4, {selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler});
    else
       [neuronHandler,  summary] = ExecuteMeasureMicroGliaActivation(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler);
    end
end
if(maxCores > 1)
    submit(job);
    wait(job);
    %job = out(2);
    out=fetchOutputs(job);
    %ToDo: Merge NeuronHandler
    TP = sum(cell2mat(out(:,2)));
    FP = sum(cell2mat(out(:,3)));
    TN = sum(cell2mat(out(:,4)));
end
% TPPerWell = cell2mat(out(:,2));
% FPPerWell = cell2mat(out(:,3));
% TNPerWell = cell2mat(out(:,4));
% TPPerWell = TPPerWell ./ (TPPerWell+FPPerWell);
% FPPerWell = FPPerWell ./ (TPPerWell+FPPerWell);
i=0;

if(maxCores > 1)
    for(j=optionHandler.StartWell:numel(wellList))
     %Get current NeuronHandler result and write it in main Neuronhandler
        i=i+1;
        selectedWell=wellList(j);
        selectedWell=selectedWell{1};    
    
        currentNeuronHandler = out(i,1);
        currentNeuronHandler=currentNeuronHandler{1};
        %Merge Neuron positions as well as neurite length matrix.
        neuronHandler.NeuronPositionsSkeletonization(selectedWell) = currentNeuronHandler.NeuronPositionsSkeletonization(selectedWell);
        neuronHandler.AreaDictionary(selectedWell) = currentNeuronHandler.AreaDictionary(selectedWell);
        neuronHandler.NeuronPositionsEdgeFill(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFill(selectedWell);
        neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
        neuronHandler.NeuronStatMatrix(selectedWell) = currentNeuronHandler.NeuronStatMatrix(selectedWell);
        neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = currentNeuronHandler.NeuronPositionsSkelDeleted(selectedWell);
        %neuronHandler.NeuriteLengthMatrix(selectedWell) = currentNeuronHandler.NeuriteLengthMatrix;
    end
end
close(h);
if(strcmp(channelName,'Oligo'))
     handles.OligoCoordinates = neuronHandler;
else
    handles.NeuronCoordinates = neuronHandler;
end
guidata(handles.figure1, handles);


 function [neuronHandler, summary] = ExecuteMeasureMicroGliaActivation(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler)    
        %Get first filter from selected Well
        selectedWellLong=selectedWell;
        if(length(selectedWell) == 2)
                  selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
        end
        imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
        nucleusImage = imread(imagePathNucleusBig);
        [sizeYorig sizeXorig] = size(nucleusImage);
        firstY=1;
        firstX=1;
        lastY=sizeYorig;
        lastX=sizeXorig;
        sizeY = sizeYorig;
        sizeX = sizeXorig;
         if(optionHandler.FilterMulPa==1)
             path = [foldername '/ConvertedCellomics/' selectedWell '.mat'];
             TP=0;
             FP=0;
             TN=0;
             if(exist(path,'file')) 
                 load(path);
                 %For every positive Filter
                 filterCount = numel(FMask.PositiveFilters);
                 for(i=1:filterCount);
                    f1 = FMask.PositiveFilters{i};
                    firstInd=find(f1,1, 'first');
                    lastInd=find(f1,1, 'last');
                    [firstY firstX] = ind2sub([sizeY sizeX], firstInd);
                    [lastY lastX] = ind2sub([sizeY sizeX], lastInd);
                    [neuronHandler summary] = MorphologyMicroglia(selectedWell,firstY,lastY, firstX, lastX,0,optionHandler,neuronHandler,NucleusM,foldername,channelName);                     
                end
             end
         else
             if(lastY - firstY > 10000 || lastX - firstX > 10000)        
                sizeXsmall = (lastX-firstX)/2;
                sizeYsmall = (lastY-firstY)/2;
                [neuronHandler summary] = MorphologyMicroglia(selectedWell,firstY,sizeYsmall,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyMicroglia(selectedWell,firstY,sizeYsmall,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyMicroglia(selectedWell,sizeYsmall+1,lastY,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyMicroglia(selectedWell,sizeYsmall+1,lastY,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
              else
                [neuronHandler summary] = MorphologyMicroglia(selectedWell,firstY,lastY, firstX, lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);  
              end
              
         end    
         
        % system('taskkill /fi "WINDOWTITLE eq FBird"');
        %Calculate only for filter 
        
function [neuronHandler summary] = MorphologyMicroglia(selectedWell, yStartPos, yEndPos, xStartPos, xEndPos, save, optionHandler, neuronHandler, NucleusM, foldername, channelName)
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end
imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
nucleusImage = imread(imagePathNucleusBig);


% Load filter if xStart is unequal 1
if(xStartPos > 1)
filterNumber = selectedWellLong(1,3);
filterNumber = str2num(filterNumber);
path = [foldername '/ConvertedCellomics/' selectedWellLong '.mat'];
if(exist(path,'file'))
        load(path);
        imageHandler.PositiveFilters = FMask.PositiveFilters;
        imageHandler.NegativeFilters = FMask.NegativeFilters;
end
FMask=0;

    %Get Filter
    filterList = imageHandler.PositiveFilters;
    filterNumber = size(filterList);
    filter = filterNumber(2);
    filterMasks = zeros(size(nucleusImage));
    for(u=1:filter)
    filterNumber = u;    
    filterHandle = filterList(filterNumber);
    filterHandle = filterHandle{1};
    filterMask = (filterHandle);
    %Remove excluded parts of picture
    if(iscell(imageHandler.NegativeFilters) > 0)
        negativeMask = imageHandler.NegativeFilters;
        for i=1:numel(negativeMask)
          workMask = negativeMask{i};
          workMask = int8(workMask -1);
          workMask = logical(workMask .* (-1));
          workMask = sparse(logical(workMask));
          filterMask = (logical((filterMask .* workMask)));
        end
    end
    filterMasks = filterMasks + filterMask;
    end
    workMask=0;
    filterMask = filterMasks;
    %Save mask:
    f{1} = selectedWellLong;
    f{2} ='Mask.png';
    newFileM = [f{1} f{2}];
    imwrite(filterMask,[foldername '/ConvertedCellomics/' newFileM]);
    %Multiply mask with image
    nucleusImage = uint8(nucleusImage).*uint8(full(filterMask));
    % With the help of the filterMask we can now determine the closest
    % sourround box coordinates
    Coordinates = regionprops(filterMask,'BoundingBox');
    xStartPos = Coordinates(1).BoundingBox(1);
    xEndPos = xStartPos + Coordinates(1).BoundingBox(3);
    yStartPos = Coordinates(1).BoundingBox(2);
    yEndPos = yStartPos + Coordinates(1).BoundingBox(4);
else
    filterMask = logical(ones(size(nucleusImage)));
    %Save mask:
    f{1} = selectedWellLong;
    f{2} ='Mask.png';
    newFileM = [f{1} f{2}];
    imwrite(filterMask,[foldername '/ConvertedCellomics/' newFileM]);
     %Multiply mask with image
    nucleusImage = uint8(nucleusImage).*uint8(full(filterMask));
    % With the help of the filterMask we can now determine the closest
    % sourround box coordinates
    Coordinates = regionprops(filterMask,'BoundingBox');
    xStartPos = Coordinates(1).BoundingBox(1);
    xEndPos = xStartPos + Coordinates(1).BoundingBox(3);
    yStartPos = Coordinates(1).BoundingBox(2);
    yEndPos = yStartPos + Coordinates(1).BoundingBox(4);
end






[sizeYorig sizeXorig] = size(nucleusImage);
%yStartPos = 1;
%yEndPos= sizeYorig;
%xStartPos = 1;
%xEndPos = sizeXorig;
%Rounding of Coordinates gives errors, if sizes are bigger than image,
%therefore we will automatically set them to max values:
if(yStartPos<1)
    yStartPos = 1;
end
if(yEndPos>sizeYorig)
    yEndPos = sizeYorig;
end
if(xStartPos<1)
    xStartPos = 1;
end
if(xEndPos>sizeXorig)
    xEndPos = sizeXorig;
end
NucleusM = NucleusM(yStartPos:yEndPos,xStartPos:xEndPos);

%[D IDX] = bwdist(full(NucleusM));
TN=0;
if(isfield(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
    NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
    NeuronManualM = NeuronManualM(yStartPos:yEndPos,xStartPos:xEndPos);
    TN = nnz(NeuronManualM);
else
    NeuronManualM = sparse(logical(zeros(sizeYorig,sizeXorig)));
    if(isa(neuronHandler.ManualNeuronPositionsSparse,'containers.Map') == 0)
        neuronHandler.ManualNeuronPositionsSparse = containers.Map();
    end
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
end
%Load neurite and nucleus image from loadcaparticlesbrain 
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end

% From here onwards starts the actual algorithm!
%These are for now set fixed parameters, but can later be included into
%Option menue


% Image locations:
imagePathBinary = [foldername '/ConvertedCellomics/' selectedWellLong 'Binary' optionHandler.FileEnding];
imagePathNucleusBinary = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
binaryimgNeurite = imread(imagePathBinary);
BinaryNucleus = imread(imagePathNucleusBinary);
if(xStartPos > 1)
binaryimgNeurite = uint8(binaryimgNeurite).*uint8(full(filterMask));
BinaryNucleus = uint8(BinaryNucleus).*uint8(full(filterMask));
end


sizeY = yEndPos - yStartPos+1;
sizeX = xEndPos - xStartPos+1;
[sizeYOrig sizeXOrig]= size(binaryimgNeurite);

binaryimgNeurite = binaryimgNeurite(yStartPos:yEndPos,xStartPos:xEndPos);
BinaryNucleus = BinaryNucleus(yStartPos:yEndPos,xStartPos:xEndPos);
Dummy = zeros(sizeYOrig,sizeXOrig);
Dummy2 = zeros(sizeYOrig,sizeXOrig);
Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = binaryimgNeurite;
Dummy2(yStartPos:yEndPos,xStartPos:xEndPos) = BinaryNucleus;
binaryimgNeurite = Dummy;
BinaryNucleus = Dummy2;
newfile = [selectedWellLong 'Binary' optionHandler.FileEnding];
imwrite(binaryimgNeurite,[foldername '/ConvertedCellomics/' newfile] );
newfile2 = [selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
imwrite(BinaryNucleus,[foldername '/ConvertedCellomics/' newfile2] );
%Check if NeuronStatMatrix exists
if(~isa(neuronHandler.NeuronStatMatrix,'containers.Map'))
    neuronHandler.NeuronStatMatrix = containers.Map();
end
if(~isa(neuronHandler.OmniNucleusPositions,'containers.Map'))
    neuronHandler.OmniNucleusPositions = containers.Map();
end
currentNeuronStat = containers.Map();


        
        NeuronM = NucleusM;
       
        Dummy = zeros(sizeYorig ,sizeXorig);
        Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = NeuronM;
        NeuronM = Dummy;
        summary=0;





neuronHandler.NeuronPositionsEdgeFill(selectedWell) = NeuronM;







% --------------------------------------------------------------------
function ExportMatrices_Callback(hObject, eventdata, handles)
% hObject    handle to ExportMatrices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
imageHandler = handles.ImageHandler;
foldername = imageHandler.Foldername;
guidata(handles.figure1, handles);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
[check] = ExportMatrices (foldername, csvHandler, neuronHandler)
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end


% --------------------------------------------------------------------
function LoadGanglia_Callback(hObject, eventdata, handles)
% hObject    handle to LoadGanglia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
neuronHandler = CSVCoordinates();
neuronHandler.ManualNeuronPositions = containers.Map();
csvHandler = CSVCoordinates();
 imageHandler = ImageHandler();
 optionHandler = Option();
 handles.OptionHandler = optionHandler;
optionHandler.FileEnding = '.png';
[wellList folder ScalingFactor csvHandler] = GangliaImport(csvHandler);

optionHandler.ScalingFactor = ScalingFactor;
 imageHandler.Foldername = folder;
set(handles.lbWell, 'String', wellList);
handles.CSVCoordinates = csvHandler;
handles.NeuronCoordinates = neuronHandler;
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles);
[token,remain] = strtok(folder, '\');
lengthRemain = length(remain);
while(length(remain)>0)
[token,remain] = strtok(remain, '\');
end
a{1} = token;
a{2} = '_Raw'
a{3}= '.mat';
FileName = [a{1} a{2} a{3}];
%Pathname is folder and a '\'
a{1} = folder;
a{2} = '\';
PathName = [a{1} a{2}];
SaveFixedNeuronPositions(FileName,PathName,handles);


% --------------------------------------------------------------------
function AnalyzeGanglia_Callback(hObject, eventdata, handles)
% hObject    handle to AnalyzeGanglia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles = guidata(handles.figure1);
optionHandler = handles.OptionHandler;
expname = inputdlg('Experiment name');
token = expname{1};
a{1} = token;
a{2} = '.xlsx';
ExcelName = [a{1} a{2}];
%foldername = uigetdir;
%Call EdgeFill Neurons and Skeleton Neurons and play with parameters in
%Option Handler
Question = questdlg('Are you analyzing macrophages?')
if(strcmp(Question,'No'))
    Macrophages = 0;
else
    Macrophages = 1;
end
imageHandler = handles.ImageHandler;
foldername = imageHandler.Foldername;
ExecuteAnalyzeGanglia_Callback(handles,1,'Neurite');
guidata(handles.figure1, handles);
csvHandler = handles.CSVCoordinates;
neuronHandler = handles.NeuronCoordinates;
[summary] = MeasureGanglia (foldername, csvHandler, ExcelName, Macrophages)
[token,remain] = strtok(foldername, '\');
lengthRemain = length(remain);
while(length(remain)>0)
[token,remain] = strtok(remain, '\');
end
a{1} = token;
a{2} = '_Analyzed'
a{3}= '.mat';
FileName = [a{1} a{2} a{3}];
%Pathname is folder and a '\'
a{1} = foldername;
a{2} = '\';
PathName = [a{1} a{2}];
SaveFixedNeuronPositions(FileName,PathName,handles);
catch errorObj
    % If there is a problem, we display the error message
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
end

function ExecuteAnalyzeGanglia_Callback(handles,onlyCompFill,channelName)
h = waitbar(0.5,'Neuronal morpholgy analyzed. This could take a few hours. Please be patient.', 'WindowStyle', 'modal');
handles = guidata(handles.figure1);
imageHandler = handles.ImageHandler;
csvHandler = handles.CSVCoordinates;
optionHandler = handles.OptionHandler;

if(strcmp(channelName,'Oligo'))
    competingHandler = handles.NeuronCoordinates;
    if(isfield(handles,'OligoCoordinates'))
       neuronHandler = handles.OligoCoordinates;
    else
       neuronHandler = CSVCoordinates();
    end
else
    neuronHandler = handles.NeuronCoordinates;
    competingHandler = 0;
end
wellList = get(handles.lbWell, 'string');
foldername = imageHandler.Foldername;
[sizeY sizeX] = size(imageHandler.NucleusImage);
%Select STARTWELL

if(numel(neuronHandler.NeuronPositionsSkeletonization) ==0)
    neuronHandler.NeuronPositionsSkeletonization = containers.Map();
end
if(numel(neuronHandler.AreaDictionary) ==0)
    neuronHandler.AreaDictionary = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsSkelDeleted) ==0)
    neuronHandler.NeuronPositionsSkelDeleted = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFill) ==0)
    neuronHandler.NeuronPositionsEdgeFill = containers.Map();
end
if(numel(neuronHandler.NeuronPositionsEdgeFillSecond) ==0)
    neuronHandler.NeuronPositionsEdgeFillSecond = containers.Map();
end
if(numel(neuronHandler.NeuronStatMatrix) ==0)
    neuronHandler.NeuronStatMatrix = containers.Map();
end

%Create Matlab job
%Leave one Core out to ensure that GUI isn't freezing
maxCores =optionHandler.NumberOfCores;
if(maxCores > 1)
    c = parcluster();
    c.NumWorkers = maxCores;
    job = createJob(c,'JobData',optionHandler);
end

for(j=optionHandler.StartWell:numel(wellList))
    selectedWell=wellList(j);
    selectedWell=selectedWell{1};
    selectedWellLong=selectedWell;
    if(length(selectedWell) == 2)
              selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
    end
    NucleusM = csvHandler.CellPosMatrix(selectedWell);    
    if(size(NucleusM,2) > sizeX || size(NucleusM,1) > sizeY)
        NucleusM = NucleusM(1:sizeY,1:sizeX);
        csvHandler.CellPosMatrix(selectedWell) = NucleusM;
        handles.CSVCoordinates = csvHandler;
        guidata(handles.figure1, handles);
    end
    if(maxCores > 1)
        %Check if watershed images are available. If not create them before
        %processing as ImageJ operation can't be done in parallel.
        %Check if there is already a watershed image. If not, create it.
        imagePathWatershed = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
        if(~exist(imagePathWatershed,'file'))
            WatershedNuclei(handles);
        end
        createTask(job, @ExecuteMeasureAnalyzeGanglia_Callback, 4, {selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler});
    else
       [neuronHandler,  summary] = ExecuteMeasureAnalyzeGanglia_Callback(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler);
    end
end
if(maxCores > 1)
    submit(job);
    wait(job);
    %job = out(2);
    out=fetchOutputs(job);
    %ToDo: Merge NeuronHandler
    TP = sum(cell2mat(out(:,2)));
    FP = sum(cell2mat(out(:,3)));
    TN = sum(cell2mat(out(:,4)));
end
% TPPerWell = cell2mat(out(:,2));
% FPPerWell = cell2mat(out(:,3));
% TNPerWell = cell2mat(out(:,4));
% TPPerWell = TPPerWell ./ (TPPerWell+FPPerWell);
% FPPerWell = FPPerWell ./ (TPPerWell+FPPerWell);
i=0;

if(maxCores > 1)
    for(j=optionHandler.StartWell:numel(wellList))
     %Get current NeuronHandler result and write it in main Neuronhandler
        i=i+1;
        selectedWell=wellList(j);
        selectedWell=selectedWell{1};    
    
        currentNeuronHandler = out(i,1);
        currentNeuronHandler=currentNeuronHandler{1};
        %Merge Neuron positions as well as neurite length matrix.
        neuronHandler.NeuronPositionsSkeletonization(selectedWell) = currentNeuronHandler.NeuronPositionsSkeletonization(selectedWell);
        neuronHandler.AreaDictionary(selectedWell) = currentNeuronHandler.AreaDictionary(selectedWell);
        neuronHandler.NeuronPositionsEdgeFill(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFill(selectedWell);
        neuronHandler.NeuronPositionsEdgeFillSecond(selectedWell) = currentNeuronHandler.NeuronPositionsEdgeFillSecond(selectedWell);
        neuronHandler.NeuronStatMatrix(selectedWell) = currentNeuronHandler.NeuronStatMatrix(selectedWell);
        neuronHandler.NeuronPositionsSkelDeleted(selectedWell) = currentNeuronHandler.NeuronPositionsSkelDeleted(selectedWell);
        %neuronHandler.NeuriteLengthMatrix(selectedWell) = currentNeuronHandler.NeuriteLengthMatrix;
    end
end
close(h);
if(strcmp(channelName,'Oligo'))
     handles.OligoCoordinates = neuronHandler;
else
    handles.NeuronCoordinates = neuronHandler;
end
guidata(handles.figure1, handles);


 function [neuronHandler, summary] = ExecuteMeasureAnalyzeGanglia_Callback(selectedWell,sizeY,sizeX,optionHandler,neuronHandler,NucleusM,foldername,onlyCompFill,channelName,competingHandler)    
        %Get first filter from selected Well
        selectedWellLong=selectedWell;
        if(length(selectedWell) == 2)
                  selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
        end
        imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
        nucleusImage = imread(imagePathNucleusBig);
        [sizeYorig sizeXorig] = size(nucleusImage);
        firstY=1;
        firstX=1;
        lastY=sizeYorig;
        lastX=sizeXorig;
        sizeY = sizeYorig;
        sizeX = sizeXorig;
         if(optionHandler.FilterMulPa==1)
             path = [foldername '/ConvertedCellomics/' selectedWell '.mat'];
             TP=0;
             FP=0;
             TN=0;
             if(exist(path,'file')) 
                 load(path);
                 %For every positive Filter
                 filterCount = numel(FMask.PositiveFilters);
                 for(i=1:filterCount);
                    f1 = FMask.PositiveFilters{i};
                    firstInd=find(f1,1, 'first');
                    lastInd=find(f1,1, 'last');
                    [firstY firstX] = ind2sub([sizeY sizeX], firstInd);
                    [lastY lastX] = ind2sub([sizeY sizeX], lastInd);
                    [neuronHandler summary] = MorphologyGanglia(selectedWell,firstY,lastY, firstX, lastX,0,optionHandler,neuronHandler,NucleusM,foldername,channelName);                     
                end
             end
         else
             if(lastY - firstY > 10000 || lastX - firstX > 10000)        
                sizeXsmall = (lastX-firstX)/2;
                sizeYsmall = (lastY-firstY)/2;
                [neuronHandler summary] = MorphologyGanglia(selectedWell,firstY,sizeYsmall,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyGanglia(selectedWell,firstY,sizeYsmall,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyGanglia(selectedWell,sizeYsmall+1,lastY,sizeXsmall+1,lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
                [neuronHandler summary] = MorphologyGanglia(selectedWell,sizeYsmall+1,lastY,firstX,sizeXsmall,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);
              else
                [neuronHandler summary] = MorphologyGanglia(selectedWell,firstY,lastY, firstX, lastX,1,optionHandler,neuronHandler,NucleusM,foldername,channelName);  
              end
              
         end    
         
        % system('taskkill /fi "WINDOWTITLE eq FBird"');
        %Calculate only for filter 
        
function [neuronHandler summary] = MorphologyGanglia(selectedWell, yStartPos, yEndPos, xStartPos, xEndPos, save, optionHandler, neuronHandler, NucleusM, foldername, channelName)
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end
imagePathNucleusBig = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
nucleusImage = imread(imagePathNucleusBig);


% Load filter if xStart is unequal 1
if(xStartPos > 1)
filterNumber = selectedWellLong(1,3);
filterNumber = str2num(filterNumber);
path = [foldername '/ConvertedCellomics/' selectedWellLong '.mat'];
if(exist(path,'file'))
        load(path);
        imageHandler.PositiveFilters = FMask.PositiveFilters;
        imageHandler.NegativeFilters = FMask.NegativeFilters;
end
FMask=0;

    %Get Filter
    filterList = imageHandler.PositiveFilters;
    filterNumber = size(filterList);
    filter = filterNumber(2);
    filterMasks = zeros(size(nucleusImage));
    for(u=1:filter)
    filterNumber = u;    
    filterHandle = filterList(filterNumber);
    filterHandle = filterHandle{1};
    filterMask = (filterHandle);
    %Remove excluded parts of picture
    if(iscell(imageHandler.NegativeFilters) > 0)
        negativeMask = imageHandler.NegativeFilters;
        for i=1:numel(negativeMask)
          workMask = negativeMask{i};
          workMask = int8(workMask -1);
          workMask = logical(workMask .* (-1));
          workMask = sparse(logical(workMask));
          filterMask = (logical((filterMask .* workMask)));
        end
    end
    filterMasks = filterMasks + filterMask;
    end
    workMask=0;
    filterMask = filterMasks;
    %Save mask:
    f{1} = selectedWellLong;
    f{2} ='Mask.png';
    newFileM = [f{1} f{2}];
    imwrite(filterMask,[foldername '/ConvertedCellomics/' newFileM]);
    %Multiply mask with image
    nucleusImage = uint8(nucleusImage).*uint8(full(filterMask));
    % With the help of the filterMask we can now determine the closest
    % sourround box coordinates
    Coordinates = regionprops(filterMask,'BoundingBox');
    xStartPos = Coordinates(1).BoundingBox(1);
    xEndPos = xStartPos + Coordinates(1).BoundingBox(3);
    yStartPos = Coordinates(1).BoundingBox(2);
    yEndPos = yStartPos + Coordinates(1).BoundingBox(4);
end






[sizeYorig sizeXorig] = size(nucleusImage);
%yStartPos = 1;
%yEndPos= sizeYorig;
%xStartPos = 1;
%xEndPos = sizeXorig;
%Rounding of Coordinates gives errors, if sizes are bigger than image,
%therefore we will automatically set them to max values:
if(yStartPos<1)
    yStartPos = 1;
end
if(yEndPos>sizeYorig)
    yEndPos = sizeYorig;
end
if(xStartPos<1)
    xStartPos = 1;
end
if(xEndPos>sizeXorig)
    xEndPos = sizeXorig;
end
NucleusM = NucleusM(yStartPos:yEndPos,xStartPos:xEndPos);

%[D IDX] = bwdist(full(NucleusM));
TN=0;
if(isfield(neuronHandler.ManualNeuronPositionsSparse,selectedWell))
    NeuronManualM = neuronHandler.ManualNeuronPositionsSparse(selectedWell);
    NeuronManualM = NeuronManualM(yStartPos:yEndPos,xStartPos:xEndPos);
    TN = nnz(NeuronManualM);
else
    NeuronManualM = sparse(logical(zeros(sizeYorig,sizeXorig)));
    if(isa(neuronHandler.ManualNeuronPositionsSparse,'containers.Map') == 0)
        neuronHandler.ManualNeuronPositionsSparse = containers.Map();
    end
    neuronHandler.ManualNeuronPositionsSparse(selectedWell) = NeuronManualM;
end
%Load neurite and nucleus image from loadcaparticlesbrain 
selectedWellLong=selectedWell;
if(length(selectedWell) == 2)
          selectedWellLong = [selectedWell(1) '0' selectedWell(2)];
end

% From here onwards starts the actual algorithm!
%These are for now set fixed parameters, but can later be included into
%Option menue


% Image locations:
imagePathNeurite = [foldername '/ConvertedCellomics/' selectedWellLong 'NeuriteBig' optionHandler.FileEnding];
imagePathOligo = [foldername '/ConvertedCellomics/' selectedWellLong 'OligoBig' optionHandler.FileEnding];
imagePathNucleus = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBig' optionHandler.FileEnding];
imagePathNucleusBinary = [foldername '/ConvertedCellomics/' selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
imgNeurite = imread(imagePathNeurite);
imgOligo = imread(imagePathOligo);
imgNucleus = imread(imagePathNucleus);
BinaryNucleus = imread(imagePathNucleusBinary);
if(xStartPos > 1)
imgNeurite = uint8(imgNeurite).*uint8(full(filterMask));
imgNucleus = uint16(imgNucleus).*uint16(full(filterMask));
BinaryNucleus = uint8(BinaryNucleus).*uint8(full(filterMask));
imgOligo  = uint8(imgOligo ).*uint8(full(filterMask));
end


sizeY = yEndPos - yStartPos+1;
sizeX = xEndPos - xStartPos+1;
[sizeYOrig sizeXOrig]= size(imgNeurite);

imgNeurite = imgNeurite(yStartPos:yEndPos,xStartPos:xEndPos);
imgNucleus = imgNucleus(yStartPos:yEndPos,xStartPos:xEndPos);
BinaryNucleus = BinaryNucleus(yStartPos:yEndPos,xStartPos:xEndPos);
imgOligo = imgOligo(yStartPos:yEndPos,xStartPos:xEndPos);
Dummy = uint8(zeros(sizeYOrig,sizeXOrig));
Dummy2 = uint16(zeros(sizeYOrig,sizeXOrig));
Dummy3 = zeros(sizeYOrig,sizeXOrig);
Dummy4 = uint8(zeros(sizeYOrig,sizeXOrig));
Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = imgNeurite;
Dummy2(yStartPos:yEndPos,xStartPos:xEndPos) = imgNucleus;
Dummy3(yStartPos:yEndPos,xStartPos:xEndPos) = BinaryNucleus;
Dummy4(yStartPos:yEndPos,xStartPos:xEndPos) = imgOligo;
imgNeurite = Dummy;
imgNeuriteSmall = imresize(imgNeurite,0.1);
imgNucleus = Dummy2;
BinaryNucleus = Dummy3;
imgNucleusSmall = imresize(imgNucleus,0.1);
imgOligo = Dummy4;
imgOligoSmall = imresize(imgOligo,0.1);
newfile = [selectedWellLong 'NeuriteBig' optionHandler.FileEnding];
imwrite(imgNeurite,[foldername '/ConvertedCellomics/' newfile] );
newfile3 = [selectedWellLong 'NeuriteSmall' optionHandler.FileEnding];
imwrite(imgNeuriteSmall,[foldername '/ConvertedCellomics/' newfile3] );
newfile2 = [selectedWellLong 'NucleusBigWatershed' optionHandler.FileEnding];
imwrite(BinaryNucleus,[foldername '/ConvertedCellomics/' newfile2] );
newfile4 = [selectedWellLong 'NucleusBig' optionHandler.FileEnding];
imwrite(imgNucleus,[foldername '/ConvertedCellomics/' newfile4] );
newfile5 = [selectedWellLong 'NucleusSmall' optionHandler.FileEnding];
imwrite(imgNucleusSmall,[foldername '/ConvertedCellomics/' newfile5] );
newfile6 = [selectedWellLong 'OligoBig' optionHandler.FileEnding];
imwrite(imgOligo,[foldername '/ConvertedCellomics/' newfile6] );
newfile7 = [selectedWellLong 'OligoSmall' optionHandler.FileEnding];
imwrite(imgOligoSmall,[foldername '/ConvertedCellomics/' newfile7] );




%Check if NeuronStatMatrix exists
if(~isa(neuronHandler.NeuronStatMatrix,'containers.Map'))
    neuronHandler.NeuronStatMatrix = containers.Map();
end
if(~isa(neuronHandler.OmniNucleusPositions,'containers.Map'))
    neuronHandler.OmniNucleusPositions = containers.Map();
end
currentNeuronStat = containers.Map();


        
        NeuronM = NucleusM;
       
        Dummy = zeros(sizeYorig ,sizeXorig);
        Dummy(yStartPos:yEndPos,xStartPos:xEndPos) = NeuronM;
        NeuronM = Dummy;
        summary=0;





neuronHandler.NeuronPositionsEdgeFill(selectedWell) = NeuronM;


% --------------------------------------------------------------------
function SCG_Morphology_Callback(hObject, eventdata, handles)
% hObject    handle to SCG_Morphology (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
neuronHandler = CSVCoordinates();
neuronHandler.ManualNeuronPositions = containers.Map();
csvHandler = CSVCoordinates();
 imageHandler = ImageHandler();
 optionHandler = Option();
 optionHandler.MinSizeNeuriteArea = 25000;
 handles.OptionHandler = optionHandler;
optionHandler.FileEnding = '.png';
[wellList folder ScalingFactor csvHandler] = SCG(csvHandler);
optionHandler.ScalingFactor = ScalingFactor;
 imageHandler.Foldername = folder;
set(handles.lbWell, 'String', wellList);
handles.CSVCoordinates = csvHandler;
handles.NeuronCoordinates = neuronHandler;
handles.ImageHandler = imageHandler;
guidata(handles.figure1, handles);
