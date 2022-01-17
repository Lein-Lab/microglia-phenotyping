function varargout = ClickTool(varargin)
% CLICKTOOL MATLAB code for ClickTool.fig
%      CLICKTOOL, by itself, creates a new CLICKTOOL or raises the existing
%      singleton*.
%
%      H = CLICKTOOL returns the handle to a new CLICKTOOL or the handle to
%      the existing singleton*.
%
%      CLICKTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLICKTOOL.M with the given input arguments.
%
%      CLICKTOOL('Property','Value',...) creates a new CLICKTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ClickTool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ClickTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ClickTool

% Last Modified by GUIDE v2.5 15-Feb-2013 10:49:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ClickTool_OpeningFcn, ...
                   'gui_OutputFcn',  @ClickTool_OutputFcn, ...
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


% --- Executes just before ClickTool is made visible.
function ClickTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ClickTool (see VARARGIN)

% Choose default command line output for ClickTool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ClickTool wait for user response (see UIRESUME)
% uiwait(handles.clickToolFigure);


% --- Outputs from this function are returned to the command line.
function varargout = ClickTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pbFinish.
function pbFinish_Callback(hObject, eventdata, handles)
% hObject    handle to pbFinish (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pbDeleteLastClick.
function pbDeleteLastClick_Callback(hObject, eventdata, handles)
% hObject    handle to pbDeleteLastClick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function clickToolFigure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clickToolFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Set Axes from guidata
% Get Image from ImageHandler and Extract Current Filter Region minus all
% Negative Filter regions if applicable
%handles = guidata(handles.clickToolFigure);
%imageHandler = handles.ImageHandler;
h = findall(0,'type','figure','Name','figure1')
imageHandler = h.imageHandler
%csvHandler = handles.CSVCoordinates;
%neuronHandler = handles.NeuronCoordinates;
positiveFilterList = imageHandler.PositiveFilters;
negativeFilterList = imageHandler.NegativeFilters;
lbFilterString = get(handles.lbFilter,'String');
%For every positive entry in Filter List: Create entry in Listbox
for i=1:numel(positiveFilterList)
    fulltext = strcat({lbFilterString,'|',i})
end
set(handles.lbFilter, 'String', fulltext);
currentImage = imfuse(imageHandler.NucleusImage,imageHandler.NeuriteImage);



% --- Executes on selection change in lbFilter.
function lbFilter_Callback(hObject, eventdata, handles)
% hObject    handle to lbFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lbFilter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lbFilter


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
