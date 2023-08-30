function varargout = specifyCond(varargin)
% SPECIFYCOND MATLAB code for specifyCond.fig
%      SPECIFYCOND, by itself, creates a new SPECIFYCOND or raises the existing
%      singleton*.
%
%      H = SPECIFYCOND returns the handle to a new SPECIFYCOND or the handle to
%      the existing singleton*.
%
%      SPECIFYCOND('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPECIFYCOND.M with the given input arguments.
%
%      SPECIFYCOND('Property','Value',...) creates a new SPECIFYCOND or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before specifyCond_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to specifyCond_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help specifyCond

% Last Modified by GUIDE v2.5 02-Feb-2022 22:47:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @specifyCond_OpeningFcn, ...
                   'gui_OutputFcn',  @specifyCond_OutputFcn, ...
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


% --- Executes just before specifyCond is made visible.
function specifyCond_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to specifyCond (see VARARGIN)

% Choose default command line output for specifyCond
handles.output = hObject;
if(nargin > 3)
    for index = 1:2:(nargin-3)
        if nargin-3==index, break, end
        switch lower(varargin{index})
            case 'tissues'; handles.tId = varargin{index+1};
        end
    end
end
def = {'white', 0.126;
    'gray', 0.276;
    'csf', 1.65;
    'bone', 0.01;
    'skin', 0.465;
    'air', 2.5e-14};
tbl = cell(length(handles.tId),4); 
tbl(:,1) = mat2cell(handles.tId,ones(length(handles.tId),1));
tbl(:,2) = mat2cell(zeros(length(handles.tId),1),ones(length(handles.tId),1));
if ~exist(fullfile(fileparts(mfilename('fullpath')),'custom_conductivities.mat'),'file')
    if length(handles.tId) >= 6; tbl(1:6,3:4) = def; end
else
    c = load(fullfile(fileparts(mfilename('fullpath')),'custom_conductivities.mat'));
    if length(handles.tId) == length(c.cond)
        tbl(:,2:4) = c.cond(:,2:4);
    else
        tbl(1:6,3:4) = def;
    end
end
set(handles.table,'Data',tbl);
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = specifyCond_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Outputs from this function are returned to the command line.
function varargout = table_CellEditCallback(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in confirm.
function confirm_Callback(hObject, eventdata, handles)
% hObject    handle to confirm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cond = get(handles.table,'Data');
if sum(cellfun(@isempty,cond(:,1))) > 0; disp(['Please Specify Name for ALL ' length(handles.tId) ' Tissues']); return; end
if sum(cellfun(@isempty,cond(:,2))) > 0; disp(['Please Specify Conductivity for ALL ' length(handles.tId) ' Tissues']); return; end
save(fullfile(fileparts(mfilename('fullpath')),'custom_conductivities.mat'),'cond')
close(gcf)
