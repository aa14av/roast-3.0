function varargout = file_select(varargin)
% ====================================================================================
% Created by: Alejandro Albizu for the Center of Cognitive Aging and Memory
% University of Florida
% Email: aa14av@gmail.com
% Created: 10/16/18
% Updated: 4/1/20
% ====================================================================================
% FILE_SELECT MATLAB code for file_select.fig
%      FILE_SELECT, by itself, creates a new FILE_SELECT or raises the existing
%      singleton*.
%
%      H = FILE_SELECT returns the handle to a new FILE_SELECT or the handle to
%      the existing singleton*.
%
%      FILE_SELECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FILE_SELECT.M with the given input arguments.
%
%      FILE_SELECT('Property','Value',...) creates a new FILE_SELECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before file_select_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to file_select_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help file_select

% Last Modified by GUIDE v2.5 01-Apr-2020 14:35:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @file_select_OpeningFcn, ...
    'gui_OutputFcn',  @file_select_OutputFcn, ...
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

% --- Executes just before file_select is made visible.
function file_select_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to file_select (see VARARGIN)
if(nargin > 3)
    for index = 1:2:(nargin-3)
        if nargin-3==index, break, end
        switch lower(varargin{index})
            case 'dir'; handles.filedir = varargin{index+1};
            case 'filetype'; handles.filetype = varargin{index+1};
            case 'filter'; handles.filter = varargin{index+1};
        end
    end
    % Choose default command line output for file_select
    handles.output = hObject;
    handles.rootDir = pwd;
    
    % Current Directory
    if isfield(handles,'filedir')
        handles.cdir = handles.filedir;
    else
        handles.cdir = pwd;
    end
    
    % Set PWD_string
    set(handles.pwd_string,'String',handles.cdir);
    
    %Set UP Dir
    up = {fileparts(handles.cdir);fileparts(fileparts(handles.cdir));fileparts(fileparts(fileparts(handles.cdir)))};
    set(handles.Up_Dir_String,'String',up);
    
    % SETUP VARARGIN
    if isfield(handles,'filetype')
        handles.file_ext = handles.filetype;
    else
        handles.file_ext = '.*';
    end
    
    % Set files
    %     row = 1;
    %     for s = handles.subs
    %         for d = handles.days
    %             for n = handles.scans
    %                 if row == 1
    %                     handles.file_srch = dir(fullfile(handles.cdir,strcat(handles.task,'-',handles.subs(s),'-',handles.days(d),handles.file_ext)));
    %                     row = row + 1;
    %                 elseif row > 1
    %                     handles.file_srch = [handles.file_srch; dir(fullfile(handles.cdir,strcat(handles.task,'-',handles.subs(s),'-',handles.days(d),handles.file_ext)))];
    %                 end
    %             end
    %         end
    %     end
    files = dir(fullfile(handles.cdir));
    for repr = 1:length(files)
        pwd_all_files{repr,1} = char(files(repr).name(:)');
        pwd_file_idx(repr,1) = files(repr).isdir(:) & ~strcmp(files(repr).name,'.');
    end
    if ~isfield(handles,'file_srch')
        dirfiles = dir(fullfile(handles.cdir,['*' handles.file_ext]));
        handles.pwd_files = {dirfiles.name};
    else
        if isempty(handles.file_srch)
            dirfiles = dir(fullfile(handles.cdir,['*' handles.filetype]));
            handles.pwd_files = {dirfiles.name};
        else
            for loop = 1:length(handles.file_srch)
                handles.pwd_files(loop,:) = {handles.file_srch(loop).name};
            end
        end
    end
    set(handles.filebox,'String',handles.pwd_files);
%     disp(handles.pwd_files)
    
    %dirbox
    pwd_dir = pwd_all_files(pwd_file_idx == 1);
    set(handles.dirbox,'String',pwd_dir);
    handles.files_out = {};
    set(handles.statbox,'String',' No Files Selected');
    set(handles.sel_box,'String',' ');
else
    % Choose default command line output for file_select
    handles.output = hObject;
    % Root Directory
    handles.rootDir = pwd;
    % Current Directory
    handles.cdir = pwd;
    % Set PWD_string
    set(handles.pwd_string,'String',handles.cdir);
    %Set UP Dir
    up = {fileparts(handles.cdir);fileparts(fileparts(handles.cdir));fileparts(fileparts(fileparts(handles.cdir)))};
    set(handles.Up_Dir_String,'String',up);
    
    % Set files
    files = dir(fullfile(handles.cdir));
    for repr = 1:length(files)
        pwd_all_files{repr,1} = char(files(repr).name(:)');
        pwd_file_idx(repr,1) = files(repr).isdir(:) & ~strcmp(files(repr).name,'.');
    end
    
    handles.pwd_files = pwd_all_files(pwd_file_idx == 0);
    set(handles.filebox,'String',handles.pwd_files);
    
    pwd_dir = pwd_all_files(pwd_file_idx == 1);
    if length(pwd_dir) < 1
        pwd_dir = {'.','..'};
    end
    handles.files_out = {};
    set(handles.dirbox,'String',pwd_dir);
    set(handles.statbox,'String',' No Files Selected');
    set(handles.sel_box,'String',' ');
    set(handles.edit_filter,'String','*.*');
    % Update handles structure
end
handles.finish = 0;
% UIWAIT makes file_select wait for user response (see UIRESUME)
guidata(hObject, handles);
uiwait(handles.file_select)

function pwd_string_Callback(hObject, eventdata, handles)
% hObject    handle to pwd_string (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
changedir = get(hObject,'String');
% Current Directory
handles.cdir = changedir;
% Set PWD_string
set(handles.pwd_string,'String',handles.cdir);
%Set UP Dir
up = {fileparts(handles.cdir);fileparts(fileparts(handles.cdir));fileparts(fileparts(fileparts(handles.cdir)))};
set(handles.Up_Dir_String,'String',up);

% Set files
files = dir(fullfile(handles.cdir));
for repr = 1:length(files)
    pwd_all_files{repr,1} = char(files(repr).name(:)');
    pwd_file_idx(repr,1) = files(repr).isdir(:) & ~strcmp(files(repr).name,'.');
end
%filebox
if isfield(handles,'file_srch')
    % Set files
    row = 1;
    for s = handles.subs
        for d = handles.days
            for n = handles.scans
                if row == 1
                    handles.file_srch = dir(fullfile(handles.cdir,strcat(handles.task,'-',handles.subs(s),'-',handles.days(d),handles.file_ext)));
                    row = row + 1;
                elseif row > 1
                    handles.file_srch = [handles.file_srch; dir(fullfile(handles.cdir,strcat(handles.task,'-',handles.subs(s),'-',handles.days(d),handles.file_ext)))];
                end
            end
        end
    end
    if isempty(handles.file_srch)
        dirfiles = dir(fullfile(handles.cdir,['*' handles.filetype]));
        handles.pwd_files = {dirfiles.name};
    else
        for loop = 1:length(handles.file_srch)
            handles.pwd_files(loop,:) = {handles.file_srch(loop).name};
        end
    end
else
    handles.pwd_files = pwd_all_files(pwd_file_idx == 0);
end
set(handles.filebox,'String',handles.pwd_files);
%dirbox
pwd_dir = pwd_all_files(pwd_file_idx == 1);
set(handles.dirbox,'String',pwd_dir);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pwd_string_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pwd_string (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in Up_Dir_String.
function Up_Dir_String_Callback(hObject, eventdata, handles)
% hObject    handle to Up_Dir_String (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String')); %returns Up_Dir_String contents as cell array
handles.cdir = fullfile(contents{get(hObject,'Value')}); %returns selected item from Up_Dir_String

% Set PWD_string
set(handles.pwd_string,'String',handles.cdir);

%Set UP Dir
up = {fileparts(handles.cdir);fileparts(fileparts(handles.cdir));fileparts(fileparts(fileparts(handles.cdir)))};
set(handles.Up_Dir_String,'String',up);

% Set files
files = dir(fullfile(handles.cdir));
if isempty(files)
    disp('There''s Nothing Here!')
else
    for repr = 1:length(files)
        pwd_all_files{repr,1} = char(files(repr).name(:)');
        pwd_file_idx(repr,1) = files(repr).isdir(:) & ~strcmp(files(repr).name,'.');
    end
    
    %filebox
    if isfield(handles,'file_srch')
        % Set files
        row = 1;
        for s = handles.subs
            for d = handles.days
                for n = handles.scans
                    if row == 1
                        handles.file_srch = dir(fullfile(handles.cdir,strcat(handles.task,'-',handles.subs(s),'-',handles.days(d),handles.file_ext)));
                        row = row + 1;
                    elseif row > 1
                        handles.file_srch = [handles.file_srch; dir(fullfile(handles.cdir,strcat(handles.task,'-',handles.subs(s),'-',handles.days(d),handles.file_ext)))];
                    end
                end
            end
        end
        if isempty(handles.file_srch)
            dirfiles = dir(fullfile(handles.cdir,['*' handles.filetype]));
            handles.pwd_files = {dirfiles.name};
        else
            for loop = 1:length(handles.file_srch)
                handles.pwd_files(loop,:) = {handles.file_srch(loop).name};
            end
        end
    else
        handles.pwd_files = pwd_all_files(pwd_file_idx == 0);
    end
    set(handles.filebox,'String',handles.pwd_files);
    
    pwd_dir = pwd_all_files(pwd_file_idx == 1);
    set(handles.dirbox,'String',pwd_dir);
    set(handles.Up_Dir_String,'Value',1);
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Up_Dir_String_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Up_Dir_String (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in filebox.
function filebox_Callback(hObject, eventdata, handles)
% hObject    handle to filebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String')); % returns filebox contents as cell array
filebox_rm = contents{get(hObject,'Value')}; % returns selected item from filebox\

% file output
if ~isfield(handles,'files_out')
    handles.files_out = {fullfile(handles.cdir,filebox_rm)};
else
    handles.files_out = [handles.files_out; {fullfile(handles.cdir,filebox_rm)}];
end

% Remove selected from filebox
if strcmp(get(handles.sel_box,'String'),' ')
    handles.file_seld = contents(get(hObject,'Value'));
    handles.fullfile_seld = fullfile(handles.cdir,handles.file_seld);
else
    handles.file_seld = [handles.file_seld; contents{get(hObject,'Value')}];
    handles.fullfile_seld = [handles.fullfile_seld; fullfile(handles.cdir,handles.file_seld)];
end
handles.pwd_files(strcmp(handles.pwd_files,filebox_rm)) = [];

set(handles.filebox,'String',handles.pwd_files)
set(handles.sel_box,'String',handles.file_seld)
set(handles.statbox,'String',strcat(num2str(length(handles.file_seld)),' Files Selected'))
set(handles.filebox,'Value',1);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function filebox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in dirbox.
function dirbox_Callback(hObject, eventdata, handles)
% hObject    handle to dirbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String')); %returns dirbox contents as cell array
seld_dir = contents{get(hObject,'Value')}; %returns selected item from Up_Dir_String

if strcmp(seld_dir,'.')
    handles.cdir = handles.cdir;
elseif strcmp(seld_dir,'..')
    handles.cdir = fileparts(handles.cdir);
else
    handles.cdir = fullfile(handles.cdir,seld_dir);
end

% Set PWD_string
set(handles.pwd_string,'String',handles.cdir);

%Set UP Dir
up = {fileparts(handles.cdir);fileparts(fileparts(handles.cdir));fileparts(fileparts(fileparts(handles.cdir)))};
set(handles.Up_Dir_String,'String',up);

% Set files
files = dir(fullfile(handles.cdir));

if isempty(files)
    disp('There''s Nothing Here!')
else
    for repr = 1:length(files)
        pwd_all_files{repr,1} = char(files(repr).name(:)');
        pwd_file_idx(repr,1) = files(repr).isdir(:) & ~strcmp(files(repr).name,'.');
    end
    
    % filebox
    if isfield(handles,'file_srch')
        % Set files
        row = 1;
        for s = handles.subs
            for d = handles.days
                for n = handles.scans
                    if row == 1
                        handles.file_srch = dir(fullfile(handles.cdir,strcat(handles.task,'-',handles.subs(s),'-',handles.days(d),handles.file_ext)));
                        row = row + 1;
                    elseif row > 1
                        handles.file_srch = [handles.file_srch; dir(fullfile(handles.cdir,strcat(handles.task,'-',handles.subs(s),'-',handles.days(d),handles.file_ext)))];
                    end
                end
            end
        end
        if isempty(handles.file_srch)
            dirfiles = dir(fullfile(handles.cdir,['*' handles.filetype]));
            handles.pwd_files = {dirfiles.name};
        else
            for loop = 1:length(handles.file_srch)
                handles.pwd_files(loop,:) = {handles.file_srch(loop).name};
            end
        end
    else
        dirfiles = dir(fullfile(handles.cdir,['*' handles.filetype]));
        handles.pwd_files = {dirfiles.name};
    end
    set(handles.filebox,'String',handles.pwd_files);
    
    % dirbox
    pwd_dir = pwd_all_files(pwd_file_idx == 1);
    set(handles.dirbox,'String',pwd_dir);
    set(handles.dirbox,'Value',1);
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function dirbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dirbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
help
% Hint: get(hObject,'Value') returns toggle state of help

% --- Executes on button press in filter.
function filter_Callback(hObject, eventdata, handles)
% hObject    handle to filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
files = dir(fullfile(handles.cdir,handles.filter));
if isempty(files)
    disp('There''s Nothing Here!')
else
    for repr = 1:length(files)
        pwd_all_files{repr,1} = char(files(repr).name(:)');
        pwd_file_idx(repr,1) = files(repr).isdir(:) & ~strcmp(files(repr).name,'.');
    end
    
    handles.pwd_files = pwd_all_files(pwd_file_idx == 0);
    set(handles.filebox,'String',handles.pwd_files);
end
guidata(hObject, handles);

function edit_filter_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.filter = get(hObject,'String'); % returns contents of edit_filter as text
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_filter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in sel_box.
function sel_box_Callback(hObject, eventdata, handles)
% hObject    handle to sel_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sel_box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sel_box
contents = cellstr(get(hObject,'String')); %returns sel_box contents as cell array
file_rm = contents{get(hObject,'Value')};

% Remove clicked files from sel_box
handles.file_seld(strcmp(handles.file_seld,file_rm)) = [];

% file output
handles.files_out(get(hObject,'Value')) = [];

handles.pwd_files = sort([file_rm; handles.pwd_files]);

set(handles.filebox,'String',handles.pwd_files)
set(handles.sel_box,'String',handles.file_seld)
set(handles.statbox,'String',strcat(num2str(length(handles.file_seld)),' Files Selected'))
set(handles.sel_box,'Value',1);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sel_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sel_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of reset
% Current Directory
handles.cdir = pwd;
% Set PWD_string
set(handles.pwd_string,'String',handles.cdir);
%Set UP Dir
up = {fileparts(handles.cdir);fileparts(fileparts(handles.cdir));fileparts(fileparts(fileparts(handles.cdir)))};
set(handles.Up_Dir_String,'String',up);
% Set files
files = dir(fullfile(handles.cdir));
for repr = 1:length(files)
    pwd_all_files{repr,1} = char(files(repr).name(:)');
    pwd_file_idx(repr,1) = files(repr).isdir(:) & ~strcmp(files(repr).name,'.');
end
%filebox
if isfield(handles,'file_srch')
    % Set files
    row = 1;
    for s = handles.subs
        for d = handles.days
            for n = handles.scans
                if row == 1
                    handles.file_srch = dir(fullfile(handles.cdir,strcat(handles.task,'-',handles.subs(s),'-',handles.days(d),handles.file_ext)));
                    row = row + 1;
                elseif row > 1
                    handles.file_srch = [handles.file_srch; dir(fullfile(handles.cdir,strcat(handles.task,'-',handles.subs(s),'-',handles.days(d),handles.file_ext)))];
                end
            end
        end
    end
    if isempty(handles.file_srch)
        dirfiles = dir(fullfile(handles.cdir,['*' handles.filetype]));
        handles.pwd_files = {dirfiles.name};
    else
        for loop = 1:length(handles.file_srch)
            handles.pwd_files(loop,:) = {handles.file_srch(loop).name};
        end
    end
else
    handles.pwd_files = pwd_all_files(pwd_file_idx == 0);
end
set(handles.filebox,'String',handles.pwd_files);
pwd_dir = pwd_all_files(pwd_file_idx == 1);
if length(pwd_dir) < 1
    pwd_dir = {'.','..'};
end
set(handles.dirbox,'String',pwd_dir);
set(handles.statbox,'String',' No Files Selected');
set(handles.sel_box,'String',' ');
set(handles.edit_filter,'String','*.*');
set(handles.filebox,'Value',1); set(handles.dirbox,'Value',1); set(handles.sel_box,'Value',1);set(handles.Up_Dir_String,'Value',1);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function statbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to statbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of done

disp(strcat(num2str(length(handles.files_out)),' Files Selected !'))
if isequal(get(handles.file_select, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(handles.file_select);
else
    % The GUI is no longer waiting, just close it
    delete(handles.file_select);
end


% --- Outputs from this function are returned to the command line.
function varargout = file_select_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.files_out;
delete(handles.file_select)

% --------------------------------------------------------------------
function select_all_Callback(hObject, eventdata, handles)
% hObject    handle to select_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'file_seld')
    if isempty(handles.file_seld)
        handles.file_seld = handles.pwd_files;
    else
        handles.file_seld = sort(unique([handles.file_seld; handles.pwd_files]));
    end
else
    handles.file_seld = unique(handles.pwd_files);
end

if ~isfield(handles,'files_out')
    handles.files_out = sort(unique(cellfun(@(x) fullfile(handles.cdir,x),handles.file_seld,'UniformOutput',0)));
else
    if isempty(handles.files_out)
        handles.files_out = cellfun(@(x) fullfile(handles.cdir,x),handles.file_seld,'UniformOutput',0);
    else
        handles.files_out = sort(unique([handles.files_out; cellfun(@(x) fullfile(handles.cdir,x),handles.file_seld,'UniformOutput',0)]));
    end
end

% Remove selected from filebox
handles.pwd_files = [];

set(handles.filebox,'String',handles.pwd_files)
set(handles.sel_box,'String',handles.file_seld)
set(handles.statbox,'String',[num2str(length(handles.file_seld)) ' ' 'Files Selected'])
set(handles.filebox,'Value',1);
guidata(hObject, handles);


% --------------------------------------------------------------------
function deselect_all_Callback(hObject, eventdata, handles)
% hObject    handle to deselect_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.pwd_files)
    handles.pwd_files = unique(handles.file_seld);
else
    handles.pwd_files = sort(unique([handles.pwd_files; handles.file_seld]));
end

% Remove selected from filebox
handles.file_seld = [];

% Output files
handles.files_out = [];

set(handles.filebox,'String',handles.pwd_files)
set(handles.sel_box,'String',handles.file_seld)
set(handles.statbox,'String',[num2str(length(handles.file_seld)) ' ' 'Files Selected'])
set(handles.filebox,'Value',1);
guidata(hObject, handles);
