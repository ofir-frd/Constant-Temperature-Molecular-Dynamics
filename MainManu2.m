function varargout = MainManu2(varargin)
% MAINMANU2 MATLAB code for MainManu2.fig
%      MAINMANU2, by itself, creates a new MAINMANU2 or raises the existing
%      singleton*.
%
%      H = MAINMANU2 returns the handle to a new MAINMANU2 or the handle to
%      the existing singleton*.
%
%      MAINMANU2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAINMANU2.M with the given input arguments.
%
%      MAINMANU2('Property','Value',...) creates a new MAINMANU2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MainManu2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MainManu2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MainManu2

% Last Modified by GUIDE v2.5 18-Jun-2014 20:12:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MainManu2_OpeningFcn, ...
                   'gui_OutputFcn',  @MainManu2_OutputFcn, ...
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


% --- Executes just before MainManu2 is made visible.
function MainManu2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MainManu2 (see VARARGIN)

% Choose default command line output for MainManu2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MainManu2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MainManu2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function NumParticles_Callback(hObject, eventdata, handles)
% hObject    handle to NumParticles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumParticles as text
%        str2double(get(hObject,'String')) returns contents of NumParticles as a double


% --- Executes during object creation, after setting all properties.
function NumParticles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumParticles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RunButton.
function RunButton_Callback(hObject, eventdata, handles)
% hObject    handle to RunButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%read input data and examine if its acceptable
btnValParticles = str2double(get(handles.NumParticles,'String'));
btnValInitialTemp = str2double(get(handles.ValInitialTemp,'String'));
btnValProbability = str2double(get(handles.ValProbability,'String'));
btnMvalue = str2double(get(handles.Mvalue,'String'));

btnDimentionscontents = get(handles.NumofDimentions,'String'); 
btnDimentionsValue = str2double(btnDimentionscontents{get(handles.NumofDimentions,'Value')});

if   isnan(btnValParticles) ||  btnValParticles<=0 || isnan(btnValInitialTemp) ||  btnValInitialTemp<0 || isnan(btnValProbability) ||  btnValProbability<=0 || isnan(btnMvalue) ||  btnMvalue<=0

     msgbox('Please review your numbers.','Error: Incorrect Input.','error');
    return;

    % run task as selected:
    else

    switch get(get(handles.RadioPanel,'SelectedObject'),'Tag')
        case 'RadioA',  Project_a(btnDimentionsValue,btnValParticles,btnValInitialTemp);
        case 'RadioB',  Project_b(btnDimentionsValue,btnValParticles,btnValInitialTemp,btnValProbability);
        case 'RadioC',  Project_c(btnDimentionsValue,btnValParticles,btnValInitialTemp);
        case 'RadioD',  Project_b(btnDimentionsValue,btnValParticles,btnValInitialTemp,btnValProbability);
        case 'RadioE',  Project_e(btnDimentionsValue,btnValParticles,btnValInitialTemp,btnMvalue);
    end
end

function ValInitialTemp_Callback(hObject, eventdata, handles)
% hObject    handle to ValInitialTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ValInitialTemp as text
%        str2double(get(hObject,'String')) returns contents of ValInitialTemp as a double


% --- Executes during object creation, after setting all properties.
function ValInitialTemp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ValInitialTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CloseButton.
function CloseButton_Callback(hObject, eventdata, handles)
% hObject    handle to CloseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all;


function ValProbability_Callback(hObject, eventdata, handles)
% hObject    handle to ValProbability (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ValProbability as text
%        str2double(get(hObject,'String')) returns contents of ValProbability as a double


% --- Executes during object creation, after setting all properties.
function ValProbability_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ValProbability (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in NumofDimentions.
function NumofDimentions_Callback(hObject, eventdata, handles)
% hObject    handle to NumofDimentions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns NumofDimentions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from NumofDimentions


% --- Executes during object creation, after setting all properties.
function NumofDimentions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumofDimentions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


set(handles.ValInitialTemp, 'String', get(handles.slider1, 'Value'));



% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function Mvalue_Callback(hObject, eventdata, handles)
% hObject    handle to Mvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mvalue as text
%        str2double(get(hObject,'String')) returns contents of Mvalue as a double


% --- Executes during object creation, after setting all properties.
function Mvalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
