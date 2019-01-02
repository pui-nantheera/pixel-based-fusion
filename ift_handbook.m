function varargout = ift_handbook(varargin)
% IFT_HANDBOOK M-file for ift_handbook.fig
%      IFT_HANDBOOK, by itself, creates a new IFT_HANDBOOK or raises the existing
%      singleton*.
%
%      H = IFT_HANDBOOK returns the handle to a new IFT_HANDBOOK or the handle to
%      the existing singleton*.
%
%      IFT_HANDBOOK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IFT_HANDBOOK.M with the given input arguments.
%
%      IFT_HANDBOOK('Property','Value',...) creates a new IFT_HANDBOOK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ift_handbook_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ift_handbook_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help ift_handbook

% Last Modified by GUIDE v2.5 12-Jul-2004 11:47:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ift_handbook_OpeningFcn, ...
                   'gui_OutputFcn',  @ift_handbook_OutputFcn, ...
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


% --- Executes just before ift_handbook is made visible.
function ift_handbook_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ift_handbook (see VARARGIN)

% Choose default command line output for ift_handbook
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ift_handbook wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ift_handbook_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
