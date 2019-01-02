function varargout = ift_help(varargin)
%ift_help.m, v 0.1
%===========================================================================
%
%                     IMAGE FUSION TOOLBOX HELP            
%
%===========================================================================
%                      Eduardo Fernandez Canga
%
%                       University of Bristol 
%
%                     Copyright (c) 2003 - 2004
%===========================================================================

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end

end

% --------------------------------------------------------------------
% Close itself
% --------------------------------------------------------------------
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)
close(handles.figure1)
