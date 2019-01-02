function varargout = ift_hist(varargin)
%ift_hist.m, v 0.1
%=========================================================================
%
%                       IMAGE FUSION TOOLBOX
%
%=========================================================================
%
%        This subroutine handles the operations performed through the
%        History window. (Load, Save, Run and Clear)
%
%=========================================================================
%
%                      Eduardo Fernandez Canga
%
%                       University of Bristol 
%
%                     Copyright (c) 2003 - 2004
%=========================================================================

try
    [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
catch
    errordlg(lasterr,'OPEN ERROR','modal')
    disp(lasterr)
end

% --------------------------------------------------------------------
%  Add Line to History
% --------------------------------------------------------------------
function Add_Hist(h, hs, hist_line)
a=get(hs.History,'String');
a{size(a,1)+1,1}=[hist_line];
set(hs.History,'String',a);
set(hs.History,'Value',size(a,1))


% --------------------------------------------------------------------
%  Clear History
% --------------------------------------------------------------------
function Clear_Hist(h, hs)
set(hs.History,'String','');

% --------------------------------------------------------------------
%  Load History
% --------------------------------------------------------------------
function Load_Hist(h, hs)
try
    pathname = hs.histpath;
    [filename, pathname] = uigetfile([pathname '*.txt'], 'Load History');
    if pathname==0
        return
    end

    file=[pathname,filename];
    fid=fopen([pathname,filename],'r');
    if fid < 0
        error('Error: cant open file');
        return
    end
    r=1;
    while r
        a=fgets(fid);
        if a==-1
            break
        end
        history{r}=a;
        r=r+1;
    end
    if fclose(fid)
        error('Error: cant close file')
    end
    set(hs.History,'Value',1);
    set(hs.History,'String',deblank(history));

    guidata(h,hs) 
catch
    errordlg(lasterr,'LOAD ERROR','modal')
    return
end


% --------------------------------------------------------------------
%  Save History
% --------------------------------------------------------------------
function Save_Hist(h, hs)
try
    pathname = hs.histpath;
    [filename, pathname] = uiputfile([pathname '*.txt'], 'Save History');
    if pathname==0
        return
    end
    
    file=[pathname,filename];
    fid=fopen([pathname,filename],'w');
    if fid < 0
        error('Error: cant open file');
        return
    end
    history=get(hs.History,'String');
    fprintf(fid,['%s\r\n'],history{1:size(history,1)});
    if fclose(fid)
        error('Error: cant close file')
    end

    hs.histpath=pathname;
    guidata(h,hs) 
catch
    errordlg(lasterr,'SAVE ERROR','modal')
    return
end

% --------------------------------------------------------------------
%  Read History. Reads one line of the history when it is pressed
% --------------------------------------------------------------------
function History_Callback(h, hs)
global ims imout
set(h,'Enable','off');
pause(0.2)
s=get(h,'String');
v=get(h,'Value');
try
    [hs.operation hs.outs] = interpreter(s{v},hs.operation, hs.outs);
    
    guidata(h,hs)
    ift('Change_Op_Callback',h, hs,0)
    ift('Refresh_All',h, hs)
    ift('change_plot',hs)
catch
    set(h,'Enable','on');
    errordlg(lasterr,'HISTORY ERROR','modal')
end
set(h,'Enable','on');

% --------------------------------------------------------------------
%  Run History 
% --------------------------------------------------------------------
function Run_Hist(h, hs)
set(hs.History,'Enable','off');
pause(0.2)
t=get(hs.History,'String');
p=0;
for i=1:size(t,1)
    try
        tic
        [s,ims,imout]=...
            interpreter(t{i},[hs.operation hs.outs],ims,imout);
        p=toc + p;
        hs.operation=s(1); hs.outs=s(2);
        guidata(h,hs)
        Change_Op_Callback(h, hs,0)
        Refresh_All(h, hs)
        change_plot(hs)
    catch
        errordlg(lasterr,'HISTORY ERROR','modal')
    end    
end
set(hs.History,'Enable','on');
titlestr = ['HIST RUN comp_time =' num2str(p) ' sec;'];
Add_Hist(h,hs,titlestr)