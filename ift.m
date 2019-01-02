function varargout = ift(varargin)
%ift.m, v 0.33
%=========================================================================
%
%                       IMAGE FUSION TOOLBOX            
%
%=========================================================================
%                      Eduardo Fernandez Canga
%
%                       University of Bristol 
%
%                     Copyright (c) 2003 - 2004
%=========================================================================
global ims imout time segmap segwei overlay

if nargin == 0  % LAUNCH GUI

	fig = openfig('ift.fig','new');

	% Generate a structure of handles to pass to callbacks, and store it. 
	hs = guihandles(fig);

    ims=0;        % input images
    imout=0;      % output images
    time=[];
    segmap=0;     % segmentation map, for region fusion 
    segwei=[-1 -1 -1];    % segmentation weight, for region fusion
    overlay=0;
    
    hs.select=0;  % selected images (only the numbers of the images)
   
    hs.outs=1;    % flag for ploting: 1 == one output image per input image 
                    %                        (registration, processing,...)
                    %                 0 == only one output image (fusion)

    hs.imgpath='';   % path for images
    hs.histpath='';  % path for the history
    hs.operation=0;  % operation selected: 0 == none
                     %                     1 == registration
                     %                     2 == processing
                     %                     3 == fusion
                     %                     4 == metrics
    
    % group similar handles into arrays of handles
    % axes handles
    for i=1:3
        command=['hs.Axes_im' num2str(i)];
        hs.Axes(i)=eval(command);
        command=['hs.Axes_Out' num2str(i)];
        hs.AxesOut(i)=eval(command);
        command=['hs.slider' num2str(i)];
        hs.Sliders(i)=eval(command);
        command=['hs.ImNr' num2str(i)];
        hs.ImNr(i)=eval(command);
    end
    
    % fusion controls' handles
    hs.fusion(1)  = hs.fus_met;
    hs.fusion(2)  = hs.fus_bot;
    hs.fusion(3)  = hs.levels;
    hs.fusion(4)  = hs.levels_text;
    hs.fusion(5)  = hs.lisq_coef;
    hs.fusion(6)  = hs.lisq_coef_text;
    hs.fusion(7)  = hs.low;
    hs.fusion(8)  = hs.low_text;
    hs.fusion(9)  = hs.high;
    hs.fusion(10) = hs.high_text;
    hs.fusion(11) = hs.sfbs;
    hs.fusion(12) = hs.sfbs_text;
    hs.fusion(13) = hs.sfthres;
    hs.fusion(14) = hs.sfthres_text;
    hs.fusion(15) = hs.area;
    hs.fusion(16) = hs.area_text;
    hs.fusion(17) = hs.segprio;
    hs.fusion(18) = hs.segprio_text;
    hs.fusion(19) = hs.seghigh;
    hs.fusion(20) = hs.Set_Weight;
    
    % process controls' handles
    hs.proc(1)  = hs.process;    
    hs.proc(2)  = hs.process_app;    
    hs.proc(3)  = hs.dist_met;
    hs.proc(4)  = hs.noise_var;
    hs.proc(5)  = hs.noise_mean;
    hs.proc(6)  = hs.gausvar_text;
    hs.proc(7)  = hs.gausmean_text;
    hs.proc(8)  = hs.speckvar_text;
    hs.proc(9)  = hs.saltdens_text;
    hs.proc(10) = hs.seglevels_text;
    hs.proc(11) = hs.seglevels;
    hs.proc(12) = hs.segtype;
    hs.proc(13) = hs.segtype_text;
    hs.proc(14) = hs.Seg_Over;
    hs.proc(15) = hs.Seg_Map;
    
    
    % metrics controls' handles
    hs.metr(1)  = hs.Metr_met;
    hs.metr(2)  = hs.Metr_app;
    hs.metr(3)  = hs.Metr_Param1;
    hs.metr(4)  = hs.Metr_Param1_text;
    hs.metr(5)  = hs.Metr_Param2;
    hs.metr(6)  = hs.Metr_Param2_text;
    hs.metr(7)  = hs.Metr_Param3;
    hs.metr(8)  = hs.Metr_Param3_text;
    hs.metr(9)  = hs.Metr_Param4;
    hs.metr(10) = hs.Metr_Param4_text;
    hs.metr(11) = hs.Metr_Param5;
    hs.metr(12) = hs.Metr_Param5_text;
    hs.metr(13) = hs.Metr_Param6;
    hs.metr(14) = hs.Metr_Param6_text;
    hs.metr(15) = hs.Metr_Results;

    % save handles
    guidata(fig, hs);    

    % update plot
    change_plot(hs)
    Refresh_All(gca,hs)

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

    try
        [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
    catch
        errordlg(lasterr,'OPEN ERROR','modal')        
        disp(lasterr)
    end

end

% --------------------------------------------------------------------
%  Change Main Operation (Registration, Processing, Fusion or Metrics)
% --------------------------------------------------------------------
function Change_Op_Callback(h,hs,r)
global ims imout
% 'r' is the number corresponding with the operation selected
% see 'switch' down for more details

if r == 0 % called by interpreter
    switch hs.operation
        case 1, Activate_Registration_Callback(h,hs)
        case 2, Activate_Processing_Callback(h,hs)
        case 3, Activate_Fusion_Callback(h,hs)
        case 4, Activate_Metrics_Callback(h,hs)
    end
    return
end

if r~=hs.operation && hs.operation && any(imout(:))
    % test if: change of operation
    %          operation at all
    %          there are output images
    
    % when changing operation outputs and inputs are changed 
    txt=sprintf(['Warning: You are going to change operation.\n',...
            'Which input would you like to preserve?\n']);
    button = questdlg(txt,'Warning','Current Inputs',...'Both',...
        'Current Outputs','Current Outputs');
    if strcmp(button,'Current Outputs')
        ims=imout;
        guidata(h,hs)
        Refresh_All(h,hs)
        change_plot(hs)
        
        ift_hist('Add_Hist',h,hs, 'UPDA')
        
    elseif strcmp(button,'Both')
        p=size(ims,5)+(1:size(imout,5));
        ims(:,:,:,:,p)=imout;
        guidata(h,hs)
        Refresh_All(h,hs)
        change_plot(hs)
        
        ift_hist('Add_Hist',h,hs, 'ADDB')
        
    end
end

if r~=hs.operation
    hs.operation=r;
    guidata(h,hs)
    switch hs.operation
        case 1, imout=0; 
            Activate_Registration_Callback(h,hs)
        case 2, imout=0; 
            Activate_Processing_Callback(h,hs)
        case 3, imout=0; 
            Activate_Fusion_Callback(h,hs)
        case 4, Activate_Metrics_Callback(h,hs)
    end
end

% --------------------------------------------------------------------
%  Registration Menu Selected
% --------------------------------------------------------------------
function Activate_Registration_Callback(h, hs)
errordlg('Function Not implemented yet','REGISTRATION ERROR','modal')

% --------------------------------------------------------------------
%  Processing Menu Selected
% --------------------------------------------------------------------
function Activate_Processing_Callback(h,hs)
% Handles objects' visivility
set(hs.FusPanel,'Visible','off')
set(hs.MetPanel,'Visible','off')
set(hs.ProcPanel,'Visible','on')
set(hs.proc,'visible','off')
set(hs.proc(1:2),'visible','on')

p=get(hs.process,'Value');
hs.outs=1; % 1 == an output per input image
           % 0 == only one output
switch p
    case 1, % Noise
        set(hs.dist_met,'visible','on')
        val = get(hs.dist_met,'Value');
        v(1:4)={'off'};
        if val~=4
            v(val)={'on'};
            v(4)={'on'};
        end
        set(hs.noise_var,'visible',char(v(4)))
        set(hs.noise_mean,'visible',char(v(1)))
        set(hs.gausvar_text,'visible',char(v(1)))
        set(hs.gausmean_text,'visible',char(v(1)))
        set(hs.saltdens_text,'visible',char(v(2)))
        set(hs.speckvar_text,'visible',char(v(3)))
       
    case 2, % Segmentation
        set(hs.proc(10:15),'visible','on')
        
    case 3, % Colouring
    case 4, % Motion Flow
        
end

guidata(h,hs) 
change_plot(hs)
Refresh_All(h,hs)

% --------------------------------------------------------------------
%  Display Overlay Segmentation
% --------------------------------------------------------------------
function Seg_Over(h, hs)
global imout segmap overlay

v=get(h,'Value');
set(hs.Seg_Map,'Value',~v)
hs.outs = 1;        
if any(segmap)
    if v
        imout = overlay;        
    else
        imout = segmap/max(segmap(:));
        imout = imout(:,:,[1,1,1],:,[1,1]);
        if size(imout,5) == 1
            hs.outs = 0;
        end
    end
end
guidata(h,hs)
change_plot(hs)
Refresh_All(h,hs)

% --------------------------------------------------------------------
%  Display Segmentation Map
% --------------------------------------------------------------------
function Seg_Map(h, hs)
global imout segmap overlay

v=get(h,'Value');
set(hs.Seg_Over,'Value',~v)
hs.outs = 1;        
if any(segmap)
    if v
        imout = segmap/max(segmap(:));
        imout = imout(:,:,[1,1,1],:,[1,1]);
        if size(imout,5) == 1
            hs.outs = 0;
        end
    else
        imout = overlay;        
    end
end
guidata(h,hs)
change_plot(hs)
Refresh_All(h,hs)

% --------------------------------------------------------------------
%  Apply Process
% --------------------------------------------------------------------
function Process_Callback(h, hs)
global ims imout segmap overlay
% Check images selected
param{1,1}='sele';
if get(hs.App_Sel,'Value') 
    param{1,2}=hs.select;
else 
    param{1,2}=1:size(ims,5);
end

try
    if ~any(ims)
        errordlg('Please Select Images To Process','PROCESS ERROR','modal')
        return
    end
    set(hs.work_text,'visible','on')
    pause(0.01) 
    
    % NOTE: Any change here, will affect the function proc.m and
    %       interpreter.m
    p=get(hs.process,'Value');
    
    if p==1 % Noise
        p=10*p+get(hs.dist_met,'Value');
        
        % Reading parameters and storing in a cellarray
        param{end+1,1}='vari';
        param{end,2}=str2double(char(get(hs.noise_var,'String'))); 

        param{end+1,1}='dens';
        param{end,2}=str2double(char(get(hs.noise_var,'String'))); 

        param{end+1,1}='medi';
        param{end,2}=str2double(char(get(hs.noise_mean,'String')));

        imout=ims;
        tic
        [imout, titlestr] = ...
            proc(ims,p,param);  % check proc.m for more information
        r=toc;
        
    elseif p==2 % Segmentation
        p=10*p+1;

        % Reading parameters and storing in a cellarray
        param{end+1,1}='lvl';
        param{end,2}=get(hs.seglevels,'Value')+1;
        
        param{end+1,1}='stype';
        param{end,2}=get(hs.segtype,'Value');

        tic
        [overlay, titlestr, segmap] = ...
            proc(ims,p,param); % check proc.m for more information   
        r=toc;
        if get(hs.Seg_Over,'Value')
            imout=overlay;
        else
            imout=segmap(:,:,[1 1 1],:,:)/max(segmap(:));
        end
    elseif p==3 % Colouring
        c = uisetcolor([1 1 1],'Choose a Color');
        param{end+1,1} = 'c';
        param{end,2} = c;
        tic        
        p=10*p+1;
        [imout, titlestr] = ...
            proc(ims,p,param); % check proc.m for more information 
        r=toc;
    elseif p==4 % Motion Flow
        tic        
        p=10*p+1;
        [hs.vectors, titlestr] = ...
            proc(ims,p,param); % check proc.m for more information 
        r=toc;
    end
catch
    set(hs.work_text,'visible','off')
    errordlg(lasterr,'PROCESS ERROR','modal')
    return
end

titlestr=[titlestr ' % ' num2str(r) ' sec'];

ift_hist('Add_Hist',h,hs,titlestr)

Refresh_All(h,hs)
guidata(h,hs) 
set(hs.work_text,'visible','off')

% --------------------------------------------------------------------
%  Fusion Method Selection
% --------------------------------------------------------------------
function Activate_Fusion_Callback(h, hs)
global imout
set(hs.MetPanel,'Visible','off')
set(hs.ProcPanel,'Visible','off')
set(hs.FusPanel,'Visible','on')
%set(hs.fusion(1:2),'visible','on');

hs.outs=0; % 0 == only one output image
                % 1 == one output image per input

if size(imout,5)>1
    imout=0;
end

meth = get(hs.fus_met,'Value'); % value of selected method
vis ={'off' 'on'};


s(5) = meth == 18;                           % segmentation priority 
s(6) = s(5) & get(hs.seghigh,'Value') == 2;
s(1) = (meth >= 5 & meth ~= 17 & meth ~=16); % decomp leveles, high and low
s(2) = ~s(5) & s(1) & get(hs.high,'Value')>1;% area 
s(3) = meth == 17;                           % block size and threshold
s(4) = meth == 15;                           % coefficients
s(7) = ~s(5)&s(4);   %lisq coef

s=s+1;

set(hs.levels_text,   'visible', vis{ s(1) })
set(hs.levels,        'visible', vis{ s(1) })
set(hs.high_text,     'visible', vis{ s(1) })
set(hs.high,          'visible', vis{ s(1)-s(5)+1})
set(hs.low_text,      'visible', vis{ s(1) })
set(hs.low,           'visible', vis{ s(1) })
set(hs.area_text,     'visible', vis{ s(2) })
set(hs.area,          'visible', vis{ s(2) })
set(hs.sfbs_text,     'visible', vis{ s(3) })
set(hs.sfbs,          'visible', vis{ s(3) })
set(hs.sfthres_text,  'visible', vis{ s(3) })
set(hs.sfthres,       'visible', vis{ s(3) })
set(hs.lisq_coef_text,'visible', vis{ s(7) })
set(hs.lisq_coef,     'visible', vis{ s(7) })
set(hs.segprio_text,  'visible', vis{ s(5) })
set(hs.segprio,       'visible', vis{ s(5) })
set(hs.seghigh,       'visible', vis{ s(5) })
set(hs.Set_Weight,    'visible', vis{ s(6) })
guidata(h,hs)
change_plot(hs)
Refresh_All(h,hs)

% --------------------------------------------------------------------
%  Fusion Button
% --------------------------------------------------------------------
function Fusion_Callback(h, hs)
global ims imout segmap segwei
% Check images selected
param{1,1}='sele';
if get(hs.App_Sel,'Value') 
    param{1,2}=hs.select;
else 
    param{1,2}=1:size(ims,5);
end

if ~any(ims)
    errordlg('Please Select Images To Fusion','FUSION ERROR','modal')
else
    v=get(hs.fus_met,'value');
    try
        set(hs.work_text,'visible','on')
        pause(0.01)        
   % **********************************************************************
   % NOTE: Any change here, will affect the function fuse.m and interpreter
   %***********************************************************************

        imout=0;
        % Reading parameters and storing in a cellarray
        param{end+1,1}='lvl';
        param{end,2}=get(hs.levels,'Value')+1;

        param{end+1,1}='low';
        param{end,2}=str2double(char(get(hs.low,'String')));        

        r=(get(hs.area,'String'));
        param{end+1,1}='high';
        param{end,2}=[get(hs.high,'Value'),...
            str2double(char(r(get(hs.area,'Value'))))];        

        param{end+1,1}='lisq_coef';
        param{end,2}=get(hs.lisq_coef,'Value');

        param{end+1,1}='bs';
        param{end,2}=str2double(char(get(hs.sfbs,'String')));

        param{end+1,1}='thrs';
        param{end,2}=str2double(char(get(hs.sfthres,'String')));

        param{end+1,1}='seghigh';
        param{end,2}=get(hs.seghigh,'Value');        

        param{end+1,1}='segprio';
        param{end,2}=get(hs.segprio,'Value');        

        param{end+1,1}='segmap';
        param{end,2}=segmap; 

        param{end+1,1}='segwei';
        param{end,2}=segwei;
        
        param{end+1,1}='hand_hist';
        param{end,2}=hs.History;
                
        tic
        [imout,titlestr]=fuse(ims(:,:,:,:,param{1,2}),v,param);
        r=toc;
        titlestr=[titlestr ' % ' num2str(r) ' sec'];
        
    catch
        set(hs.work_text,'visible','off')
        errordlg(lasterr,'FUSION ERROR','modal')
        return
    end
    set(hs.work_text,'visible','off')
    if numel(imout) ~= 1
        Refresh_All(h, hs)
        
        ift_hist('Add_Hist',h,hs, titlestr)
        
        guidata(h,hs) 
    end
end

% --------------------------------------------------------------------
%  Set Weight for Regions
% --------------------------------------------------------------------
function Set_Weight(h,hs)
global ims segmap segwei overlay

if segmap == 0
    errordlg('No Segmentation Map','HISTORY ERROR','modal')
    return
end
point = pclick;
r=hs.Axes==point(1);
img=str2double(char(get(hs.ImNr(r),'String')));
img=img-size(ims,5)/2;
%p=str2double(char(get(hs.frame_num,'String')));

region = segmap(round(point(3)),round(point(2)));

pos1 = find(segwei(:,1) == img);
pos2 = find(segwei(pos1,2) == region);

if isempty(pos2)
    txt = {''};
else
    txt = mat2cell(num2str( segwei(pos1(pos2)+2*size(segwei,1))));
end

wei = str2double(cell2mat(inputdlg('Insert Weight',...
              'Segmentation Weighting',1,txt)));
if segwei(1,1:2) == [-1 -1]
    segwei = [img region wei point(2) point(3)];
elseif ~strcmp(txt,'')
    segwei(pos1(pos2)+2*size(segwei,1)) = wei;
else
    segwei(end+1,1:5) = [img region wei point(2) point(3)];
end
guidata(h,hs)
r=text(point(2),point(3),num2str(wei));
set(r,'Color','red')




% --------------------------------------------------------------------
%  Metrics Menu Selection
% --------------------------------------------------------------------
function Activate_Metrics_Callback(h, hs)
set(hs.MetPanel,'Visible','on')
set(hs.ProcPanel,'Visible','off')
set(hs.FusPanel,'Visible','off')
set(hs.metr,'visible','on')

meth = get(hs.Metr_met,'Value'); % value of selected method
vis ={'off' 'on'};

s(1) = (meth >= 4);
s(2) = (meth > 5);

if meth == 4
    set(hs.Metr_Param1_text, 'String', 'BlockSize');
    set(hs.Metr_Param1, 'String', 8);
    set(hs.Metr_Param2_text, 'String', 'EdgeWeight');
    set(hs.Metr_Param2, 'String', 0.2);
elseif meth == 5
    set(hs.Metr_Param1_text, 'String', 'Levels');
    set(hs.Metr_Param2_text, 'String', 'Wavf');
    set(hs.Metr_Param3_text, 'String', 'Recons');
    set(hs.Metr_Param1, 'String', 6);
    set(hs.Metr_Param2, 'String', 'bior6.8');
    set(hs.Metr_Param3, 'String', 0);
else
    set(hs.Metr_Param1_text, 'String', 'Param 1:');
    set(hs.Metr_Param2_text, 'String', 'Param 2:');
    set(hs.Metr_Param3_text, 'String', 'Param 3:');
    set(hs.Metr_Param4_text, 'String', 'Param 4:');
    set(hs.Metr_Param5_text, 'String', 'Param 5:');
    set(hs.Metr_Param6_text, 'String', 'Param 6:');
    set(hs.Metr_Param1, 'String', 0);
    set(hs.Metr_Param2, 'String', 0);
    set(hs.Metr_Param3, 'String', 0);
    set(hs.Metr_Param4, 'String', 0);
    set(hs.Metr_Param5, 'String', 0);
    set(hs.Metr_Param6, 'String', 0);
end
s=s+1;

set(hs.Metr_Param1_text, 'visible', vis{ s(1) })
set(hs.Metr_Param1,      'visible', vis{ s(1) })
set(hs.Metr_Param2_text, 'visible', vis{ s(1) })
set(hs.Metr_Param2,      'visible', vis{ s(1) })
set(hs.Metr_Param3_text, 'visible', vis{ s(2) })
set(hs.Metr_Param3,      'visible', vis{ s(2) })
set(hs.Metr_Param4_text, 'visible', vis{ s(2) })
set(hs.Metr_Param4,      'visible', vis{ s(2) })
set(hs.Metr_Param5_text, 'visible', vis{ s(2) })
set(hs.Metr_Param5,      'visible', vis{ s(2) })
set(hs.Metr_Param6_text, 'visible', vis{ s(2) })
set(hs.Metr_Param6,      'visible', vis{ s(2) })

% --------------------------------------------------------------------
%  Apply Metric Button
% --------------------------------------------------------------------
function Metrics_Callback(h, hs)
global ims imout
param{1,1}='sele';
if get(hs.App_Sel,'Value') 
    param{1,2}=hs.select;
else 
    param{1,2}=1:size(ims,5);
end

if ~any(ims)
    errordlg('Please Select Images To Compute','METRICS ERROR','modal')
else
    v=get(hs.Metr_met,'value');
    try
        set(hs.work_text,'visible','on')
        pause(0.01)
        % NOTE: Any change here, will affect the function metric.m and
        %       interpreter.m
        hs.metric_out=0;
        
        % You can play around whit this parameters
        param{end+1,1}='p1';
        param{end,2}=str2double(char(get(hs.Metr_Param1,'String')));
        param{end+1,1}='p2';
        param{end,2}=str2double(char(get(hs.Metr_Param2,'String')));
        param{end+1,1}='p3';
        param{end,2}=str2double(char(get(hs.Metr_Param3,'String')));
        param{end+1,1}='p4';
        param{end,2}=str2double(char(get(hs.Metr_Param4,'String')));
        param{end+1,1}='p5';
        param{end,2}=str2double(char(get(hs.Metr_Param5,'String')));
        param{end+1,1}='p6';
        param{end,2}=str2double(char(get(hs.Metr_Param6,'String')));
        
        % This ones are for metrics already implemented
        param{end+1,1}='bs';
        param{end,2}=str2double(char(get(hs.Metr_Param1,'String')));
        param{end+1,1}='cedge';
        param{end,2}=str2double(char(get(hs.Metr_Param2,'String')));
        param{end+1,1}='lvl';
        param{end,2}=str2double(char(get(hs.Metr_Param1,'String')));
        param{end+1,1}='wavf';
        param{end,2}=char(get(hs.Metr_Param2,'String'));
        param{end+1,1}='recon';
        param{end,2}=str2double(char(get(hs.Metr_Param3,'String')));
        tic
        [result_str,titlestr]=metric(ims,imout,v,param);
        r=toc;
        titlestr=[titlestr ' % ' num2str(r) ' sec'];
        set(hs.Metr_Results,'String',result_str);  
    catch
        set(hs.work_text,'visible','off')
        errordlg(lasterr,'METRICS ERROR','modal')
        return
    end
end


set(hs.work_text,'visible','off')
if numel(imout)~=1
    Refresh_All(h, hs)
    
    ift_hist('Add_Hist',h,hs, titlestr)
    
    guidata(h,hs)
end

% --------------------------------------------------------------------
%  Save Metrics Results
% --------------------------------------------------------------------
function Save_Metric_List(h,hs)
try
    pathname = hs.histpath;
    [filename, pathname] = uiputfile([pathname '*.txt'],...
        'Save Metric Results');
    fid=fopen([pathname,filename],'w');
    if fid < 0
        error('Error: cant open file');
        return
    end
    history=get(hs.Metr_Results,'String');
    fprintf(fid,'%s\r\n',history{1:size(history,1)});
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
%  Delete_All: Callback
% --------------------------------------------------------------------
function Delete_All(h, hs)
global ims imout time

button = questdlg('Are you sure you want to clear the input images?',...
    'Warning','Yes','No','No');
if strcmp(button,'Yes')
    ims=0;
    imout=0;
    time=[];
    guidata(h,hs)
    Refresh_All(h,hs)
    change_plot(hs)
    
    ift_hist('Add_Hist',h,hs, 'DELE Input 0')
    
end

% --------------------------------------------------------------------
%  Delete Single Input (through context menu)
% --------------------------------------------------------------------
function Del_Input(h, hs)
global ims imout time

% if ~any(ims(:)) % no image to delete
%     return
% end

v = get(gcbf,'CurrentAxes');  % gets handle of selected image
r = hs.Axes==v;               % compares with all handles of input image
v = str2double(char(get(hs.ImNr(r),'String')));  % reads number of image from 
                                                 % its block number

c=size(ims,5); % number of input sensors
% remove selected image
if c<2 
    ims=0;
    time=[];
else
    p=1:c;
    p=p(p~=(v));
    ims=ims(:,:,:,:,p);
    if ~isempty(time)
        time=time(:,p);
    end
end

% set sliders' values
m=get(hs.Sliders(r),'Max');

slider_step(1) = 1/max((m-1),1);
slider_step(2) = 1/max((m-1),1);
set(hs.Sliders,'sliderstep',slider_step,'max',max(1,m-1))

if m>=2
    for i=2:3
        r=get(hs.Sliders(i),'Value');
        if r>=m
            set(hs.Sliders(i),'Value',r-1);
        end
    end
end
% update images
Refresh_All(h,hs)
guidata(h,hs)

Select_Input(h, hs, v) % rearange selection

ift_hist('Add_Hist',h,hs, ['DELE Input ' num2str(v)])

change_plot(hs)

% --------------------------------------------------------------------
%  View Output through context menu
% --------------------------------------------------------------------
function View_Output(h, hs)
global imout

v=get(gcbf,'CurrentAxes');
r=hs.AxesOut==v;
n=str2double(char(get(hs.ImNr(r),'String')));
p=str2double(char(get(hs.frame_num,'String')));

if numel(imout)==1
    errordlg('There are no output images','VIEW ERROR','modal')
    return
end
figure,image(imout(:,:,:,p,n))
title(['Output Image ' num2str(n)])
axis image

% --------------------------------------------------------------------
%  View Input through context menu
% --------------------------------------------------------------------
function View_Input(h, hs)
global ims

v=get(gcbf,'CurrentAxes');
r=hs.Axes==v;
n=str2double(char(get(hs.ImNr(r),'String')));
p=str2double(char(get(hs.frame_num,'String')));

if numel(ims)==1
    errordlg('There are no input images','VIEW ERROR','modal')
    return
end
figure,image(ims(:,:,:,p,n))
title(['Input Image ' num2str(n)])
axis image

% --------------------------------------------------------------------
%  Input Selection
% --------------------------------------------------------------------
function Select_Input(h, hs, n)
if nargin < 3
    v=get(gcbf,'CurrentAxes');
    r=hs.Axes==v;
    n=str2double(char(get(hs.ImNr(r),'String')));
    if hs.select==0            % if first select
        hs.select=n;           % just select
    elseif any(hs.select == n) % if already selected
        hs.select=hs.select(hs.select ~= n); % remove from selection
    else
        hs.select=[hs.select n]; % add to selection
    end
elseif any(hs.select >= n)
    hs.select=[hs.select(find(hs.select<n)) hs.select(find(hs.select>n))-1 ];
end

% activate or desactivate apply buttons
if isempty(hs.select) || hs.select==0
    hs.select=0;
    set(hs.App_All,'Enable','off')
    set(hs.App_Sel,'Enable','off')
    set(hs.App_All,'Value',1)
    set(hs.App_Sel,'Value',0)
else
    set(hs.App_All,'Enable','on')
    set(hs.App_Sel,'Enable','on')
end
        
guidata(h,hs)
Refresh_All(h,hs)

% --------------------------------------------------------------------
%  Executes on any slider movement. It is usefull to refresh images
% --------------------------------------------------------------------
function Refresh_All(h, hs, m)
global ims imout

% Change/Refresh picture
if nargin < 3 
    m=1:3; % update all
end

f=size(ims,4);
n=size(ims,5);

% checks that the number of the frame is in range
v=min(f-1,round(get(hs.frame_sel,'Value')));
set(hs.frame_sel,'Value',v)    % set frame slider value
set(hs.frame_num,'String',v+1) % set frame number

fr=v+1; % frame number

% v will be a vector with the values of the three sensors displayed
if n < 4
    v = 1:3;
else
    v = min(n,round(cell2mat(get(hs.Sliders,'Value')))'+1); 
end

for i=m
   set(hs.Sliders(i),'Value',v(i)-1)
   set(hs.ImNr(i),'String',v(i))
   % M1 input image
   if size(ims,5)<(v(i)) || numel(ims)==1
       M1=180/255*ones(34,41,3);
   else 
       M1=ims(:,:,:,fr,v(i));
   end
   
   % M2 output image
   if ~hs.outs && numel(imout) ~= 1
       M2=imout(:,:,:,fr,1);
   elseif size(imout,5)<(v(i)) || numel(imout)==1
       M2=180/255*ones(34,41,3);
   else
       M2=imout(:,:,:,fr,v(i));
   end
   
   set(gcf,'CurrentAxes',hs.AxesOut(i))
   image(M2,'UIContextMenu',hs.Mod_Output);
   axis image
      
   set(gcf,'CurrentAxes',hs.Axes(i))
   %if selected, frame is red, otherwise frame is blue
   if any(hs.select == v(i)) 
       set(hs.Axes(i),'XColor',[1 0 0],'YColor',[1 0 0],'LineWidth',1)
   else
       set(hs.Axes(i),'XColor',[197 251 254]/255,'YColor',[197 251 254]/255,'LineWidth',.5)
   end
   image(M1,'UIContextMenu',hs.Mod_Input);
   axis image

   guidata(h,hs)
end

% Activate or deactivate context menus
if any(ims(:))
    set(hs.Del_Input,'Enable','on');   set(hs.Delete_All,'Enable','on')
    set(hs.View_Input,'Enable','on');  set(hs.Select_Input,'Enable','on')
    set(hs.Save_Video_In,'Enable','on');
else
    set(hs.Del_Input,'Enable','off');  set(hs.Delete_All,'Enable','off')
    set(hs.View_Input,'Enable','off'); set(hs.Select_Input,'Enable','off')
    set(hs.Save_Video_In,'Enable','off');
end
if any(imout(:))
    set(hs.update,'Enable','on');      set(hs.Save_Output,'Enable','on')
    set(hs.Save_Video_Out,'Enable','on');  set(hs.View_Output,'Enable','on')
else
    set(hs.update,'Enable','off');     set(hs.Save_Output,'Enable','off')
    set(hs.Save_Video_Out,'Enable','off'); set(hs.View_Output,'Enable','off')
end


% --------------------------------------------------------------------
%  Executes on button press in zoom.
% --------------------------------------------------------------------
function zoom_Callback(hs)
r=get(hs.zoom,'Value');
if r
    zoom on
else
    zoom off
end

% --------------------------------------------------------------------
%  Executes on button press in zoom all.
% --------------------------------------------------------------------
function zoom_all_Callback(h, hs)
set(hs.Axes,'Ylim',ylim,'Xlim',xlim)
set(hs.AxesOut,'Ylim',ylim,'Xlim',xlim)

% --------------------------------------------------------------------
%  Executes on button press in grid.
% --------------------------------------------------------------------
function grid_Callback(hs)
r=get(hs.grid,'Value');
if r
   for i=1:3
        set(gcf,'CurrentAxes',hs.Axes(i))
        set(hs.Axes(i),'XtickMode','auto','YtickMode','auto')
        grid on
        set(gcf,'CurrentAxes',hs.AxesOut(i))
        set(hs.AxesOut(i),'XtickMode','auto','YtickMode','auto')
        grid on
    end
else
    for i=1:3
        a=[];
        set(gcf,'CurrentAxes',hs.Axes(i))
        set(hs.Axes(i),'Xtick',a,'Ytick',a)
        grid off
        set(gcf,'CurrentAxes',hs.AxesOut(i))
        set(hs.AxesOut(i),'Xtick',a,'Ytick',a)
        grid off
    end
end

% --------------------------------------------------------------------
%  Plot Positions. Changes the Axe's Position 
% --------------------------------------------------------------------
function change_plot(hs)
change_position(hs) %external m file
grid_Callback(hs)
zoom_Callback(hs)

% --------------------------------------------------------------------
%  Update. Move Output Images to Input Images
% --------------------------------------------------------------------
function Update_Callback(h, hs)
global ims imout
txt=sprintf('This operation will clear input images.\n\n',...
    'Do you want to proceed?');
button = questdlg(txt,'Warning','Yes','No','No');
if strcmp(button,'Yes')
    ims=0;
    ims=imout;
    imout=0;
    
    guidata(h,hs)
    Refresh_All(h,hs)
    change_plot(hs)
    
    ift_hist('Add_Hist',h,hs, 'UPDA')
    
    guidata(h,hs)
end

% --------------------------------------------------------------------
% Plays Sequence of Images
% --------------------------------------------------------------------
function frame_play(h, hs)
global ims
for i=1:size(ims,4)
    set(hs.frame_sel,'Value',i-1)
    Refresh_All(h,hs)
    pause(0.3)
end

% --- Executes on button press in tests.
function tests_Callback(h, hs)
global ims imout
keyboard

%
% to use in the future
%
function video(h,hs)
mov=avifile('testa.avi','FPS',2,'COMPRESSION','NONE','QUALITY',100);
[b a c d e] = size(ims);
frame=100*ones(e*b+10*(e+1),2*a+30);
for i=1:d
    for j=1:e
        frame((j-1)*(10+b)+10+(1:b),10+(1:a))=ims(:,:,:,i,j);
    end
    frame(round((e*b-b+10*e)/2)+(1:b),20+((a+1):2*a))=imout(:,:,:,i,1);
    frame=round(frame);
    frame=min(max(frame,0),255);
    f(i)=im2frame(round(frame)+1,gray(256));
    mov=addframe(mov,f(i));
end
close(mov);

% --------------------------------------------------------------------
function Exit(h, hs)
close(hs.Interface)

