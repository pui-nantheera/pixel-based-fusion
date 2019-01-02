function change_position(handles)
%change_position.m, v 0.1 2004/01/29
%===========================================================================
%               Eduardo Fernandez Canga - University of Bristol
%
%                        Copyright (c) 2004
%===========================================================================
%     function change_position(handles)
%
%      Part of the IFT interface. Sets the positions of axes, sliders and
%    number boxes regarding the number of input and output images.

global ims imout
% ***** Positions *********
% *** 2 images
A2_x=0.0243; A2_y=0.5343; % axes position
A2_s=[.36 .45];           % axes size

S2_x=0.0106; S2_y=0.5690; % sliders position
S2_s=[0.0107 0.3857];     % sliders size

N2_x=0.0048; N2_y=0.5343; % number position
N2_s=[0.0186 0.02];       % number size

D2=0.46;                  % distance

% *** 3 images
A3_x=0.0293; A3_y=0.677; %0.6643;
A3_s=[.3193 .32];

S3_x=0.0156; S3_y=0.6990;
S3_s=[0.0107 0.2857];

N3_y=0.6643;

D3=0.324;  %0.3272

% *** frame selector
Fs_x = 0.05; Fs_y=.25;
Fs_s = [.22 .015];

Fn_x = 0.28; Fn_y=.247;
Fn_s = [0.0186 0.02];

Fp_x = 0.31; Fp_y=.247;
Fp_s = [0.02 0.02];


%**********************************************
cin=size(ims,5);     % number of input sensors
cout=size(imout,5);  % number of output modalities

slider_step(1) = 1/((max(cin-1,1)));
slider_step(2) = 1/(max(cin-1,1));
set(handles.Sliders,'sliderstep',slider_step,'max',max(cin-1,1))

frames=size(ims,4);  % number of frames

slider_step(1) = 1/((max(frames-1,1)));
slider_step(2) = 1/(max(frames-1,1));
set(handles.frame_sel,'sliderstep',slider_step,'max',max(frames-1,1))
if frames < 2
    set(handles.frame_sel,'Visible','off');
    set(handles.frame_num,'Visible','off');    
    set(handles.frame_play,'Visible','off');        
else
    set(handles.frame_sel,'Visible','on');
    set(handles.frame_num,'Visible','on');    
    set(handles.frame_play,'Visible','on');    
end    

if cin<4
    set(handles.Sliders,'Visible','off');
    set(handles.ImNr,'Visible','off');    
else
    set(handles.Sliders,'Visible','on');
    set(handles.ImNr,'Visible','on');        
end


% Set Input Axes possitions regarding number of input images
if cin == 1 | cin > 3
    for i=1:3
        position_a=(i==1)*[A2_x (A2_y-(i-1)*D2) A2_s]+eps+[0 -.25 0 0];
        position_s=(i==1)*[S2_x (S2_y-(i-1)*D2) S2_s]+eps+[0 -.25 0 0];
        position_n=(i==1)*[N2_x (N2_y-(i-1)*D2) N2_s]+eps+[0 -.25 0 0];

        set(handles.Axes(i),'Position',position_a)
        set(handles.Sliders(i),'Position',position_s)
        set(handles.ImNr(i),'Position',position_n)
    end
    position_fs = [Fs_x Fs_y Fs_s];    
    position_fn = [Fn_x Fn_y Fn_s];    
    position_fp = [Fp_x Fp_y Fp_s];    
    set(handles.frame_sel,'Position',position_fs)    
    set(handles.frame_num,'Position',position_fn)    
    set(handles.frame_play,'Position',position_fp)        
elseif cin == 2
    for i=1:3
        position_a=(i~=3)*[A2_x (A2_y-(i-1)*D2) A2_s]+eps;
        position_s=(i~=3)*[S2_x (S2_y-(i-1)*D2) S2_s]+eps;
        position_n=(i~=3)*[N2_x (N2_y-(i-1)*D2) N2_s]+eps;
        set(handles.Axes(i),'Position',position_a)
        set(handles.Sliders(i),'Position',position_s)
        set(handles.ImNr(i),'Position',position_n)
    end
    position_fs = [Fs_x Fs_y - .2 Fs_s];    
    position_fn = [Fn_x Fn_y - .2 Fn_s];    
    position_fp = [Fp_x Fp_y - .2 Fp_s];    
    set(handles.frame_sel,'Position',position_fs)    
    set(handles.frame_num,'Position',position_fn)    
    set(handles.frame_play,'Position',position_fp)        
elseif cin==3
    for i=1:3
        position_a=[A3_x (A3_y-(i-1)*D3) A3_s];
        position_s=[S3_x (S3_y-(i-1)*D3) S3_s];
        position_n=[N2_x (N3_y-(i-1)*D3) N2_s];
        set(handles.Axes(i),'Position',position_a)
        set(handles.Sliders(i),'Position',position_s)
        set(handles.ImNr(i),'Position',position_n)
    end

    position_fs = [Fs_x Fs_y - .24 Fs_s];    
    position_fn = [Fn_x Fn_y - .24 Fn_s];    
    position_fp = [Fp_x Fp_y - .24 Fp_s];    
    set(handles.frame_sel,'Position',position_fs)    
    set(handles.frame_num,'Position',position_fn)    
    set(handles.frame_play,'Position',position_fp)        
end

% Set Output axes possitions regarding input images and operation
if (cin == 1) | (handles.outs == 0) | (cin > 3)
    for i=1:3
        position_a=(i<2)*[A2_x (A2_y-(i-1)*D2) A2_s]+eps+[0.37 -.25 0 0];
        set(handles.AxesOut(i),'Position',position_a)
    end
elseif (cin == 2) & handles.outs
    for i=1:3
        position_a=(i<3)*[A2_x (A2_y-(i-1)*D2) A2_s]+eps+[0.37 0 0 0];            
        set(handles.AxesOut(i),'Position',position_a)
    end
elseif (cin == 3) & handles.outs
    for i=1:3
        position_a=[A3_x (A3_y-(i-1)*D3) A3_s]+[0.37 0 0 0];
        set(handles.AxesOut(i),'Position',position_a)
    end
end