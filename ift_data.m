function varargout = ift_data(varargin)
%ift_data.m, v 0.1
%=========================================================================
%
%                       IMAGE FUSION TOOLBOX
%
%=========================================================================
%
%        This subroutine handles the IFT load/save data operations.
%
%=========================================================================
%
%                      Eduardo Fernandez Canga
%
%                       University of Bristol 
%
%                     Copyright (c) 2003 - 2005
%=========================================================================

try 
    [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
catch
    errordlg(lasterr,'OPEN ERROR','modal')
    disp(lasterr)
end

% --------------------------------------------------------------------
%  Load_Image: Callback
% --------------------------------------------------------------------
function Load_Image(h, hs)
global ims imout time

% select and read the images
if length(hs.imgpath)
    [files,pathname] = uigetfile('*.*','Select Input Images',hs.imgpath );
else    
    [files,pathname] = uigetfile('*.*','Select Input Images');
end
files = {files};    
if pathname==0
    return
end
hs.imgpath=pathname;

for i=1:size(files,2)
    filename = files{i};
    M1=0;
    if filename~=0
        M1=LoadIm([pathname filename]);
    else
        return
    end

    if numel(M1)==1
        continue
    end

    % add it to the 3D matrix
    if (~any(ims(:)))
        c=0;
        ims=M1;
    else
        %check size
        if any(size(M1) ~= size(ims(:,:,:,1,1))) && any(ims(:))
            button = questdlg(sprintf(['Wrong size: ' filename '\n\n'...
             'This image has different size than the existing ones.\n\n'...
             'Do you want to discard it or start a new set of images?']),...
             'Warning','New Set','Discard','Discard');
            if strcmp(button,'New Set')
                c=0;
                ims=M1;
                
                ift_hist('Add_Hist',h,hs, 'DELE Input 0')
                
            else
                continue
            end
        else
            c=size(ims,5);
            ims(:,:,:,:,c+1)=M1;
        end
    end
    hs.imgpath=pathname;
    
    ift_hist('Add_Hist',h,hs, ['LOAD ' pathname,filename])
    
    guidata(h,hs)
    %update image
    if c<3
        set(hs.Sliders(c+1),'Value',c)
        ift('Refresh_All',h, hs,c+1)
    end
    ift('change_plot',hs)
end

% --------------------------------------------------------------------
%  Load_Sequence: Callback
% --------------------------------------------------------------------
function Load_Sequence(h, hs)
global ims imout time

% select and read the images
if length(hs.imgpath)
    [files,pathname] = uigetfile('*.*','Select Input Images',hs.imgpath );
else    
    [files,pathname] = uigetfile('*.*','Select Input Images');
end

if pathname==0
    return
end
hs.imgpath=pathname;

p = (numel(ims)>1)*size(ims,5); % last modality/sensor number

for i=1:size(files,2)
    if strcmp(files{i}((end-2):end),'avi')
        done = 1;
        y=avi2vol([pathname files{i}]);
        if (~any(ims(:)))
            ims=y;
            time=(1:size(y,4))';
        else 
            ims(:,:,:,:,p+i)=y;
            time(:,p+i)=(1:size(y,4))';
        end
        hs.imgpath=pathname;
        
        ift_hist('Add_Hist',h,hs, ['LOAV seq = ' num2str(p+i) '; '  pathname,files{1}])
       
        guidata(h,hs)
    else
        done = 0;
        break
    end
    ift('change_plot',hs)
    ift('Refresh_All',h,hs)
end

if done
    return
end

for i=1:size(files,2)
    filename = files{i};
    M1=0;
    if filename~=0
        M1=LoadIm([pathname filename]);
    else
        return
    end
    
    if numel(M1)==1
        continue
    end
    
    % add it to the 5D matrix
    if (~any(ims(:)))
        ims=M1;
    else
        %check size
        if any(size(M1) ~= size(ims(:,:,:,1,1))) && any(ims(:))
            error('Wrong sequence: images have different size');    
        else
            ims(:,:,:,i,p+1)=M1;
        end
    end
    hs.imgpath=pathname;
    
    ift_hist('Add_Hist',h,hs, ['LOAS seq = ' num2str(p+1) '; '  pathname,filename])
    
    guidata(h,hs)
end
time(:,p+1)=(1:size(ims(:,:,:,:,p+1),4))';

ift('change_plot',hs)
ift('Refresh_All',h,hs)

% --------------------------------------------------------------------
%  Save Output/Input Video : Callback
% --------------------------------------------------------------------
function Save_Video(h,hs,io)
global imout ims
% io = 1 input images
% io = 0 output images

v = get(gcbf,'CurrentAxes');  % gets handle of selected image
if io
    r = hs.Axes==v;    % compares with all handles of input image
else
    r = hs.AxesOut==v; % compares with all handles of output image
end
v = str2double(char(get(hs.ImNr(r),'String')));  % reads number of image from
                                                 % its block number

[filename, pathname] = uiputfile([hs.imgpath '*.avi'],...
    'Save Video');
if filename
    if io
        vol2avi(ims(:,:,:,:,v),[pathname filename]);
    else
        vol2avi(imout(:,:,:,:,v),[pathname filename]);
    end
end

% --------------------------------------------------------------------
%  Save Output Image
% --------------------------------------------------------------------
function Save_Image(h, hs, save_all)
global imout

if save_all % save_all = 1 if the function is called from the 'File' menu
            % save_all = 0 if the function is called by right click in 
            % specific image
    n=1:size(imout,5);
else
    v=get(gcbf,'CurrentAxes');
    r=hs.AxesOut==v;
    n=str2double(char(get(hs.ImNr(r),'String')));
end
try
    if numel(imout)==1
        errordlg('There are no output images','SAVE ERROR','modal')
        return
    end
    for i=n
        pathname = hs.imgpath;
        if size(imout,4) == 1 % single images
            [filename, pathname] = uiputfile([pathname '*.*'], ...
                ['Save Output Image ' num2str(i)]);
            if filename==0
                return
            end
            file=sprintf('%s%s',pathname,filename);
            imwrite((imout(:,:,:,:,i)),file);
        else % sequences
            [filename, pathname] = uiputfile([pathname '*.*'], ...
                ['Save Output Sequence ' num2str(i)]);
            if filename==0
                return
            end
            for j=1:size(imout,4)
                file=[pathname filename(1:find(filename=='.')-1) ...
                    sprintf('%.2d',j) filename(find(filename=='.'):end)];                
                imwrite((imout(:,:,:,j,i)),file);
            end
        end
    end
    guidata(h,hs) 
    
catch
    errordlg(lasterr,'SAVE ERROR','modal')
end