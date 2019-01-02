function mymetric2(h, handles)
%mymetric1.m, v 0.1 2003/12/07

v1='off';
v2='off';
vf='off';
    
if handles.reference & ...
    all(size(handles.im1)==size(handles.im2)) & ...
    all(size(handles.im1)==size(handles.im))
    
    if (isfield(handles,'im1'))
        RMSE=rmse(handles.im,handles.im1);
        RMSE=['RMSE: ',num2str(RMSE)];
        set(handles.text_rmse1,'string',RMSE);
        v1='on';
    end
    if (isfield(handles,'im2'))
        RMSE=rmse(handles.im,handles.im2);
        RMSE=['RMSE: ',num2str(RMSE)];
        set(handles.text_rmse2,'string',RMSE);
        v2='on'
    end
    if (isfield(handles,'fusion')) & all(size(handles.im)==size(handles.fusion))
        RMSE=rmse(handles.im,handles.fusion);
        RMSE=['RMSE: ',num2str(RMSE)];
        set(handles.text_rmsef,'string',RMSE);
        vf='on';
        
        q=petmetric(handles.im1b,handles.im2b,handles.fusion);
        pet=['Quality: ',num2str(q)];
        set(handles.petmetric,'string',pet);
    end
elseif (isfield(handles,'fusion')) & ...
        all(size(handles.im1)==size(handles.fusion)) & ...
        all(size(handles.im2)==size(handles.fusion))
    
    RMSE=rmse(handles.fusion,handles.im1);
    RMSE=['RMSE: ',num2str(RMSE)];
    set(handles.text_rmse1,'string',RMSE);
    RMSE=rmse(handles.fusion,handles.im2);
    RMSE=['RMSE: ',num2str(RMSE)];
    set(handles.text_rmse2,'string',RMSE);
    v1='on';
    v2='on';
    %q=petmetric(handles.im1b,handles.im2b,handles.fusion);
    %pet=['Quality: ',num2str(q)];
    %set(handles.petmetric,'string',pet);
end

set(handles.text_rmse1,'visible',v1)
set(handles.text_rmse1b,'visible',v1)
set(handles.text_rmse2,'visible',v2)
set(handles.text_rmse2b,'visible',v2)
set(handles.text_rmsef,'visible',vf)
set(handles.text_rmsefb,'visible',vf)
