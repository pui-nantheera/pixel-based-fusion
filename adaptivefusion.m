function [inputmetrics,fusemetrics,history,v,param]=adaptivefusion(im,imfa,framenr,inputmetrics,fusemetrics,history,param)

    % compute metrics from input frames
    
    % analyse imfa - fusion history (1,i-1)
    
    % select fusion method and create params
if mod(framenr,2)
    v = 1;
else
    v=6;    
    paramst=cell2struct(param(:,2),param(:,1));
    paramst.lvl=2;
    paramst.high=[1 0];
    paramst.low=0;
    param(:,2)=struct2cell(paramst);
end

history{end,1}=v;
history{end,2}=param;
history{end+1,1}=[];