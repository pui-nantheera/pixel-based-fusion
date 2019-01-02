function y=avi2vol(filename,start,finish)
%avi2vol.m, v 0.1 2004/02/18
%===========================================================================
%               Eduardo Fernandez Canga - University of Bristol
%
%                        Copyright (c) 2004
%===========================================================================
%     function y = avi2vol (filename)

mov=aviread(filename);
mov_inf=aviinfo(filename);
if nargin < 2
    finish = mov_inf.NumFrames;
    start = 1;
elseif nargin < 3
    finish = mov_inf.NumFrames;
end
if finish > mov_inf.NumFrames
    finish = mov_inf.NumFrames;
end
    

if strcmp(mov_inf.ImageType,'truecolor')
    for i=start:finish
        y(:,:,:,i)=double((mov(i).cdata))/255;
    end
elseif strcmp(mov_inf.ImageType,'indexed')
    for i=start:finish
        y(:,:,:,i)=double(ind2rgb(mov(i).cdata,mov(i).colormap));
    end
elseif strcmp(mov_inf.ImageType,'grayscale')
    disp('avi2vol.m... if this works remove this line')
    if isfield(r,'MaxValue')
        p=r.MaxValue;
    elseif isfield(r,'BitDepth')
        p=2^r.BitDepth - 1;
    end
    for i=start:finish
        y(:,:,:,i)=double(mov(i).cdata)/p;
    end        
end    
        