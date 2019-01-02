function y=vol2avi(x, filename, rate)
%vol2avi.m, v 0.1 2003/12/07
%===========================================================================
%               Eduardo Fernandez Canga - University of Bristol
%
%                        Copyright (c) 2003
%===========================================================================
%     function y = vol2avi(x,filename)
if nargin < 3
    rate = 3;
end

%x=min(max(round(x),0),255);

[a b c d e] = size(x);
for j=1:e   
    if e>1
        file=[filename(1:end-3) num2str(j) '.avi'];
    else
        file=filename;
    end
    mov=avifile(file,'FPS',rate,'COMPRESSION','NONE','QUALITY',100);
    for i=1:d
        f=im2frame(x(:,:,:,i,e));
        mov=addframe(mov,f);
    end
    mov=close(mov);
end
y=0;