function [imageout]=imblur(imagein)
%imblur.m,v 1.0 2001/10/29 17:35:40$
%===========================================================================
%               Eduardo Fernandez Canga - University of Bath
%
%                        Copyright (c) 2001
%===========================================================================
%
%            [imageout] = imblur(imagein)
%
%	Blurs an area selected by the user
%

[y,x]=size(imagein);
xx=num2str(x);
yy=num2str(y);
r=1;
while r
    x1=[' Enter x1: (1 - ',xx,') '];
    x2=[' Enter x2: (1 - ',xx,') '];
    y1=[' Enter y1: (1 - ',yy,') '];
    y2=[' Enter y2: (1 - ',yy,') '];
    prompt = { x1, x2, y1, y2};
    title = ' Select Area '; 
    lines = 1; 
    def = {num2str(round(x/4)),num2str(round(x/2)),...
            num2str(round(y/4)),num2str(round(y/2))};
    answer = inputdlg(prompt,title,lines,def,'on');
    x1=str2num(answer{1});
    x2=str2num(answer{2});
    y1=str2num(answer{3});
    y2=str2num(answer{4});
    if x1<=x & x2<=x & y1<=y & y2<=y
        r=0;
    end
end

if x1<x2
    incrx=1;
else
    incrx=-1;
end
if y1<y2
    incry=1;
else
    incry=-1;
end
imageout=imagein;
imageout(y1:incry:y2,x1:incrx:x2)=...
    conv2(impad(imagein(y1:incry:y2,x1:incrx:x2),6),gmask(2),'valid');
