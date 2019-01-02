function [fusion] = wavelet(im1,im2,nscales,consistency,wavfile,reflect)
%
%wavelet.m, v 1.1 2002/02/06 18:08:40
%===========================================================================
%               Eduardo Fernandez Canga - University of Bath
%
%                        Copyright (c) 2002
%===========================================================================
%
%        [fusion] = wavelet(im1,im2,nscales,consistency,wavfile,reflect)
%
%                   im1: input image
%                   im2: input image
%               nscales: number of scales (default 4)
%           consistency: apply consistency
%                           1=yes (default)
%                           0=no 
%               wavfile: File containing the coef (default 'bi97.wvf')
%               reflect: Edge handling: 
%                           0=wrap-around (default)
%                           1=try reflection 
%
%===========================================================================
if nargin<6, reflect=0;end
if nargin<5, wavfile='bi97.wvf';end
if nargin<4, consistency=1;end
if nargin<3, nscales=4;end

if any (size(im1)~=size(im2))
    error('Error: Different Size Images')
end

r=2^nscales;
sx=size(im1,1);
x=mod(sx,r);
if x
    im1=rowpad(im1,r-x);
    im2=rowpad(im2,r-x);
end
sy=size(im1,2);
y=mod(sy,r);
if y
    im1=colpad(im1,r-y);
    im2=colpad(im2,r-y);    
end

% Forward transform 
%
u=wt2dscl(im1, wavfile, reflect, nscales);
v=wt2dscl(im2, wavfile, reflect, nscales);

% Fusion 
%
decision=abs(v)>abs(u);

if consistency==1
    mask=ones(3,3)/9;
    decision=impad(double(decision),1);
    decision=round(conv2(decision,mask,'valid'));
end

w = v.*decision + u.*(~decision);


% Reverse transform
%   
fus=iwt2dscl(w, wavfile, reflect, nscales);

fusion=fus(1:sx,1:sy);