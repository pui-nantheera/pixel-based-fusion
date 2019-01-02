function [fusion] = cemif(im,op,r,consistency)
%
%cemif.m, v 1.2 2003/11/14 18:01:40
%===========================================================================
%               Eduardo Fernandez Canga - University of Bristol
%
%                        Copyright (c) 2003
%===========================================================================
%
%                [fusion] = cemif(im,op,r,consistency)
%
%                im: input images (3D)
%                op: selection of background fussion method
%                     0 - both backgrounds 
%                     1 - first image background  (default)
%                     2 - second image background
%                 r: averaging mask size (default 11)
%       consistency: apply consistency
%                       1=yes (default)
%                       0=no 
%
%===========================================================================

if nargin < 4, consistency=1;end
if nargin < 3, r=11;end
if nargin < 2, op=1;end

%Calculating backgrounds
back=imaver(im,r);

%Calculating foregrounds
fore=im-back;

%Fusion of foregrounds
%decision=abs(fore1)>abs(fore2);
%if consistency==1               % apply consistency check if it is active
%    mask=ones(3,3)/9;           % use 3x3 window mask
%    decision=impad(double(decision),1);
%    decision=round(conv2(decision,mask,'valid'));
%end
%fore=fore1.*decision+fore2.*(~decision);
fore=selc(fore,1);
%Fusion of backgrounds

switch 1
case op==0 % Convine all backgrounds
    m=mean(mean(im));
    back=sum(back,3)-sum(m,3)/size(m,3);
case op>0 & op<=size(im,3) & op==round(op) % Selected a single background
    back=back(:,:,op);
otherwise
    error('Invalid input for background fussion method');
end

fusion=back+fore;