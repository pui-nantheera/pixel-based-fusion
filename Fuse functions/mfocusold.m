function [imageout] = mfocus_1(im,bs,TH)
%
%mfocus.m, v 1.3 2003/11/13 16.52
%===========================================================================
%               Eduardo Fernandez Canga - University of Bristol
%
%                        Copyright (c) 2003
%===========================================================================
%
%           [imageout] = mfocus(im,bs,TH)
%
%                 im: input images
%                 bs: block size
%                 TH: threshold
%
%===========================================================================

if size(im,3)>2
    error('Method available only for two input images')
    return
end

if ( (bs < 1) | (bs > min(size(im,1),size(im,2)) ) )
   error(sprintf('Wrong blocksize'));
end

sx=size(im,1);           % check row size of the input images
x=mod(sx,bs);
if x                     % if they are NOT multiple of bs then
    im=rowpad(im,bs-x);  % add rows to make the size multiple of bs
end

sy=size(im,2);           % check col size of the input images
y=mod(sy,bs);
if y                     % if they are NOT multiple of bs then
    im=colpad(im,bs-y);  % add cols to make the size multiple of bs
end

[r,c,n] = size(im);

imr=zeros(bs,bs,n,r*c/(bs*bs));

% transform images in a 'stack of blocks'
k=0;
for i=1:bs:r
    for j=1:bs:c
        k=k+1;
        imr(1:bs,1:bs,:,k)=im(i-1+(1:bs),j-1+(1:bs),:);
    end
end

p=spfreq(imr);

% undo transformation
p=reshape([p(1:2:end) p(2:2:end)],c/bs,r/bs,n);
p=permute(p,[2 1 3]);
% build decision map
decision = -1 * (p(:,:,2) > (p(:,:,1) + TH)) + (p(:,:,1) > (p(:,:,2) + TH));

%if consistency==1               % apply consistency check if it is active
%    mask=ones(3,3)/9;           % use 3x3 window mask
%    decision=impad(decision,1);
%    decision=round(conv2(decision,mask,'valid'));
%end

% expand the block decision mask to a pixel decision mask
% dec=kron(decision,ones(bs));
[maxi,maxj] = size(decision);
dec=zeros(maxi*bs,maxj*bs);
for i = 1:bs
   for j = 1:bs
       dec(i:bs:bs*maxi,j:bs:bs*maxj)=decision;      
  end
end



imageout=im(:,:,1).*(dec==1)+im(:,:,2).*(dec==-1)+(im(:,:,1)+im(:,:,2))/2.*(dec==0);
imageout=imageout(1:sx,1:sy);