function FUSEDIMG = mfocus(ORIGIMGS,BLOCKSIZE,TH)
% mfocus Image fusion method based on spatial frequency
% 
%   FUSEDIMG = mfocus(ORIGIMGS,BLOCKSIZE,TH)
%   ORIGIMGS is a 3-D matrix or a Cell array of 2-D matrix
%       If ORIGIMGS is a 3-D matrix, the 3rd dimension is the stack of original images
%       If ORIGIMGS is a Cell array, each Cell is the original images
%   BLOCKSIZE is block size for spatial correlation
%   TH is a threshold for image selection
%
%   FUSEDIMG is a 2-D matrix. It is a grayscale fused image.

% v 2.0 2004/01/19 Eduardo Fernandez Canga - University of Bristol
% -------------------------------------------------------------------------
% if input images are cell, convert to mat
if iscell(ORIGIMGS)
    norig = length(ORIGIMGS);
    [height, width] = size(ORIGIMGS{1});
    temp = ORIGIMGS;
    clear ORIGIMGS
    ORIGIMGS = cell2mat(temp);
    ORIGIMGS = reshape(ORIGIMGS,height,width,norig);
    clear temp
end


if ( (BLOCKSIZE < 1) || (BLOCKSIZE > min(size(ORIGIMGS,1),size(ORIGIMGS,2)) ) )
   error('Wrong blocksize');
end

sx=size(ORIGIMGS,1);           % check row size of the input images
x=mod(sx,BLOCKSIZE);
if x                     % if they are NOT multiple of BLOCKSIZE then
    ORIGIMGS=rowpad(ORIGIMGS,BLOCKSIZE-x);  % add rows to make the size multiple of BLOCKSIZE
end

sy=size(ORIGIMGS,2);           % check col size of the input images
y=mod(sy,BLOCKSIZE);
if y                     % if they are NOT multiple of BLOCKSIZE then
    ORIGIMGS=colpad(ORIGIMGS,BLOCKSIZE-y);  % add cols to make the size multiple of BLOCKSIZE
end
% dimension
[r,c,n] = size(ORIGIMGS);

% transform images in a 'stack of blocks'
% ---------------------------------------
imr=zeros(BLOCKSIZE,BLOCKSIZE,n,r*c/(BLOCKSIZE*BLOCKSIZE));
k=0;
for i=1:BLOCKSIZE:r
    for j=1:BLOCKSIZE:c
        k=k+1;
        imr(1:BLOCKSIZE,1:BLOCKSIZE,:,k)=ORIGIMGS(i-1+(1:BLOCKSIZE),j-1+(1:BLOCKSIZE),:);
    end
end

% Calculates the spacial frequency
% ---------------------------------------
p=spfreq(imr);

% undo transformation
% ---------------------------------------
p = p(:);
p = reshape(p,n,c/BLOCKSIZE,r/BLOCKSIZE);
p = permute(p,[3 2 1]);

% build decision map
% ---------------------------------------
dec=max(p,[],3);
dec2=repmat(dec, [1 1 size(p,3)]);

% thresholding
dec = p >= (dec2-TH);
% expand to original image size
weight = zeros(r,c,n);
for k = 1:n
    curMat = dec(:,:,k);
    curMat = kron(curMat,ones(BLOCKSIZE));
    weight(:,:,k) = reshape(curMat,[r c]);
end

% fuse images
% ---------------------------------------
indz = repmat((sum(weight,3))==0, [1 1 n]);
weight(indz) = 1;
weight = weight./repmat(sum(weight,3), [1 1 n]);
FUSEDIMG = sum(ORIGIMGS.*weight,3);
FUSEDIMG=FUSEDIMG(1:sx,1:sy);