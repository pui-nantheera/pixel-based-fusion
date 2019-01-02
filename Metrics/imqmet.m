function [QI, QI_map] = imqmet(ORIGIMGS, FUSEDIMG, BLOCKSIZE, CEDGE)
% imqmet - Piella's Quality Index for fusion assessment
% 
%   [QI, QI_map] = imqmet(ORIGIMGS, FUSEDIMG, BLOCKSIZE, CEDGE)
%   ORIGIMGS is a 3-D matrix or a Cell array of 2-D matrix
%       If ORIGIMGS is a 3-D matrix, the 3rd dimension is the stack of original images
%       If ORIGIMGS is a Cell array, each Cell is the grayscale original images
%   FUSEDIMG is a 2-D matrix. It is a grayscale fused image.
%   BLOCKSIZE is a block size for calculating spatial similarity (default = 8)
%   CEDGE is edge coefficient (default = 0.2)
%   
%   QI is quality index
%   QI_map is quality Map of the Fused Image.

% v 0.2 2004/04/19 Eduardo Fernandez Canga University of Bristol
% -------------------------------------------------------------------------

% if input images are cell, convert to mat
if iscell(ORIGIMGS)
    norig = length(ORIGIMGS);
    [height, width, depth] = size(ORIGIMGS{1});
    if (depth>1)
        for k = 1:length(ORIGIMGS)
            ORIGIMGS{k} = rgb2gray(ORIGIMGS{k});
        end
    end
    temp = ORIGIMGS;
    clear ORIGIMGS
    ORIGIMGS = cell2mat(temp);
    ORIGIMGS = reshape(ORIGIMGS,height,width,norig);
    clear temp
end

if size(FUSEDIMG,3) > 1
    FUSEDIMG = rgb2gray(FUSEDIMG);
end

% check input
if nargin < 3
    BLOCKSIZE = 8;
end
if nargin < 4
    CEDGE = 0.2;
end
if (nargin == 1 || nargin > 4)
   QI = -Inf;
   QI_map = -1*ones(size(FUSEDIMG));
   return;
end

% input images' dimension
[r c n] = size(ORIGIMGS);

% check fused image
if [r c] ~= size(FUSEDIMG)
   QI = -Inf;
   QI_map = -1*ones(size(FUSEDIMG));
   return;
end

% Compute image fusion QI index with window size = blocksize
[~, QI_map, s] = imqind(ORIGIMGS, FUSEDIMG, BLOCKSIZE);

% salient information (s) and weight (lambda)
s_sum=sum(s,3);
lambda=s./(s_sum(:,:,ones(1,n))+eps);
% overall salience (w)
pw=1;
w=max(s,[],3);
sumw=sum(w(:));
if sumw ~=0
    w=w./sum(w(:));
    w=w.^pw;
end

QI_map=w.*sum(lambda.*QI_map,3);
QI = sum(QI_map(:));

if CEDGE ~=0
    % find `edge images' and their QI
    eims = sobel(ORIGIMGS);
    eimf = sobel(FUSEDIMG);
    [~, equality_map, s] = imqind(eims, eimf, BLOCKSIZE);

    s_sum=sum(s,3);
    lambda=s./(s_sum(:,:,ones(1,n))+eps);
    
    pw=1;
    w=max(s,[],3);
    sumw=sum(sum(w));
    if sumw ~=0
        w=w./sum(sum(w));
        w=w.^pw;
    end

    equality_map=w.*sum(lambda.*equality_map,3);
    equality = sum(equality_map(:));
    
    % same measure is computed with the `edge images'
    % for CEDGE = 0:.2:1    
    QI = QI ^ (1-CEDGE) * (equality ^ CEDGE);
end