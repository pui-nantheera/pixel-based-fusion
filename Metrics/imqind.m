function [QI, QI_map, S] = imqind(ORIGIMGS, FUSEDIMG, BLOCKSIZE)
% imqind - A universal image quality index.
%
%   [QI, QI_map, S] = imqind(ORIGIMGS, FUSEDIMG, BLOCKSIZE)
%   ORIGIMGS is a 3-D matrix or a Cell array of 2-D matrix
%       If ORIGIMGS is a 3-D matrix, the 3rd dimension is the stack of original images
%       If ORIGIMGS is a Cell array, each Cell is the grayscale original images
%   FUSEDIMG is a 2-D matrix. It is a grayscall fused image.
%   BLOCKSIZE is a block size for calculating spatial similarity (default = 8)
%   
%   QI is quality index (range [-1,1])
%   QI_map is quality Map of the Fused Image.
%   S is Salient information
%
%   Ref: Zhou Wang and A.C. Bovik. A universal image QI index. Signal Processing Letters, IEEE, 9(3):81{84, Mar.2002.

% -------------------------------------------------------------------------

% check input
if nargin < 3
    BLOCKSIZE = 8;
end

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

% number of original frames
n = size(ORIGIMGS,3);

% window 
N = BLOCKSIZE.^2;
% for summation in the window
sum2_filter = ones(BLOCKSIZE);

% -------------------------------------------------------------------------
% calculate IFQI
% -------------------------------------------------------------------------
imi_sq = ORIGIMGS .^2;
imf_sq = FUSEDIMG .^2;
imif   = ORIGIMGS .* FUSEDIMG(:,:,ones(1,n)); %every input image is multiplied by the fused image
% summation in window
imf_sum    = filter2(sum2_filter, FUSEDIMG, 'valid');
imf_sq_sum = filter2(sum2_filter, imf_sq, 'valid');
for i = 1 : n
    imi_sum(:,:,i)    = filter2(sum2_filter, ORIGIMGS(:,:,i),    'valid'); 
    imi_sq_sum(:,:,i) = filter2(sum2_filter, imi_sq(:,:,i), 'valid');
    imif_sum(:,:,i)   = filter2(sum2_filter, imif(:,:,i),   'valid');
end

imif_sum_mul    = imi_sum.*imf_sum(:,:,ones(1,n));
imif_sq_sum_mul = imi_sum.^2 + imf_sum(:,:,ones(1,n)).^2;

numerator    = 4*(N*imif_sum - imif_sum_mul).*imif_sum_mul;
denominator1 = N*(imi_sq_sum + imf_sq_sum(:,:,ones(1,n))) - imif_sq_sum_mul;
denominator  = denominator1.*imif_sq_sum_mul;

% Quality Map of the Fused Image
QI_map = ones(size(denominator));
index = (denominator1 == 0) & (imif_sq_sum_mul ~= 0);
QI_map(index) = 2*imif_sum_mul(index)./(imif_sq_sum_mul(index)+eps);
index = (denominator ~= 0);
QI_map(index) = numerator(index)./(denominator(index)+eps);
% Image fusion QI index
QI = mean2(QI_map);

if nargout > 2
    % salient information
    S = imi_sq_sum/N - (imi_sum.*imi_sum)/N^2;
end