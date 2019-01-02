function FUSEDIMG = fuse_pca(ORIGIMGS)
% fuse_pca - image fusion with PCA method
% 
%   FUSEDIMG = fuse_pca(ORIGIMGS)
%   ORIGIMGS is a 3-D matrix or a Cell array of 2-D matrix
%       If ORIGIMGS is a 3-D matrix, the 3rd dimension is the stack of original grayscale images
%       If ORIGIMGS is a Cell array, each Cell is the original grayscale images
%
%   FUSEDIMG is a 2-D matrix. It is a grayscale fused image.

%   v 1.0 16.08.99  Oliver Rockinger
%   v 2.0 14.08.12  N. Anantrasirichai (University of Bristol)
%                   Modified for N input images
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

% dimension
[height, width, totalframe] = size(ORIGIMGS);
Mreshape = reshape(ORIGIMGS, height*width, totalframe);

% apply PCA
[COEFF,~,latent] = princomp(Mreshape);

% find maximum of eigenvalues of the covariance matrix of data
[~, ind] = max(latent);

% find weight
w(1,1,:) = COEFF(:,ind)/sum(COEFF(:,ind));
w = repmat(w, [height width 1]);

% fusing
FUSEDIMG = sum(w.*ORIGIMGS,3);