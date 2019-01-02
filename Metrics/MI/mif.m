function [QI,QInorm]=mif(ORIGIMGS,FUSEDIMG)
% mif - Mutual information for fusion assessment
% 
%   [QI,QInorm]=mif(ORIGIMGS,FUSEDIMG)
%   ORIGIMGS is a 3-D matrix or a Cell array of 2-D matrix
%       If ORIGIMGS is a 3-D matrix, the 3rd dimension is the stack of original images
%       If ORIGIMGS is a Cell array, each Cell is the original images
%   FUSEDIMG is a 2-D matrix. It is a grayscale fused image.
%   
%   QI is mutual information
%   QInorm is normalized mutual information (see An Overlap invariant entropy
%       measure of 3D medical image alignment, by C. Sudholme et al for more information)

% v1.0 2004/11/14 Eduardo Fernandez Canga University of Bristol
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

if any(size(ORIGIMGS(:,:,1)) ~= size(FUSEDIMG))
    error('Wrong image size')
end

QI = 0;
QInorm = 0;
for i = 1:size(ORIGIMGS,3)
    % AL:
    if nargout > 1
        [mitemp,minormtemp]=MI_norm(ORIGIMGS(:,:,i),FUSEDIMG);
        QInorm = QInorm + minormtemp;
    else
        mitemp=MI_norm(ORIGIMGS(:,:,i),FUSEDIMG);
    end
    QI = QI     + mitemp;
end
QI     = QI     / size(ORIGIMGS,3);
QInorm = QInorm / size(ORIGIMGS,3);