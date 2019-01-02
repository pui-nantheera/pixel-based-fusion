function I = rgb2grayn(X)
%rgb2grayn.m, v 0.1 2004/05/01
%===========================================================================
%               Eduardo Fernandez Canga - University of Bristol
%
%                        Copyright (c) 2004
%===========================================================================
%     function I = rgb2grayn(X)
%
%  Inputs:   X - N-Dimension matrix of rgb set of images. RGB components must
%                be in the first three dimension
%
% Outputs:   I - N-Dimension matrix of gray imges.
%
%===========================================================================
%
%  Example:
%      
%        X = imread('your_image1.jpg');
%        X(:,:,:,1,2) = imread('your_image2.jpg');
%        % size(X) would be  [R  C  3  1  2] where R and C are the number
%        % of rows and columns in the images
%        I = rgb2grayn(X)
%        % size(I) would be  [R  C  1  1  2] the third dimension is now
%        % size 1 and the images are grayscale images.
%
%   Function developed for the IFT.
%
%   See also IFT (Add_Img_Callback), LOADIM.
origSize=size(X);

% Determine if input includes a 3-D array 
threeD = (size(X,3)==3);

% Calculate transformation matrix
T = inv([1.0 0.956 0.621; 1.0 -0.272 -0.647; 1.0 -1.106 1.703]);
coef = T(1,:)';

if threeD
    %RGB
    % Shape input matrix so that it is a n x 3 array and initialize output
    % matrix  
    X = permute (X,[1,2,4:ndims(X),3]);
    
    X = reshape(X(:),[prod(origSize)/3,3]);
    sizeOutput = [origSize(1), origSize(2), 1, origSize(4:end)];
    
    % Do transformation
    if isa(X, 'double')
        I = X*coef;
        I = min(max(I,0),1);
    else
        error('Wrong input: X must be double')
    end
else
    error('Wrong input: input must be a RGB set of images')
end

I = reshape(I,sizeOutput); 