function FUSEDIMG = fuse_avg(ORIGIMGS, OPTIONS)
% fuse_avg - Fusion by averaging or selection
% 
%   FUSEDIMG = fuse_avg(ORIGIMGS, OPTIONS)
%   ORIGIMGS is a 3-D matrix or a Cell array of 2-D matrix
%       If ORIGIMGS is a 3-D matrix, the 3rd dimension is the stack of original grayscale images
%       If ORIGIMGS is a Cell array, each Cell is the original grayscale images
%   OPTIONS switch for selection type
%    	OPTIONS == 0: average (default)
%    	OPTIONS == 1: select first input
%      	OPTIONS == 2: select second input
%    	OPTIONS == n: select N th input
%
%   FUSEDIMG is a 2-D matrix. It is a grayscale fused image.

%   v 1.0 16.08.99  Oliver Rockinger
%   v 2.0 12.11.03  Eduardo Fernandez Canga (University of Bristol)
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

% set default
if nargin < 2 || isempty(OPTIONS)
    OPTIONS = 0;
end
% read number of images
dimens= size(ORIGIMGS,3);

switch 1
    case OPTIONS==0,
        % average all images
        FUSEDIMG = mean(ORIGIMGS,3);
        
    case OPTIONS>0 & OPTIONS<=dimens & OPTIONS==round(OPTIONS),
        % select one image
        FUSEDIMG = ORIGIMGS(:,:,OPTIONS);
        
    otherwise, error('unknown option');
end;
