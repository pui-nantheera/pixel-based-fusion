function QI = petmetric(ORIGIMGS,FUSEDIMG,SIGMOIDS,SIGMOIDA)
% petmetric - Petrovic and Xydeas Metric
% 
%   QI = petmetric(ORIGIMGS,FUSEDIMG,[SIGMOIDS],[SIGMOIDA])
%   ORIGIMGS is a 3-D matrix or a Cell array of 2-D matrix
%       If ORIGIMGS is a 3-D matrix, the 3rd dimension is the stack of original grayscale images
%       If ORIGIMGS is a Cell array, each Cell is the original grayscale images
%   FUSEDIMG is a 2-D matrix. It is a grayscale fused image.
%   SIGMOIDS - [sigma k] values that defines the sigmoid for edge strenght (optional)
%              default: sigma = 0.7, k = 11
%   SIGMOIDA - [sigma k] values that defines the sigmoid for edge angle (optional)
%              default: sigma = 0.8, k = 24
%
%   QI is quality score

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

% number of original images
n = size(ORIGIMGS,3);

% Calculating non-linear parameters
if nargin < 3 
    sigmag = .7;
    kg = 11;
else
    sigmag = SIGMOIDS(1);
    kg = SIGMOIDS(2);
end
% parameter for edge strength
gamag=(1+exp(-kg*(1-sigmag)));

if nargin < 4
    sigmaa=.8;
    ka=24;
else
    sigmaa = SIGMOIDA(1);
    ka = SIGMOIDA(2);
end
% parameter for edge orientation
gamaa=(1+exp(-ka*(1-sigmaa)));

% edge information
[gi,ai]=sobel(ORIGIMGS);
[gf,af]=sobel(FUSEDIMG);

% relative strength and orientation change values
Gain=(min(gi,gf(:,:,ones(1,n)))+eps)./(max(gi,gf(:,:,ones(1,n)))+eps);
Gain(isnan(Gain))=0;
Alpha=abs(abs(ai-af(:,:,ones(1,n)))-(pi/2))/(pi/2);

% edge strength and orientation preservation values
QGain=gamag./(1+exp(-kg*(Gain-sigmag)));
QAlpha=gamaa./(1+exp(-ka*(Alpha-sigmaa)));

% overall edge information preservation values
QI=sqrt(QGain.*QAlpha);

% normalised weighted performance metric
QI=QI.*gi;
QI=sum(QI(:))/sum(gi(:));