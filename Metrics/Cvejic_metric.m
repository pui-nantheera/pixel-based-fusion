function QI=Cvejic_metric(ORIGIMGS,FUSEDIMG,BLOCKSIZE)
% Cvejic_metric - A similarity metric for assessment of image fusion algorithms by Cvejic et al. (University of Bristol)
% 
%   QI = Cvejic_metric(ORIGIMGS,FUSEDIMG,BLOCKSIZE)
%   ORIGIMGS is a 3-D matrix or a Cell array of 2-D matrix
%       If ORIGIMGS is a 3-D matrix, the 3rd dimension is the stack of original images
%       If ORIGIMGS is a Cell array, each Cell is the original images
%   FUSEDIMG is a 2-D matrix. It is a grayscale fused image.
%   BLOCKSIZE is a block size (integer) used for calculating spatial similarity (default BLOCKSIZE = 4)
%   QI is quality score (0-1)
%
%   Ref: Nedeljko Cvejic, Artur Loza, David Bull, and Nishan Canagarajah. A similarity metric for assessment of image
%   fusion algorithms. International Journal of Signal Processing, 2(3):178-182, 2005.

% v1.0 2005 Cvejic University of Bristol
% v1.1 17.09.12 Pui University of Bristol: Comments and speed
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

% input images
ORIGIMGS=double(ORIGIMGS);
FUSEDIMG=double(FUSEDIMG);
numOrigImg = size(ORIGIMGS,3);

% check inputs
if nargin < 3
    BLOCKSIZE=4;
end

% numbers of blocks
s=size(ORIGIMGS);
a=floor(s(1)/BLOCKSIZE);
b=floor(s(2)/BLOCKSIZE);

% similarity in spatial domain between the input image and the fused image
% -------------------------------------------------------------------------
simw = zeros(numOrigImg,a*b);
down = zeros(numOrigImg,a*b) + eps;
for numf = 1:numOrigImg
    count = 1;
    % for each block
    for i=1:a,
        for j=1:b,
            % current block
            x=ORIGIMGS((i-1)*BLOCKSIZE+1:i*BLOCKSIZE,(j-1)*BLOCKSIZE+1:j*BLOCKSIZE, numf);
            z=FUSEDIMG((i-1)*BLOCKSIZE+1:i*BLOCKSIZE,(j-1)*BLOCKSIZE+1:j*BLOCKSIZE);
            x=reshape(x',1,BLOCKSIZE^2);
            z=reshape(z',1,BLOCKSIZE^2);
            mean_x=mean(x);
            mean_z=mean(z);
            % find covarience
            sigma_xz=(1/(BLOCKSIZE^2-1))*sum((x-mean_x).*(z-mean_z));
            down(:,count) = down(:,count) + sigma_xz;
            simw(numf,count) = sigma_xz;
            count = count + 1;
        end
    end
end
simw = min(1, max(0, simw./down));


% find Universal Image Quality Index (Wang and Bovik): between x and z
% -------------------------------------------------------------------------
QX = zeros(numOrigImg,a*b);
for numf = 1:numOrigImg
    count = 1;
    for i=1:a,
        for j=1:b,
            x=ORIGIMGS((i-1)*BLOCKSIZE+1:i*BLOCKSIZE,(j-1)*BLOCKSIZE+1:j*BLOCKSIZE , numf);
            z=FUSEDIMG((i-1)*BLOCKSIZE+1:i*BLOCKSIZE,(j-1)*BLOCKSIZE+1:j*BLOCKSIZE);
            x=reshape(x',1,BLOCKSIZE^2);
            z=reshape(z',1,BLOCKSIZE^2);
            mean_x=mean(x);
            mean_z=mean(z);
            var_x=std(x)^2;
            var_z=std(z)^2;
            sigma_xz=(1/(BLOCKSIZE^2-1))*sum((x-mean_x).*(z-mean_z));
            up=4*sigma_xz*mean_x*mean_z;
            down=(var_x+var_z)*(mean(x)^2+mean_z^2)+eps;
            QX(numf,count) = up/down;
            count = count + 1;
        end
    end
end

% quality index from weighting average
Q= sum(QX.*simw);
QI=mean(Q);