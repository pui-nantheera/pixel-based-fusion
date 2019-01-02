function FUSEDCOEF = fuse_dwb(ORIGIMGS, LEVELS, HIGH, LOW)
%   fuse_dwb - Image fusion with Discrete Wavelet Transform
%
%   FUSEDIMG = fuse_dwb(ORIGIMGS, LEVELS, HIGH, LOW)
%   ORIGIMGS is a 3-D matrix: input images
%   LEVELS is maximum decomposition level   (default = 4)
%   HIGH is coefficient selection highpass  (default = 1) (see selc.m) 
%   LOW is coefficient selection base image (default = 0) (see selb.m) 
%
%   FUSEDIMG is a 2-D matrix. It is a grayscale fused image.
%
%   Examples:   I(:,:,1) = im2double(imread('UNcamp_2N_i_20.bmp'));
%               I(:,:,2) = im2double(imread('UNcamp_2N_v_20.bmp'));
%               levels = 3; % decomposition levels
%               high = 1;   % select max absolute coefficients
%               low = 0;    % average base level
%               F = fuse_dwb(I, levels, high, low);
%               figure; imshow(F)

%   v 1.0 16.08.99  Oliver Rockinger
%   v 2.0 12.11.03  Eduardo Fernandez Canga (University of Bristol)
%                   Modified for N input images
% -------------------------------------------------------------------------
%%  check inputs and put default values
if nargin < 2 || isempty(LEVELS)
    LEVELS = 4;
end
if nargin < 3 || isempty(HIGH)
    HIGH = 1;
end
if nargin < 4 || isempty(LOW)
    LOW = 0;
end

% define filters, padd with zeros due to phase distortions
h1 = [-1  2  6  2 -1  0  0]   / (4*sqrt(2));
g1 = [ 0  0 -2  4 -2  0  0]   / (4*sqrt(2));
h2 = [ 0  0  0  2  4  2  0]     / (4*sqrt(2));
g2 = [ 0 -1 -2  6 -2 -1  0] / (4*sqrt(2));

% cells for selected images
E = cell(3,LEVELS);
I0 = ORIGIMGS;

% loop over decomposition depth -> analysis
rows = zeros(1,LEVELS);
cols = zeros(1,LEVELS);
for i = 1:LEVELS 
  % calculate and store actual image size 
  [r c n]  = size(ORIGIMGS); 
  rows(i) = r; cols(i)  = c;

  % images
  clear Z A1 A2 A3 A4
  for t = 1:n  
      Z(:,:,t)  = dec(conv2(es(ORIGIMGS(:,:,t), 7, 1), g1, 'valid'),2); %16
      A1(:,:,t) = dec(conv2(es(Z(:,:,t), 7, 2), g1','valid'),1); %8
      A2(:,:,t) = dec(conv2(es(Z(:,:,t), 7, 2), h1','valid'),1); %8
      Z(:,:,t)  = dec(conv2(es(ORIGIMGS(:,:,t), 7, 1), h1, 'valid'),2); %16
      A3(:,:,t) = dec(conv2(es(Z(:,:,t), 7, 2), g1','valid'),1); %8
      A4(:,:,t) = dec(conv2(es(Z(:,:,t), 7, 2), h1','valid'),1); %8
  end
  % select coefficients and store them
  E(1,i) = {selc(A1, HIGH)};%4
  E(2,i) = {selc(A2, HIGH)};%4
  E(3,i) = {selc(A3, HIGH)};%4
 
  % copy input image for next decomposition stage
  ORIGIMGS = A4;  

end;

% select base coefficients of last decompostion stage
if LOW < 0
  A4 = selb(A4,LOW,I0);
else
  A4 = selb(A4,LOW);
end

% loop over decomposition depth -> synthesis
for i = LEVELS:-1:1
  % undecimate and interpolate (rows)
  A4 = conv2(es(undec(A4    ,1), 3, 2), h2', 'valid'); %3  
  A3 = conv2(es(undec(E{3,i},1), 3, 2), g2', 'valid'); %3 
  A2 = conv2(es(undec(E{2,i},1), 3, 2), h2', 'valid'); %3 
  A1 = conv2(es(undec(E{1,i},1), 3, 2), g2', 'valid'); %3 

  % undecimate and interpolate (columns)  
  A4 = conv2(es(undec(A4+A3,2), 3, 1), h2, 'valid'); %5
  A2 = conv2(es(undec(A2+A1,2), 3, 1), g2, 'valid'); %5 

  % add images and select valid part 
  A4 = A4 + A2;
  A4 = A4(5:5+rows(i)-1,5:5+cols(i)-1);
end;

% copy image
FUSEDCOEF = A4;