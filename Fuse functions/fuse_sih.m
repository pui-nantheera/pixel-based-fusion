function FUSEDIMG = fuse_sih(ORIGIMGS, LEVELS, HIGH, LOW)
%   fuse_sih - Image fusion with Shift-Invariant Discrete Wavelet Transform
%
%   FUSEDIMG = fuse_sih(ORIGIMGS, LEVELS, HIGH, LOW)
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
%               F = fuse_sih(I, levels, high, low);
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

% cells for selected images
E = cell(3,LEVELS);

% loop over decomposition depth -> analysis
rows = zeros(1,LEVELS);
cols = zeros(1,LEVELS);
for i = 1:LEVELS 
  % calculate and store actual image size 
  [r c n]  = size(ORIGIMGS); 
  rows(i) = r; cols(i)  = c;
  
  % define actual filters (inserting zeros between coefficients)
  h1 = [zeros(1,floor(2^(i-2))), 0.5, zeros(1,floor(2^(i-1)-1)), 0.5, zeros(1,max([floor(2^(i-2)),1]))];
  g1 = [zeros(1,floor(2^(i-2))), 0.5, zeros(1,floor(2^(i-1)-1)), -0.5, zeros(1,max([floor(2^(i-2)),1]))];
  fh = floor(length(h1)/2);

  clear Z A1 A2 A3 A4
  for t = 1:n
      Z(:,:,t)  = conv2(es(ORIGIMGS(:,:,t), fh, 1), g1, 'valid');
      A1(:,:,t) = conv2(es(Z(:,:,t), fh, 2), g1','valid');
      A2(:,:,t) = conv2(es(Z(:,:,t), fh, 2), h1','valid');
      Z(:,:,t)  = conv2(es(ORIGIMGS(:,:,t), fh, 1), h1, 'valid');
      A3(:,:,t) = conv2(es(Z(:,:,t), fh, 2), g1','valid');
      A4(:,:,t) = conv2(es(Z(:,:,t), fh, 2), h1','valid');
  end
  % select coefficients and store them
  E(1,i) = {selc(A1, HIGH)};
  E(2,i) = {selc(A2, HIGH)};
  E(3,i) = {selc(A3, HIGH)};
 
 	% copy input image for next decomposition stage
  ORIGIMGS = A4;  
end;

% select base coefficients of last decompostion stage
A4 = selb(ORIGIMGS,LOW);

% loop over decomposition depth -> synthesis
for i = LEVELS:-1:1
  % define actual filters (inserting zeros between coefficients)
  h2 = fliplr([zeros(1,floor(2^(i-2))), 0.5, zeros(1,floor(2^(i-1)-1)), 0.5, zeros(1,max([floor(2^(i-2)),1]))]);
  g2 = fliplr([zeros(1,floor(2^(i-2))), 0.5, zeros(1,floor(2^(i-1)-1)), -0.5, zeros(1,max([floor(2^(i-2)),1]))]);
  fh = floor(length(h2)/2);
  
  % filter (rows)
  A4 = conv2(es(A4, fh, 2), h2', 'valid');   
  A3 = conv2(es(E{3,i}, fh, 2), g2', 'valid'); 
  A2 = conv2(es(E{2,i}, fh, 2), h2', 'valid'); 
  A1 = conv2(es(E{1,i}, fh, 2), g2', 'valid'); 

  % filter (columns)  
  A4 = conv2(es(A4+A3, fh, 1), h2, 'valid');  
  A2 = conv2(es(A2+A1, fh, 1), g2, 'valid');  

  % add images 
  A4 = A4 + A2;
end;

% copy image
FUSEDIMG = A4;