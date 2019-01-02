function FUSEDIMG = fuse_gra(ORIGIMGS, LEVELS, HIGH, LOW)
%   fuse_gra - Image fusion with Gradient Pyramid 
%
%   FUSEDIMG = fuse_gra(ORIGIMGS, LEVELS, HIGH, LOW)
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
%               F = fuse_gra(I, levels, high, low);
%               figure; imshow(F)

%   v 1.0 16.08.99  Oliver Rockinger
%   v 2.0 12.11.03  Eduardo Fernandez Canga (University of Bristol)
%                   Modified for N input images
%   v 2.01 25.09.12 Pui (University of Bristol)
%                   Debugged  dimension mismatch
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
% define filters 
w  = [1 4 6 4 1] / 16;
v  = [1 2 1] / 4;
d1 = [1 -1];
d2 = [0 -1; 1 0] / sqrt(2);
d3 = [-1 1];
d4 = [-1 0; 0 1] / sqrt(2);
% compute derivatives
d1e = conv2(d1,d1); 
d1e = [zeros(1,3); d1e; zeros(1,3)];
d2e = conv2(d2,d2);
d3e = d1e';
d4e = conv2(d4,d4);

% cells for selected images
E = cell(1,LEVELS);

% loop over decomposition depth -> analysis
rows = zeros(1,LEVELS);
cols = zeros(1,LEVELS);
for i = 1:LEVELS 
  % calculate and store actual image size 
  [r c n] = size(ORIGIMGS);
  rows(i) = r; cols(i)  = c;
  
  % check if image expansion necessary 
  if (floor(r/2) ~= r/2), ew(1) = 1; else ew(1) = 0; end;
  if (floor(c/2) ~= c/2), ew(2) = 1; else ew(2) = 0; end;

  % perform expansion if necessary
  if (any(ew))
  	ORIGIMGS = adb(ORIGIMGS,ew);
    [r c n]  = size(ORIGIMGS); 
  end	
 
  % perform filtering 
  G = zeros(r,c,n);
  clear Z
  for t = 1 : n
      G(:,:,t) = conv2(conv2(es2(ORIGIMGS(:,:,t),2), w, 'valid'),w', 'valid');
      Z(:,:,t) = es2(ORIGIMGS(:,:,t)+conv2(conv2(es2(ORIGIMGS(:,:,t), 1), v, 'valid'), v', 'valid'), 1);
  end
%   G = convn(convn(es2(ORIGIMGS,2), w, 'valid'),w', 'valid'); % slower
%   Z = es2(ORIGIMGS+convn(convn(es2(ORIGIMGS, 1), v, 'valid'), v', 'valid'), 1); slower

  
  % compute directional derivatives
  B  = zeros(r,c);
  D = zeros(r,c,n);
  for t = 1 : n
      D(:,:,t) = conv2(Z(:,:,t), d1e, 'valid');
  end
  
  B  = B + selc(D, HIGH);
  for t = 1 : n
      D(:,:,t) = conv2(Z(:,:,t), d2e, 'valid');
  end
  
  B  = B + selc(D, HIGH);
  for t = 1 : n
      D(:,:,t) = conv2(Z(:,:,t), d3e, 'valid');
  end

  B  = B + selc(D, HIGH);

  for t = 1 : n
      D(:,:,t) = conv2(Z(:,:,t), d4e, 'valid');
  end
  B  = B + selc(D, HIGH);
  
  % store coefficients
  E(i) = {-B/8};

  % decimate
  ORIGIMGS = dec2(G);
end;

% select base coefficients of last decompostion stage
ORIGIMGS = selb(ORIGIMGS,LOW);

% loop over decomposition depth -> synthesis
for i = LEVELS:-1:1
  % undecimate and interpolate 
  MT = conv2(conv2(es2(undec2(ORIGIMGS), 2), 2*w, 'valid'), 2*w', 'valid');
  % add coefficients
  ORIGIMGS = MT + E{i};
  % select valid image region 
  ORIGIMGS = ORIGIMGS(1:rows(i),1:cols(i));
end;
% copy image
FUSEDIMG = ORIGIMGS;
