function FUSEDIMG = fuse_mod(ORIGIMGS, LEVELS, HIGH, LOW)
%   fuse_mod - Image fusion with morphological difference pyramid
%
%   FUSEDIMG = fuse_mod(ORIGIMGS, LEVELS, HIGH, LOW)
%   ORIGIMGS is a 3-D matrix: input images
%   LEVELS is maximum decomposition level  (default = 4)
%   HIGH is coefficient selection highpass (default = 1) (see selc.m) 
%   LOW is coefficient selection base image (default = 0)(see selb.m) 
%
%   FUSEDIMG is a 2-D matrix. It is a grayscale fused image.
%
%   Examples:   I(:,:,1) = im2double(imread('UNcamp_2N_i_20.bmp'));
%               I(:,:,2) = im2double(imread('UNcamp_2N_v_20.bmp'));
%               levels = 3; % decomposition levels
%               high = 1;   % select max absolute coefficients
%               low = 0;    % average base level
%               F = fuse_mod(I, levels, high, low);
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

  % check and store new size
  [r c n] = size(ORIGIMGS);

  clear O
  % gray scale opening
  for im=1:n
	O(:,:,im) = ordfilt2(ordfilt2(es2(ORIGIMGS(:,:,im), 3), 1, ones(5)), 25, ones(5));
  end
  
  % gray scale closing 
  for im=1:n
	O(:,:,im) = ordfilt2(ordfilt2(O(:,:,im), 25, ones(5)), 1, ones(5));
  end
	
  % select valid image region 
  O2 = O(4:r+3,4:c+3,:);
  
  % decimate
  Z = dec2(O2);

  
  % decimate, undecimate and dilate
  for im=1:n
      O(:,:,im) = ordfilt2(es2(undec2(dec2(O2(:,:,im))), 3), 49, ones(7));
  end

  % select valid image region 
  O = O(4:r+3,4:c+3,:);

  % select coefficients and store them
  E(i) = {selc(ORIGIMGS-O, HIGH)};
  
  % copy tmp images
  ORIGIMGS = Z;
end;  

% select base coefficients of last decompostion stage
ORIGIMGS = selb(ORIGIMGS,LOW);

% loop over decomposition depth -> synthesis
for i = LEVELS:-1:1
  % dilate 
  ORIGIMGS = ordfilt2(es2(undec2(ORIGIMGS), 3), 49, ones(7));
  % select valid image region 
  ORIGIMGS = ORIGIMGS(4:ceil(rows(i)/2)*2+3,4:ceil(cols(i)/2)*2+3);
  % add coefficients
  ORIGIMGS = ORIGIMGS + E{i};
  % select valid image region 
  ORIGIMGS = ORIGIMGS(1:rows(i),1:cols(i));
end;

FUSEDIMG = ORIGIMGS;
