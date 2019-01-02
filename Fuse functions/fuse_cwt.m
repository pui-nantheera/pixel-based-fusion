function FUSEDCOEF = fuse_cwt(ORIGIMGS, LEVELS, HIGH, LOW, BIORT, QSHIFT, GAIN, FINETUNE)
%   fuse_cwt - Image fusion with Dual-tree complex wavelet transfrom
%
%   FUSEDIMG = fuse_cwt(ORIGIMGS, LEVELS, HIGH, LOW, BIORT, QSHIFT, GAIN, FINETUNE)
%   ORIGIMGS is a 3-D matrix: input images
%   LEVELS is maximum decomposition level   (default = 4)
%   HIGH is coefficient selection highpass  (default = 1) (see selcComplex.m) 
%   LOW is coefficient selection base image (default = 0) (see selb.m) 
%   BIORT  - 'antonini'   => Antonini 9,7 tap filters. (default)
%            'legall'     => LeGall 5,3 tap filters.
%            'near_sym_a' => Near-Symmetric 5,7 tap filters.
%            'near_sym_b' => Near-Symmetric 13,19 tap filters.
%   QSHIFT -   'qshift_06' => Quarter Sample Shift Orthogonal (Q-Shift)
%                              10,10 tap filters,
%                              (only 6,6 non-zero taps). (default)
%              'qshift_a' =>  Q-shift 10,10 tap filters,
%                              (with 10,10 non-zero taps, unlike qshift_06).
%              'qshift_b' => Q-Shift 14,14 tap filters.
%              'qshift_c' => Q-Shift 16,16 tap filters.
%              'qshift_d' => Q-Shift 18,18 tap filters.
%   GAIN  is for boosting high-pass coefficients (sharpening: gain > 1)
%              GAIN = [(gain for 1st level) (gain for 2nd level) ...] 
%   FINETUNE is aplying phase of DT-CWT to finely adjust misalignment 
%              between input images (default = 0 not apply)
%
%   FUSEDIMG is a 2-D matrix. It is a grayscale fused image.
%
%   Examples:   I(:,:,1) = im2double(imread('UNcamp_2N_i_20.bmp'));
%               I(:,:,2) = im2double(imread('UNcamp_2N_v_20.bmp'));
%               levels = 3; % decomposition levels
%               high = 1;   % select max absolute coefficients
%               low = 0;    % average base level
%               % Using default biort and qshift, 
%               F = fuse_cwt(I, levels, high, low);
%               figure; imshow(F)
%               F = fuse_cwt(I, levels, high, low, 'near_sym_a', 'qshift_a');
%               figure; imshow(F)
%               % Using default biort and qshift, and apply sharpening with
%               % gain for first level = 2 and second level = 1.6
%               F = fuse_cwt(I, levels, high, low, [], [], [2 1.6]);

%   v 1.0 16.08.99  Oliver Rockinger
%   v 2.0 12.11.03  Eduardo Fernandez Canga (University of Bristol)
%                   Modified for N input images
%   v 2.1 20.09.12  Nantheera Anantrasirichai (University of Bristol)
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
if (nargin < 5) || isempty(BIORT)
    BIORT = 'antonini';
end
if (nargin < 6) || isempty(QSHIFT)
    QSHIFT = 'qshift_06';
end
if nargin < 7 || isempty(GAIN)
    GAIN = ones(1,LEVELS); % no highpass boost
elseif length(GAIN)~=LEVELS
    temp = GAIN;
    GAIN = ones(1,LEVELS);
    GAIN(1:min(length(temp),length(GAIN))) = temp(1:min(length(temp),length(GAIN)));
end
GAIN = ones(6,1)*GAIN;

if nargin < 8 || isempty(FINETUNE)
    FINETUNE = 0;
end

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
n = size(ORIGIMGS,3);

% DT-CWT transformation
for i=1:n
    [Yl(:,:,i),Yh(:,:,:,i)] = dtwavexfm2(ORIGIMGS(:,:,i),LEVELS,BIORT,QSHIFT); 
end

% fusing hight-pass coefficients
Outh = cell(LEVELS, 1); 
for i1 = 1:LEVELS 

    % select coefficients and store them 
    temp=permute(cell2mat(Yh(i1,:,:,:)),[1 2 4 3]);
    % for each band
    for i = 1:6
        Outh{i1}(:,:,i) = selcComplex(temp(:,:,:,i), HIGH); 
        
        if FINETUNE
            Zh = sum(temp(:,:,:,i),3)./abs(sum(temp(:,:,:,i),3));
            Outh{i1}(:,:,i) = Zh.*abs(Outh{i1}(:,:,i));
        end
    end;
    
end;

% fusing low-pass coefficients
Outl = selb(Yl,LOW); 

% reverse DT-CWT
FUSEDCOEF = dtwaveifm2(Outl,Outh,BIORT,QSHIFT, GAIN); 
% crop if wrong size
if any(size(FUSEDCOEF)~=size(ORIGIMGS(:,:,1)))
    FUSEDCOEF = FUSEDCOEF(1:size(ORIGIMGS,1),1:size(ORIGIMGS,2));
end