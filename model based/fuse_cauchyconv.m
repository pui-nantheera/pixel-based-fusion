function FUSEDCOEF = fuse_cauchyconv(ORIGIMGS, LEVELS, LOW, BIORT, QSHIFT, GAIN, FINETUNE)
%   fuse_cauchyconv - Image fusion with Dual-tree complex wavelet transfrom
%   and Cauchy probability
%
%   FUSEDCOEF = fuse_cauchyconv(ORIGIMGS, LEVELS, WINDOWSIZE, LOW)
%   ORIGIMGS is a 3-D matrix: input images
%   LEVELS is maximum decomposition level   (default = 4)
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
%               low = 0;    % average base level
%               % Using default biort and qshift, 
%               F = fuse_cauchyconv(I, levels, high, low);
%               figure; imshow(F)

%   v 1.0 09.07.07  Tao wan(University of Bristol)
%   v 2.0 02.04.10  Mayank Agrawal
%                   Modified for use for remote sensing and surveillence images
%   v 2.1 20.09.12  Nantheera Anantrasirichai (University of Bristol)
%                   Modified to use DT-CWT by Nick Kingsbury
% -------------------------------------------------------------------------
%%  check inputs and put default values

% Number of Scale
if nargin < 2
    LEVELS = 4;
end
% 
if nargin < 3
    LOW = 0;
end
if (nargin < 4) || isempty(BIORT)
    BIORT = 'antonini';
end
if (nargin < 5) || isempty(QSHIFT)
    QSHIFT = 'qshift_06';
end
if nargin < 6
    GAIN = ones(6,LEVELS); % no highpass boost
elseif length(GAIN)~=LEVELS
    temp = GAIN;
    GAIN = ones(1,LEVELS);
    GAIN(1:min(length(temp),length(GAIN))) = temp(1:min(length(temp),length(GAIN)));
    GAIN = ones(6,1)*GAIN;
end
if nargin < 7
    FINETUNE = 0;
end

% check image inputs
if size(ORIGIMGS,3)>2
    error('Sorry this fusion method supports only 2 images!!');
end

%%
% place input
imageA_double=ORIGIMGS(:,:,1);
imageB_double= ORIGIMGS(:,:,2);


% symmetric extension
L = length(imageA_double); % length of original images.
N = L+2^(LEVELS);     % length after extension.

% DT-CWT transform
[Yl(:,:,1),W1] = dtwavexfm2(imageA_double,LEVELS,BIORT,QSHIFT); 
[Yl(:,:,2),W2] = dtwavexfm2(imageB_double,LEVELS,BIORT,QSHIFT); 

for scale = 1:LEVELS
    for band = 1:6
        
        % find weight and combine coefficients
        COEFS = cat(3, W1{scale}(:,:,band), W2{scale}(:,:,band));
        Outh{scale}(:,:,band) = selcModel(COEFS, 1);
        
        if FINETUNE
            sumZh = W1{scale}(:,:,band) + W2{scale}(:,:,band);
            Zh = sumZh./abs(sumZh);
            Outh{scale}(:,:,band) = Zh.*abs(Outh{scale}(:,:,band));
        end
    end
end


% Approximation fusion
Outl = selb(Yl,LOW); 

% reverse DT-CWT
FUSEDCOEF = dtwaveifm2(Outl,Outh,BIORT,QSHIFT, GAIN); 