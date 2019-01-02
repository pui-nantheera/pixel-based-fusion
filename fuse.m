function [FUSEDIMG,FUSEDIMGcolour,CTIME] = fuse (varargin)
% fuse - Image fusion tools
% 
%   [FUSEDIMG,FUSEDIMGcolour,CTIME] = fuse(ORIGIMGS, METHOD, PARAM);
%   Inputs: 
%   ORIGIMGS is a 3-D matrix or a Cell array of 2-D matrix
%       If ORIGIMGS is a 3-D matrix, the 3rd dimension is the stack of original images (grayscale only)
%       If ORIGIMGS is a Cell array, each Cell is the original images. If
%       one of the input images are colour, ORIGIMGS must be a Cell array
%   METHOD is a fusion method name (string) or a method index (integer) (see below).
%   PARAM is parameter lists for METHOD (see below).
%   Outputs:
%   FUSEDIMG is a fused Image (grayscale).
%   FUSEDIMGcolour is a fused Image (RGB if applicable).
%   CTIME is computational time of fusion process.
%
%   METHOD and PARAM:
%   1: Average (no PARAM)       (fuse_avg.m for details)
%   2: PCA (no PARAM)           (fuse_pca.m for details)
%   3: Maximum (no PARAM)       (selc.m for details)
%   4: Minimum (no PARAM)       (selc.m for details)
%   5: Spatial Frequency             (fuse_spafrq.m for details)
%      PARAM:   'blockSize' is a block size for spatial correlation (default = 8)
%               'threshold' for the blocks to include in the fusion (default = 1)
%   6: Laplacian Pyramid        (fuse_lap.m for details)
%      PARAM:   'levels', 'low', 'high' (see below)
%   7: FSD Pyramid  (Filter-Substract-Decimate Laplacian Pyramid)  (fuse_fsd.m for details)
%      PARAM:   'levels', 'low', 'high' (see below)
%   8: Ratio Pyramid            (fuse_rat.m for details)
%      PARAM:   'levels', 'low', 'high' (see below)
%   9: Contrast Pyramid         (fuse_con.m for details)
%      PARAM:   'levels', 'low', 'high' (see below)
%   10: Gradient Pyramid        (fuse_gra.m for details)
%      PARAM:   'levels', 'low', 'high' (see below)
%   11: Morphological Pyramid         (fuse_mor.m for details)
%      PARAM:   'levels', 'low', 'high' (see below)
%   12: DWT (Discrete Wavelet Transform)    (fuse_dwb.m for details)
%      PARAM:   'levels', 'low', 'high' (see below)
%   13: LISQ (Lifting Scheme on Quincunx Grids) (fuse_lisq.m for details)
%      PARAM:   'levels', 'low', 'high' (see below)
%   14: SIDWT (Shift-Invariant Discrete Wavelet Transform)  (fuse_sih.m for details)
%      PARAM:   'levels', 'low', 'high' (see below)
%   15: CWT (Complex Wavelet Transform)     (fuse_cwt.m for details)
%      PARAM:   'levels', 'low', 'high' (see below)
%               'biort' is a filter for decomposing first level 
%                   'antonini'   => Antonini 9,7 tap filters. (default)
%                   'legall'     => LeGall 5,3 tap filters.
%                   'near_sym_a' => Near-Symmetric 5,7 tap filters.
%                   'near_sym_b' => Near-Symmetric 13,19 tap filters.
%               'qshift' is a filter for decomposing higher level
%                   'qshift_06' => Quarter Sample Shift Orthogonal (Q-Shift)
%                                   10,10 tap filters,
%                                   (only 6,6 non-zero taps). (default)
%                   'qshift_a' =>  Q-shift 10,10 tap filters,
%                                   (with 10,10 non-zero taps, unlike qshift_06).
%                   'qshift_b' => Q-Shift 14,14 tap filters.
%                   'qshift_c' => Q-Shift 16,16 tap filters.
%                   'qshift_d' => Q-Shift 18,18 tap filters.
%   16: Model              (fuse_model.m for details)
%      PARAM:   'method' is the model-based method for fusing high-pass coeffients
%                   METHOD == 1: Weigted average scheme
%                   METHOD == 2: FLOM algorithm
%                   METHOD == 3: FLOM-Cauchy algorithm 
%                   METHOD == 4: Meridian distribution
%                   METHOD == 5: Cauchy probability (very slow. Recommend:
%                                use 'Stats' with 'method' = 2)
%   17: Stats (Fusion and Denoising using statistic model)  (fuse_dtcwt.m for details)
%      PARAM:   'method' is the model-based method for fusing high-pass coeffients
%                   METHOD == 1: Laplacian
%                   METHOD == 2: Cauchy probability
%                   METHOD == 3: Generalised Gaussian (default)
%                   METHOD == 4: Alpha-Stable
%
%   PARAM details:
%   'level' is the number of decomposition levels (default = 4)
%   'low' is a method for fusing low-pass coefficients (default = 0)
%   'high' is a method for fusing hig-pass coefficients (selc.m for details)
%       It is an column or row array. The first element is selection type.
%       HIGH(1) == 1: choose max(abs) (default)
%       HIGH(1) == 2: salience / match measure with threshold == .75 (as proposed by Burt et al)      
%       HIGH(1) == 3: choose max with consistency check (as proposed by Li et al)     
%       HIGH(1) == 4: simple choose max
%       The second element is for HIGH(1) == 2 and 3. It is the window size (default = 3)
%
%   Examples:   I{1} = im2double(imread('image1.png'));
%               I{2} = im2double(imread('image2.png'));
%
%               fusedImgContrast = fuse(I, 'CONTRAST PYRAMID','levels',3,'high',2,'low',0);
%               figure; imshow(fusedImgContrast); title('Fused image using Contrast pyramid');
%
%               [fusedImgDTCWT,~,ctime] = fuse(I, 'CWT','levels',4,,'high',1,'low',0,'biort','near_sym_a','qshift','qshift_d');
%               figure; imshow(fusedImgDTCWT); title(['Fused image using DT-CWT, processing time ',num2str(ctime),' sec']);

% v 1.0 2004 Eduardo Fernandez Canga - University of Bristol
% v 1.1 2012 Pui Added comments, improved reading inputs


%% check input arg
% =========================================================================
if nargin < 2
    error('Input images and fusion method are needed!');
end
if rem(nargin,2)
    error('Wrong input arguments!')
end

%% check input images
% =========================================================================
inputImgs = varargin{1};
% number of original images
totalImgs = length(inputImgs);
% convert to input mat
origImgsC = []; % for store colour channels
if iscell(varargin{1})
    ORIGIMGS = zeros(size(inputImgs{1},1), size(inputImgs{1},2), totalImgs);
    
    for k = 1:totalImgs
        curImg = inputImgs{k};
        if size(curImg,3)==3
            % other colour conversion may be used, e.g. rgb2ycbcr
            % if colours are ignore, rgb2gray can be used. Also change the
            % function at the end of the process for converting back to rgb
            % We also provide rgb2grayn.m in this toolbox
            curImg = rgb2hsv(curImg);
            % shift intensity to the first component.
            % If rgb2ycbcr is used, you don't need this line.
            curImg = curImg(:,:,[3 1 2]);
            origImgsC{k} =  curImg(:,:,2:3);
        end
        ORIGIMGS(:,:,k) = curImg(:,:,1);
    end
else
    ORIGIMGS = inputImgs;
end
clear inputImgs
%% check method
% =========================================================================
method = varargin{2};
switch upper(method)
    case {1,'AVERAGE'}
        v = 1;
        fusefcn  = @fuse_avg;
        titlestr = 'FUSE AVR ';
    case {2,'PCA'}
        v = 2;
        fusefcn  = @fuse_pca;
        titlestr = 'FUSE PCA ';
    case {3,'MAXIMUM'}
        v = 3;
        fusefcn  = @selc;
        fuseinput{2} = 4;
        titlestr = 'FUSE MAX ';
    case {4,'MINIMUM'}
        v = 4;
        fusefcn  = @selc;
        fuseinput{2} = 5;
        titlestr = 'FUSE MIN ';
    case {5,upper('Spatial Frequency'), upper('SpatFreq')}
        v = 5;
        fusefcn  = @fuse_spafrq;
        titlestr=['FUSE MFO'];
    case {6,'LAPLACIAN PYRAMID','LAPLACIAN'}
        v = 6;
        fusefcn  = @fuse_lap;
        titlestr = 'FUSE LAP';
    case {7,'FSD PYRAMID','FSD'}
        v = 7;
        fusefcn  = @fuse_fsd;
        titlestr = 'FUSE FSD';
    case {8,'RATIO PYRAMID','RATIO'}
        v = 8;
        fusefcn  = @fuse_rat;
        titlestr = 'FUSE RAT';
    case {9,'CONTRAST PYRAMID','CONTRAST'}
        v = 9;
        fusefcn  = @fuse_con;
        titlestr = 'FUSE CON';
    case {10,'GRADIENT PYRAMID','GRADIENT'}
        v = 10;
        fusefcn  = @fuse_gra;
        titlestr = 'FUSE GRA';
    case {11,upper('Morphological Pyramid'),upper('Morphological')}
        v = 11;
        fusefcn  = @fuse_mod;
        titlestr = 'FUSE MOD';
    case {12,'DWT',upper('Discrete Wavelet Transform')}
        v = 12;
        fusefcn  = @fuse_dwb;
        titlestr = 'FUSE DWT';
    case {13,'LISQ'}
        v = 13;
        fusefcn  = @fuse_LISQ;
        fuseinput{5} = 4;
        titlestr = 'FUSE LIS';
    case {14,'SIDWT'}
        v = 14;
        fusefcn  = @fuse_sih;
        titlestr = 'FUSE SIH';
    case {15,'CWT',upper('Complex Wavelet')}
        v = 15;
        fusefcn  = @fuse_cwt;
        biort = 'antonini';
        qshift = 'qshift_06';
        titlestr = 'FUSE CWT';
    case {16,upper('Model')}
        v = 16;
        fusefcn  = @fuse_model;
        biort = 'antonini';
        qshift = 'qshift_06';
        titlestr = 'FUSE Model based CWT';
    case {17,upper('Stats')}
        v = 17;
        fusefcn  = @fuse_dtcwt;
        biort = 'antonini';
        qshift = 'qshift_06';
        titlestr = 'FUSE Model based DTCWT';
    otherwise
        error('Wrong method!!');
end % end of switch



%% check parameters
% =========================================================================

% paste default values
high = [1 3];   % select maximum magnitude of high-pass coefficients
low  = 0;       % average low-pass coefficients
levels = 4;     % decomposition levels
blockSize = 8;  % blocksize and threshold for Spatial Frequency
threshold = 1;
% read parameters if define
if nargin > 2
    for k = 3:2:nargin
        if isnumeric(varargin{k+1})
            eval([varargin{k},'= [',num2str(varargin{k+1}),'];']);
        elseif (v==14)
            eval([varargin{k},'=''',varargin{k+1},'''']);
        else
            warning(['Wrong input value for "',varargin{k},'" - default is used']);
        end
    end
end
% repaste if define
if (v==5)
    % parameters for mfocus
    fuseinput{2} = blockSize;
    fuseinput{3} = threshold;
elseif (v > 5)
    % parameters for transform-domain fusion methods
    fuseinput{2} = levels;
    fuseinput{3} = high;
    fuseinput{4} = low;
    if (v>=15) % for DT-CWT
        fuseinput{5} = biort;
        fuseinput{6} = qshift;
    end
    if (v>=16) % for model-based
        fuseinput{3} = method;
    end
end

% add some more info for display at the end for some methods
if (v==5)
    titlestr = [titlestr,' block size=',num2str(blockSize),'; threshold=', num2str(threshold),';'];
elseif (v > 5)
    titlestr = [titlestr,' levels=',num2str(levels),'; high=[',num2str(high),']; low=',num2str(low),';'];
end
if (v>=15) % for DT-CWT
    titlestr = [titlestr,'  biort=''' fuseinput{5} '''; ' 'qshift=''' fuseinput{6} ''';'];
elseif (v==13)
    titlestr = [titlestr,' lisq_coef=' num2str(fuseinput{5}),';'];
end

%% FUSION PROCESS
% =========================================================================
fuseinput{1} = ORIGIMGS;
tic % check processing time
FUSEDIMG = fusefcn(fuseinput{:});
CTIME = toc;

%% OUTPUT
% =========================================================================

% combine colour channel
FUSEDIMGcolour = FUSEDIMG;
if ~isempty(origImgsC)
    % average is used here
    sumColour = zeros(size(ORIGIMGS,1),size(ORIGIMGS,2),2);
    count = 0;
    for k = length(origImgsC)
        if sum(origImgsC{k}(:)) > 0
            count = count + 1;
            sumColour = sumColour + origImgsC{k};
        end
    end
    sumColour = sumColour/count;
    if sum(sumColour(:))>0
        FUSEDIMGcolour = cat(3,sumColour, FUSEDIMG);
        FUSEDIMGcolour = hsv2rgb(FUSEDIMGcolour);
    end
end
disp(titlestr);
disp(['Running time = ',num2str(CTIME),' seconds']);