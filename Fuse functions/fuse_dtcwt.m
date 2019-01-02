function FUSEDCOEF = fuse_dtcwt(ORIGIMGS, LEVELS, HIGH, LOW, BIORT, QSHIFT, VARIATE, DENOISE, GAIN, FINETUNE)
% fuse_dtcwt - Dual-Tree Complex Wavelet Fusion and Denoising using statistic model
%
%   FUSEDCOEF = fuse_dtcwt(ORIGIMGS, LEVELS, HIGH, LOW, BIORT, QSHIFT, VARIATE, DENOISE, GAIN, FINETUNE)
%   ORIGIMGS is a 3-D matrix: input images
%   LEVELS is maximum decomposition level   (default = 4)
%   HIGH is coefficient selection highpass according to defined model
%           HIGH == 'LAP' or 1: Laplacian
%           HIGH == 'CAU' or 2: Cauchy
%           HIGH == 'GEN' or 3: Generalised Gaussian (default)
%           HIGH == 'ALP' or 4: Alpha-Stable
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
%   VARIATE is variable relationship
%           VARIATE = 1: Univariate
%           VARIATE = 2: Bivariate (default)
%   DENOISE applies denosing in wavelet tranform domain
%           DENOISE = 0: No denoising (default)
%           DENOISE = 0: Denoising 
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
%               high = 'ALP';   % Alpha-Stable distribution to model highpass coefficients
%               % Using default biort and qshift, 
%               F = fuse_dtcwt(I, levels, high, low);
%               figure; imshow(F)
%               var = 1; %univariate
%               gain = [2 1.4 1.1 1]; % Add sharpening
%               F = fuse_dtcwt(I, levels, high, low, [], [], var, [], gain);
%               figure; imshow(F)

%   v 1.0 19.03.09  Artur Loza (University of Bristol)
%   v 1.1 20.09.12  Nantheera Anantrasirichai (University of Bristol)
%                   Modified to use more options for low-pass coefficient fusion
%                   Modified to use DT-CWT by Nick Kingsbury
%                   Added GAIN and FINETUNE
% -------------------------------------------------------------------------

%%  check inputs and put default values
if nargin < 2 || isempty(LEVELS)
    LEVELS = 4;
end
% Set the HIGH and the corresponding filter
if nargin < 3
    error('Fusion method is required!!');
elseif ~ischar(HIGH)
    switch HIGH
        case 1, HIGH = 'LAP';
        case 2, HIGH = 'CAU';
        case 3, HIGH = 'GEN';
        case 4, HIGH = 'ALP';
    end
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
if nargin < 7 || isempty(VARIATE)
    VARIATE = 2;
end
if nargin < 8 || isempty(DENOISE)
    DENOISE = 0;
end
if nargin < 9 || isempty(GAIN)
    GAIN = ones(6,LEVELS); % no highpass boost
elseif length(GAIN)~=LEVELS
    temp = GAIN;
    GAIN = ones(1,LEVELS);
    GAIN(1:min(length(temp),length(GAIN))) = temp(1:min(length(temp),length(GAIN)));
    GAIN = ones(6,1)*GAIN;
end
if nargin < 10 || isempty(FINETUNE)
    FINETUNE = 0;
end

%%
% parameters
N = 3;
sas = 0;
totalImages = size(ORIGIMGS,3);
convsize = 'same';%'valid';
th = 0.75; % threshold

fmeth = upper(deblank(HIGH(1:min(3,length(HIGH)))));
switch fmeth   %         win b den ma bi
    case 'LAP', params = {N  1  0  0  0}; % univar Laplacian (local)
    case 'BLA', params = {N  1  0  0  1}; % bivar  Laplacian (local)
    case 'ULS', params = {N  1  1  0  0}; % univar Laplacian shrinkage  (local)
    case 'BLS', params = {N  1  1  0  1}; % bivar  Laplacian shrinkage  (local)
    case 'GEN', params = {N nan 0  0  0}; % GGD (with B estimation)
    case 'CAU', params = {N  1  0  1  0}; sas = 1; % univar Cauchy (local?)
    case 'BCA', params = {N  1  0  1  1}; sas = 1; % bivar  Cauchy (global?)
    case 'UCS', params = {N  1  1  1  0}; sas = 1; % univar Cauchy shrinkage (global)
    case 'BCS', params = {N  1  1  1  1}; sas = 1; % bivar  Cauchy shrinkage (global)
    case 'ALP', params = {N nan 0  1  0}; sas = 1; % SAS (with estimation)
end

params{5} = VARIATE-1;
params{3} = DENOISE;

%win b  den         ma     bi
[win,b,denois_fuse,sign_ma,parent] = deal(params{:});

% fix cov for multiple inputs
if totalImages > 2 && ~strcmp(fmeth,'WAN')
    thr = inf;
end % or use mtple inputs only w/WAN

% to bypass fusion set this:  
if totalImages == 1 && denois_fuse % if denoising-only
    thr = inf;
end

% determine the method for b
B{1} = b(1);
B{2} = b(end);
B{3} = b(end);
if isnan(b),
    b_est = 1;
else
    b_est = 0;
end

% averaging window
if win == 2.5,
    winfilt = [0 1 0; 1 1 1; 0 1 0];
elseif win == 1
    winfilt = 1;
else   %   winfilt = GausWin2(1.1,win);%   winfilt = winfilt * winfilt';
    winfilt = ones(win);
end
winfilt = winfilt/sum(winfilt(:));


%% DT-CWT --------------------------------------------------------

% Forward dual-tree DWT, either FSfarras or AntonB function can be used to combute the stage 1 filters
clear Xl Xh
for imi = 1 : totalImages
    [Xl(:,:,imi),Xh{imi}] = dtwavexfm2(ORIGIMGS(:,:,imi),LEVELS,BIORT,QSHIFT); 
    % b/c we need one more level, for bi-variate HIGH
    if parent
        [~,Xh{imi}] = dtwavexfm2(ORIGIMGS(:,:,imi),LEVELS+1,BIORT,QSHIFT); 
    end
    % Noise variance estimation using robust median estimator..
    tmp = real(Xh{imi}{1});% Selesnic uses {1}{1}{1}{1} but this seems to work better
    Nsig(1,1,imi) = median(abs(tmp(:)))/0.6745; % Nsig2 = Nsig^2; % (1:2)
end


%% ANALYSE COEFFS 

Outh = Xh{1}(1:LEVELS);
for scale = 1 : LEVELS,
    clear Y* S* A 
    for dir = 1 : 6 % orientation
        for re_im = 1 : 2 % real\imaginary
            
            % ESTIMATE shape parameter
            % ---------------------------------------------------------
            for imi = 1 : totalImages % # of inputs
                
                % complex coefficients
                Yj00(:,:,imi) = Xh{imi}{scale}(:,:,dir);
                
                if parent
%                     Yjp1(:,:,imi) = expand(Xh{imi}{scale+1}(:,:,dir));
                    [r,c] = size(Yj00(:,:,imi));
                    temp = expand(Xh{imi}{scale+1}(:,:,dir));
                    Yjp1(:,:,imi) = temp(1:r,1:c);
                else
                    Yjp1(:,:,imi)  = zeros(size(Yj00(:,:,imi)));
                end
                
                if re_im == 1 % real part
                    Yj00(:,:,imi) = real(Yj00(:,:,imi));
                    Yjp1(:,:,imi) = real(Yjp1(:,:,imi));
                else % im part
                    Yj00(:,:,imi) = imag(Yj00(:,:,imi));
                    Yjp1(:,:,imi) = imag(Yjp1(:,:,imi));
                end
                
                % ESTIMATE shape parameter
                if sas
                    B{imi} = ShapeParSAS(Yj00(:,:,imi),b_est,B{imi});
                else  % GGD
                    B{imi} = ShapeParGGD(Yj00(:,:,imi),b_est,B{imi});
                end
                
            end % imi
            
            % COMPUTE SALIENCY & MATCH
            % ---------------------------------------------------------
            b = min([B{:}]);
            
            if sas
                b = max(1,b);
                for imi = 1 : totalImages
                    S2(:,:,imi) = SaliencySAS(Yj00(:,:,imi),Yjp1(:,:,imi),winfilt,0,b);
                end
                
                if totalImages > 1
                    MA = MatchMeasureSAS(Yj00,winfilt,b);
                else % if denoising-only
                    MA = 0;
                end
                
            else % GGD
                
                % this function computes everything in "laplacian" way (analogously), fusion and denoising together
                for imi = 1 : totalImages
                    S2(:,:,imi) = SaliencyGGD(Yj00(:,:,imi),Yjp1(:,:,imi),winfilt,0,b,fmeth,convsize);
                end
                
                if totalImages > 1
                    MA = MatchMeasureGGD(Yj00,Yjp1,S2,winfilt,b,fmeth,sign_ma,convsize);
                else % if denoising-only
                    MA = 0;
                end
                
            end
            
            % DENOISE (or not)
            % ---------------------------------------------------------
            if denois_fuse
                
                if sas
                    
                    for imi = 1 : totalImages
                        if re_im==1, % estimate from real only
                            K = 2:4; % 1:5 for global?
                            gamm = CauchynoisyecfMomFit(Yj00(:,:,imi),Yjp1(:,:,imi),Nsig(1,1,imi),winfilt,K,parent);
                        end
                        Yj00(:,:,imi) = Cauchyshrink(Yj00(:,:,imi),Yjp1(:,:,imi),gamm,Nsig(1,1,imi),parent);
                    end
                    
                else % GGD
                    
                    A = DenoiseGGD(Nsig,Yj00,Yjp1,S2); % shrinkage function
                    Yj00 = Yj00 .* A; % apply shrinkage to coeffs
                    S2 = S2 .* A.^2;  % consequently, update saliency with shrinkage
                    MA = MA .* prod(A,3); % update matching with shrinkage
                end
            end
            
            % FUSION
            % ---------------------------------------------------------
            if totalImages > 1
                m_weight = MA > th;
                m12 = S2(:,:,1) > S2(:,:,end);
                w_min = 0.5 * ( 1 - (1-MA)/(1-th) );
                w_max = 1-w_min;
                
                Y_fused  = ...
                    ~ m_weight .* ( m12.*Yj00(:,:,1) + ~m12.*Yj00(:,:,totalImages) ) ... selection
                    + m_weight .* ( m12.*(Yj00(:,:,  1).*w_max + Yj00(:,:,totalImages).*w_min) ...  weighting S1>S2
                    +             ~m12.*(Yj00(:,:,totalImages).*w_max + Yj00(:,:,1   ).*w_min) ); % weighting S2>S1
                
            else % if denoising-only
                Y_fused = Yj00(:,:,1);
            end
            Z_fused{re_im} = Y_fused;
        end 
        
        % pass the fused coefficents on
        Outh{scale}(:,:,dir) = Z_fused{1} + 1i*Z_fused{2};
        
        if FINETUNE
            clear temp
            for imi = 1 : totalImages 
                temp(:,:,imi) = Xh{imi}{scale}(:,:,dir);
            end
            Zh = sum(temp,3)./abs(sum(temp,3));
            Outh{scale}(:,:,dir) = Zh.*abs(Outh{scale}(:,:,dir));
        end
    end 
end 


%% APPROXIMATION FUSION - LOWPASS COEFFICIENTS
Outl = selb(Xl,LOW); 


%% INVERSE TRANSFORM 
FUSEDCOEF = dtwaveifm2(Outl,Outh,BIORT,QSHIFT,GAIN);
% crop if wrong size
if any(size(FUSEDCOEF)~=size(ORIGIMGS(:,:,1)))
    FUSEDCOEF = FUSEDCOEF(1:size(ORIGIMGS,1),1:size(ORIGIMGS,2));
end

% -------------------------------------------------------------------------
%                                FUNCTIONS 
% -------------------------------------------------------------------------

function [S2 B] = SaliencySAS(Yj00,Yjp1,winfilt,b_est,B)

Yj00 = abs(Yj00); % remove abs() in following equations
Yjp1 = abs(Yjp1);

% compute saliency (dispersion of SaS components)
k1 = conv2a(log(Yj00+eps),winfilt,'same');

if b_est % estimate SaS alpha (or 'b') - globally
    k1a = mean2(  log(Yj00+eps) ); % =0, for Y=-eps!
    k2a = mean2( (log(Yj00+eps)-k1a).^2 );
    B = min(1./sqrt(6*k2a/pi^2-0.5),2);
end

S2 = exp(B.*k1+psi(1)*(1-B)).^(1./B);


% -------------------------------------------------------------------------
function B = ShapeParSAS(Yj00,b_est,B)

if ~b_est, return, end

k1a = mean2(  log(abs(Yj00+eps)) );
k2a = mean2( (log(abs(Yj00+eps))-k1a).^2 );
B = min(1./sqrt(6*k2a/pi^2-0.5),2);


% -------------------------------------------------------------------------
function a = GGDcoeff(fmeth,b)

% a shortcut for coefficient resulting from definition or ML of GGD when
% estimating variance

switch fmeth
    case 'LAP', a = 2; % b = 1, univariate
    case {'BLA' 'SHU' 'SHR'}, a = 3/4; % Selesnic bivariate HIGH
    case {'GGD' 'SHM'}, a = gamma(3/b)/gamma(1/b) * b^(2/b); % GGD
    otherwise, a = 1;
end

% -------------------------------------------------------------------------
function [S2 B] = SaliencyGGD(Yj00,Yjp1,winfilt,b_est,B,fmeth,convsize)

% compute saliency (dispersion of GGD components)

if b_est % estimate GLOBALLY
    B = laplacpar3(Yj00);
end

a = GGDcoeff(fmeth,B); % normalising constant

R = sqrt(abs(Yj00).^2 + abs(Yjp1).^2);
S2 = a * conv2a( R.^B, winfilt, convsize ).^(2/B); % variance

% -------------------------------------------------------------------------
function B = ShapeParGGD(Yj00,b_est, B)

if ~b_est, return, end

B = laplacpar3(Yj00);

% -------------------------------------------------------------------------
function MA = MatchMeasureSAS(Yj00,winfilt,B)

x1 = Yj00(:,:,1);
x2 = Yj00(:,:,end);

if length(B) > 1
    p = max(1,min(B{:}));
else
    p = B;
end

% symmetric coefficient of covariation
Nom1 = conv2a(x1.*sign(x2).*abs(x2).^(p-1),winfilt,'same');
Nom2 = conv2a(x2.*sign(x1).*abs(x1).^(p-1),winfilt,'same');

Den1 = conv2a(abs(x2).^p,winfilt,'same');
Den2 = conv2a(abs(x1).^p,winfilt,'same');
MA = Nom1.*Nom2./(Den1.*Den2+eps);

% -------------------------------------------------------------------------
function MA = MatchMeasureGGD(Yj00,Yjp1,S2,winfilt,b,fmeth,sign_ma,convsize)

if strcmp(fmeth,'AVE'),
    MA = 1;
    return,
end

if length(b) > 1,
    b = mean([b{:}]); % sqrt(B{1}.*B{2}); % but this mean value is not used for S2
end

if sign_ma, % get 'possible' sign
    signY = sign(Yj00); % covariance, with sign (e.g. =0 for orthogonal)
else
    signY = 1;
end

R = signY .* sqrt(abs(Yj00).^2 + abs(Yjp1).^2);

% (1) this is generalised cov
MA = conv2a( prod(R,3).^(b/2), winfilt, convsize ).^(2/b);

% (2) this is normal cov w/o or w/sign
% MA = conv2a( prod(R,3), winfilt, convsize );

% (3) this normal covariance (univariate)
% MA = conv2a( prod(Yj00,3), winfilt, convsize );

a = GGDcoeff(fmeth,b); % normalising constant
MA = a * MA; % b/c S2 = a * ORIGIMGS;

MA = (2 * MA + eps) ./ (sum(S2,3) + eps); % abs(MA)?
%         [min(MA(:)) max(MA(:))]

MA = max(min(MA,1),-1); % threshold the matching measure % w = w / max(w(:));

% -------------------------------------------------------------------------
function A = DenoiseGGD(Nsig,Yj00,Yjp1,S2)

% compute denosing parameters
Nsig2 = repmat(Nsig.^2,[size(S2,1) size(S2,2) 1]);

% get the signal variance
s2 = S2-Nsig2;
s2 = s2 .* (s2 > 0); % else A = ones(size(S2));

% shrinkage function
R = sqrt(abs(Yj00).^2 + abs(Yjp1).^2); % =|Yj00| if parent not used
T = sqrt(3)*Nsig2 ./ sqrt(s2+eps); % threshold
A = (R - T)./(R+eps);
A = A .* (A > 0); % this is shrinkage function

% -------------------------------------------------------------------------
function  [p,s]=laplacpar3(x,p)

% Log estimation of the parameters 's' and 'p' of a generalized Laplacian (gL)
% by Artur Loza 23/06/2006

Y = log( abs(x(:)) + eps);
%mu_y = mean( Y );
N = length(Y);

if nargin == 1 || isempty(p) || p == 0
    
    ord = 2;
    
    % Define the function to be inverted
    p = .1 :.001: 2;
    
    var_p = psi( ord-1 , 1./p ) ./ p.^ord;
    
    % make it converge to 0 (experi-mental)
    if N <= 64
        warning('small sample size var = var - 1')
        var_p = var_p - 1;
    end
    
    % compute the inverse of var_x in the input point, thus resulting
    % the second parameter, p.
    var_y = var( Y );
    
    p = interp1(var_p, p, var_y, 'spline');
    
end

if p < 0, warning('p negative, setting to 1'), p = 1; end

if nargout > 1
    mu_y = mean( Y );
    s = exp( mu_y - psi(1./p)/p );
end

% -------------------------------------------------------------------------
% CauchynoisyecfMomFit
function gama=CauchynoisyecfMomFit(data,datap,Nsig,winfilt,K,parent)%,lev_noise
% Moments estimate of Cauchy parameter
% K is a vector containing the two points where the ECF is to be estimated

convsize = 'same';%'valid';
% this part is 'sasecf'
f(1,1,:) = (pi*K(:)/256)';

f0 = repmat(f(1,1,:),[size(data) 1]);
f = f0;
data = data(:,:,ones(1,size(f,3)));

exp_jfd_c = exp(1i * data .* f);

if parent
    datap = datap(:,:,ones(1,size(f,3)));
    exp_jfd_p = exp(1i * datap .* f);
    exp_jfd = (exp_jfd_c + exp_jfd_p)/2;
else
    exp_jfd = exp_jfd_c;
end

% this part is 'CauchynoisyecfMomFit'
FNs = conv2a(exp_jfd,winfilt,convsize); % FNs = mean(exp(j*f*data'),2);
PNs = abs(FNs);
gama = mean((-log(PNs.^2)-Nsig^2*f0.^2)./(2*f0),3);

if all(ismember(1:5,K))
    thr  = eps;
elseif all(ismember(2:4,K))
    thr = 1;
else
    error('K?')
end

gama = max(gama,thr); 

% -------------------------------------------------------------------------
function [w1] = Cauchyshrink(y1,y2,disp,Nsig,parent) % bi
% Bivariate Cauchy Shrinkage Function
% Usage :
%      [w1] = biCauchyshrink(y1,y2,disp,Nsig)
% INPUT :
%      y1 - a noisy coefficient value
%      y2 - the corresponding parent value
%      disp  - signal dispersion
%      Nsig - noise standard deviation
% OUTPUT :
%      w1 - the denoised coefficient
if parent
    p = (disp.^2+3*Nsig^2)./(1+y2.^2./y1.^2) - y1.^2./3;
    q = -2*y1.^3/27 + (y1./3).*(disp.^2+3*Nsig^2)./(1+y2.^2./y1.^2) - disp.^2.*y1./(1+y2.^2./y1.^2);
else
    p =  disp.^2 + 2*Nsig^2  - y1.^2/3;
    q = -2*y1.^3/27 - 2/3*disp.^2.*y1 + 2/3*Nsig^2*y1;
end


DD = p.^3/27 + q.^2/4;
w1 = y1/3 + (abs(-q/2 + sqrt(DD)).^(1./3.)).*sign(-q/2 + sqrt(DD)) + ...
    (abs(-q/2 - sqrt(DD)).^(1./3.)).*sign(-q/2 - sqrt(DD));

% -------------------------------------------------------------------------
function [Y, first, last] = GetCMat(ORIGIMGS,s)

sx(1)  = size(ORIGIMGS,1);
sx(2)  = size(ORIGIMGS,2);

if any(isinf(s))
    dim = find(s == inf);
    s(dim) = sx(dim);
end

s = min([s; sx]);
ds = (sx-s)/2;

first = 1  + floor(ds);
last  = sx - ceil (ds);

Y = ORIGIMGS( first(1):last(1) , first(2):last(2),:);
% -------------------------------------------------------------------------
function y = conv2a(x, winfilt, convsize)

% a shortcut to conv2: conv2 is skipped if filtering window is larger
% than image to be filtered (useful when a global estimate is computed)

if all(size(winfilt) < [size(x,1) size(x,2)]) % convolve
    
    y = convn(x, winfilt, convsize);
    
else % multiply and convolve the central parts
    
    winfilt = GetCMat(winfilt,size(x));
    
    y = sum2(x .* winfilt); % assuming the window is symmetrical W(n)=W(-n)
    
end