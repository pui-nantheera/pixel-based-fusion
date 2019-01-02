function Y = selcModel(COEFS, METHOD, WINDOW, TH)
% selcModel - coefficinet selection for highpass components 
% 
%   Y = selcModel(COEFS, METHOD)
%   COEFS is a 3-D matrix
%   METHOD switch for selection type
%          METHOD == 1: Weigted average scheme
%          METHOD == 2: FLOM algorithm
%          METHOD == 3: FLOM-Cauchy algorithm (default)
%          METHOD == 4: Meridian distribution
%          METHOD == 5: Cauchy probability
%   WINDOW is WINDOW size (default = 7);
%          METHOD(3): 
%
%   Y is a 2-D matrix. It is a combined coefficients.

%   v 1.0 Dr. Alin Achim and Wan Tao
%   v 2.0 02.04.10  Mayank Agrawal
%                   Modified for use for remote sensing and surveillence images
%   v 2.1 20.09.12  Nantheera Anantrasirichai (University of Bristol)
%                   Comments and put into cases
% -------------------------------------------------------------------------
% check image inputs
if size(COEFS,3)>2
    error('Sorry this fusion method supports only 2 images!!');
end
% if input images are cell, convert to mat
if iscell(COEFS)
    norig = length(COEFS);
    [height, width] = size(COEFS{1});
    temp = COEFS;
    clear COEFS
    COEFS = cell2mat(temp);
    COEFS = reshape(COEFS,height,width,norig);
    clear temp
end
% check inputs
if nargin < 2
    METHOD = 3;
end
if nargin < 3
    WINDOW = 7;
end
if nargin < 4
    TH = 0.75;
end

% Complex coefficients for image 1
%Real part
Y1_coef_real = real(COEFS(:,:,1));
% imaginary part
Y1_coef_imag = imag(COEFS(:,:,1));
% The corresponding coefficients from image 2
%Real part
Y2_coef_real = real(COEFS(:,:,2));
% imaginary part
Y2_coef_imag = imag(COEFS(:,:,2));

% switch to method to find weight
% -------------------------------------------------------------------------
switch(METHOD(1))
        
    case 1
        % Weigted average scheme
        % disp(' Weigted average...');
        
        % Real part:
        % compute salience
        S1 = conv2(symextend(Y1_coef_real.^2, floor(WINDOW/2)), ones(WINDOW), 'valid');
        S2 = conv2(symextend(Y2_coef_real.^2, floor(WINDOW/2)), ones(WINDOW), 'valid');
        % compute match
        MA = conv2(symextend(Y1_coef_real.*Y2_coef_real, floor(WINDOW/2)), ones(WINDOW), 'valid');
        MA = 2 * MA ./ (S1 + S2 + eps);
        % selection
        m1 = MA > TH; m2 = S1 > S2;
        w1 = (0.5 - 0.5*(1-MA) / (1-TH));
        
        Y_coef_real  = (~m1) .* ((m2.*Y1_coef_real) + ((~m2).*Y2_coef_real));
        Y_coef_real  = Y_coef_real + (m1 .* ((m2.*Y1_coef_real.*(1-w1))+((m2).*Y2_coef_real.*w1) + ...
            ((~m2).*Y2_coef_real.*(1-w1))+((~m2).*Y1_coef_real.*w1)));
        
        % Imaginary part:
        % compute salience
        S1 = conv2(symextend(Y1_coef_imag.^2, floor(WINDOW/2)), ones(WINDOW), 'valid');
        S2 = conv2(symextend(Y2_coef_imag.^2, floor(WINDOW/2)), ones(WINDOW), 'valid');
        % compute match
        MA = conv2(symextend(Y1_coef_imag.*Y2_coef_imag, floor(WINDOW/2)), ones(WINDOW), 'valid');
        MA = 2 * MA ./ (S1 + S2 + eps);
        % selection
        m1 = MA > TH; m2 = S1 > S2;
        w1 = (0.5 - 0.5*(1-MA) / (1-TH));
        Y_coef_imag  = (~m1) .* ((m2.*Y1_coef_imag) + ((~m2).*Y2_coef_imag));
        Y_coef_imag  = Y_coef_imag + (m1 .* ((m2.*Y1_coef_imag.*(1-w1))+((m2).*Y2_coef_imag.*w1) + ...
            ((~m2).*Y2_coef_imag.*(1-w1))+((~m2).*Y1_coef_imag.*w1)));
        
    case 2
        % FLOM Algorithm
        % disp(' FLOM Algorithm...');
        % Achim A et al. 'Complex wavelet domain image fusion based on
        % fractional lower order moment', IEEE conf. information fusion,
        % 2005
        
        % Real part:
        % compute salience (dispersion of SaS components)
        [alpha1_real,S1] = sas_salience(WINDOW,Y1_coef_real);
        [alpha2_real,S2] = sas_salience(WINDOW,Y2_coef_real);
        % compute match (symmetric coefficient of covariation)
        p=1;
        MA = sas_match(WINDOW,Y1_coef_real,Y2_coef_real,p);
        % selection
        m1 = MA > TH; m2 = S1 > S2;
        w1 = (0.5 - 0.5*(1-MA) / (1-TH));
        Y_coef_real  = (~m1) .* ((m2.*Y1_coef_real) + ((~m2).*Y2_coef_real));
        Y_coef_real  = Y_coef_real + (m1 .* ((m2.*Y1_coef_real.*(1-w1))+((m2).*Y2_coef_real.*w1) + ...
            ((~m2).*Y2_coef_real.*(1-w1))+((~m2).*Y1_coef_real.*w1)));
        
        % Imaginary part:
        % compute salience
        [alpha1_imag,S1] = sas_salience(WINDOW,Y1_coef_imag);
        [alpha2_imag,S2] = sas_salience(WINDOW,Y2_coef_imag);
        % compute match
        MA = sas_match(WINDOW,Y1_coef_imag,Y2_coef_imag,p);
        % selection
        m1 = MA > TH; m2 = S1 > S2;
        w1 = (0.5 - 0.5*(1-MA) / (1-TH));
        Y_coef_imag  = (~m1) .* ((m2.*Y1_coef_imag) + ((~m2).*Y2_coef_imag));
        Y_coef_imag  = Y_coef_imag + (m1 .* ((m2.*Y1_coef_imag.*(1-w1))+((m2).*Y2_coef_imag.*w1) + ...
            ((~m2).*Y2_coef_imag.*(1-w1))+((~m2).*Y1_coef_imag.*w1)));
        
    case 3
        % using DTCWT for the FLOM-Cauchy Algorithm
        % disp(' FLOM-Cauchy Algorithm...');
        
        % Real part:
        % compute salience (dispersion of SaS components)
        [alpha1_real,S1] = sas_salience_flomcauchy(WINDOW,Y1_coef_real);
        [alpha2_real,S2] = sas_salience_flomcauchy(WINDOW,Y2_coef_real);
        % compute match (symmetric coefficient of covariation)
        p=1;
        MA = sas_match(WINDOW,Y1_coef_real,Y2_coef_real,p);
        % selection
        m1 = MA > TH; m2 = S1 > S2;
        w1 = (0.5 - 0.5*(1-MA) / (1-TH));
        Y_coef_real  = (~m1) .* ((m2.*Y1_coef_real) + ((~m2).*Y2_coef_real));
        Y_coef_real  = Y_coef_real + (m1 .* ((m2.*Y1_coef_real.*(1-w1))+((m2).*Y2_coef_real.*w1) + ...
            ((~m2).*Y2_coef_real.*(1-w1))+((~m2).*Y1_coef_real.*w1)));
        
        % Imaginary part:
        % compute salience
        [alpha1_imag,S1] = sas_salience_flomcauchy(WINDOW,Y1_coef_imag);
        [alpha2_imag,S2] = sas_salience_flomcauchy(WINDOW,Y2_coef_imag);
        % compute match
        MA = sas_match(WINDOW,Y1_coef_imag,Y2_coef_imag,p);
        % selection
        m1 = MA > TH; m2 = S1 > S2;
        w1 = (0.5 - 0.5*(1-MA) / (1-TH));
        Y_coef_imag  = (~m1) .* ((m2.*Y1_coef_imag) + ((~m2).*Y2_coef_imag));
        Y_coef_imag  = Y_coef_imag + (m1 .* ((m2.*Y1_coef_imag.*(1-w1))+((m2).*Y2_coef_imag.*w1) + ...
            ((~m2).*Y2_coef_imag.*(1-w1))+((~m2).*Y1_coef_imag.*w1)));
        
    case 4
        % meridian distribution
        weights_real = weight_computation_Meridian(WINDOW,Y1_coef_real,Y2_coef_real);
        weights_imag = weight_computation_Meridian(WINDOW,Y1_coef_imag,Y2_coef_imag);
        Y_coef_real = Y1_coef_real .* weights_real(:,:,1) + Y2_coef_real .* weights_real(:,:,2);
        Y_coef_imag = Y1_coef_imag .* weights_imag(:,:,1) + Y2_coef_imag .* weights_imag(:,:,2);
        
    case 5
        % Cauchy probability
        % disp('Computing Cauchy probability...');
        weights_real = weight_computation(WINDOW,Y1_coef_real,Y2_coef_real);
        weights_imag = weight_computation(WINDOW,Y1_coef_imag,Y2_coef_imag);
        Y_coef_real = Y1_coef_real .* weights_real(:,:,1) + Y2_coef_real .* weights_real(:,:,2);
        Y_coef_imag = Y1_coef_imag .* weights_imag(:,:,1) + Y2_coef_imag .* weights_imag(:,:,2);
        
    otherwise,
        error('unkown option');
end

% combine coefficients
% -------------------------------------------------------------------------
Y = real(Y_coef_real) + 1i*real(Y_coef_imag);
 



