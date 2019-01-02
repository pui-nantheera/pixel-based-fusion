function SNR_DWT = DWTSNR(I0,I1,J,wavf,recon)

%DWTSNR.m, 2004/06/03 Artur Loza
%
%===========================================================================
% 
%     function SNR_DWT = DWTSNR(I0,I1[,J][,wavf][,recon])
%
% Distortion measure computed as SNR in the wavelet domain 
%
%  Inputs:  I0,I1 - input images (orignal first); uint8 or double
%           J     - number of levels (optional, default floor(min(log2(rc))))
%           wavf  - wavelet (as in Wavelet Toolbox, optional, default 'bior6.8')
%           recon - use reconstructed coefficients (optional, default 0)
%
% Outputs:  SNR_DWT - distortion measure
%
%===========================================================================
%
% Example: I0=imread('pout.tif');
%          I1 = imdistor(I0,[30 8 4],'grid','PSNR');
%          SNR_DWT = DWTSNR(I0,I1)
%
%===========================================================================    

rc= size(I0);
k = 1:rc(1); % spatial indices 
l = 1:rc(2);

if nargin < 3, 
  J = floor(min(log2(rc))); % max number of decomp levels 
end

if nargin < 4, 
  wavf = 'bior6.8'; % 'bior1.1' 'bior3.7' 'bior3.9' 'bior5.5' 'bior6.8'
end

if nargin < 5, 
  recon = 0; % use reconstructed coefficients
end

p = 2; % norm for the SNR_DWT equation
s =.5; % masking

I(:,:,1) = double(I0); 
I(:,:,2) = double(I1);


for K = 1 : 2 % orignal first, then distorted
  
  [C,S] = wavedec2(I(:,:,K),J,wavf); 
  
  if recon, cA = wrcoef2('a',C,S,wavf,J); % reconstruct 
  else cA = appcoef2(C,S,wavf,J); % extract the level J approximation coefficients 
  end
  
  for j = J :-1: 1
    
    n = 2^j; 
    kn = ceil(k(2:2:end)/n); % spatial indices at level j (in the paper floor?)
    ln = ceil(l(2:2:end)/n); % even b/c dyaddown keeps even samples
    
    if recon
      
      cH = wrcoef2('h',C,S,wavf,j); % reconstruct 
      cV = wrcoef2('v',C,S,wavf,j); 
      cD = wrcoef2('d',C,S,wavf,j); 
      
    else
      
      [cH,cV,cD] = detcoef2('all',C,S,j); % extract the j detail coefficients 
      
      cH = wkeep(cH,ceil(rc/n)); % take the central part
      cV = wkeep(cV,ceil(rc/n));
      cD = wkeep(cD,ceil(rc/n));
      
      cH = cH(kn,ln); % replicate spatially corresponding coeffs
      cV = cV(kn,ln); 
      cD = cD(kn,ln);
      
    end    
    
    c(:,j,K) = [cH(:); cV(:); cD(:)];
    
  end % j
  
end % K

j = J:-1:1;  
jsp2 = repmat(2.^(-j*s*p),[size(c,1) 1]);

SNR_DWT = ...
  20 * log10( ( ...
  sum( max( jsp2.*abs(c(:,:,1)         ).^p ,[],2) ) / ...
  sum( max( jsp2.*abs(c(:,:,1)-c(:,:,2)).^p ,[],2) ) )^(1/p));