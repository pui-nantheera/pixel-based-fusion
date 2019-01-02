function SNR_DWT = DWTSNR_loop(I0,I1,J,wavf,recon)

%DWTSNR.m, 2004/06/03 Artur Loza
%
%===========================================================================
% 
%     function SNR_DWT = DWTSNR_loop(I0,I1[,J][,wavf][,recon])
% helper function that computes DWTSNR for a collection of images
% SNR_DWT = DWTSNR(I0,I1), I0 ref. image(s) R x C x N0, I1: R x C x N1
% images will be compared N0-->N1 (corresponding ones, and not not all 
% the combinations)
%
% Example: I0(:,:,1)=imread('pout.tif'); 
%          I1(:,:,1) = imnoise(I0,'salt & pepper', 0.02);
%          I1(:,:,2) = imnoise(I0,'gaussian',0,0.01); 
%          I0(:,:,2) = I0;
%          SNR_DWT = DWTSNR_loop(I0(:,:,1)  ,I1(:,:,1)  ,6,'bior6.8')
%          SNR_DWT = DWTSNR_loop(I0(:,:,1)  ,I1(:,:,1:2),6,'bior6.8')
%          SNR_DWT = DWTSNR_loop(I0(:,:,1:2),I1(:,:,1)  ,6,'bior6.8')
%          SNR_DWT = DWTSNR_loop(I0(:,:,1:2),I1(:,:,1:2),6,'bior6.8')
%          figure, 
%          subplot(131), imshow(I0(:,:,1)), 
%          subplot(132), imshow(I1(:,:,1)), 
%          subplot(133), imshow(I1(:,:,2))

if nargin < 3, 
  J = floor(min(log2(rc))); % max number of decomp levels 
end

if nargin < 4, 
  wavf = 'bior6.8'; % 'bior1.1' 'bior3.7' 'bior3.9' 'bior5.5' 'bior6.8'
end

if nargin < 5, 
  recon = 0; % use reconstructed coefficients
end

if nargin < 6
  c = []; % coefficients has not been supplied
end

nI0 = size(I0,3); 
nI1 = size(I1,3);
kmax = max(nI0,nI1); 
ki0 = nI0 > 1; % volume increments
ki1 = nI1 > 1; % 0 for one image, 1 for multiple

for k = 1 : kmax
  
  k0 = k*ki0 + ~ki0; % image index
  k1 = k*ki1 + ~ki1;
  
  % compute SNR
  [SNR_DWT(k),c] = DWTSNR(I0(:,:,k0),I1(:,:,k1),J,wavf,recon,c);
  
  if ki0 > 0, c(1,1,1) = nan; end % nan -- coeffs have not been computed
  if ki1 > 0, c(1,1,2) = nan; end
  
end