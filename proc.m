function [imout,titlestr,impr2] = proc(ims,operation,param)
%proc.m, v 0.1 2004/02/17
%===========================================================================
%               Eduardo Fernandez Canga - University of Bristol
%
%                        Copyright (c) 2004
%===========================================================================
%     function [impr1, impr2, titlestr] = proc(ims,op,{param})
%
%  Inputs:  ims - Input Images to be fused
%            op - Fuse Method (see below)
%         param - Cell Array with Input Parameters needed (see below)
%                 First column contains the names of the parameters
%                 Second column contains the values of the parameters
%
% Outputs: imout - Output Images
%          impr2 - Secondary Output
%          titlestr - String with information about the performed process
%
%===========================================================================
%     
%
%    operation:
%      op = 11 -> Gaussian Noise 
%                Input parameters needed: var med
%
%      op = 12 -> Salt and Pepper Noise 
%                Input parameters needed: den
%
%      op = 13 -> Speckle Noise
%                Input parameters needed: var
%
%      op = 14 -> Poisson Noise (no input arguments)
%
%      op = 21 -> Segmentation
%                Input parameters needed: stype lvl
%                     stype: 1 or 'Joint' -> Joint segmentation
%                     stype: 2 or 'Unimode' -> Unimode segmentation
%
%      op = 31 -> Colouring
%                Input parameters needed: c (color [r g b])
%
%      op = 41 -> Motion Flow
%
%===========================================================================

imout=0;
impr2=0;

% read parameters
for i = 1:size(param,1)
    command = param{i,1};
    command = [command '= (param{i,2});'];
    eval(command);
end

% Remove the comment from the next line to gain access to the command
% window. Then, type param if you want to see the parameters you are
% sending to this function
% keyboard

imout = ims;

switch operation
case 11, %gaussian noise
    ims = checkgray(ims);
    imn = imnoise(ims(:,:,:,:,sele),'gaussian', medi ,vari);    
%    imn = double(imnoise(uint8(255*ims(:,:,:,:,sele)),'gaussian', medi ,vari))/255;
    if size(imn,3) == 1
        imout(:,:,:,:,sele)= imn(:,:,[1 1 1],:,:);
    else
        imout(:,:,:,:,sele)= imn(:,:,:,:,:);
    end
    titlestr=['DIST GAUS medi=' num2str(medi) '; vari=' num2str(vari) '; '];
    
case 12, %salt & pepper
    ims = checkgray(ims);
    imn = imnoise(ims(:,:,:,:,sele),'salt & pepper', dens);
%    imn = double(imnoise(uint8(255*ims(:,:,:,:,sele)),'salt & pepper', dens))/255;    
    if size(imn,3) == 1
        imout(:,:,:,:,sele)= imn(:,:,[1 1 1],:,:);
    else
        imout(:,:,:,:,sele)= imn(:,:,:,:,:);
    end
    titlestr=['DIST SALT dens=' num2str(dens) '; ']; 
    
case 13, %speckle
    ims = checkgray(ims);
    imn = imnoise(ims(:,:,:,:,sele),'speckle', vari);           
%    imn = double(imnoise(uint8(255*ims(:,:,:,:,sele)),'speckle', vari))/255;           
    if size(imn,3) == 1
        imout(:,:,:,:,sele)= imn(:,:,[1 1 1],:,:);
    else
        imout(:,:,:,:,sele)= imn(:,:,:,:,:);
    end
    titlestr=['DIST SPEC vari=' num2str(vari) '; '];
    
case 14, %poisson
    ims = checkgray(ims);
    imn = double(imnoise(uint8(255*ims(:,:,:,:,sele)),'poisson'))/255;           
    if size(imn,3) == 1
        imout(:,:,:,:,sele)= imn(:,:,[1 1 1],:,:);
    else
        imout(:,:,:,:,sele)= imn(:,:,:,:,:);
    end
    titlestr='DIST POIS ';
    
case 21, %segmentation only gray and a single frame
    ims = rgb2grayn(ims);
    if size(ims,4)>1
        txt=sprintf('\nSegmentation has not been implemented for sequences yet.');
        error(txt);
    end
    ims=squeeze(ims);
    imout=ims;

    [imout(:,:,sele),intolay,impr2,intmap] = Segment_Rob(ims(:,:,sele), stype, lvl);
    
    titlestr=['SEGM stype=' num2str(stype) '; lvl=' num2str(lvl) ';'];
    
    imout=permute(imout,[1 2 4 5 3]);
    impr2=permute(impr2,[1 2 4 5 3]);
    
case 31,
    imout(:,:,1,1:end,sele) = (c(1)*ims(:,:,1,:,sele));
    imout(:,:,2,1:end,sele) = (c(2)*ims(:,:,2,:,sele));
    imout(:,:,3,1:end,sele) = (c(3)*ims(:,:,3,:,sele));

    titlestr = ['COLO c=[' num2str(c) '];'];
case 41,
    error('Method not implemented yet')
end    
titlestr = [titlestr ' sele=' num2str(sele) ';'];

if size(imout,3) == 1;
    imout = imout(:,:,[1 1 1],:,:);
end



%==========================================================================
function ims = checkgray(ims)
if all(all(ims(:,:,1) == ims(:,:,2) & ims(:,:,1) == ims(:,:,3))) %grayscale
    ims = ims(:,:,1,:,:);
end
