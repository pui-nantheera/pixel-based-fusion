function [imageout] = imaver(imagein,template)
%imaver.m, v 1.1 2002/01/22 13:25:40$
%===========================================================================
%               Eduardo Fernandez Canga - University of Bath
%
%                        Copyright (c) 2002
%===========================================================================
%
%                  [imageout] = imaver(imagein,template)
%
%               imagein : input image
%               template: averaging mask size
%
%===========================================================================

if (mod(template,2)==0)
    error(sprintf('Template size must be odd'));
end
tmpsiz=(template-1)/2;
imageaux=impad(imagein,tmpsiz);            %padd the input image to perform 
mask(1:template,1:template)=1/template^2;   %convoultion with the averaging
imageout=convn(imageaux,mask,'valid');      %mask, select valid area.
