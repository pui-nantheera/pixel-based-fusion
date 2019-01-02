function [mask] = gmask(sigma)
%gmask.m,v 1.1 2002/02/06 18:08:40
%===========================================================================
%                   Eduardo Fernandez - University of Bath
%
%                        Copyright (c) 2002
%===========================================================================
%
%		  [mask] = gmask(sigma)                   
% 

g = 3*sigma;
for i = -g : g
    for j = -g : g
        mask (i+g+1,j+g+1)= sqrt(i^2 + j^2);
    end
end
mask=exp(-mask/sigma);
mask=mask/sum(sum(mask));
   