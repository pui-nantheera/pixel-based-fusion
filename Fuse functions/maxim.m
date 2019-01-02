function [maxmat] = maxim(mat1,mat2,th)
%maxim.m,v 1.11 2002/03/14 11:07:40
%===========================================================================
%               Eduardo Fernandez Canga - University of Bath
%
%                        Copyright (c) 2002
%===========================================================================
%
%             [maxmat] = maxim(mat1,mat2,[th])     
%
%	mat1, mat2: input images
%		th: threshold (default: 0)
%
%	    maxmat: output image, for each point it takes the maximum 
%		absoulte value of both images or the average value if
%		the difference is smaller than the threshold
%
%===========================================================================
%
%   

if nargin<3, th=0;end

decision1=abs(mat1)>(abs(mat2)+th);
decision2=(abs(mat1)+th)<abs(mat2);
both=(ones(size(decision1))-decision1-decision2)/2;

maxmat=mat1.*decision1+mat2.*decision2+(mat1+mat2).*both;