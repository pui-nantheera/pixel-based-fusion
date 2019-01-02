function [SF] = spfreq (mat)
%
%spfreq.m,v 1.11 2002/03/04 17:37:40
%===========================================================================
%               Eduardo Fernandez Canga - University of Bath
%
%                        Copyright (c) 2002
%===========================================================================
%
%       	[SF] = spfreq (mat)              
%
%	Calculates the spacial frequency of the input matrix
%
%===========================================================================
%
%   

[maxi,maxj]=size(mat);

rff=[1,-1];
RF=sum(sum((convn(mat,rff,'valid')).^2))/(maxi*maxj);
cff=[1;-1];
CF=sum(sum((convn(mat,cff,'valid')).^2))/(maxi*maxj);

SF=sqrt(RF+CF);