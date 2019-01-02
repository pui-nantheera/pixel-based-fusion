function [imageout] = colpad(imagein,numcols)
%colpad.m,v 1.0 2001/10/29 17:35:40$
%===========================================================================
%               Eduardo Fernandez Canga - University of Bath
%
%                        Copyright (c) 2001
%===========================================================================
%
%            [imageout] = colpad(imagein,numcols)

[r,c,n]=size(imagein);

imageaux=impad(imagein,numcols);
imageout=imageaux(1+numcols:r+numcols,1+numcols:c+2*numcols,:);