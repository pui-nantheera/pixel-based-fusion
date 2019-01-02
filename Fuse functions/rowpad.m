function [imageout] = rowpad(imagein,numrows)
%rowpad.m,v 1.0 2001/10/29 17:35:40
%===========================================================================
%               Eduardo Fernandez Canga - University of Bath
%
%                        Copyright (c) 2001
%===========================================================================
%
%            [imageout] = rowpad(imagein,numrows)

[maxi,maxj,n]=size(imagein);

imageaux=impad(imagein,numrows);
imageout=imageaux(1+numrows:maxi+2*numrows,1+numrows:maxj+numrows,:);