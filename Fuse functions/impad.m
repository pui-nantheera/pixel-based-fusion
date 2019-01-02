function [imageout] = impad(imagein,sizepad)
%impad.m,v 1.0 2001/10/29 17:35:40$
%===========================================================================
%               Eduardo Fernandez Canga - University of Bath
%
%                        Copyright (c) 2001
%===========================================================================
%
%            [imageout] = impad(imagein,sizepad)
%
%	Image padding using median filter
%

[r,c,n]=size(imagein);


for i=1:n
    % Fill the corners first
    imageout(1:(sizepad+1),1:sizepad+1,i)=median(median(imagein(1:2,1:2,i)));
    imageout(1:(sizepad+1),(0:sizepad)+c+sizepad,i)=...
        median(median(imagein(1:2,(c-1):c,i)));
    imageout((0:sizepad)+r+sizepad,1:(sizepad+1),i)=...
        median(median(imagein((r-1):r,1:2,i)));
    imageout((0:sizepad)+r+sizepad,(0:sizepad)+c+sizepad,i)=...
        median(median(imagein((r-1):r,(c-1):c,i)));

    % Then the original image in the centre
    imageout((1:r)+sizepad,(1:c)+sizepad,i)=imagein(:,:,i);

    % Finaly the sides
    % left - right
    for j=2:(r-1)
       imageout(j+sizepad,1:sizepad,i)=median(median(imagein(j+(-1:1),1:3,i)));
       imageout(j+sizepad,c+sizepad+(1:sizepad),i)=...
           median(median(imagein(j+(-1:1),c-(0:2),i)));
    end

    % up - down
    for j=2:(c-1)
       imageout(1:sizepad,j+sizepad,i)=median(median(imagein(1:3,j+(-1:1),i)));
       imageout(r+sizepad+(1:sizepad),j+sizepad,i)=...
           median(median(imagein(r-(0:2),j+(-1:1),i)));
   end
end