function imksize(filename,option,method)
%imksize.m, v 0.1 2004/01/28
%===========================================================================
%               Eduardo Fernandez Canga - University of Bristol
%
%                        Copyright (c) 2004
%===========================================================================
%   function imksize(filename,option,[method])
%       
%       reads image in 'filename' and scale it into a 1024x1024 image using: 
%           option == 0 -> zero padding
%           option ~= 0 -> streching by interpolation
%               method == 1 -> nearest neighbour (default)
%               method == 2 -> bilinear
%               method == 3 -> bicubic


if nargin < 3 & option ~= 0
    method = 1;
end

[im,map]=imread(filename);
if isind(im) & ~isempty(map)  
   im = ind2rgb(im,map);
   im = uint8(256*im);
end


filename = filename(1:find(filename == '.')-1);

if option
    switch method
        case 1, im=imresize(im,[1024 1024]);
            filename = [filename '_sn.bmp'];
        case 2, im=imresize(im,[1024 1024],'bilinear');
            filename = [filename '_sl.bmp'];            
        case 3, im=imresize(im,[1024 1024],'bicubic');
            filename = [filename '_sc.bmp'];            
    end
else
    [r c x] = size(im);
    z = uint8(zeros(1024,1024,x));    
    for i=1:x
        z(floor((1024-r)/2)+(1:r),floor((1024-c)/2)+(1:c),i)=im(:,:,i);
    end
    im=z;

    clear [z r c x]   
    filename = [filename '_zp.bmp'];                
end


imwrite(im,filename)