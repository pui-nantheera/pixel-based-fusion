function Im = LoadIm(FileName)
%LoadIm.m, v 0.1 2004/02/17
%===========================================================================
%               Eduardo Fernandez Canga - University of Bristol
%
%                        Copyright (c) 2004
%===========================================================================
%     function Im = LoadIm (filename)
%
% Loads image stored in FileName and stores in Im (rgb double)


[Im,ma] = imread(FileName);
imfinfo(FileName);
r=imfinfo(FileName);
if strcmp(r.ColorType,'truecolor')
    Im = double(Im)/255;
elseif strcmp(r.ColorType,'indexed')
    Im =(ind2rgb(Im,ma));
elseif strcmp(r.ColorType ,'grayscale')
        if isfield(r,'MaxValue')
            p=r.MaxValue;
        elseif isfield(r,'BitDepth')
            p=2^r.BitDepth - 1;
        end


    Im = double(Im(:,:,[1 1 1]))/p;
end
