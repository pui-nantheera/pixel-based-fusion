function out=impstrip(imponescale,crop);

imponescale = imponescale(crop+1:end-crop,crop+1:end-crop,:);
sz = size(imponescale);
temp = reshape(imponescale,sz(1),6*sz(2));
out = [real(temp);imag(temp)];