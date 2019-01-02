function Y = dec2(X);
%Y = dec2(X) downsampling of a matrix by 2
%
%    X - input matrix
%
%    Y - output matrix

%    (Oliver Rockinger 16.08.99)

[r c n] = size(X);
Y     = X(1:2:r, 1:2:c,:);
