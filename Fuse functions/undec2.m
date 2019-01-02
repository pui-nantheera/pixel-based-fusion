function Y = undec2(X)
%Y = undec2(X) upsampling of a matrix by 2
%
%    X - input matrix
%
%    Y - output matrix

%    (Oliver Rockinger 16.08.99)

[r c n] = size(X);
Y     = zeros(2*r, 2*c,n); 

Y(1:2:2*r,1:2:2*c,:) = X;







