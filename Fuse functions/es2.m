function Y = es2(X, n)
%Y = ES2(X, n) symmetric extension of a matrix on all borders
%    
%    X  - input matrix
%    n  - number of rows/columns to extend
%
%    Y  - extended matrix

%    (Oliver Rockinger 16.08.99)
%       Adapted for 3D matrix by Eduardo Fernandez Canga


[r c p] = size(X);
Y                        = zeros(r+2*n, c+2*n,p);
Y(n+1:n+r,n:-1:1,:)        = X(:,2:1:n+1,:); 
Y(n+1:n+r,n+1:1:n+c,:)     = X;
Y(n+1:n+r,n+c+1:1:c+2*n,:) = X(:,c-1:-1:c-n,:);
Y(n:-1:1,n+1:c+n,:)        = X(2:1:n+1,:,:); 
Y(n+r+1:1:r+2*n,n+1:c+n,:) = X(r-1:-1:r-n,:,:); 
 



