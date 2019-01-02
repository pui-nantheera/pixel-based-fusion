function Y = es(X, n, wo)
%Y = ES2(X, n) symmetric extension of a matrix on selected borders
%    
%    X  - input matrix
%    n  - number of rows/columns to extend
%    wo - where to extend
%       wo == 1: left and right
%				wo == 2: up and bottom
%				wo == 3: all borders
%
%    Y  - extended matrix

%    (Oliver Rockinger 16.08.99)

[r c p] = size(X);

if wo == 1
  Y = zeros(r, c+2*n,p);
  Y(:,n:-1:1,:)  = X(:,2:1:n+1,:); 
  Y(:,n+1:1:n+c,:) = X(:,:,:);
  Y(:,n+c+1:1:c+2*n,:) = X(:,c-1:-1:c-n,:); 
end;

if wo == 2
  Y = zeros(r+2*n, c,p);
  Y(n:-1:1,:,:)  = X(2:1:n+1,:,:); 
  Y(n+1:n+r,:,:) = X(:,:,:);
  Y(n+r+1:1:r+2*n,:,:) = X(r-1:-1:r-n,:,:); 
end;

if wo == 3
  Y = zeros(r+2*n, c+2*n,p);
  Y(n+1:n+r,n:-1:1,:)  = X(:,2:1:n+1,:); 
  Y(n+1:n+r,n+1:1:n+c,:) = X;
  Y(n+1:n+r,n+c+1:1:c+2*n,:) = X(:,c-1:-1:c-n,:);
  Y(n:-1:1,n+1:c+n,:)  = X(2:1:n+1,:,:); 
  Y(n+1:n+r,n+1:c+n,:) = X;
  Y(n+r+1:1:r+2*n,n+1:c+n,:) = X(c-1:-1:c-n,:,:); 
end; 