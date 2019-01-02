function Y = adb(X, bd)
%Y = ADB(X, bd) adds rows and/or columns by duplication on the 
%    lower resp. right side of the input matrix
%    
%    X     - input matrix
%    bd(1) - number of rows to add 
%    bd(2) - number of columns to add
%  
%    Y     - extended matrix

%    (Oliver Rockinger 16.08.99)

[r c n] = size(X);

% copy interior
Y            = zeros(r+bd(1),c+bd(2),n);
Y(1:r,1:c,:) = X;

% add rows
if (bd(1) > 0)
  Y(r+1:r+bd(1),1:c,:) = X(r-1:-1:r-bd(1),1:c,:);
end;
% add columns
if (bd(2) > 0)
  Y(1:r,c+1:c+bd(2),:) = X(1:r,c-1:-1:c-bd(2),:);
end;
% add corner
if (bd(1) > 0 & bd(2) > 0)
  Y(r+1:r+bd(1),c+1:c+bd(2)) = X(r-1:-1:r-bd(1),c-1:-1:c-bd(2));
end;  