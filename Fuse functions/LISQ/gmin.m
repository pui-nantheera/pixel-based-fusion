function g = gmin(A, B)
%------------------------------------------------------------------------------
% UNUSED OBSOLETE??
%------------------------------------------------------------------------------
if ~all(size(A) == size(B))
  error(' gmin - dimensions do not match ')
end
D=(A<B);
g= D.*A+(~D).*B;
%------------------------------------------------------------------------------
