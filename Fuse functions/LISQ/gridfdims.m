function [hx, hy, lx, ly] = gridfdims(F)
%------------------------------------------------------------------------------
%
% <Paul.de.Zeeuw@cwi.nl> dd001127
%
% This function determines the dimensions of gridfunction F.
%
% [hx, hy] = gridfdims(F) yields meshwidths in x- and y-dir.
%
% [hx, hy, lx, ly] = gridfdims(F) yields meshwidths in x- and y-dir
%                    and dimensions of domain of F.
%
%         x m
%     o---->
%     |
%   y |            [n,m] = size(F)
%     v
%   n
%
% See also: mupq
%------------------------------------------------------------------------------
[n, m] = size(F);
if nargout == 2
  if m<n
    hx = 1.0/m;
    hy = hx;
  else
    hy = 1.0/n;
    hx = hy;
  end
elseif nargout == 4
  if m<n
    lx = 1.0;
    ly = (n*lx)/m;
    hx = lx/m;
    hy = hx;
  else
    ly = 1.0;
    lx = (m*ly)/n;
    hy = ly/n;
    hx = hy;
  end
else
  error(' gridfdims - wrong number of output arguments ')
end
%------------------------------------------------------------------------------
