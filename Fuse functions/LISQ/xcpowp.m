function xcp = xcpowp(F, p, c)
%------------------------------------------------------------------------------
%
% <Paul.de.Zeeuw@cwi.nl> dd001127
%
% This function computes the gridfunction xcp of size F with values (x-c)^p at
% its gridpoints. Piecewise constant approximation is assumed.
%
%   Orientation
%
%         x
%     o---->
%     |
%   y |
%     v
%    
%   
%------------------------------------------------------------------------------
if isempty(F)
  error(' xcpowp - gridfunction is empty ')
else
  [n, m] = size(F);
  if p == 0
    xcp = ones(n,m);
  else  
    [hx, dummy1, lx, dummy2] = gridfdims(F);
    xfirst =  0 + hx/2 - c;
    xlast  = lx - hx/2 - c;
    lineofxcp = linspace(xfirst, xlast, m);
    xcp = lineofxcp(ones(1,n),:);   % Tony's trick
    if p ~= 1
      xcp = xcp.^p;
    end
  end
end  
%------------------------------------------------------------------------------
