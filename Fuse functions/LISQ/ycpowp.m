function ycp = ycpowp(F, p, c)
%------------------------------------------------------------------------------
%
% <Paul.de.Zeeuw@cwi.nl> dd001127
%
% This function computes the gridfunction ycp of size F with values (y-c)^p at 
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
  error(' ycpowp - gridfunction is empty ')
else
  [n, m] = size(F);
  if p == 0
    ycp = ones(n,m);
  else   
    [dummy1, hy, dummy2, ly] = gridfdims(F);
    yfirst =  0 + hy/2 - c;
    ylast  = ly - hy/2 - c;
    lineofycp = linspace(yfirst, ylast, n)';
    ycp = lineofycp(:,ones(m,1)); 
    if p ~= 1
      ycp = ycp.^p;
    end
  end
end  
%------------------------------------------------------------------------------
