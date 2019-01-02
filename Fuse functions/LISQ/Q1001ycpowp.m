function [F10ycp, F01ycp] = Q1001ycpowp(F10, F01, p, c)
%------------------------------------------------------------------------------
%
% <Paul.de.Zeeuw@cwi.nl> dd001212
%
% This function computes the gridfunction {F10ycp U F01ycp} of size {F10 U F01}
% with values (x-c)^p at its gridpoints. Piecewise constant approximation is 
% assumed.
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
if isempty(F10) | isempty(F01)
  error(' Q1001ycpowp - at least one colour not present (empty) ')
else
  [n10, m10] = size(F10);
  [n01, m01] = size(F01);  
  if p == 0
    F10ycp = ones(n10, m10);
    F01ycp = ones(n01, m01);    
  else  
    [hx, hy] = Q1001gridfdims(F10, F01);
%    
    yfirst =  0 + hy/2 - c;
    ylast  =  yfirst + (m01 - 1) * 2 * hy;
    lineofycp = linspace(yfirst, ylast, n01)';
    F01ycp = lineofycp(:,ones(m01,1)); 
    if p ~= 1
      F01ycp = F01ycp.^p;
    end
%    
    yfirst =  0 + 3*hy/2 - c;
    ylast  =  yfirst + (m10 - 1) * 2 * hy;
    lineofycp = linspace(yfirst, ylast, n10)';
    F10ycp = lineofycp(:,ones(m10,1)); 
    if p ~= 1
      F10ycp = F10ycp.^p;
    end 
  end  
end  
%------------------------------------------------------------------------------
