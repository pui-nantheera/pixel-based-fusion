function [c, mass] = masscenter(F)
%------------------------------------------------------------------------------
%
% <Paul.de.Zeeuw@cwi.nl> dd001211
%
% This function computes the center of mass of gridfunction F (seen as a  
% density distribution) .
%
% Beware if F assumes not mere positive values, does it make sense?
%
% See also: m00, m10, m01, mu10, mu01
%------------------------------------------------------------------------------
denomina = m00(F);
if denomina == 0
  error(' masscenter - demoninator vanishes ')
else
  cx = m10(F)/denomina;
  cy = m01(F)/denomina;  
end
% whether cx, cy < 0 etc. could be checked here (see above warning).
if nargout == 1
  c = [cx cy];
elseif nargout == 2
  c = [cx cy];
  mass = denomina;
else
  error(' masscenter - wrong number of output arguments ')
end
%------------------------------------------------------------------------------
