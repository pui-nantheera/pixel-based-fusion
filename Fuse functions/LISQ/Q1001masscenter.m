function [c, mass] = Q1001masscenter(F10, F01)
%------------------------------------------------------------------------------
%
% <Paul.de.Zeeuw@cwi.nl> dd001212
%
% This function computes the center of mass of gridfunction {F10 U F01}
% (seen as a density distribution).
%
% Beware if {F10 U F01} assumes not mere positive values, does it make sense?
%
% See also: m00Q1001, m10Q1001, m01Q1001, masscenter
%------------------------------------------------------------------------------
denomina = m00Q1001(F10, F01);
if denomina == 0
  error(' Q1001masscenter - demoninator vanishes ')
else
  cx = m10Q1001(F10, F01)/denomina;
  cy = m01Q1001(F10, F01)/denomina;  
end
% whether cx, cy < 0 etc. could be checked here (see above warning).
if nargout == 1
  c = [cx cy];
elseif nargout == 2
  c = [cx cy];
  mass = denomina;
else
  error(' Q1001masscenter - wrong number of output arguments ')
end
%------------------------------------------------------------------------------
