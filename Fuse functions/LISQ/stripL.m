function S = stripL(A)
%------------------------------------------------------------------------------
%
% Design and implementation by:
% Dr. Paul M. de Zeeuw   <Paul.de.Zeeuw@cwi.nl>   http://www.cwi.nl/~pauldz/
% Last Revision: June 6, 1999.
% Copyright 1999-2002 Stichting CWI, Amsterdam
%------------------------------------------------------------------------------
[n, m] = size(A);
Y=A.';
S=[Y(2:m,:)].';
%------------------------------------------------------------------------------
