function P = putcolor11(A11, sizeP)
%------------------------------------------------------------------------------
%
% Gridfunction A11 is upsampled.
%
%
% Design and implementation by:
% Dr. Paul M. de Zeeuw   <Paul.de.Zeeuw@cwi.nl>   http://www.cwi.nl/~pauldz/
% Last Revision: June 6, 1999.
% Copyright 1999-2002 Stichting CWI, Amsterdam
%------------------------------------------------------------------------------
[n, m] = size(A11);
if nargin == 2
  nP = sizeP(1);
  mP = sizeP(2);
  if nP < 2*n
    error(' putcolor11 - 1st dimension of P too small ')
  end
  if mP < 2*m
    error(' putcolor11 - 2nd dimension of P too small ')
  end
elseif nargin == 1
  nP = 2*n+1;
  mP = 2*m+1;
else
  error(' putcolor11 - wrong number of arguments ')
end
P=reshape(linspace(0,0,nP*mP),nP,mP);
P(2:2:nP, 2:2:mP)=A11;
%------------------------------------------------------------------------------
