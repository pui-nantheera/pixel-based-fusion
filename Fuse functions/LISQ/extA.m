function E =extA(A, c)
%------------------------------------------------------------------------------
%
% Design and implementation by:
% Dr. Paul M. de Zeeuw   <Paul.de.Zeeuw@cwi.nl>   http://www.cwi.nl/~pauldz/
% Last Revision: June 6, 1999.
% Copyright 1999-2002 Stichting CWI, Amsterdam
%------------------------------------------------------------------------------
E=extU(extD(extR(extL(A, c), c), c), c);
%------------------------------------------------------------------------------
