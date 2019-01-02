function X = QLiftRec2MinMin(C, S)
%-----------------------------------------------------------------------------
% QLiftRec2MinMin 
% Multilevel 2-D reconstruction by inverting the lifting scheme and using
% quincunx grids.
%
% The MinMin scheme has been proposed by Heijmans and Goutsias, see e.g.
%    H.J.A.M. Heijmans, J. Goutsias,
%    Multiresolution signal decomposition schemes.
%    Part 2: morphological wavelets.
%    CWI Report PNA-R9905, Amsterdam, 1999.
%    http://www.cwi.nl/ftp/CWIreports/PNA/PNA-R9905.pdf
%
% Calls for: 
%            retrieveQ1001, retrieveR,
%            getcolor01, getcolor10, getcolor00, getcolor11,
%            putcolor01, putcolor10, putcolor00, putcolor11.        
% See also: QLiftDec2MaxMin, QLiftRec2 
%
% Design and implementation by:
% Dr. Paul M. de Zeeuw   <Paul.de.Zeeuw@cwi.nl>   http://www.cwi.nl/~pauldz/
% Last Revision: February 3, 2003.
% Copyright 1999-2003 Stichting CWI, Amsterdam
%-----------------------------------------------------------------------------
% Note: argument list of this function might be extended with an argument that
%       points to another level than 1, nargin could be checked for this.
% Firstly, check input data
%
if isempty(C)
  error(' QLiftRec2MinMin - empty decomposition ');
else
  if isempty(S)
    error(' QLiftRec2MinMin - empty bookkeeping ');
  end
end
%
N = numoflevs(S); 
if mod(N, 2) == 1
  error(' QLiftRec2MinMin - only an even number of levels is accepted ');
end
%
%Secondly, start reconstruction
%
for lev=N:-2:1
%
%  The Inverse Scheme proceeds from rectangular grid to quincunx grid.
   if lev >= N
     APPROX00 = retrieveR(lev, 'a', C, S);     % "even slots"
   end
   DETAIL11 = retrieveR(lev, 'd', C, S);       % "odd slots"
   minO = min( min(min(APPROX00)), min(min(DETAIL11)) );
   maxO = max( max(max(APPROX00)), max(max(DETAIL11)) );
   cmin = minO-(maxO-minO);
   cmax = maxO+(maxO-minO);
%
%  Stage: undo update
   sizeA00 = size(APPROX00);
   Q0011A00 = APPROX00 - ...
              min(zeros(sizeA00), synA00Qmin(DETAIL11, sizeA00, cmax));
   clear APPROX00 sizeA00;      
%
%  Stage: undo predict
   Q0011A11 = DETAIL11 + synA11Qmin(Q0011A00, size(DETAIL11), cmax);
   clear DETAIL11; 
%     
%  Merge
%  The union Q0011 of Q0011A00 & Q0011A11 now contains the approximation on 
%  the next scale (with index lev-1).
%
%  The Inverse Scheme proceeds from quincunx grid to rectangular grid.
%
%  The "even slots" are in the colours 00 and 11,
%  the "odd slots"  are in the colours 10 and 01.
%
%  Stage: undo update
   [DETAIL10, DETAIL01] = retrieveQ1001(lev-1, 'd', C, S);
   A00 = Q0011A00 - ...
         min(zeros(size(Q0011A00)), synA00min(DETAIL10, DETAIL01, cmax));
   clear Q0011A00; 
   A11 = Q0011A11 - ...
         min(zeros(size(Q0011A11)), synA11min(DETAIL10, DETAIL01, cmax));
   clear Q0011A11;         
%
%  Stage: undo predict  
   A10 = DETAIL10 + synA10min(A11, A00, cmax);
   A01 = DETAIL01 + synA01min(A11, A00, cmax);
   sizeR = size(DETAIL10) + size(DETAIL01);
   clear DETAIL10 DETAIL01;
%
%  Merge
   APPROX00 = putcolor00(A00, sizeR) + putcolor11(A11, sizeR) + ...
              putcolor10(A10, sizeR) + putcolor01(A01, sizeR);
%  APPROX00 is used in the next iteration of the "for"-loop!
   clear A00 A11 A10 A01;
%   
%  Warning: not yet verified whether sizeQ and sizeR will be consistent for all
%           possible griddimensions of original image. 
%
   if (lev-1) == 1
     X = APPROX00; clear APPROX00;
   end
end
%-----------------------------------------------------------------------------
