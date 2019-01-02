function AA00 = Q2R(EA11, A00)
%------------------------------------------------------------------------------
% UNUSED OBSOLETE ??
%------------------------------------------------------------------------------
[nEA11, mEA11]=size(EA11);
[nA00, mA00]  =size(A00);
if nEA11 ~= (nA00+1)
  error(' Q2R - EA11 and A00 do not match ')
end
if mEA11 ~= (mA00+1)
  error(' Q2R - EA11 and A00 do not match ')
end
%
ALD=stripL(stripD(EA11));
AA00=gmin(A00,ALD);       % N.B.
clear ALD;
%
AUR=stripU(stripR(EA11));
AA00=gmin(AA00,AUR);
clear AUR;
%
ALU=stripL(stripU(EA11));
AA00=gmin(AA00,ALU);
clear ALU;
%
ADR=stripD(stripR(EA11));
AA00=gmin(AA00,ADR);
clear ADR;
%
%------------------------------------------------------------------------------
