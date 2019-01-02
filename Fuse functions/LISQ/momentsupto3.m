function [mass, mus, orthos] = momentsupto3(F,N)
%------------------------------------------------------------------------------
%
% <Paul.de.Zeeuw@cwi.nl> dd021220
%
% Computation of central moments (upto 3rd order) is called for, followed by 
% call for similitude invariants, orthogonal invariants and similitude with 
% orthogonal invariants combined.
%
% Input:
% F         = Gridfunction defined on a rectangular grid
%
% Output:
% mass      = mass of gridfunction
% mus       = vector of length 7, containing the 2nd & 3rd order moments 
%             mu20 mu11 mu02 mu30 mu21 mu12 mu03
% sims      = Central moments made invariant to similitude transforms
% orthos    = Central moments made orthogonally invariant
% simorthos = Invariants w.r.t. both similitude and orthogonal transformations
%
% See pages 180, 181, 185 in:
% Ming-Kuei Hu  Visual Pattern Recognition by Moment Invariants,
% IRE Transactions on Information Theory, pp. 179--187 (1962).
%
% See also: masscenter, mupq, HUinvariants, Q1001momentsupto3
%------------------------------------------------------------------------------
% Mass and center
[c, mass] = masscenter(F);
%
% Central moments
mu20 = mupq(F, 2, 0, c);
mu11 = mupq(F, 1, 1, c);
mu02 = mupq(F, 0, 2, c);
%
mu30 = mupq(F, 3, 0, c);
mu21 = mupq(F, 2, 1, c);
mu12 = mupq(F, 1, 2, c);
mu03 = mupq(F, 0, 3, c);
%
mus = [mu20 mu11 mu02 mu30 mu21 mu12 mu03];
%
% Invariants
[sims,orthos,simorthos]= HUinvariants(mass, mus);
%
%------------------------------------------------------------------------------
