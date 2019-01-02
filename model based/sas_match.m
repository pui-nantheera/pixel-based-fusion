function [MA] = sas_match(windowsize,x1,x2,p)

% SAS_MATCH - Local estimates of symmetric coefficient of covariation
% for a SaS vector [x1 x2].
% [alpha,gama] = SAS_MATCH(windowsize,x1,x2,p) returns the symmetric
% covariation coefficient estimated between x1 and x2

% Reference:
% Garel B, d'Estampes L and Tjostheim D: "Revealing some unexpected
% dependence properties of linear combinations of stable random variables
% using symmetric covariations," Communications in Statistics - Theory and
% Methods, Vol. 33, No 4, pp. 769-786, 2004

%   Copyright (c) Alin Achim 16/11/2004. 
%   Alin.Achim@bris.ac.uk

windowfilt = ones(1,windowsize)/windowsize;

Nom1 = conv2(windowfilt,windowfilt,x1.*sign(x2).*x2.^(p-1),'same');
Nom2 = conv2(windowfilt,windowfilt,x2.*sign(x1).*x1.^(p-1),'same');
Den1 = conv2(windowfilt,windowfilt,(abs(x2)).^p,'same');
Den2 = conv2(windowfilt,windowfilt,(abs(x1)).^p,'same');
MA = (Nom1./(Den1+eps)).*(Nom2./(Den2+eps));