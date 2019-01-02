function [alpha,gama] = sas_salience_flomcauchy(windowsize,data)

% SAS_SALIENCE - Local estimates for Flom-Cauchy when alpha=1
% distribution function parameters from array DATA.
% The estimation is based on the method of second order cumulants.
% [alpha,gama] = sas_salience_flomcauchy(windowsize,data) returns the parameters
% ALPHA and GAMMA estimated in neighborhoods around every point in 'DATA'

%   Author: Alin Achim 15/11/2004.
%   Alin.Achim@bris.ac.uk
%   Revised by Author:21/02/2006
% Modified: Mayank Agarawal
%   ma6879@bris.ac.uk
%Last Modified: 26.02.2010

windowfilt = ones(1,windowsize)/windowsize;
k1 = conv2(windowfilt,windowfilt,log(abs(data+eps)),'same');
k2 = conv2(windowfilt,windowfilt,(log(abs(data+eps))-k1).^2,'same');
% alpha = min(1./sqrt(6*k2/pi^2-0.5),2);
%Cauchy case
alpha=1;
gama = exp(alpha.*k1+psi(1)*(1-alpha)).^(1./alpha);
% gama = alpha.^(1./alpha).*gama;% To switch from S0 parametrization to
%S2 parametrization. That is to make gama in the Gaussian case be equal
%to the standard deviation. If S0 was fine for you, then just comment
%out this last line!