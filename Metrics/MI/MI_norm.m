function [h,h_norm]=MI_norm(image_1,image_2)

% Takes a pair of images and returns the mutual information Ixy using joint entropy function JOINT_H.m

h = joint_h_C(image_1,image_2); % Joint histogram
h = h/sum(h(:));
hy=sum(h);  % sum of the rows of normalized joint histogram
hx=sum(h,2); % sum of columns of normalized joint histogram


Hy= -sum(hy.*(log2(hy+(hy==0))));
Hx= -sum(hx.*(log2(hx+(hx==0))));
h_xy = -sum(sum(h.*(log2(h+(h==0)))));

% % Renyi:
% alph = 0.5;
% Hy = log2(sum(hy.^alph))/(1-alph);
% Hx = log2(sum(hx.^alph))/(1-alph);
% h_xy = log2(sum(h(:).^alph))/(1-alph);

h = Hx + Hy - h_xy;
% AL:
if nargout > 1
  h_norm = (Hx + Hy)/h_xy;
end
% Mutual information computed from eq. (6) of: J. Pluim, J. Maintz and M.
% Viergever "Mutual information based registration of medical images : a
% survey"




