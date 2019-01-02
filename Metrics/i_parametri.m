function [a_j,a_o,varargout]=i_parametri(a,varargin);

% Racunanje ivicnih parametara, jacine i orijentacije
% ulazne slike a, za svaki piksel
% velicina moze da definise da filtriranje vraca samo sredisnji deo slike bez sirenja
%    [a_j,a_o,[sx,sy]]=i_parametri(a,[velicina])
% V P, Februar 2000: 

% Calcualting of the edge parameters, strength and
% orientation of the input image, while for each pixel size can be defined in a way that the 
% filtering returns only the middle part of the image without widening 

% Definsemo marginu: 
% Defining margin
marg=0.000034;

% Definisemo velicinu rezultata: 
% Defining the size of the result
if length(varargin)>0 & strcmp(varargin{1},'valid')
  velicina='valid';
else
  velicina='same';
end

% Racunamo sobel vrednosti: 
% Calculating Sobel values
[sx_a,sy_a]=sobel_xy(a,velicina);

% Racunamo ivicne parametre: 
% Calc. edge parameters
a_j=sqrt(sx_a.^2+sy_a.^2);
a_o=atan((sy_a+marg)./(sx_a+marg));

% Ukoliko vracamo istu velicinu nuliramo ivice u prvim i poslednjim
% redovima i kolonama: If we are returning the same size, we are putting
% the first and the last rows and columns to zero value.
if strcmp(velicina,'same')
  [v,s]=size(a_j);
  a_j(1,:)=zeros(1,s); a_j(:,1)=zeros(v,1); a_j(v,:)=zeros(1,s); a_j(:,s)=zeros(v,1); 
end

% Ukoliko je trazeno vracamo x i y gradijente: If required, returning x and
% y gradients
if nargout>2,   varargout(1)={sx_a};    end
if nargout>3,   varargout(2)={sy_a};    end