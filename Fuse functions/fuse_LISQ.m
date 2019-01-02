function FUSEDIMG = fuse_LISQ(ORIGIMGS, LEVELS, HIGH, LOW, FILTERNO)
%   fuse_LISQ - Image fusion with Lifting Scheme on Quincunx Grids
%
%   FUSEDIMG = fuse_LISQ(ORIGIMGS, LEVELS, HIGH, LOW, FILTERNO)
%   ORIGIMGS is a 3-D matrix: input images
%   LEVELS is maximum decomposition level (default=4)
%   HIGH is coefficient selection highpass (default=1)  (see selc.m) 
%   LOW is coefficient selection base image (default=0) (see selb.m) 
%   FILTERNO is the index (1-7) of the filter number of the lifting scheme
%           1, filtername='maxmin' (default)
%           2, filtername='minmin';
%           3, filtername='maxmax';
%           4, filtername='Neville2';
%           5, filtername='Neville4';
%           6, filtername='Neville6';
%           7, filtername='Neville8';
%
%   FUSEDIMG is a 2-D matrix. It is a grayscale fused image.
%
%   Examples:   I(:,:,1) = im2double(imread('UNcamp_2N_i_20.bmp'));
%               I(:,:,2) = im2double(imread('UNcamp_2N_v_20.bmp'));
%               levels = 3;  % decomposition levels
%               high = 1;    % select max absolute coefficients
%               low = 0;     % average base level
%               filterno = 4 % Using Neville2 for filter, 
%               F = fuse_LISQ(I, levels, high, low, 4);
%               figure; imshow(F)

% Dr. Paul ORIGIMGS. de Zeeuw   <Paul.de.Zeeuw@cwi.nl>   http://www.cwi.nl/~pauldz/
% Report: http://www.cwi.nl/ftp/CWIreports/PNA/PNA-R0224.pdf
% Copyright 1999-2003 Stichting CWI, Amsterdam
%------------------------------------------------------------------------------
% We present an example of an application of LISQ to image fusion.
% Two similar, though different, images are fused into one image that is meant
% to unify the information included in both originals.
% It is not claimed to represent the present state of the art in image fusion.
%--------------------------------------------------------------------------
%Added to the Image Fusion Tool box (v5) April 2003 
%by John Lewis
%University of Bristol

%%  check inputs and put default values
if nargin < 2 || isempty(LEVELS)
    LEVELS = 4;
end
if nargin < 3 || isempty(HIGH)
    HIGH = 1;
end
if nargin < 4 || isempty(LOW)
    LOW = 0;
end
if nargin < 5 || isempty(FILTERNO)
    FILTERNO = 1;
end

n=size(ORIGIMGS,3);

switch (FILTERNO) 
    case 1, filtername='maxmin';
    case 2, filtername='minmin';
    case 3, filtername='maxmax';
    case 4, filtername='Neville2';
    case 5, filtername='Neville4';
    case 6, filtername='Neville6';
    case 7, filtername='Neville8';
end;

% Multilevel 2-D decomposition by the lifting scheme
% --------------------------------------------------
for i=1:n
    [Cd(:,:,i),Sd(:,:,i)] = QLiftDec2(ORIGIMGS(:,:,i),LEVELS,filtername);
end

% Create C from CA and CB (only the details)
% --------------------------------------------------
C = Cd(:,:,1);  S = Sd(:,:,1);
for level = 1:LEVELS
   rectgrids = mod(level, 2)+1;
   
   % extracts gridfunction
   if rectgrids == 2
       for i=1:n
           [F1(:,:,i), F2(:,:,i)] = retrieveQ1001(level, 'd', Cd(:,:,i), Sd(:,:,i));
       end
   else
       for i=1:n
           [F1(:,:,i)] = retrieveR(level, 'd', Cd(:,:,i), Sd(:,:,i));
       end
   end
   for no = 1:rectgrids
      if no == 1
        fc = F1;
      else
        fc = F2;
      end

      fcab = selc(fc, HIGH);
       
      if rectgrids == 2
        if no == 1
          [first, last] = whatcoef2QL(level, '10', 'd', S);
        else
          [first, last] = whatcoef2QL(level, '01', 'd', S);
        end
      else
        [first, last] = whatcoef2QL(level, 'none', 'd', S);
      end
      C(first:last) = fcab;
   end
   clear ('F1','F2')
end

% Update C from CA and CB (only the approximation)
% --------------------------------------------------
clear fc
for i=1:n
    fc(:,:,i)=retrieveR(level, 'a', Cd(:,:,i), Sd(:,:,i));
end


fcab = selb(fc, LOW);


[first, last] = whatcoef2QL(level, 'none', 'a', S);
C(first:last) = fcab;

% Reconstruction
% --------------------------------------------------
FUSEDIMG = QLiftRec2(C,S,filtername);