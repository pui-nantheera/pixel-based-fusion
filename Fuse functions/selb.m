function Y = selb(COEFS, OPTIONS, I)
% selb - coefficient selection for base image
% 
%   FUSEDCOEF = selc(COEFS, OPTIONS)
%   COEFS is a 3-D matrix
%   OPTIONS switch for selection type
%          OPTIONS ==  0: average
%          OPTIONS ==  1: select first coefficients
%          OPTIONS ==  2: select second coefficients
%          OPTIONS ==  n: select N th coefficients
%          OPTIONS == -1: PCA 
%               for this option, I is require. It is 3D matrix of the original images.
%          OPTIONS == -2: max is selected
%          OPTIONS == -3: min is selected
%          OPTIONS == -4: variance is used
%          OPTIONS == -5: weight to midtone (gray)
%          OPTIONS == -6: min is selected and shadow areas are adjusted (for IR/VI)
%          OPTIONS == -7: curve adjust and min is selected (for IR/VI)
%
%   FUSEDCOEF is a 2-D matrix. It is a grayscale combined coefficients.

%   v 1.0 16.08.99  Oliver Rockinger
%   v 2.0 12.11.03  Eduardo Fernandez Canga (University of Bristol)
%                   Modified for N input images
% -------------------------------------------------------------------------

dimens= size(COEFS,3);

switch 1
    case OPTIONS==0, Y = mean(COEFS,3);
    case OPTIONS>0 & OPTIONS<=dimens & OPTIONS==round(OPTIONS),
        Y = COEFS(:,:,OPTIONS);
        % AL:
    case OPTIONS==-1, % PCA#
        if nargin >= 3
            X0 = reshape(I, size(I,1)*size(I,2),[]);
            X1 = reshape(COEFS, size(COEFS,1)*size(COEFS,2),[]);
            
            [EVecs, EVals] = eig(cov(X0));
            EVals = diag(EVals);
            
            Y = X1*EVals/sum(EVals);
            Y = reshape(Y, size(COEFS,1),[]);
        else
            Y = fuse_pca(COEFS);
        end
        
    case OPTIONS==-2, % max
        Y = max(abs(COEFS),[],3);
        
     case OPTIONS==-3, % min
        Y = min(abs(COEFS),[],3);
        
    case OPTIONS==-4, % ~vars
        
        X1 = reshape(COEFS, size(COEFS,1)*size(COEFS,2),[]);
        
        EVals = var(X1)';
        
        Y = X1*EVals/sum(EVals);
        
        Y = reshape(Y, size(COEFS,1),[]);
        
    case OPTIONS==-5, % weight gray (midtone)
        midValue = max(abs(COEFS(:)))/2;
        w = midValue - abs(midValue - abs(COEFS)) + eps;
        w = w./repmat(sum(w,3), [1 1 size(COEFS,3)]);
        Y = sum(w.*COEFS,3);
        
    case OPTIONS==-6 % for IR1
        minY = min(abs(COEFS),[],3);
        avgY = mean(abs(COEFS),3);
        th = max(abs(COEFS(:)))*0.2;
        w = th - minY;
        w = w./max(w(:));
        w(w<0) = 0;
        Y = w.*avgY*0.8 + (1-w).*minY;
        
    case OPTIONS==-7 % for IR2
        maxValue = max(abs(COEFS(:)));
        midValue = maxValue/2;
        % mark upper value
        markUp = abs(COEFS)>midValue;
        % parameters
        c = 0.6;
        m = (1-c)/0.5;
        g = 0.01;
        b = 50;
        a = 10;
        xx = abs(COEFS)./maxValue;
        xx = max(eps, min(1-eps, xx));
        lowerPart = (b - a.*(log(1./xx -1))).*g.*(m.*xx+c);
        lowerPart = max(0,min(1,lowerPart));
        xx = 1 - abs(COEFS)./maxValue;
        xx = max(eps, min(1-eps, xx));
        upperPart = 1 - (b - a.*(log(1./xx -1))).*g.*(m.*xx+c);
        upperPart = max(0,min(1,upperPart));
        adjX = lowerPart;
        adjX(markUp) = upperPart(markUp);
        adjX = adjX.*maxValue;
        % weight average
        Y = min(adjX,[],3);
        Y = Y.*sum(COEFS,3)./abs(sum(COEFS,3));
        
        
        
    otherwise, error('unknown option');
end;
