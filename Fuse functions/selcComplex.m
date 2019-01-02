function Y = selcComplex(COEFS, OPTIONS)
% selcComplex - coefficinet selection for highpass components
%
%   FUSEDCOEF = selc(COEFS, OPTIONS)
%   COEFS is a 3-D matrix
%   OPTIONS switch for selection type
%          OPTIONS(1) == 1: choose max(abs)
%          OPTIONS(1) == 2: salience / match measure with threshold == .75 (as proposed by Burt et al)
%               OPTIONS(2): area
%          OPTIONS(1) == 3: choose max with consistency check (as proposed by Li et al)
%               OPTIONS(2): area
%
%   FUSEDCOEF is a 2-D matrix. It is a grayscale combined coefficients.

%   v 1.0 16.08.99  Oliver Rockinger
%   v 2.0 12.11.03  Eduardo Fernandez Canga (University of Bristol)
%                   Modified for N input images
% -------------------------------------------------------------------------


[r c n] =size (COEFS);
% switch to method
switch(OPTIONS(1))
    case 1,
        % choose max(abs)
        map=max(abs(COEFS),[],3);
        map=map(:,:,ones(n,1));
        map=map==abs(COEFS);
        Y=sum(map.*COEFS,3)./sum(map,3);
        
    case 2, %not adapted for N input coeff yet
        % Burts method
        if n>2
            perror('Method not implemented for more than 2 input images');
            return
        end
        um = OPTIONS(2); th = .75;
        M1=COEFS(:,:,1);
        M2=COEFS(:,:,2);
        % compute salience
        S1 = conv2(es2(abs(M1).*abs(M1), floor(um/2)), ones(um), 'valid');
        S2 = conv2(es2(abs(M2).*abs(M2), floor(um/2)), ones(um), 'valid');
        % compute match
        MA = conv2(es2(abs(M1).*abs(M2), floor(um/2)), ones(um), 'valid');
        MA = 2 * MA ./ (S1 + S2 + eps);
        % selection
        m1 = MA > th; m2 = S1 > S2;
        w1 = (0.5 - 0.5*(1-MA) / (1-th));
        Y  = (~m1) .* ((m2.*M1) + ((~m2).*M2));
        Y  = Y + (m1 .* ((m2.*M1.*(1-w1))+((m2).*M2.*w1) + ((~m2).*M2.*(1-w1))+((~m2).*M1.*w1)));
        
    case 3,
        % Lis method
        if n>2
            perror('Method not implemented for more than 2 input images');
            return
        end
        um = OPTIONS(2);
        M1=COEFS(:,:,1);
        M2=COEFS(:,:,2);
        % first step
        % local maximum filters
        A1 = ordfilt2(abs(es2(M1, floor(um/2))), um*um, ones(um));
        A2 = ordfilt2(abs(es2(M2, floor(um/2))), um*um, ones(um));
        % second step
        mm = (conv2((A1 > A2), ones(um), 'valid')) > floor(um*um/2);
        Y  = (mm.*M1) + ((~mm).*M2);
        
    otherwise,
        error('unkown option');
end