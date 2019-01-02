function [y, Nsigfused, Nsig] = FUSEcontr_enh(X1, X2, FMETHOD, LEVELS, ENHMETHOD)


% LEVELS = number of scales
% X = input image
% FMETHOD = fuse method (set to GGD)
% ENHMETHOD = enhancement method
% Nsig  = noise signal
% Nsigfused = 
% SETTINGS --------------------------------------------------
if nargin < 3, LEVELS = floor(log2(min([size(X1,1) size(X1,2)]))-4); end

% [y Nsigfused Nsig] = deal(0);

N = 3;% 2.5; 7;%% LEVELS = 5;
sas = 0;
totalImgs = size(X1,3);
convsize = 'same';%'valid';%
dwt_opt = 'cplx'; RI = 2; % dwt_opt = 'real'; RI = 1;
FMETHOD = upper(deblank(FMETHOD));
th = 0.75; % threshold

switch FMETHOD %       win b den ma bi
  case 'MAX', params = {1  1  0  0  0}; th = inf;
  case 'AVE', params = {1  1  0  0  0}; th = 0;
  case 'WA' , params = {N  2  0  1  0};
  case 'WAN' , params = {N  2  0  1  0};
  case 'LAP', params = {N  1  0  0  0}; % univar Laplacian (local)
  case 'BLA', params = {N  1  0  0  1}; % bivar  Laplacian (local)
  case 'ULS', params = {N  1  1  0  0}; % univar Laplacian shrinkage  (local)
  case 'BLS', params = {N  1  1  0  1}; % bivar  Laplacian shrinkage  (local)
  case 'GGD', params = {N nan 0  0  0}; % GGD (with B estimation)
  case 'CAU', params = {N  1  0  1  0}; sas = 1; % univar Cauchy (local?)
  case 'BCA', params = {N  1  0  1  1}; sas = 1; % bivar  Cauchy (global?)
  case 'UCS', params = {N  1  1  1  0}; sas = 1; % univar Cauchy shrinkage (global)
  case 'BCS', params = {N  1  1  1  1}; sas = 1; % bivar  Cauchy shrinkage (global)
  case 'SAS', params = {N nan 0  1  0}; sas = 1; % SAS (with estimation)?
  case 'SHM', params = {1  1  1  0  1}; th = inf; % MAX and shrinkage

end

%win th     b  den        ma     bi
[win,b,denois_fuse,sign_ma,parent] = deal(params{:});

% !!! imposing other than default (below, if any)!
% here...
% IDEA: for multiple inputs only selection mode, by thr=inf;
% until I fix cov for multiple inputs
if totalImgs > 2 && ~strcmp(FMETHOD,'WAN')
  th = inf;
end % or use mtple inputs only w/WAN

% to bypass fusion set this:  ----------------------------
if totalImgs == 1
  th = inf; %totalImgs = 1;% x2 = x1;
end

% determine the method for b
if numel(b) == 1,
  B(1:totalImgs) = {b};
else
  B{1:totalImgs} = b;
end
if isnan(b),
  b_est = 1;
else
  b_est = 0;
end

% averaging window
if win == 2.5,
  winfilt = [0 1 0; 1 1 1; 0 1 0];
elseif win == 1
  winfilt = 1;
else   %   winfilt = GausWin2(1.1,win);%   winfilt = winfilt * winfilt';
  winfilt = ones(win);
end
winfilt = winfilt/sum(winfilt(:));


%% DT-CWT --------------------------------------------------------
% Forward dual-tree DWT, either FSfarras or AntonB function can be used to combute the stage 1 filters
load nor_dualtree % run normaliz_coefcalc_dual_tree to generate this mat file.
[Faf, Fsf] = AntonB; % [Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;

% if parent, J1 = LEVELS+1; else J1 = LEVELS; end
Nsig = zeros(1,totalImgs);
for imi = 1 : totalImgs
  %   W12{imi} = DTree2D(X(:,:,imi),LEVELS,J1,Faf,af,nor,dwt_opt); % symextend inside - duplicates pixels.
  [W121{imi} L1{imi}] = DTree2Djp1(X1(:,:,imi),LEVELS,Faf,af,nor,dwt_opt);
  [W122{imi} L2{imi}] = DTree2Djp1(X2(:,:,imi),LEVELS,Faf,af,nor,dwt_opt);
  % Noise variance estimation using robust median estimator..
  tmp = W121{imi}{1}{1}{1}{1}; %me: {3} Selesnic uses {1}{1}{1}{1} but this seems to work better
  Nsig(imi) = median(abs(tmp(:)))/0.6745; % Nsig2 = Nsig^2; % (1:2) Nsig(1,1,imi)
end
W = cell(size(W121{1}));


%% ANALYSE COEFFS -----------------------------------------------
for scale = 1 : LEVELS, clear Y* S* A % 3 2 1 1

  % pre-locate the matrices
  szmat= GetCplxCoeffs(W12,1,1,scale,1,1,dwt_opt,1,parent);
  [Yj001 Yjp11] = deal(zeros([size(szmat) totalImgs]));
  [Yj002 Yjp12] = deal(zeros([size(szmat) totalImgs]));
  AA = deal(zeros([size(szmat) 3]));
  [S2 S2b] = deal(zeros(size(szmat)));

  for dir = 1 : 2, clear AA % orientation
    for dir1 = 1 : 3 % orientation
      for re_im = 1 : RI % real\imaginary
        for imi = 1 : totalImgs % # of inputs

          % complex coefficients
          [Yj001(:,:,imi),Yjp11(:,:,imi)] = GetCplxCoeffs( ...
            W121,imi,re_im,scale,dir,dir1,dwt_opt,RI,parent);

          [Yj002(:,:,imi),Yjp12(:,:,imi)] = GetCplxCoeffs( ...
            W122,imi,re_im,scale,dir,dir1,dwt_opt,RI,parent);


%% ESTIMATE shape parameter ----------------------------------------

          % diagnostics: Yj00(:,:,imi) = Yj00(:,:,1);%ones(size());

          if sas

            B{imi} = ShapeParSAS(Yj001(:,:,imi),b_est,B{imi});
            %[S2(:,:,imi) B{imi}] = SaliencySAS(Yj00(:,:,imi),winfilt,b_est,B{imi});

          else  % GGD

            B{imi} = ShapeParGGD(Yj001(:,:,imi),b_est,B{imi});
            %[S2(:,:,imi) B{imi}] = SaliencyGGD(Yj00(:,:,imi),Yjp1(:,:,imi),winfilt,b_est,B{imi},FMETHOD,convsize);

          end

        end % imi

%% ====================================================================
%% COMPUTE SALIENCY --------------------------------------------------

        b = min([B{:}]);%mean([B{:}]); %         if b_est, b = mean([B{:}]); end

        if sas

          b = max(1,b);
          for imi = 1 : totalImgs
            %             S2(:,:,imi) = SaliencySAS(Yj00(:,:,imi),             0,winfilt,0,b);
            %             S2b(:,:,imi) = SaliencySAS(Yj00(:,:,imi),Yjp1(:,:,imi),winfilt,0,b);
            S2(:,:,imi) = SaliencySAS(Yj001(:,:,imi),Yjp11(:,:,imi),winfilt,0,b);
          end

          %           MA = MatchMeasureSAS(Yj00,  0  ,winfilt,b);
          %           MAb = MatchMeasureSAS(Yj00,Yjp1,winfilt,b);
          MA = MatchMeasureSAS(Yj001,Yjp11,winfilt,b); % MAb = MA;
          S2b = S2;
          S2g = S2b;

        else % GGD

          % Yjp10 = Yjp1;

          for imi = 1 : totalImgs
            %             S2(:,:,imi) = SaliencyGGD(Yj00(:,:,imi),    0         ,winfilt,0,b,FMETHOD,convsize);

            mu_log_x = conv2a( log(abs(Yj001(:,:,imi)) + eps) + ...
              log(abs(Yjp1(:,:,imi)) + eps), winfilt, convsize ); % variance
            S2(:,:,imi) = exp( mu_log_x - psi(1./b)/b ); % this is s not \sigma^2, use only for decision map

            S2b(:,:,imi) = SaliencyGGD(Yj001(:,:,imi),Yjp1(:,:,imi),winfilt,0,b,FMETHOD,convsize);
            %             S2b(:,:,imi) = SaliencyGGD(Yj00(:,:,imi),Yjp1(:,:,imi),winfilt,0,b,FMETHOD,convsize);
            %             S2(:,:,imi) = SaliencyGGD(Yj00(:,:,imi),Yjp1(:,:,imi),winfilt,0,b,FMETHOD,convsize);
          end %  Yjp1(:,:,imi)
          S2g = S2b;

          %           [MA, ma] = MatchMeasureGGD(Yj00, 0  ,S2 ,winfilt,b,FMETHOD,sign_ma,convsize);
          %           MAb      = MatchMeasureGGD(Yj00,Yjp1,S2b,winfilt,b,FMETHOD,sign_ma,convsize);
          %         MA      = MatchMeasureGGD(Yj00,Yjp1,S2b,winfilt,b,FMETHOD,sign_ma,convsize);
          %           MA    = MatchMeasureGGD(Yj00,Yjp1,S2,winfilt,b,FMETHOD,sign_ma,convsize);

          %           for imi = 1 : totalImgs
          %             S3(:,:,imi) = SaliencyGGD(Yj00(:,:,imi),Yjp1(:,:,imi),winfilt,0,b,FMETHOD,convsize);
          %           end % Yjp1(:,:,imi)
          %           MA = MatchMeasureGGD(Yj00,Yjp1,S2,winfilt,b,FMETHOD,sign_ma,convsize);

        end

%% DENOISE (or not) -----------------------------------------------------
        T1 = 0;
        if denois_fuse

          if sas

            for imi = 1 : totalImgs
              if re_im==1, % estimate from real only
                K = 2:4; % 1:5 for global?
                gamm = CauchynoisyecfMomFit(Yj001(:,:,imi),Yjp1(:,:,imi),Nsig(imi),winfilt,K,parent);
                S2g(:,:,imi) = gamm;
              end

              %CauchynoisyecfMomFit(Yj00(:,:,imi),Yjp1(:,:,imi),Nsig(imi),winfilt,K,parent);
              Yj00(:,:,imi) = biCauchyshrink(Yj001(:,:,imi),Yjp11(:,:,imi),gamm,Nsig(imi),parent);

              %           for imi = 1 : totalImgs
              %             S2(:,:,imi) = SaliencySAS(Yj00(:,:,imi), 0 ,winfilt,0,b);
              %           end %  Yjp1(:,:,imi)
              %
              %           MA = MatchMeasureSAS(Yj00,winfilt,b);

            end



          else % GGD
            % (1) LAP s2, (2) Gauss s2 (3) ML est. according to the model
            %             Yjp1 = Yjp10;
            %             % compute LAP variance
            %             % s2 = 2 * conv2a( abs(Yj00), winfilt, convsize ).^2;
            %             s2 = conv2a( abs(Yj00), winfilt, convsize ).^2;
            %             S2 = s2;
            %             before 1* used, by mistake, better results

            %             compute GAUS variance
            %             s2 = conv2a( abs(Yj00).^2, winfilt, convsize );
            s2 = S2b(:,:,imi); % ML estimate

            [A, s2, Nsig2, T1] = DenoiseGGD(Nsig,Yj001,Yjp11,s2); % Ai shrinkage function
            %            S2 = s2; % using shrinkage s2
            S2g = s2;

%% shrinkage
            Yj001 = Yj001 .* A; % apply shrinkage to coeffs
            %             Yj00(~A) = 0; % only apply zeroing

%%% ==============================================================

            %             S2 = S2 .* A.^2;  % consequently, update saliency with shrinkage
            %             MA = MA .* prod(A,3); % update matching with shrinkage
            %             MA = (2 * ma * prod(A,3) + eps) ./ (sum(S2.* A,3) + eps); % abs(MA)?
            %             MA = max(min(MA,1),-1);
            %             S2 = s2; % keep the order if using shrinkage s2
            %             for imi = 1 : totalImgs
            %               S2(:,:,imi) = SaliencyGGD(Yj00(:,:,imi),  0   ,winfilt,0,b,FMETHOD,convsize);
            %             end
            %             MA = MatchMeasureGGD(Yj00, 0 ,S2,winfilt,b,FMETHOD,sign_ma,convsize); % Yjp1
            %             S2 = S3;
            % %             but MA,S2 ~ sqrt(abs(Yj00).^2 + abs(Yjp1).^2)
            %             MA = MA .* prod(A,3); % update matching with shrinkage
            %             MA = 1;

          end

          %         else % if denois_fuse
          %           T1 = 0;
        end

%%% COMPUTE MATCH 2 MA = MatchMeasureGGD(Yj00,Yjp1,S2,winfilt,b,FMETHOD,sign_ma,convsize);

%% CONTRAST ENHANCEMENT
        % s2 = conv2a( abs(Yj00).^2, winfilt, convsize );
        s2 = S2g;
%         Yj00_1 = Yj00; % store non-enhanced coeffs for visualisation

        switch upper(ENHMETHOD)
          
          case {'HEQ'}
          
          Amax = max(Yj001(:)); 
          Yj001 = Amax*histeq(Yj001/Amax);
            
          case {'VBE'}%,'AHE'}

            % detail enhancement
%             s2 = Yj00.^2;  % using coeff instead of var
            Adet = exp( -s2 ./ max(s2(:)) ) - exp(-1) + 1;
            AA(:,:,dir1) = Adet; % store enh. func. for use with approx

            A = Adet.^2;  % experimental (different than in the paper)
%             A = Adet; 

            % limit boost of the high value coeffs - optional
            clip_ind = abs(Yj001) .* A >= .75*max(abs(Yj001(:)));   
            A( clip_ind ) = 1;
            
            % clip at the max coeff
            clip_ind = abs(Yj001) .* A >= max(abs(Yj001(:)));    % A( clip_ind ) = 1;
            A( clip_ind ) = max(abs(Yj001(:)))./abs(Yj001(clip_ind));

            % apply enhacement
            Yj002 = A .* Yj002;

            %%
          case {'SH2','AHE'} % quite extreme contrast enhacement

            % detail enhancement
            Adet = exp( -s2 ./ max(s2(:)) ) - exp(-1) + 1;
            AA(:,:,dir1) = Adet; % store enh. func. for use with approx

            % approx enhancement
            Yj00a = L1{imi}{scale}{re_im}{dir}; % using coeff instead of var
            s2a = Yj00a.^2; %  s2a = conv2a( abs(Yj00a).^2, winfilt, convsize );
            Aapp = exp( -s2a ./ max(s2a(:)) ) - exp(-1) + 1;

            % combined
            % sqrt results in less enhancement but also less saturation
            % sqrt can be removed if saturation is not a concern
%            A = sqrt(Aapp .* Adet); % very strong contrast enh. if this is used
            A = (Aapp .* Adet); % even stronger contrast enh. if this is used
            % this might be equivalent to Adet^2            

            % limit boost of the high value coeffs - optional
            clip_ind = abs(Yj002) .* A >= .75*max(abs(Yj002(:)));   
            A( clip_ind ) = 1;

            % clip at the max coeffs
            clip_ind = abs(Yj002) .* A >= max(abs(Yj002(:)));    % A( clip_ind ) = 1;
            A( clip_ind ) = max(abs(Yj002(:)))./abs(Yj002(clip_ind));

            % apply enhacement
            Yj002 = A .* Yj002;

          case 'MPL' % just multiplying subbands - look below
          
          case 'GAG'
            A = 1;

          case 'LCE' 
            epps = .01; gam = .85;
            s = sqrt(conv2a( abs(Yj00).^2, winfilt, convsize ));
            A = ( (1-epps)*s/max(s(:)) + epps ).^(gam-1);
            Yj00 = A .* Yj00;
        end

        % other enhancement options
        switch upper(ENHMETHOD)
%           case {'VBE' 'SH2' 'AHE' 'LCE'}
%             enh_meth1 = 'NIL';
          case 'GAG'
            enh_meth1  = 'GAG';
             %enh_meth1  = 'NIL';
            Y_fused = WaveletContrastEnh(enh_meth1,Yj00,T1,A); % max(T1,.05)
          otherwise
            Y_fused = Yj002;
%             enh_meth1 = 'NIL'; %T1 = 0; A = 1;
        end


%% FUSION -----------------------------------------------
        % pass the fused coefficents on
        W2 = SetCplxCoeffs(W2,Y_fused,re_im,scale,dir,dir1,dwt_opt,RI);

      end % for re_im
    end % for dir
  end % for dir1
end % for scale


%% APPROXIMATION FUSION ------------------------------------------------
for dir = 1 : 2
  for re_im = 1 : RI

    %     if strcmp(dwt_opt,'real')
    %       W{LEVELS+1}{dir} = .5*(W12{1}{end}{dir} +  ...
    %         imadjust(W12{totalImgs}{end}{dir},[],[],.5));
    %     else
    %       W{LEVELS+1}{re_im}{dir} = .5*(W12{1}{end}{re_im}{dir} + ...
    %         imadjust(W12{totalImgs}{end}{re_im}{dir},[],[],.5));
    %     end

    if strcmp(dwt_opt,'real')
      W2{LEVELS+1}{dir} = mean(W12{1}{end}{dir} + W12{totalImgs}{end}{dir});
    else
      W2{LEVELS+1}{re_im}{dir} = zeros(size(W12{1}{end}{re_im}{dir}));
      for imi = 1 : totalImgs

        Yj00a1 = W121{imi}{end}{re_im}{dir};
        Yj00a2 = W122{imi}{end}{re_im}{dir};

        switch upper(ENHMETHOD)

%%
          case 'VBE'
            
            % approx enhancement
            s2a = Yj00a.^2;
            Aapp = exp( -s2a ./ max(s2a(:)) ) - exp(-1) + 1;
            
            A = Aapp;

            % clip at the max coeffs
            clip_ind = abs(Yj00a) .* A >= max(abs(Yj00a(:)));    % A( clip_ind ) = 1;
            A( clip_ind ) = max(abs(Yj00a(:)))./abs(Yj00a(clip_ind));

            % apply enhancement
            Yj00a  = WaveletContrastEnh('LGA',Yj00a,0,A);
            W{LEVELS+1}{re_im}{dir} = Yj00a; %nrm_Yj00a*Yj00a/norm(Yj00a(:));


          case 'SH2'

            % approx enhancement
            s2a = Yj00a1.^2;
            Aapp = exp( -s2a ./ max(s2a(:)) ) - exp(-1) + 1;

            % detail enhancement
            Adet = mean(AA, 3);
            
            % combined enhacement  sqrt removed
            A = (Adet .* Aapp);

            % limit boost of the high value coeffs - no difference

            % clip at the max coeffs
            clip_ind = abs(Yj00a2) .* A >= max(abs(Yj00a2(:)));    % A( clip_ind ) = 1;
            A( clip_ind ) = max(abs(Yj00a2(:)))./abs(Yj00a2(clip_ind));

            % apply enhancement             %             nrm_Yj00a = norm(Yj00a(:));
            Yj00a2  = WaveletContrastEnh('LGA',Yj00a2,0,A);

            W2{LEVELS+1}{re_im}{dir} = Yj00a2; %nrm_Yj00a*Yj00a/norm(Yj00a(:));

%%
          case 'AHE'
     
%             % approx enhancement
%             s2a = Yj00a.^2;
%             Aapp = exp( -s2a ./ max(s2a(:)) ) - exp(-1) + 1;
%             A = Aapp;
%             % clip at the max coeffs
%             clip_ind = abs(Yj00a) .* A >= max(abs(Yj00a(:)));    % A( clip_ind ) = 1;
%             A( clip_ind ) = max(abs(Yj00a(:)))./abs(Yj00a(clip_ind));
%             % apply enhancement           
%             Yj00a = A .* Yj00a;
            
            Yj00a0 = Yj00a;
            Amax = max(Yj00a0(:)); % simple HE Yj00a = Amax*histeq(Yj00a0/Amax);
            
            Amean = mean(Yj00a0(:));
            
%             round(100*Amax/Amean)/10000
            
            Yj00a = Amax*adapthisteq(Yj00a0/Amax, ...
              'NumTiles', max( [ round(size(Yj00a)/8); 2 2 ] ) ...
              ,...
              'NBins',64 ...
              , ... round(Amax/Amean)/100/2
              'ClipLimit', .03 ... .001 gives higher contrast (deeper)  , ...   'Distribution','' ... rayleigh exponential uniform
              ); % rayleigh boosts a lot but not very natural, exponential good but a bit dark, uniform good
            %, ...  'Alpha',.5 ... the smaller the brighter .8 s = 1/(sqrt(2)*mean(abs(Yj00a(:)/Amax)))
            
            % optional
            Yj00a  = WaveletContrastEnh('LGA',Yj00a,0,1);
            
            W{LEVELS+1}{re_im}{dir} = Yj00a;

%             [xsort, isort] = sort(Yj00a0(:));            
            
          case 'MPL'

            % approx enhancement
            s2a = Yj00a.^2;
            Aapp = exp( -s2a ./ max(s2a(:)) ) - exp(-1) + 1;

            % detail enhancement
            Adet = mean(AA, 3);
            
            % combined enhacement  sqrt removed
            A = (Adet .* Aapp);

            clip_ind = abs(Yj00a) .* A >= max(abs(Yj00a(:)));    % A( clip_ind ) = 1;
            A( clip_ind ) = max(abs(Yj00a(:)))./abs(Yj00a(clip_ind));

%             A = .8*mean2(A) * ones(size(A));
            A = .6*mean2(A) * ones(size(A));

            % apply enhancement             %             nrm_Yj00a = norm(Yj00a(:));
            Yj00a  = WaveletContrastEnh('LGA',Yj00a,0,A);

            W{LEVELS+1}{re_im}{dir} = Yj00a; %nrm_Yj00a*Yj00a/norm(Yj00a(:));

          
          otherwise
            W2{LEVELS+1}{re_im}{dir} = Yj00a2;
        end
      end
    end

    % Y_fusedA = W{LEVELS+1}{re_im}{dir};
    % pokaz(Yj00,Yjp1,Y_fusedA(:,:,[1 1]),[],[],[]),

%% VIS --------------------------------------------------

    if 0 %all( [re_im dir dir1 scale]==[2 1 1 2] )

%% show what's happening to the coeffs

      figure('Position',[20 10 1240 920])
      subplot(221), implotnf(Yj00a), title('approx. coeffs')
      subplot(222), implotnf(A), title('enh func')
      subplot(223), implotnf(W{LEVELS+1}{re_im}{dir}./(Yj00a+eps)), title('enh func')
      subplot(224), implotnf(W{LEVELS+1}{re_im}{dir}), title('enhanced coeffs')
      figure('Position',[20 10 1240 920])
      subplot(231), implotnf(abs(Yj00_1)), title('abs(coeffs)')
      subplot(232), implotnf(s2), title('local s^2')
      subplot(233), implotnf(abs(Y_fused)), title('abs(enhanced coeffs)')
      subplot(234), implotnf(Adet), title('detail enhancement function')
      subplot(235), implotnf(Aapp), title('approximation enhancement function')
      subplot(236), implotnf(A), title('combined enhancement function')
    end

  end
end

% Noise variance estimation using robust median estimator..
tmp = W2{1}{1}{1}{3}; % Selesnic uses {1}{1}{1}{1} but this seems to work better
Nsigfused = median(abs(tmp(:)))/0.6745; % Nsig2 = Nsig^2; % (1:2)


%% INVERSE TRANSFORM --------------------------------------------------
if strcmp(dwt_opt,'real')
  y = idualtree2D(W2, LEVELS, Fsf, sf);
else
  W2 = unnormcoef(W2,LEVELS,nor);
  y = icplxdual2D(W2, LEVELS, Fsf, sf);
end

% extract the image
y = GetCMat(y,[size(X1,1) size(X1,2)]);

% END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% alcplxdual2D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [w, L] = alcplxdual2D(x, LEVELS, Faf, af)

% Dual-Tree Complex 2D Discrete Wavelet Transform
%
% USAGE:
%   w = cplxdual2D(x, LEVELS, Faf, af)
% INPUT:
%   x - 2-D array
%   LEVELS - number of stages
%   Faf{i}: first stage filters for tree i
%   af{i}:  filters for remaining stages on tree i
% OUTPUT:
%   w{j}{i}{d1}{d2} - wavelet coefficients
%       j = 1..LEVELS (scale)
%       i = 1 (real part); i = 2 (imag part)
%       d1 = 1,2; d2 = 1,2,3 (orientations)
%   w{LEVELS+1}{m}{n} - lowpass coefficients
%       d1 = 1,2; d2 = 1,2
% EXAMPLE:
%   x = rand(256);
%   LEVELS = 5;
%   [Faf, Fsf] = FSfarras;
%   [af, sf] = dualfilt1;
%   w = cplxdual2D(x, LEVELS, Faf, af);
%   y = icplxdual2D(w, LEVELS, Fsf, sf);
%   err = x - y;
%   max(max(abs(err)))
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

% modified by Artur Loza to get lowpass coeffs for j = 2 : LEVELS, 24/8/2006

% normalization
x = x/2;

for m = 1:2
  for n = 1:2
    j = 1;
    [lo w{j}{m}{n}] = afb2D(x, Faf{m}, Faf{n});
    L{j}{m}{n} = lo; % modified here
    for j = 2:LEVELS
      [lo w{j}{m}{n}] = afb2D(lo, af{m}, af{n});
      L{j}{m}{n} = lo; % modified here
    end
    w{LEVELS+1}{m}{n} = lo;
  end
end

for j = 1:LEVELS
  for m = 1:3
    [w{j}{1}{1}{m} w{j}{2}{2}{m}] = pm(w{j}{1}{1}{m},w{j}{2}{2}{m});
    [w{j}{1}{2}{m} w{j}{2}{1}{m}] = pm(w{j}{1}{2}{m},w{j}{2}{1}{m});
  end
end

%% aldualtree2D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [w, L] = aldualtree2D(x, LEVELS, Faf, af)

% 2D Dual-Tree Discrete Wavelet Transform
%
% USAGE:
%   w = dualtree2D(x, LEVELS, Faf, af)
% INPUT:
%   x - M by N array
%   LEVELS - number of stages
%   Faf - first stage filters
%   af - filters for remaining stages
% OUPUT:
%   w{j}{d1}{d2} - DWT coefficients
%       j = 1..LEVELS, k = 1..2, d = 1..3
%   w{LEVELS+1}{k} - lowpass coefficients
%       k = 1..2
% % EXAMPLE:
%   x = rand(256);
%   LEVELS = 3;
%   [Faf, Fsf] = FSfarras;
%   [af, sf] = dualfilt1;
%   w = dualtree2D(x, LEVELS, Faf, af);
%   y = idualtree2D(w, LEVELS, Fsf, sf);
%   err = x - y;
%   max(max(abs(err)))
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

% modified by Artur Loza to get lowpass coeffs for j = 1 : LEVELS, 08-Dec-2006

% normalization
x = x/sqrt(2);

% Tree 1
[lo w{1}{1}] = afb2D(x, Faf{1});      % stage 1
for j = 2:LEVELS
  [lo w{j}{1}] = afb2D(lo, af{1});  % remaining stages
  L{j}{1} = lo; % modified here
end
w{LEVELS+1}{1} = lo;                       % lowpass subband

% Tree 2
[lo w{1}{2}] = afb2D(x, Faf{2});      % stage 1
for j = 2:LEVELS
  [lo w{j}{2}] = afb2D(lo, af{2});  % remaining stages
  L{j}{2} = lo; % modified here
end
w{LEVELS+1}{2} = lo;                       % lowpass subband

% sum and difference
for j = 1:LEVELS
  for m = 1:3
    A = w{j}{1}{m};
    B = w{j}{2}{m};
    w{j}{1}{m} = (A+B)/sqrt(2);
    w{j}{2}{m} = (A-B)/sqrt(2);
  end
end

%% DTree2D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = DTree2D(x,J0,J1,Faf,af,nor,opt)

% a shortcut for wavelet analysis
x = symextend0(x,2^(J1-1)); % symmetric extension

if opt == 'real'
  [W L] = aldualtree2D(x, J1, Faf, af);
elseif opt == 'cplx'
  [W L] = alcplxdual2D(x, J1, Faf, af);
  W = normcoef(W,J1,nor);
end

if J1 == J0+1
  W{end} = L{J0}; % replace approx LEVELS+1 with approx LEVELS, keep detail LEVELS+1
end


%% DTree2Djp1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W, L] = DTree2Djp1(x,LEVELS,Faf,af,nor,opt)

% a shortcut for wavelet analysis
% we want to have a decompositon for LEVELS levels
% but with LEVELS+1 details available:
% W{1:LEVELS} details 1..LEVELS
% W{LEVELS+1} details LEVELS+1
% W{LEVELS+2} approx  LEVELS

J1 = LEVELS+1;

x = symextend0(x,2^(J1-1)); % symmetric extension

if opt == 'real'
  [W L] = aldualtree2D(x, J1, Faf, af);
elseif opt == 'cplx'
  [W L] = alcplxdual2D(x, J1, Faf, af);
  W = normcoef(W,J1,nor);
end

W{J1+1} = L{LEVELS}; % replace approx LEVELS+1 with approx LEVELS, keep detail LEVELS+1


%% SaliencySAS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [S2 B] = SaliencySAS(Yj00,winfilt,b_est,B)
function [S2, B] = SaliencySAS(Yj00,Yjp1,winfilt,b_est,B)

Yj00 = abs(Yj00); % remove abs() in following equations
Yjp1 = abs(Yjp1);

% disabled
% Yjp1 = 0;
%
% % combine child and parent coeffs
% if any(Yjp1 > 0)
%   Yj00 = Yj00 + Yjp1;
% end

% compute saliency (dispersion of SaS components)

% k1 = conv2a(log(Yj00+eps),winfilt,'same');
k1 = conv2a(log(Yj00+Yjp1+eps),winfilt,'same');
% should be /2 b/c 2x more coeffs but the window is 1/N
% but it doesn't matter when S2y compared S2x

% k1 = conv2a(log(abs(Yj00+eps)),winfilt,'same'); % abs(Y+eps) may still be =0, for Y=-eps!

if b_est % estimate SaS alpha (or 'b') - globally
  k1a = mean2(  log(Yj00+eps) ); % =0, for Y=-eps!
  k2a = mean2( (log(Yj00+eps)-k1a).^2 );
  B = min(1./sqrt(6*k2a/pi^2-0.5),2);
  % is it correct?? or shall I take square?? (also below): may not matter b/c no influence on MA
end

S2 = exp(B.*k1+psi(1)*(1-B)).^(1./B); % FLOM's gamma


%% ShapeParSAS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = ShapeParSAS(Yj00,b_est,B)

if ~b_est, return, end

k1a = mean2(  log(abs(Yj00+eps)) );
k2a = mean2( (log(abs(Yj00+eps))-k1a).^2 );
B = min(1./sqrt(6*k2a/pi^2-0.5),2);
% is it correct?? or shall I take square?? (also below): may not matter b/c no influence on MA


%% GGDcoeff%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function a = GGDcoeff(FMETHOD,b)

% a shortcut for coefficient resulting from definition or ML of GGD when
% estimating variance


% why was it commented before???
switch FMETHOD
  case {'LAP' 'ULS'}, a = 2; % b = 1, univariate
  case {'BLA' 'BLS'}, a = 3/4; % Selesnic bivariate model (sqrt(3)/2)^2
  case 'GGD', a = gamma(3/b)/gamma(1/b) * b^(2/b); % GGD
  otherwise,
    a = 1;
end


%% SaliencyGGD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S2, B] = SaliencyGGD(Yj00,Yjp1,winfilt,b_est,B,FMETHOD,convsize)

% compute saliency (dispersion of GGD components)

if b_est % estimate GLOBALLY
  B = laplacpar3(Yj00);
end

a = GGDcoeff(FMETHOD,B); % normalising constant

% saliency = mean(|X|^b) (MLE)                 Yj00(:,:,end) = Yj00(:,:,1);
% R = Yj00, WA & parent not used
% or abs(Yj00) if sign not used
R = sqrt(abs(Yj00).^2 + abs(Yjp1).^2);
S2 = a * conv2a( R.^B, winfilt, convsize ).^(2/B); % variance

% variance = g(3/b)/g(1/b) * s^2/b
%           S2(:,:,imi) = gamma(3./B{imi})./gamma(1./B{imi}) .* ( B{imi}.*S(:,:,imi) ).^(2./B{imi});
% S2(:,:,imi) = s.^2;           %           S2(:,:,imi) = max(S2(:,:,imi)-Nsig(imi).^2,eps);

%% ShapeParGGD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = ShapeParGGD(Yj00,b_est, B)


if ~b_est, return, end

B = laplacpar3(Yj00);


%% MatchMeasureSAS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MA = MatchMeasureSAS(Yj00,Yjp1,winfilt,B)

x1 = Yj00(:,:,1);
x2 = Yj00(:,:,end);

% x1 = abs(Yj00(:,:,1));
% x2 = abs(Yj00(:,:,end));

x1p = Yjp1(:,:,1);
x2p = Yjp1(:,:,end);

if length(B) > 1
  p = max(1,min(B{:}));
else
  p = B;
end

% p = 2;

% symmetric coefficient of covariation
Nom1 = conv2a(x1.*sign(x2).*abs(x2).^(p-1),winfilt,'same');
Nom2 = conv2a(x2.*sign(x1).*abs(x1).^(p-1),winfilt,'same');

Nom1p = conv2a(x1p.*sign(x2p).*abs(x2p).^(p-1),winfilt,'same');
Nom2p = conv2a(x2p.*sign(x1p).*abs(x1p).^(p-1),winfilt,'same');

Den1 = conv2a(abs(x2).^p,winfilt,'same');
Den2 = conv2a(abs(x1).^p,winfilt,'same');

Den1p = conv2a(abs(x2p).^p,winfilt,'same');
Den2p = conv2a(abs(x1p).^p,winfilt,'same');

MA = ( (Nom1+Nom1p) .* (Nom2+Nom2p)+eps ) ./ ...
  ( (Den1+Den1p) .* (Den2+Den2p) +eps ); % all terms should be /2, but here makes no difference
%
%MA = (Nom1.*Nom2+eps)./(Den1.*Den2+eps);
% was: MA = Nom1.*Nom2./(Den1.*Den2+eps);

%% MatchMeasureGGD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MA,ma] = MatchMeasureGGD(Yj00,Yjp1,S2,winfilt,b,FMETHOD,sign_ma,convsize)

if strcmp(FMETHOD,'AVE'),
  MA = 1;
  return,
end

if length(b) > 1,
  b = mean([b{:}]); % sqrt(B{1}.*B{2}); % but this mean value is not used for S2
end

if sign_ma, % get 'possible' sign
  %   MA = conv2a(prod(Yj00,3)+prod(Yjp1,3), winfilt, convsize);
  %   signY = sign(MA); % what about 0?
  %  signY = sign(prod(Yj00,3)); % covariance, with sign (e.g. =0 for orthogonal)
  signY = sign(Yj00); % covariance, with sign (e.g. =0 for orthogonal)
else
  signY = 1;
end

totalImgs = size(Yj00,3);

R = signY .* sqrt(abs(Yj00).^2 + abs(Yjp1).^2); % = Yj00, WA & parent not used
% or abs(Yj00) if sign not used

if strcmp(FMETHOD,'WAN')
  MA = totalImgs/2 * conv2a( prod(R,3).^(2/totalImgs), winfilt, convsize ); % b = 2
else
  % (1) this is generalised cov
  MA = conv2a( prod(R,3).^(b/2), winfilt, convsize ).^(2/b);
end

% (2) this is normal cov w/o or w/sign
% MA = conv2a( prod(R,3), winfilt, convsize );

% (3) this normal covariance (univariate)
%MA = conv2a( prod(Yj00,3), winfilt, convsize );

a = GGDcoeff(FMETHOD,b); % normalising constant
MA = a * MA; % b/c S2 = a * X;
ma = MA;

MA = (2 * MA + eps) ./ (sum(S2,3) + eps); % abs(MA)?
%         [min(MA(:)) max(MA(:))]

MA = max(min(MA,1),-1); % threshold the matching measure % w = w / max(w(:));

%         if denois_fuse == 1 % use signal variance
%           MA = MA-(parent+1)*prod(Nsig,3);
%           MA = MA .* (MA > 0);
%         end


%             MA = MA*a*winfilt(1)*(1+0);
%           MA = gamma(3./b)./gamma(1./b) .* ( b.*MA ).^(2./b);

% here MA is computed as for any B
%           b = sqrt(B{1}.*B{2}); % mean([B{1} B{2}]);
%         BM = mean([B{:}]);
%         signY = sign(Yj00(:,:,1) .* Yj00(:,:,2));
%         MA = LocalSpatialStats2(signY.*abs(Yj00(:,:,1) .* Yj00(:,:,2)).^(BM/2),0,0,Np,1); % abs
%         MA = gamma(3./BM)./gamma(1./BM) .* ( BM.*MA
%         ).^(2./BM);


%% DenoiseGGD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, s2, Nsig2, T] = DenoiseGGD(Nsig,Yj00,Yjp1,S2)

% , Ai

% compute denosing parameters
Nsig2 = repmat(Nsig.^2,[size(S2,1) size(S2,2) 1]);

% get the signal variance
s2 = S2-Nsig2;
s2 = s2 .* (s2 > 0); % else A = ones(size(S2));

% shrinkage function
R = sqrt(abs(Yj00).^2 + abs(Yjp1).^2); % =|Yj00| if parent not used
T = sqrt(3)*Nsig2 ./ sqrt(s2+eps); % threshold
A = (R - T)./(R+eps);
A = A .* (A > 0); % this is shrinkage function

% % R1 = sqrt(abs(Yj00).^2); % =|Yj00| if parent not used
% Ai = (R + T)./(R+eps);
% Ai = Ai .* (A > 0); % this is shrinkage function


%% laplacpar3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [p,s]=laplacpar3(x,p)

% Log estimation of the parameters 's' and 'p' of a generalized Laplacian (gL)

% by Artur Loza 23/06/2006

Y = log( abs(x(:)) + eps);
%mu_y = mean( Y );
N = length(Y);

if nargin == 1 || isempty(p) || p == 0

  ord = 2;

  % Define the function to be inverted
  p = .1 :.001: 2;

  var_p = psi( ord-1 , 1./p ) ./ p.^ord;

  % make it converge to 0 (experi-mental)
  if N <= 64
    warning('small sample size var = var - 1')
    var_p = var_p - 1;
  end

  % compute the inverse of var_x in the input point, thus resulting
  % the second parameter, p.
  var_y = var( Y ); % or mean( (Y - mu_y).^ord ); %

  p = interp1(var_p, p, var_y, 'spline');

end

if p < 0, warning('p negative, setting to 1'), p = 1; end

if nargout > 1
  mu_y = mean( Y );
  s = exp( mu_y - psi(1./p)/p );
end

% p = .001 :.001: 2.5; no difference
%k1 = mean( Y ); k2 = mean( (Y-k1).^2 );var_y = k2; % + eps var( Y )
% var_x = var( Y , 1); % + eps


%% CauchynoisyecfMomFit%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gama=CauchynoisyecfMomFit(data,datap,Nsig,winfilt,K,parent)%,lev_noise

% Moments estimate of Cauchy parameter
% K is a vector containing the two points where the ECF is to be estimated

%[t_m,fi_m]=sasecf(K,data); fi_m=abs(fi_m);
%gama=(-log(fi_m))/t_m;

% if all(datap==0)
%   data=data(:);
% else
%   data = [data(:); datap(:)]; % use parent
% end
%
% [f,FNs] = sasecf(K,data);
% PNs = abs(FNs);
%
% gama = mean((-log(PNs.^2)-sigma^2*f.^2)./(2*f));
%
% gama = max(gama,eps);

% slightly different results when symextend0/floor(size(winfilt,1)/2)) and
% 'same' is used --> different boundary distortions and the matrix is
% shifted by 1 pixel!

convsize = 'same';%'valid';
% this part is 'sasecf'
f(1,1,:) = (pi*K(:)/256)';

f0 = repmat(f(1,1,:),[size(data) 1]);
% data = symextend(data,(size(winfilt,1)+1)/2);
f = f0;%repmat(f(1,1,:),[size(data) 1]);%
data = data(:,:,ones(1,size(f,3)));

exp_jfd_c = exp(i * data .* f);

if parent
  %   datap = symextend(datap,(size(winfilt,1)+1)/2);
  datap = datap(:,:,ones(1,size(f,3)));
  exp_jfd_p = exp(i * datap .* f);
  exp_jfd = (exp_jfd_c + exp_jfd_p)/2;
else
  exp_jfd = exp_jfd_c;
end

% this part is 'CauchynoisyecfMomFit'
FNs = conv2a(exp_jfd,winfilt,convsize); % FNs = mean(exp(j*f*data'),2);
PNs = abs(FNs);
gama = mean((-log(PNs.^2)-Nsig^2*f0.^2)./(2*f0),3);

if all(ismember(1:5,K))
  thr  = eps;
elseif all(ismember(2:4,K))
  thr = 1;
else
  error('K?')
end

gama = max(gama,thr); % K=1:5; --> disper = max(gama,eps);

%% biCauchyshrink%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [w1] = biCauchyshrink(y1,y2,disp,Nsig,parent) % bi
% Bivariate Cauchy Shrinkage Function
% Usage :
%      [w1] = biCauchyshrink(y1,y2,disp,Nsig)
% INPUT :
%      y1 - a noisy coefficient value
%      y2 - the corresponding parent value
%      disp  - signal dispersion
%      Nsig - noise standard deviation
% OUTPUT :
%      w1 - the denoised coefficient

% y1=y1(1,1)
% y2=y2(1,1)
% disp
% Nsig
% A = 1+y2^2/y1^2
% B = -y1*(1+y2^2/y1^2)
% C = disp.^2+3*Nsig^2
% D = -y1.*disp.^2

if parent
  p = (disp.^2+3*Nsig^2)./(1+y2.^2./y1.^2) - y1.^2./3;
  q = -2*y1.^3/27 + (y1./3).*(disp.^2+3*Nsig^2)./(1+y2.^2./y1.^2) - disp.^2.*y1./(1+y2.^2./y1.^2);
else
  p =  disp.^2 + 2*Nsig^2  - y1.^2/3;
  q = -2*y1.^3/27 - 2/3*disp.^2.*y1 + 2/3*Nsig^2*y1;
end

% p = (disp.^2+3*Nsig^2)./(1+y2.^2./y1.^2) - y1.^2./3;
% q = -2*y1.^3/27 + (y1/3).*(disp.^2+3*Nsig^2)./(1+y2.^2./y1.^2) - disp.^2*y1./(1+y2.^2./y1.^2);

DD = p.^3/27 + q.^2/4;
w1 = y1/3 + (abs(-q/2 + sqrt(DD)).^(1./3.)).*sign(-q/2 + sqrt(DD)) + ...
  (abs(-q/2 - sqrt(DD)).^(1./3.)).*sign(-q/2 - sqrt(DD));

% w1 = y1/3 + (-q/2 + (p.^3/27 + q.^2/4).^(1/2)).^(1/3) + (-q/2 - (p.^3/27 + q.^2/4).^(1/2)).^(1/3);
% w2 = y1/3 + Eps*(-q/2 + (p.^3/27 + q.^2/4).^(1/2)).^(1/3) + Eps^2*(-q/2 - (p.^3/27 + q.^2/4).^(1/2)).^(1/3);
% w3 = y1/3 + Eps^2*(-q/2 + (p.^3/27 + q.^2/4).^(1/2)).^(1/3) + Eps*(-q/2 -
% (p.^3/27 + q.^2/4).^(1/2)).^(1/3);

%% sasecf%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t,fi]=sasecf(K,data);
% Compute empirical characteristic function (ECF) of data
% K should be a vector of points at which to compute the ECF

data=data(:);
%k = [1:K]';
K = K(:);
t = pi*K/256;
%t = pi*K/50;

%N = length(data);
%for i=1:length(t)
%    fi(i) = (1/N)*sum(exp(j*t(i)*data));
%end

fi = mean(exp(j*t*data'),2);



%% GetCMat%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y, first, last] = GetCMat(X,s)

sx(1)  = size(X,1);
sx(2)  = size(X,2);

if any(isinf(s))
  dim = find(s == inf);
  s(dim) = sx(dim);
end

% if any(s > sx)
%   dim = find(s > sx);
%   s(dim) = sx(dim);
% end

s = min([s; sx]);

ds = (sx-s)/2;

first = 1  + floor(ds);
last  = sx - ceil (ds);

Y = X( first(1):last(1) , first(2):last(2),:);

%% conv2a%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = conv2a(x, winfilt, convsize)

% a shortcut to conv2: conv2 is skipped if filtering window is larger
% than image to be filtered (useful when a global estimate is computed)

% if all(size(winfilt) < [size(x,1) size(x,2)]) % convolve

y = convn(x, winfilt, convsize);

% else % multiply and convolve the central parts
%
%   winfilt = GetCMat(winfilt,size(x));
%
%   y = sum2(x .* winfilt); % assuming the window is symmetrical W(n)=W(-n)
%
% end

%% Cauchyshrink%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [w1] = Cauchyshrink(y,disp,Nsig)  %,a,b,c,d
% % Cauchy Shrinkage Function
% % Usage :
% %      [w1] = biCauchyshrink(y,disp,Nsig)
% % INPUT :
% %      y - a noisy coefficient value
% %      disp  - signal dispersion
% %      Nsig - noise standard deviation
% % OUTPUT :
% %      w1 - the denoised coefficient
%
% %  a = 1;
% %  b = -y;
% %  c = disp^2+2*Nsig^2;
% %  d = -disp^2*y;
%
% % p  = c/a - b.*b/a/a/3.
% % q  = (2.*b.*b.*b/a/a/a - 9.*b.*c/a/a + 27.*d/a) / 27.
%
% p = disp.^2 + 2*Nsig^2 - y.^2/3;
% q = -2*y.^3/27 - 2/3*disp.^2.*y + 2/3*Nsig^2*y;
%
% DD = p.^3/27 + q.^2/4;
% w1 = y/3 + (abs(-q/2 + sqrt(DD)).^(1./3.)).*sign(-q/2 + sqrt(DD)) + ...
%   (abs(-q/2 - sqrt(DD)).^(1./3.)).*sign(-q/2 - sqrt(DD));
% % whos p q DD
% % pause
% % Eps = ((-1+i*sqrt(3))/2);
% %
%
% %         temp1 = -q/2. + sqrt(DD);
% %         temp2 = -q/2. - sqrt(DD);
% %         u = abs(temp1)^(1./3.);
% %         v = abs(temp2)^(1./3.);
% %         if (temp1 < 0.) u=-u; end
% %         if (temp2 < 0.) v=-v; end
% %         y1  = u + v;
% % temp1 = b/a/3.;
% %       w11 = y1-temp1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = symextend0(x,Nnum)

% y(:,1:Nnum) = fliplr(x(:,1:Nnum));
% use it with Nnum=floor(size(winfilt,1)/2)
y = [fliplr(x(:,1:Nnum)) x x(:,end:-1:end-Nnum+1)];
y = [flipud(y(1:Nnum,:)); y ;y(end:-1:end-Nnum+1,:)];

% AL: this version, unlike above, doesn't duplicate the first and last pixel
% use with Nnum=(size(winfilt,1)+1)/2
% y = [x(:,Nnum:-1:2)    x  x(:,end-1:-1:end-Nnum+1)];
% y = [y(  Nnum:-1:2,:); y; y(  end-1:-1:end-Nnum+1,:)];
% y = [fliplr(x(:,2:Nnum)) x x(:,end-1:-1:end-Nnum+1)];
% y = [flipud(y(2:Nnum,:)); y ;y(end-1:-1:end-Nnum+1,:)];

% AL: and this version, as above doesn't duplicate but extends by Nnum
% use it with Nnum=floor(size(winfilt,1)/2)
% y = [x(:,Nnum+1:-1:2)    x  x(:,end-1:-1:end-Nnum)];
% y = [y(  Nnum+1:-1:2,:); y; y(  end-1:-1:end-Nnum,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = symextend(x,Nnum)

% y(:,1:Nnum) = fliplr(x(:,1:Nnum));
% y = [fliplr(x(:,1:Nnum)) x x(:,end:-1:end-Nnum+1)];
% y = [flipud(y(1:Nnum,:)); y ;y(end:-1:end-Nnum+1,:)];

% AL: this version, unlike above, doesn't duplicate the first and last pixel
% use with Nnum=(size(winfilt,1)+1)/2
y = [x(:,Nnum:-1:2)    x  x(:,end-1:-1:end-Nnum+1)];
y = [y(  Nnum:-1:2,:); y; y(  end-1:-1:end-Nnum+1,:)];
% y = [fliplr(x(:,2:Nnum)) x x(:,end-1:-1:end-Nnum+1)];
% y = [flipud(y(2:Nnum,:)); y ;y(end-1:-1:end-Nnum+1,:)];

% AL: and this version, as above doesn't duplicate but extends by Nnum
% use it with Nnum=floor(size(winfilt,1)/2)
% y = [x(:,Nnum+1:-1:2)    x  x(:,end-1:-1:end-Nnum)];
% y = [y(  Nnum+1:-1:2,:); y; y(  end-1:-1:end-Nnum,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Yj00,Yjp1] = GetCplxCoeffs(W12,imi,re_im,scale,dir,dir1,dwt_opt,RI,parent)

if strcmp(dwt_opt,'real')
  Yj00 = W12{imi}{scale}{dir}{dir1};
else
  Yj00 = W12{imi}{scale}{1}{dir}{dir1} + i*W12{imi}{scale}{2}{dir}{dir1};
end

if parent
  if strcmp(dwt_opt,'real')
    Yjp1 = expand(W12{imi}{scale+1}{dir}{dir1});
  else
    Yjp1 = expand(W12{imi}{scale+1}{1}{dir}{dir1} + i*W12{imi}{scale+1}{2}{dir}{dir1});
  end
else
  Yjp1 = zeros(size(Yj00));
end

if RI == 2 % process real and imag separately
  if re_im==1, % real
    [Yj00,Yjp1] = deal(real(Yj00),real(Yjp1)); % ,Yjm1 ,real(Yjm1)
  else % imaginary
    [Yj00,Yjp1] = deal(imag(Yj00),imag(Yjp1)); % ,Yjm1 ,imag(Yjm1)
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = SetCplxCoeffs(W,Y_fused,re_im,scale,dir,dir1,dwt_opt,RI)

if strcmp(dwt_opt,'real') % real DWT
  W{scale}{dir}{dir1} = Y_fused;
else % complex DWT
  if RI == 2 % real, imag separately
    W{scale}{re_im}{dir}{dir1} = Y_fused;
  else
    W{scale}{1}{dir}{dir1} = real(Y_fused);
    W{scale}{2}{dir}{dir1} = imag(Y_fused);
  end
end

%           if strcmp(dwt_opt,'real')
%             Yj00(:,:,imi) = W12{imi}{scale}{dir}{dir1};
%           else
%             Yj00(:,:,imi) = W12{imi}{scale}{1}{dir}{dir1} + i*W12{imi}{scale}{2}{dir}{dir1};
%           end
%
%           if parent
%             if strcmp(dwt_opt,'real')
%               Yjp1(:,:,imi) = expand(W12{imi}{scale+1}{dir}{dir1});
%             else
%               Yjp1(:,:,imi) = expand(W12{imi}{scale+1}{1}{dir}{dir1} + i*W12{imi}{scale+1}{2}{dir}{dir1});
%             end
%           else
%             Yjp1(:,:,imi) = zeros(size(Yj00(:,:,imi)));
%           end
%
%           if RI == 2 % process real and imag separately
%             if re_im==1, % real
%               [Yj00(:,:,imi),Yjp1(:,:,imi)] = deal(real(Yj00(:,:,imi)),real(Yjp1(:,:,imi))); % ,Yjm1(:,:,imi) ,real(Yjm1(:,:,imi))
%             else % imaginary
%               [Yj00(:,:,imi),Yjp1(:,:,imi)] = deal(imag(Yj00(:,:,imi)),imag(Yjp1(:,:,imi))); % ,Yjm1(:,:,imi) ,imag(Yjm1(:,:,imi))
%             end
%           end

%               if parent
%                 Yj00(:,:,imi) = biCauchyshrink(Yj00(:,:,imi),Yjp1(:,:,imi),gamm,Nsig(1,1,imi));
%               else
%               end

%         % Noise variance estimation using robust median estimator..
%         if all([scale dir dir1 re_im] == 1) %  ...or --> scale ==  1 %
%           Nsigfused = median(abs(Y_fused(:)))/0.6745;
%         end

% Noise variance estimation using robust median estimator..
%            if all([scale dir dir1 re_im] == 1) % when complex amplitudes are used? Rayleigh?
%             tmp = Yj00(:,:,imi);
%             Nsig(1,1,imi) = median(abs(tmp(:)))/0.6745;
%           end

%  +             ~m_max.*(Yj00(:,:,1).*w_max + Yj00(:,:,2).*w_min) ); % old weighting S2>S1

%         end
%   m1 .* ((m2.*Y1_coef_imag.*(1-w1))+(m2.*Y2_coef_imag.*w1) + ...
%               (~m2.*Y2_coef_imag.*(1-w1))+(~m2.*Y1_coef_imag.*w1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pokaz(Yj00,Yjp1,S2,S2b,MA,MAb,FMETHOD)

thr = 1;

figure('Position',get(0,'ScreenSize')-[0 0 0 30]),

if ~isempty(Yj00)

  Yj00 = ( Yj00 .* (Yj00 >= thr) + (Yj00 < thr) ); % log

  figure('Position',get(0,'ScreenSize')-[0 0 0 30]),
  for k = 1 : 2
    subplot('position',SubPlCoord([2 4],[1 k]));
    imagesc(Yj00(:,:,k)), axis image, colormap gray
    title(['Y_j_0_0 ' NumAlign(mean2(abs(Yj00(:,:,k))),2)])
  end
  set(gca,'YColor',[1 1 0])
end

if ~isempty(Yjp1)

  Yjp1 = ( Yjp1 .* (Yjp1 >= thr) + (Yjp1 < thr) ); % log

  for k = 1 : 2
    subplot('position',SubPlCoord([2 4],[2 k]));
    imagesc(Yjp1(:,:,k)), axis image, colormap gray
    title(['Y_j_p_1 ' NumAlign(mean2(abs(Yjp1(:,:,k))),2)])
  end
  set(gca,'YColor',[1 1 0])
end

if ~isempty(S2)

  S2 = log( S2 .* (S2 >= thr) + (S2 < thr) );

  for k = 1 : 2
    subplot('position',SubPlCoord([2 4],[1 k+2]));
    imagesc(S2(:,:,k)), axis image, colormap gray
    title(['S_2 ' NumAlign(mean2(abs(S2(:,:,k))),2)])
  end
  set(gca,'YColor',[1 1 0])
end

if ~isempty(S2b)

  subplot(144)
  imagesc(S2b(:,:,1) > S2b(:,:,end)), axis image, colormap gray,
  title('\sigma_X > \sigma_Y')%,'FontSize',12,'FontName','Times')
  NoTickLabels

  S2b = log( S2b .* (S2b >= thr) + (S2b < thr) );

  for k = 1 : 2
    subplot(1,4,k);
    %    subplot('position',SubPlCoord([2 4],[2 k+2]));
    imagesc(S2b(:,:,k)), axis image, colormap gray
    title('\sigma_X')%,'FontSize',12,'FontName','Times')
    NoTickLabels
    if k == 1
      ylabel(FMETHOD,'FontName','Times')%,'FontSize'%,12,'FontName','Times')
    end
    %     set(gca,'FontSize',16)
    %     title(['S_2b ' NumAlign(mean2(abs(S2b(:,:,k))),2)])
  end
  %   set(gca,'YColor',[1 1 0])


end

if ~isempty(MA)
  %   figure('Position',get(0,'ScreenSize')-[0 0 0 30]),
  %   subplot(2,2,1)
  subplot(143)
  imagesc(MA), axis image, colormap gray, %title('MA')
  title('M_X_Y')%,'FontSize',12,'FontName','Times')
  %    set(gca,'FontSize',16)
  NoTickLabels


  %   subplot(2,2,3)
  %   surf(double(MA)), colormap gray, set(gca,'View',[23 28]),
  %   title(['MA ' NumAlign(mean2(abs(MA)),2)])
  %   axis([0 size(MA,2) 0 size(MA,1) -1 1])
end

if ~isempty(MAb)
  subplot(2,2,2)
  imagesc(MAb), axis image, colormap gray, title('MAb')

  subplot(2,2,4)
  surf(double(MAb)), colormap gray, set(gca,'View',[23 28]),  % [28 82]
  title(['MAb ' NumAlign(mean2(abs(MAb)),2)])
  axis([0 size(MAb,2) 0 size(MAb,1) -1 1])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Yj00 = LogContrastEnh(Yj00,T1)

Yj00max = max2(abs(Yj00));

Yj00 = ...
  (abs(Yj00) > T1) .* ... zero small coefficients
  sign(Yj00) .* ... keep the sign
  log(abs(Yj00)+1) .* ... log-boost the coefficients
  Yj00max/log(Yj00max+1) ; % scale to the original range
%   ...   imadjust(abs(Yj00),stretchlim(abs(Yj00),.05)); %,[.2 .6],[.4 .8], 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Yj00 = GAGContrastEnh(Yj00,T1)

% test
% Yj00 = -10:.01:10; x = Yj00;

Yj00max = max2(abs(Yj00));
Yj00 = Yj00/Yj00max;

% thresholds
% T1 = T1/Yj00max; % T2 = max(T1, .2); % T3 = max(T2, .8);
if nargin == 1
  T1 = .03; % "denoising" thr
end
T2 = .05; T3 = .4; b = .08; c = 10;

a = 1 / ( sigmoid(c*(1-b),1,0) - sigmoid(-c*(1+b),1,0) );

u1 = sign(Yj00) .* (abs(Yj00) - T2)/(T3 - T2);
u2 = a*(T3-T2).* (sigmoid(c*(u1-b),1,0)-sigmoid(-c*(u1+b),1,0));
u3 = sign(Yj00) .* T2 + u2;

Yj00(abs(Yj00) <=  T1) = 0; % zero small coefficients
%( (abs(Yj00) <= T2) + (abs(Yj00) >  T3) ) ... unchanged
Yj00( (abs(Yj00) >  T2) & (abs(Yj00) <= T3) ) =  ... % boost the coefficients in [T2,T3]
  u3((abs(Yj00) >  T2) & (abs(Yj00) <= T3));
Yj00 = Yj00 * Yj00max ; % scale to the original range
%  .* Yj00 .* ...  sign(Yj00) .* ... keep the sign (no need --> Yj00 used)

%  Yj00 = ...
%   (abs(Yj00) >  T1) .* ... zero small coefficients
%   ...
%   ( ...
%   (T2 .* ...
%   ( (abs(Yj00) >  T2) & (abs(Yj00) <= T3) ) + ... boost the coefficients in [T2,T3]
%   ( (abs(Yj00) <= T2) + (abs(Yj00) >  T3) ) ... unchanged
%   ) ...
%   ...
%   .* Yj00 .* ...  sign(Yj00) .* ... keep the sign (no need --> Yj00 used)
%   Yj00max ; % scale to the original range

% test
% figure, plot(x,Yj00)

% compute denoising threshold
%                 s2 = conv2a( abs(Yj00).^2, winfilt, convsize );  % mean2(abs(Yj00).^2); %
%   s2 = 2 * conv2a( abs(Yj00), winfilt, convsize ).^2;
%    [A, s2, nil2, T1, nil5] = DenoiseGGD(Nsig,Yj00,Yjp1,s2); % shrinkage function
% T1 = thselect(Yj00(:),'rigrsure'); % rigrsure heursure sqtwolog minimaxi
% Yapp = W12{imi}{end}{re_im}{dir};
% T1 = alwdcbm2(Yj00,[size(Yapp); scale scale; LEVELS LEVELS],2,6*numel(Yapp));%,alpha,m
% T1 = T1/10;
% gamm = CauchynoisyecfMomFit(Yj00,Yjp1,Nsig(1,1,imi),winfilt,2:4,parent);
% ww = biCauchyshrink(Yj00,Yjp1,gamm,Nsig(1,1,imi),parent);
% Ai = zeros(size(Yj00));  Ai(abs(ww)>0) = Yj00(abs(ww)>0) ./ ww(abs(ww)>0); Ai(abs(ww)==0) = 0;
%         gamma_ce = .5; epsi = .01; tc = max(S2b(:));
%         A = ( (1-epsi)* S2b./tc + epsi ).^(gamma_ce-1);
%         epsi=0;gamma_ce=.5;SS2b=0:.01:10;A=((1-epsi)*SS2b./max(SS2b(:))+epsi).^(gamma_ce-1);figure,plot(SS2b,A)
%         SS2b=0:.01:10;A=exp(-SS2b/(max(SS2b(:))))+1;figure,plot(SS2b,A),figure,plot(SS2b,A.*SS2b), axis([0 10 0 10]),axis square
%    y = LG1ContrastEnh(SS2b,0,0,1,A,1);figure,plot(SS2b,A),figure,plot(SS2b,y), axis square

%        Yapp = W12{1}{end}{re_im}{dir}; % S2b
% y_max = max(abs(Yj00(:)));

%% VIS --------------------------------------------------

% pokaz(Yj00,Yjp1,S2,S2b,MA,[]), %close
% pokaz([],[],[],S2b,MA,[],FMETHOD),

% return

% m_w(:,:,1) = m_weight;
%         M12(:,:,1) = m_max;
%
%         m_weight = MAb > th;
%         m_max = S2b(:,:,1) > S2b(:,:,end); %=diff(S2,1,3) <= 0; saliency(I1) > saliency(I2)
%         m_w(:,:,2) = m_weight;
%         M12(:,:,2) = m_max;
%         pokaz([],[],[],[],M12(:,:,1),M12(:,:,2))       % pokaz(Yj00,m_w,S2g,M12,[],[])
%         title('S_X > S_Y, u-b')
%         pokaz([],[],[],[],m_w(:,:,1),m_w(:,:,2))       % pokaz(Yj00,m_w,S2g,M12,[],[])
%         title('MA > thr, u-b')
%
%         MA = MAb;
%         S2 = S2b;
% %         MA = MA0;
% %         S2 = S20;

% Y_fused0(:,:,1) = Y_fused;
%
% m_weight = MA > th;
%         m_max = S2(:,:,1) > S2(:,:,end); %=diff(S2,1,3) <= 0; saliency(I1) > saliency(I2)
%         w_min = 0.5 * ( 1 - (1-MA)/(1-th) );
%         w_max = 1-w_min;
%         Y_fused0(:,:,2)  = ...
%           ~ m_weight .* ( m_max.*Yj00(:,:,1) + ~m_max.*Yj00(:,:,totalImgs) ) ... selection
%           + m_weight .* ( m_max.*(Yj00(:,:,  1).*w_max + Yj00(:,:,totalImgs).*w_min) ...  weighting S1>S2
%           +             ~m_max.*(Yj00(:,:,totalImgs).*w_max + Yj00(:,:,1   ).*w_min) ); % weighting S2>S1
%
%
% Y_fused0 = log( Y_fused0 .* (Y_fused0 >= 1) + (Y_fused0 < 1) );
% for k = 1 : 2
%   figure, imagesc(Y_fused0(:,:,k)), axis image, colormap gray
%   title(['fused uni:1, bi:2 --> ' n2s(k)])
% end


%% FUSION --------------------------------------------
%         %       if denois_fuse == 1 % use signal variance
%         %         Nsig2 = repmat(Nsig.^2,[size(S2,1) size(S2,2) 1]);
%         %         S2 = S2-Nsig2;    S2 = S2 .* (S2 > 0); % get the signal variance
%         %       end
%
% %         MA0 = MA;
% %         S20 = S2;
% %         MA = MAb;
% %         S2 = S2b;
%
%         m_weight = MA > th;
%         maxS2 = max(S2,[],3); % S2(:,:,1) > S2(:,:,end); %=diff(S2,1,3) <= 0; saliency(I1) > saliency(I2)
% %         m_max = S2(:,:,1) > S2(:,:,end);
%         m_max = (repmat(maxS2, [1 1 totalImgs]) == S2);
%
%         w_min = 0.5 * ( 1 - (1-MA)/(1-th) );
%         w_max = 1-w_min;
%         w_min = w_min ./ max(totalImgs-1,1); % for > 2 inputs, averages images with w_min
%
%         % !rewritten for n-inputs!
%         Y_fused  = ...
%         ~ m_weight .*  sum( m_max.*Yj00,3) + ... SELECTION
%           m_weight .* (sum( m_max.*Yj00,3).*w_max ... WEIGHTED AVERAGE
%           +            sum(~m_max.*Yj00,3).*w_min);
%
% %         Y_fused  = ...
% %           ~ m_weight .* ( m_max.*Yj00(:,:,1) + ~m_max.*Yj00(:,:,totalImgs) ) ... selection (rewritten for n-inputs)
% %           + m_weight .* ( m_max.*(Yj00(:,:,  1) .*w_max + Yj00(:,:,totalImgs).*w_min) ...  weighting S1>S2
% %           +              ~m_max.*(Yj00(:,:,totalImgs).*w_max + Yj00(:,:,1   ).*w_min) ); % weighting S2>S1
%

% case 'MPL' % just multiplying subbands - look below
% 
%             % detail enhancement
%             Adet = exp( -s2 ./ max(s2(:)) ) - exp(-1) + 1;
%             AA(:,:,dir1) = Adet; % store enh. func. for use with approx
% 
%             % approx enhancement
%             Yj00a = L{imi}{scale}{re_im}{dir}; % using coeff instead of var
%             s2a = Yj00a.^2; %  s2a = conv2a( abs(Yj00a).^2, winfilt, convsize );
%             Aapp = exp( -s2a ./ max(s2a(:)) ) - exp(-1) + 1;
% 
%              A = (Aapp .* Adet); % very strong contrast enh. if this is used
%              % this might be equivalent to Adet^2 ?
% 
%             % clip at the max coeff
%             clip_ind = abs(Yj00) .* A >= max(abs(Yj00(:)));    % A( clip_ind ) = 1;
%             A( clip_ind ) = max(abs(Yj00(:)))./abs(Yj00(clip_ind));
% 
% %             A = .8*mean2(A) * ones(size(A));
%             A = .6*mean2(A) * ones(size(A));
%             
%             % apply enhacement
%             Yj00 = A .* Yj00;
%             %             nrm_Yj00 = norm(Yj00(:));
%             %             Yj00 = nrm_Yj00*Yj00/norm(Yj00(:));
% 
%             %             max_Yj00 = max(Yj00(:));
%             %             Yj00 = A .* Yj00;
%             %             Yj00 = max_Yj00*Yj00/max(Yj00(:));
%           