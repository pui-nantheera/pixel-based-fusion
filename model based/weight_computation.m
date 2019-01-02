function [weights, m, gama] = weight_computation(windowsize,data1,data2)

% weight_computation - weight computered based on window using Cauchy probability
% D
% 
% weights = weight_computation(windowsize,data1,data2) returns the parameter
% weights estimated weights for data1 and data2 respectively

%   Copyright (c) Tao wan 09/07/2007. 
%   tao.wan@bris.ac.uk

%windowfilt = ones(1,windowsize)/windowsize;

%k1 = conv2(windowfilt,windowfilt,log(abs(data+eps)),'same');

%gama = exp(k1);

 [nl,nc] = size(data1);
 data1_ext = symextend(data1,(windowsize+1)/2);
 data2_ext = symextend(data2,(windowsize+1)/2);
 
 
 step=1;
 for i=1:step:nc
    for j=1:step:nl
         data1_window(:,:) = data1_ext(j:windowsize+j-1,i:windowsize+i-1)';
         data2_window(:,:) = data2_ext(j:windowsize+j-1,i:windowsize+i-1)';
         
         samples(1,:) = data1_window(:)';
         samples(2,:) = data2_window(:)';
         
         [weights(j,i,:),m(j,i),gama(j,i), n_iter, final_err] = Weight_PDFPrior_tao_cauchy2(samples,0.05, 0.001);
         clear data1_window data2_window samples;
    end
 end