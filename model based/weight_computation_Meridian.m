function [weights, m, gama] = weight_computation_Meridian(windowsize,data1,data2)

% weight_computation - weight computered based on window using meridian probability

% Modified from the code written by tao-wan on cauchy convolution (Context Enhancement Through Image Fusion: A Multiresolution Approach Based on Convolution of Cauchy Distributions)
% for meridian covolution.
%   Modified by Mayank Agrawal 19/03/2010. 
%   ma6879@bris.ac.uk

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
         
         [weights(j,i,:),m(j,i),gama(j,i), n_iter, final_err] = Weight_Meridian(samples,0.05, 0.001);
         clear data1_window data2_window samples;
     end;
 end;
end

            
         
         
