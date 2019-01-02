function [w, p1, p2, n_iter, final_err] = Weight_Meridian(samples, LearnRate, ErrThres)


% It computes the weights under a Merdidian assumption for

% the prior distribution, paramater estimation of meridian case by Mayank Agrawal.
%Codes were modified from the Cauchy convolution scheme done Tao-Wan.

%

% Input:

% samples   : [TxN] , matrix whose COLUMNS are the samples

% LearnRate : [1x1] , the learning rate of the steepest ascent method

% ErrThres  : [1x1] , error threshold for the termination of the iterations

%                     of the steepest ascent method


% Output:

% w         : [Tx1] , weight vector

% p1        : [1x1] ,the mean of the meridian distribution

% p2        : [1x1] , the dispersion of the Meridian distribution

% n_iter    : [1x1] , number of iterations until convergence

% final_err : [1x1] , final convergence error

% Adapted from the code written prviously by Tao wan for Cauchy convolution
% Copyright (c)Mayank Agrawal 19/03/2010. 
% ma6879@bristol.ac.uk

[T,N] = size(samples);



%%%%%%%%%%%%%%%%%%%%%%%%%
%    Steepest Ascent    %
%%%%%%%%%%%%%%%%%%%%%%%%%

%--- Initialization

e_vec = ones(T,1); 

w_prev = e_vec/T;

n = 1;

err = Inf;

        %Parameter for Meridian case...
        m = mean(log(abs((w_prev'*samples)+eps)));
        g = exp(m);
        
        m1 = mean(samples(1,:));
        m2 = mean(samples(2,:));
        

        %% Dispersion parameter for meridian distribution
        
        m1_1 = mean(log(abs(samples(1,:)+eps)));
        m2_1 = mean(log(abs(samples(2,:)+eps)));                
        g1 = exp(m1_1);
        g2 = exp(m2_1);
   

        %--- Update

        while (n<=N && err>ErrThres)

            x = samples(:,n);
   

            u = w_prev'*x;
            
            
            Q = w_prev(1)*g1 + w_prev(2)*g2;
            Q1 = (u - w_prev(1)*m1 - w_prev(2)*m2)^2 + (w_prev(1)*g1 + w_prev(2)*g2)^2+ (2*(u - w_prev(1)*m1 - w_prev(2)*m2)*(w_prev(1)*g1 + w_prev(2)*g2)) ;
            Q2 = (x(1,1)-m1)*w_prev(1) + (x(2,1)-m2)*w_prev(2);
            Q3 = 2*(x(1,1)-m1)*((2*w_prev(1)*g1)+(w_prev(2)*g2))+(2*w_prev(2)*g1*(x(2,1)-m2));
            Q4 = 2*(x(2,1)-m2)*((2*w_prev(2)*g2)+(w_prev(1)*g1))+(2*w_prev(1)*g2*(x(1,1)-m1));
            
            w_new(1,1) = w_prev(1) + LearnRate*((-g1/Q) + (2*Q2*(x(1,1)-m1) + 2*g1*g2*w_prev(2)+2*g1^2*w_prev(1)+Q3)/Q1);
            w_new(2,1) = w_prev(2) + LearnRate*((-g2/Q) + (2*Q2*(x(2,1)-m2) + 2*g1*g2*w_prev(1)+2*g2^2*w_prev(2)+Q4)/Q1);

            w_new = abs(w_new)/(e_vec'*abs(w_new));     
            n = n+1;

            err = norm(w_new-w_prev);

            w_prev = w_new;           
                                

        end

        p1 = m;
        p2 = g2;

     
w = w_new;

n_iter = n-1;

final_err = err;
end