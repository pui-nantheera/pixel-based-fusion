function [w, p1, p2, n_iter, final_err] = Weight_PDFPrior_tao_cauchy2(samples, LearnRate, ErrThres)

% [w, n_iter, final_err] = CauchyPrior(samples, m, g, LearnRate, ErrThres, pdf_prior)

% It computes the weights under a Cauchy assumption for

% the prior distribution.

%

% Input:

% samples   : [TxN] , matrix whose COLUMNS are the samples

% LearnRate : [1x1] , the learning rate of the steepest ascent method

% ErrThres  : [1x1] , error threshold for the termination of the iterations

%                     of the steepest ascent method

% pdf_prior : string, - 'cauchy', for the Cauchy prior

%

% Output:

% w         : [Tx1] , weight vector

% p1        : [1x1] , Cauchy prior: the mean of the Cauchy distribution

% p2        : [1x1] , Cauchy prior: the dispersion of the Cauchy distribution

% n_iter    : [1x1] , number of iterations until convergence

% final_err : [1x1] , final convergence error
%Authhor: Tao Wan
%Modified by Mayank Agrawal, to be only used for Cauchy case...


[T,N] = size(samples);



%%%%%%%%%%%%%%%%%%%%%%%%%
%    Steepest Ascent    %
%%%%%%%%%%%%%%%%%%%%%%%%%

%--- Initialization

e_vec = ones(T,1);

w_prev = e_vec/T;

n = 1;

err = Inf;


m = mean(log(abs((w_prev'*samples)+eps)));

m1 = mean(samples(1,:));
m2 = mean(samples(2,:));



m1_1 = mean(log(abs(samples(1,:)+eps)));
m2_1 = mean(log(abs(samples(2,:)+eps)));
g1 = exp(m1_1);
g2 = exp(m2_1);


%--- Update

while (n<=N && err>ErrThres)
    
    x = samples(:,n);
    
    
    u = w_prev'*x;
    
    
    Q = w_prev(1)*g1 + w_prev(2)*g2;
    Q1 = (u - w_prev(1)*m1 - w_prev(2)*m2)^2 + (w_prev(1)*g1 + w_prev(2)*g2)^2;
    Q2 = (x(1,1)-m1)*w_prev(1) + (x(2,1)-m2)*w_prev(2);
    
    w_new(1,1) = w_prev(1) + LearnRate*((-g1/Q) + (2*Q2*(x(1,1)-m1) + 2*g1*g2*w_prev(2)+2*g1^2*w_prev(1))/Q1);
    w_new(2,1) = w_prev(2) + LearnRate*((-g2/Q) + (2*Q2*(x(2,1)-m2) + 2*g1*g2*w_prev(1)+2*g2^2*w_prev(2))/Q1);
    
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