alpha = 2; 
xm = 10; 

k = @(alpha) 1/(alpha - 1);

sigma = @(alpha,k,xm) k^(-1 - 1/k) * xm ^(-alpha) / (alpha - 1); 


kk = k(alpha); 
sig = sigma(alpha,kk,xm);
theta = sig/kk;

pd = makedist('GeneralizedPareto','k',kk,'sigma',sig,'theta',theta)

R = 