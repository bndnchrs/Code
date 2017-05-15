%Power Function Distribution Maximum Likelihood Estimator
%Based on equation in Table 1, modified from Evans et al. 2000 with
%correction.

function exponent = mle_power(data,xmax)

exponent=(log(xmax)-mean(log(data)))^-1-1;