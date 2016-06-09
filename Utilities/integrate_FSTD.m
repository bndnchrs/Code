function integral = integrate_FSTD(psi,func,dA,doaverage)
% This code evaluates the sum
% normalize * \int \int (psi_spec(x,y,...) * func(x,y) dA)
% where normalize = 1 if doaverage = 0
% normalize = 1/\int \int psi(x,y) dA if doaverage = 1

% 
if doaverage
normalize = (squeeze(sum(sum(bsxfun(@times,psi,dA),1),2))).^(-1); 
else
    normalize = 1; 
end

% The interior of the integral
func = bsxfun(@times,func,dA);
func = bsxfun(@times,psi,func); 

integral = normalize .* squeeze(sum(sum(func,1),2)); 


end
