function integral = sum_FSTD(psi,func,doaverage)
% This code evaluates the sum
% normalize * sum_x sum_y (psi(x,y,...) * func(x,y))
% where normalize = 1 if doaverage = 0
% normalize = 1/sum_x sum_y psi(x,y) if doaverage = 1

% This code is used when evaluating an integral over FSTD.psi, since
% FSTD.psi(r(i),h(j)) is the total area between sizes r(i) and r(i+1) and
% h(j) and h(j+1). sum(FSTD.psi(:)) is the concentration!

% For calculating integrals, utilize the code integrate_FSTD, which takes a
% distribution as input, along with dx,and dy. 

% The function returns integral, for which
% size(integral) is [size(psi,3) size(psi,4) ...]; 

normalize = 1;
    
if doaverage == 1; 
    % Calculate the concentration of psi
    normalize = 1./squeeze(sum(sum(psi,1),2));
    
end

% Evaluate the interior of the integral

% Take the sum over the inside two indices
integral = squeeze(sum(sum(bsxfun(@times,psi,func),1),2));

integral = normalize .* integral; 

end
