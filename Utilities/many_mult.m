function output = many_mult(varargin)
% This function just accumulates a bunch of bsxfun calls
output = ones(size(varargin{1})); 

for i = 1:length(varargin)
    try
    output = bsxfun(@times,output,varargin{i}); 
    catch err_mult
        i
        error(err_mult); 
    end
        
end