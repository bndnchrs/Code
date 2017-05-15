function div = KLdiv(f,g)
% Computes the Kullback-Leibler Divergence for the pdfs f and g

f = f/sum(f); 
g = g/sum(g); 

f = reshape(f,[1 length(f)]); 
g = reshape(g,[1 length(g)]); 

measure = f .* log10 (f./g); 
div = nansum(measure); 

end


