function dist = CH_dist(f,g)

% Distance metric concocted by Chris Horvat
% f and g are probability distributions that integrate to 1

f = f / sum(f); 
g = g/sum(g); 

dist = nansum(abs(f-g)./g);

end