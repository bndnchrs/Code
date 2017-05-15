function dist = ET_dist(f,g)
% Distance metric concocted by Eli Tziperman
% f and g are probability distributions that integrate to 1

f = f / sum(f); 
g = g/sum(g); 

dist = sum(abs(f-g));

end


