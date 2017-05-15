function [VC_alpha,VC_xmin,VC_L,UN_alpha,UN_xmin,UN_n,VC_p,VC_gof] = calc_power_law_stuff(distribution,dR,R)

% This accepts any distribution f(r)dr and computes the maximum likelihood
% estimate of the power law slope using the code of Virkar and Clauset

addpath('power-law')

VC_gof = [];

% We form the distribution f(r)dr
VC_distribution = distribution; 

for i = 1:size(distribution,2)
    % For each timestep
    
    % Optional tracking
    if mod(i,10) == 1
        disp(sprintf('Number %d out of %d \n',i,size(distribution,2)));
    end
    % Take the distribution at each timestep
    dist_temp = VC_distribution(:,i);
    Rtemp = R; 
    % only consider the values for which it is positive
  %  Rtemp = R(dist_temp > 0);
  %  dist_temp = dist_temp(dist_temp > 0);
    
    %% Option 1 - treat the distribution as a histogram and use rldecode to get
    % A vector of datapoints
    %    AREA = (25*1000)^2; % 25 km domain
    %    dist_temp = AREA * dist_temp;
    %    dist_temp = round(dist_temp);
    %    Rtemp = Rtemp(dist_temp >= 1);
    %    dist_temp = dist_temp(dist_temp >= 1);
    %    vals = rldecode(dist_temp,Rtemp);
    
    
    %% Option 2 - use a rejection method to produce the distribution with a random
    % number generator. This seems to do a better job than does the rldecode
    % method. I don't know why.
    
    % Convert the distribution into a function from 0 to 1. Each r value is
    % assigned to a chunk of the space in (0,1). When a random number falls
    % into the space (r_i,r_i+1), it becomes an "observation" of a floe with
    % size r_i
    
    % Convert the distribution into one on the unit interval
    VC_D = cumsum(dist_temp)/sum(dist_temp);
    
    % Make many "observations"
    r = rand(1000,1);
    % Bin these observations over the categories with their assigned [0,1]
    % values
    [p,VC_data] = histc(r,VC_D);
    % Convert this histogram into a synthetic dataset
    
    vals = round(Rtemp(VC_data+1)); % our synthetic data    
    
    [VC_alpha(i), VC_xmin(i), VC_L(i)]=plfit(vals);
   %  plplot(vals, VC_xmin(i), VC_alpha(i));
    [UN_alpha(i), UN_xmin(i), UN_n(i)]=plvar(vals,'reps',200,'silent');
    [VC_p(i), VC_gof(i)] = plpva(vals,VC_xmin(i),'reps',200,'silent');
    
    
end

disp('done')



