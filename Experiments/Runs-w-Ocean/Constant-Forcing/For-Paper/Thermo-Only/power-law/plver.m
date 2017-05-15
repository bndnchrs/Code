%function
     x = [500,150,90,81,75,75,70,65,60,58,49,47,40]; 
 function   [alpha, xmin, D] = plfit(x); 
 function   [a_err, xm_err, nt_err] = plvar(x); 
 function   [p, gof] = plpva(x, xmin);
% PLVER is Matlab code used to verify that the algorithms gave 
% real / sensible results. The code will also plot the distributions and 
% the fitted power laws. 
% In all these programs, the code treats integer values as integers so 
% long as the minimum value is < 1000, or if there are < 100 integer
% observations (a criteria I just added; previously, if the minimum value 
% was > 1000, it would use the continuous approximation;
% this is fine except in the small-sample regime, when finite-size effects 
% in the continuous estimator make it less accurate).
%   
%   (description missing)
%   
%   THIS MIGHT BE IRRELEVANT -- 
%   
%    The fitting procedure works as follows:
%    1) For each possible choice of x_min, we estimate alpha via the 
%       method of maximum likelihood, and calculate the Kolmogorov-Smirnov
%       goodness-of-fit statistic D.
%    2) We then select as our estimate of x_min, the value that gives the
%       minimum value D over all values of x_min.
%
%    Note that this procedure gives no estimate of the uncertainty of the 
%    fitted parameters, nor of the validity of the fit.
%
%    Example:
%       x = (1-rand(10000,1)).^(-1/(2.5-1));
%	[alpha, xmin, D]        = plfit(x);
%	[a_err, xm_err, nt_err] = plvar(x);
%	[p, gof]                = plpva(x, xmin);
%    For more information, try 'type plver'
%
%    See also PLVAR, PLPVA, PLFIT

% Version 1.0 (May, 2007)
% Copyright (C) 2007 Aaron Clauset (Santa Fe Institute)
% Distributed under GPL 2.0
% http://www.gnu.org/copyleft/gpl.html
% PLFIT comes with ABSOLUTELY NO WARRANTY
%%%%%%%%%%%%%%%%%%%%%%%% ABOVE FROM DRW.
% -- dataset 1
%%x = [500,150,90,81,75,75,70,65,60,58,49,47,40,40];
%%[alpha, xmin, D]        = plfit(x);
%%[a_err, xm_err, nt_err] = plvar(x);
%%[p, gof]                = plpva(x, xmin);

%%fprintf('alpha   = %4.2f +- %4.2f\nxmin    = %4.2f +- %4.2f\nntail   = %4.2f +- %4.2f\nD       = %4.2f\n', ...
 %%  alpha,a_err,xmin,xm_err,sum(x>=xmin),nt_err,D);
%%fprintf('p-value = %4.2f +- 0.03\n',p);

% alpha   =  2.71 +- 0.48
% xmin    =    47 +- 7.6
% ntail   = 12.00 +- 1.79
% D       =  0.15
% p-value =  0.58 +- 0.03

x = reshape(x,numel(x),1);
q = unique(x); n = length(x);
c = hist(x,q)'./n;
c = [[q; q(end)+1] 1-[0; cumsum(c)]]; c(find(c(:,2)<10^-10),:) = [];
cf = ([xmin:q(end)]'.^-alpha)./(zeta(alpha) - sum([1:xmin-1].^-alpha));
cf = [[xmin:q(end)+1]' 1-[0; cumsum(cf)]];
cf(:,2) = cf(:,2) .* c(find(c(:,1)==xmin),2);

figure(1);
loglog(c(:,1),c(:,2),'bo','MarkerSize',8,'MarkerFaceColor',[1 1 1]); hold on;
loglog(cf(:,1),cf(:,2),'k--','LineWidth',2); hold off;
set(gca,'FontName','Times','FontSize',16);
