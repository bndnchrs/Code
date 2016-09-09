function [Aside] = calc_side_area(FSTD,OPTS)
%% Compute the side area shared between the sea ice boundary layer
% and the open water surface layer

Per = 2*pi*(FSTD.meshRmid + OPTS.r_p) .* (FSTD.meshHmid + OPTS.h_p);
Num = FSTD.psi ./ (pi*FSTD.meshRmid.^2); 

Aside = integrate_FSTD(Num,Per,FSTD.dA,0); 


end