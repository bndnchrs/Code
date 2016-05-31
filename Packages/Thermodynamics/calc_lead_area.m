function [Al,Ao] = calc_lead_area(FSTD,THERMO)

lead_reg = integrate_FSTD(FSTD.psi,(((FSTD.Rmid + THERMO.lead_width).^2) ./ FSTD.Rmid.^2)',FSTD.dA,0);
conc = integrate_FSTD(FSTD.psi,FSTD.one,FSTD.dA,0); 

Al = lead_reg - conc; 

% Al = (sum_FSTD(FSTD.psi,((FSTD.Rmid' + THERMO.lead_width).^2) ./ (FSTD.Rmid'.^2),0)) - sum_FSTD(FSTD.psi,FSTD.one,0); 

% Total open water
openwater = 1 - integrate_FSTD(FSTD.psi,FSTD.one,FSTD.dA,0); 

% Lead region is value up to the open water fraction. If lead region exceed
% open water fraction, it is simple equal to that
Al = min(Al,openwater);
% If Al is smaller than zero, it is equal to zero. This shouldn't happen. 
Al = max(0,Al); 
% The open water that can form pancakes is what is left
Ao = openwater - Al; 
% It should not be less than zero, though this can only happen if Al is
% greater than openwater, which also should not happen!
Ao = max(0,Ao); 


end