%% update_psi
% Updated 12/8/2015 - Chris Horvat

% This script updates the floes distribution using the calculated timestep 
% dt_temp. It also updates the open water fraction and
% concentration, and removes a residual for floe sizes possessing less than
% a small percentage of the ice cover.

% New FSTD
FSTD.psi = FSTD.psi + OPTS.dt_temp*FSTD.diff;

%% Take away areas with very small concentrations
% In case we overshoot due to rounding errors, we adjust these
% (this will happen once in a while)

% All values that have less than machine precision concentration
resid_adjust = FSTD.psi.*(abs(FSTD.psi) < 1e-8);
% Delete them
FSTD.psi = FSTD.psi - resid_adjust;


% Update the number distribution
FSTD.NumberDist = FSTD.psi./(pi*FSTD.meshR.^2);
% Update the spectrum
FSTD.psi_spec = FSTD.psi ./ FSTD.dA; 

% Reset the open water and concentration
FSTD.conc = sum_FSTD(FSTD.psi,FSTD.one,0); 
FSTD.openwater = 1 - FSTD.conc;

if OCEAN.DO && THERMO.DO
    
    OCEAN.T = OCEAN.T + OPTS.dt_temp * OCEAN.dTdt;
    OCEAN.S = OCEAN.S + OPTS.dt_temp * OCEAN.dSdt;
    
end