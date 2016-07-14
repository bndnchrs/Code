%% Take away areas with very small concentrations
% In case we overshoot due to rounding errors, we adjust these
% (this will happen once in a while)

% All values that have less than machine precision concentration
% We can't due this for the thickest category. The reason is that we know
% exactly how much volume is going in, and how much area. If we suddenly
% reduce the area, this means the thickness of that category is going to be
% much higher than it should be

% Therefore this has been removed temporarily
  resid_adjust = FSTD.psi.*(abs(FSTD.psi.*FSTD.dA) < 1e-9);
% Delete them
  FSTD.psi = FSTD.psi - resid_adjust;


% Update the number distribution
FSTD.NumberDist = FSTD.psi./(pi*FSTD.meshRmid.^2);

% Reset the open water and concentration
FSTD.conc = integrate_FSTD(FSTD.psi,1,FSTD.dA,0); 
FSTD.openwater = 1 - FSTD.conc;

if OCEAN.DO && THERMO.DO
    
    OCEAN.T = OCEAN.T + OPTS.dt_temp * OCEAN.dTdt;
    OCEAN.S = OCEAN.S + OPTS.dt_temp * OCEAN.dSdt;
    
end