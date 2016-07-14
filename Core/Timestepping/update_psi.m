
%% update_psi
% Updated 12/8/2015 - Chris Horvat

% This script updates the floes distribution using the calculated timestep
% dt_temp. It also updates the open water fraction and
% concentration, and removes a residual for floe sizes possessing less than
% a small percentage of the ice cover.

% New FSTD
FSTD.psi = FSTD.psi + OPTS.dt_temp*FSTD.diff;

% How much volume is in the largest floe class.
V_max = FSTD.V_max + OPTS.dt_temp*FSTD.dV_max;

% How much area is there. We will repeat these steps later, but this is
% merely to ensure the highest thickness class isn't altered by the
% residual removal we do in a second.
A_max = integrate_FSTD(FSTD.psi(:,end),1,FSTD.dA(:,end),0);

if A_max ~= 0
    H_max = V_max / A_max;
else
    H_max = FSTD.H_max_i;
end

%% Take away areas with very small concentrations
% In case we overshoot due to rounding errors, we adjust these
% (this will happen once in a while)

% There is an issue here, which is what happens when the incoming area to
% the thickest ice category is small. Then this is eliminated by this
% residual processing. However the total volume entering this category is
% not reduced in the same way, so this needs to be addressed.

% All values that have less than machine precision concentration
FSTD.resid_adjust = FSTD.psi.*(abs(FSTD.psi.*FSTD.dA) < 1e-8);
% Delete them
FSTD.psi = FSTD.psi - FSTD.resid_adjust;

% This is the new area of the highest thickness class
A_max_new = integrate_FSTD(FSTD.psi(:,end),1,FSTD.dA(:,end),0);
% We have computed H_max = V_max / A_max
% But now A_max is reduced. Therefore we need to reduce V_max as well,
% which we will do by reducing FSTD.dV_max.
FSTD.dV_max = ((A_max_new/A_max) * V_max - FSTD.V_max)/OPTS.dt_temp;

% Update the number distribution
FSTD.NumberDist = FSTD.psi./(pi*FSTD.meshRmid.^2);

% Reset the open water and concentration
FSTD.conc = integrate_FSTD(FSTD.psi,1,FSTD.dA,0);
FSTD.openwater = 1 - FSTD.conc;

if OCEAN.DO && THERMO.DO
    
    OCEAN.T = OCEAN.T + OPTS.dt_temp * OCEAN.dTdt;
    OCEAN.S = OCEAN.S + OPTS.dt_temp * OCEAN.dSdt;
    OCEAN.H_ml = OCEAN.H_ml + OPTS.dt_temp * OCEAN.w; 
    
end