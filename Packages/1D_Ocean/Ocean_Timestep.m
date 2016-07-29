%% Now compute power budget

% Heat loss of mixed layer

OCEAN.Q_ml_SW = (1 - OCEAN.alpha) * OCEAN.Io * (1 - exp(-OCEAN.kappa_w * OCEAN.H_ml)) * OCEAN.SW;

OCEAN.Q_ml_out = ...
    OCEAN.Q_surf_ml + ... % heat to surface ocean
    OCEAN.Q_mi - ... % heat to surface ice
    OCEAN.Q_ml_SW; % heat from SW

% Thermodynamic change of ice volume
FSTD.dV_ice = integrate_FSTD(THERMO.diff,FSTD.Hmid,FSTD.dA,0);

% Evaporation from latent heat flux
OCEAN.Evap =  OCEAN.Q_LH / (OCEAN.rho_a * OCEAN.L_v); 

% Salinity loss of mixed layer
OCEAN.S_ml_out = ...
    (OPTS.rho_ice/OCEAN.rho) * (OCEAN.S_i - OCEAN.S)*FSTD.dV_ice ... % Change of salinity from ice formation/melting 
    + (1 - FSTD.conc) * (OCEAN.Precip - OCEAN.Evap) * OCEAN.S; % Change of salinity from precip/evap

g = 9.81; 

OCEAN.w = (EXFORC.Hml(FSTD.i+1) - EXFORC.Hml(FSTD.i))/(OPTS.dt);

OCEAN.deltaT = (OCEAN.T - OCEAN.T_b(OCEAN.H_ml));
OCEAN.deltaS = (OCEAN.S - OCEAN.S_b(OCEAN.H_ml));

wflag = 0; 
if OCEAN.w > 0
    wflag = 1; 
end

% Turbulent Heating from below
if OCEAN.compute_turb_deep
OCEAN.w_turb = OCEAN.kappa_turb / (OCEAN.H_ml);
else
    OCEAN.w_turb = 0; 
end

% Time evolution of salinity/temp
OCEAN.Q_base_mix = OCEAN.rho * OCEAN.cp_w * ...
    (OCEAN.w_turb + wflag * OCEAN.w) * ...
    (OCEAN.T_b(OCEAN.H_ml) - OCEAN.T); 

OCEAN.S_base_mix = (OCEAN.w_turb + wflag * OCEAN.w)* ...
    (OCEAN.S_b(OCEAN.H_ml) - OCEAN.S);

OCEAN.dSdt = (1 / OCEAN.H_ml) * ...
    (-OCEAN.S_ml_out + OCEAN.S_base_mix); 

OCEAN.dTdt = 1 / (OCEAN.rho * OCEAN.cp_w * OCEAN.H_ml) * ...
    (- OCEAN.Q_ml_out + OCEAN.Q_base_mix); 


if ADVECT.DO && OCEAN.advect_oc
    
    OCEAN.dTdt = OCEAN.dTdt + OCEAN.ustar * (OCEAN.T_ml_in - OCEAN.T_ml); 
    OCEAN.dSdt = OCEAN.dSdt + OCEAN.ustar * (OCEAN.S_ml_in - OCEAN.S_ml); 

end
