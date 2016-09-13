%% Now compute power budget

% Heat loss of mixed layer. 
OCEAN.Q_ml_SW = (1-FSTD.conc) * (1 - OCEAN.alpha) ...
    * OCEAN.Io * (1 - exp(-OCEAN.kappa_w * OCEAN.H_ml)) * OCEAN.SW;

OCEAN.Q_ml_out = ...
    OCEAN.Q_surf_ml + ... % heat to surface ocean
    OCEAN.Q_mi - ... % heat to surface ice
    OCEAN.Q_ml_SW; % heat from SW

% Thermodynamic change of ice volume
OCEAN.dV_ice = integrate_FSTD(THERMO.diff,FSTD.Hmid,FSTD.dA,0);

% Evaporation from latent heat flux
OCEAN.Evap =  OCEAN.Q_LH / (OPTS.rho_water * OCEAN.L_v); 

OCEAN.S_ml_ice = (OPTS.rho_ice/OCEAN.rho) * (OCEAN.S_i - OCEAN.S)*OCEAN.dV_ice;
OCEAN.S_ml_precip = (1 - FSTD.conc) * OCEAN.Precip * OCEAN.S; 
OCEAN.S_ml_evap = (1 - FSTD.conc) * OCEAN.Evap * OCEAN.S;

% Salinity loss of mixed layer 
OCEAN.S_ml_out = ...
    OCEAN.S_ml_ice ... % Change of salinity from ice formation/melting 
    + OCEAN.S_ml_precip ... % Change of salinity from precip (more means freshening)
    - OCEAN.S_ml_evap; % Change from evap (more means saltier)

g = 9.81; 

OCEAN.w = (EXFORC.Hml(FSTD.i+1) - EXFORC.Hml(FSTD.i))/(OPTS.dt);

OCEAN.deltaT = (OCEAN.T - OCEAN.T_b(OCEAN.H_ml));
OCEAN.deltaS = (OCEAN.S - OCEAN.S_b(OCEAN.H_ml));

wflag = 0; 
if OCEAN.w > 0
    wflag = 1; 
end

% Turbulent Heating from below
rho_oc = OCEAN.EOS(OCEAN.T,OCEAN.S); 
rho_b = OCEAN.EOS(OCEAN.T_b(OCEAN.H_ml),OCEAN.S_b(OCEAN.H_ml)); 

% If the ml density is higher than the deep density, we have convection.
% We effect this by increasing the eddy diffusivity by an order of
% magnitude, from kappa_conv to kappa_turb
if rho_b < rho_oc
    OCEAN.kappa_calc = OCEAN.kappa_conv; 
else
    OCEAN.kappa_calc = OCEAN.kappa_turb; 
end

if OCEAN.compute_turb_deep
    OCEAN.w_turb = OCEAN.kappa_calc / (OCEAN.H_ml);
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
