%% Now compute power budget

% Heat loss of mixed layer
OCEAN.Q_ml_out = ...
    OCEAN.Q_surf_ml + ... % heat to surface ocean
    OCEAN.Q_mi - ... % heat to surface ice
    (1 - OCEAN.alpha) * OCEAN.Io * (1 - exp(-OCEAN.kappa_w * OCEAN.H_ml)) * OCEAN.SW; % heat from SW

% Thermodynamic change of ice volume
FSTD.dV_ice = integrate_FSTD(THERMO.diff,FSTD.Hmid,FSTD.dA,0);

% Evaporation from latent heat flux
OCEAN.Evap =  OCEAN.Q_LH / (OCEAN.rho_a * OCEAN.L_v); 

% Salinity loss of mixed layer
OCEAN.S_ml_out = ...
    (OPTS.rho_ice/OCEAN.rho) * (OCEAN.S_i - OCEAN.S)*FSTD.dV_ice ... % Change of salinity from ice formation/melting 
    + (1 - FSTD.conc) * (OCEAN.Precip - OCEAN.Evap) * OCEAN.S; % Change of salinity from precip/evap

g = 9.81; 



%% Power input from ocean buoyancy exchange
P_B = ...
    OCEAN.H_ml * g * OCEAN.alpha_T * OCEAN.Q_ml_out / (OCEAN.rho * OCEAN.cp_w) - ...
    OCEAN.H_ml * g * OCEAN.beta_S * OCEAN.S_ml_out; % Power by mixing or removal of buoyancy from domain

% Power is dissipated more effectively when the ocean is convecting
% From Tang (1991)
if P_B < 0
    OCEAN.Power_Buoy = P_B * 1; % Losing energy 
else
    OCEAN.Power_Buoy = P_B * .8; % Gaining energy
end

%% Power input from mixing by wind

OCEAN.Power_Wave = OCEAN.c1 * exp(-OCEAN.H_ml/OCEAN.H_w) * abs(OCEAN.ustar.^3); 

%% Power input from entrainment

% The jump in buoyancy across the mixed-layer base
OCEAN.deltaB = ...
    g * OCEAN.alpha_T * (OCEAN.T - OCEAN.T_b(OCEAN.H_ml)) - ... % Temperature
    g * OCEAN.beta_S  * (OCEAN.S - OCEAN.S_b(OCEAN.H_ml)); % Salinity gap

% The amount of entrainment needed to balance the two sources
OCEAN.w = ...
    (OCEAN.Power_Buoy + OCEAN.Power_Wave) * ... % The total input from buoyancy and wave
    (OCEAN.H_ml * OCEAN.deltaB + OCEAN.c_m.^2)^(-1);

OCEAN.Power_Entrain = OCEAN.w * (OCEAN.H_ml * OCEAN.deltaB + OCEAN.c_m.^2); 

%%

wflag = 0; 
if OCEAN.w > 0
    wflag = 1; 
end

% Time evolution of salinity/temp
OCEAN.dSdt = -OCEAN.S_ml_out / (OCEAN.H_ml) + ...
    wflag * OCEAN.w * (OCEAN.S_b(OCEAN.H_ml) - OCEAN.S)/OCEAN.H_ml; 

OCEAN.dTdt = -OCEAN.Q_ml_out / (OCEAN.rho * OCEAN.cp_w * OCEAN.H_ml) + ...
    wflag * OCEAN.w * (OCEAN.S_b(OCEAN.H_ml) - OCEAN.S)/OCEAN.H_ml; 
