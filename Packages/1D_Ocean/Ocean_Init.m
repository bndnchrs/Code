% Ocean Init
% this function initializes the ocean component of the FSTD model

if ~isfield(OCEAN,'Z')
    
    OCEAN.Z = 1:500;
    
end

if ~isfield(OCEAN,'H_ml')
    
    OCEAN.H_ml = EXFORC.Hml(1); % Mixed Layer Depth of the ocean
    
end

if ~isfield(OCEAN,'S')
    
    OCEAN.S = 32.5; % Mixed Layer Salinity in psu
    
end

if ~isfield(OCEAN,'lambda_rest')
    % The restoring timescale to deep ocean values
    OCEAN.lambda_rest = 7*86400;
    
end

if ~isfield(OCEAN,'compute_turb_deep')
    OCEAN.compute_turb_deep = 0; 
end

if ~isfield(OCEAN,'kappa_turb')
    % The turbulent exchange with the deep layer below. Pick length scale
    % of 25 meters, time scale of 7 days. 
    OCEAN.kappa_turb = .001;
    
end

if ~isfield(OCEAN,'Tfrz')
    
    OCEAN.Tfrz = -1.8; % Freezing Temperature of seawater
    
end


if ~isfield(OCEAN,'T')
    
    OCEAN.T = -1.79; % Mixed Layer Temp
    
end




if ~isfield(OCEAN,'T_rest')
    
    OCEAN.T_rest = OCEAN.Tfrz; % Freezing Temperature of seawater
    
end

if ~isfield(OCEAN,'ustar_oceice')
    
    OCEAN.ustar_oceice = 1e-4; % m/s. Ice-ocean roughness velocity
    
end

if ~isfield(OCEAN,'cp_w')
    
    OCEAN.cp_w = 4185; % Specific Heat Capacity of Water (J/kg K)
    
end

if ~isfield(OCEAN,'alpha_T')
    
    OCEAN.alpha_T = 6e-5; % 1/deg C Thermal expansion coeff of seawater
    
end

if ~isfield(OCEAN,'beta_S')
    
    OCEAN.beta_S = 8e-4; % 1/psu Haline contraction coeff of seawater
    
end

if ~isfield(OCEAN,'T_0')
    OCEAN.T_0 = 0; % Deg C Reference Temperature for linear EOS
end

if ~isfield(OCEAN,'S_0')
    OCEAN.S_0 = 34; % ppt Reference Salinity for linear EOS
end

if ~isfield(OCEAN,'EOS')
    % This function is the linear equation of state. This is the default.
    OCEAN.EOS = @(T,S) OPTS.rho_water * (1 - OCEAN.alpha_T * (T - OCEAN.T_0) ...
        + OCEAN.beta_S * ( S - OCEAN.S_0));
end

if ~isfield(OCEAN,'oi_hf')
    % We calculate a heat flux between the ocean and the ice using a
    % turbulent velocity exchange parameterization (with u^*). If we don't
    % want to do this, this flag should be set to zero. Default is to
    % actually have this flux
    OCEAN.oi_hf = 1;
end

if ~isfield(OCEAN,'taui')
    
    OCEAN.taui = .5*86400; % Relaxation timescale of the ice temperature
    
end

if ~isfield(OCEAN,'do_SH')
    OCEAN.do_SH = 1;
end

if ~isfield(OCEAN,'do_LH');
    OCEAN.do_LH = 0;
end


if ~isfield(OCEAN,'T_b')
    
    
    OCEAN.T_b = @(z) OCEAN.Tfrz + (z>200).*(z-200)/200; %+ 1.8*tanh((z)/50);
    
end

if ~isfield(OCEAN,'S_b')
    
    
    OCEAN.S_b = @(z) 33 + (z>20).*(1.4 + .2*(z)/500); % + 1*tanh((z)/30);
    
end

% Make sure we have external forcing fields that are required
if ~isfield(EXFORC,'QATM')
    EXFORC.QATM = zeros(1,OPTS.nt);
end

if ~isfield(EXFORC,'TATM')
    EXFORC.TATM = zeros(1,OPTS.nt);
end

if ~isfield(EXFORC,'UATM')
    EXFORC.UATM = zeros(1,OPTS.nt);
end

if ~isfield(EXFORC,'PATM')
    EXFORC.PATM = zeros(1,OPTS.nt) + 100; % kPA
end

if ~isfield(EXFORC,'QLW')
    EXFORC.QLW = zeros(1,OPTS.nt);
end

if ~isfield(EXFORC,'QSW')
    EXFORC.QSW = zeros(1,OPTS.nt);
end

if ~isfield(EXFORC,'PRECIP')
    EXFORC.PRECIP = zeros(1,OPTS.nt); % m/s
end

if ~isfield(OCEAN,'S_i')
    OCEAN.S_i = 5;
end

if ~isfield(OCEAN,'T_ml_in')
    OCEAN.T_ml_in = OCEAN.T_0; 
end

if ~isfield(OCEAN,'S_ml_in')
    OCEAN.S_ml_in = OCEAN.S_0; 
end

if ~isfield(OCEAN,'advect_oc')
    OCEAN.advect_oc = 0; 
end

% We want to use the petty (2013) mixed layer model.
OCEAN.rho_a = 1.275; %kg/m^3 - density of air
OCEAN.c_a = 1005; % J /kg K - specific heat capacity of air
OCEAN.cw = 4190; % " - specific heat capacity of water
OCEAN.CD_o = .001; % Turbulent transfer coeff for oce - atmosphere
OCEAN.CD_i = .0013; % " for ice-atmosphere
OCEAN.epsilon = .97; % emissivity of ocean
OCEAN.epsilon_i = 1; % emissivity of ice
OCEAN.alpha = .06; % albedo of ocean
OCEAN.alpha_i = .75; % albedo of ice
OCEAN.Io = .45; % Surface permissivity
OCEAN.L_v = 2.501 * 10^6; % latent heat of vaporisation
OCEAN.L_s = 2.834 * 10^6; % latent heat of sublimation
OCEAN.ch = .006; % Stanton number for transfer of heat b/w mixed layer and sea ice
OCEAN.kappa_w = .1; % Extinction coefficient of seawater
OCEAN.c1 = .8; % wind stirring transfer
OCEAN.H_w = 10; % depth of wind reach
OCEAN.c_m = .03; % background turbulence
% Formula for the partial pressure of water at temperature T
OCEAN.pv = @(T) 2.53 * 10^8 * exp(1).^(-5420 ./ (T + 273));
