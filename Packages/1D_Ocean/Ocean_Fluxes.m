function OCEAN = Ocean_Fluxes(OCEAN,OPTS,FSTD)
% This routine updates the ocean component of the model based on the model
% developed by Petty et al (2013). 

% In this routine, we calculate:
% dHdt - the time rate of change of the ocean mixed layer depth
% dT/dt - the time rate of change of the ocean mixed layer temperature
% dS/dt - the time rate of change of the ocean mixed layer salinity
% T_s - the ocean surface layer temperature
OCEAN.rho = OCEAN.EOS(OCEAN.T,OCEAN.S);

turbvel = .1; 

OCEAN.tau = (OCEAN.U_a + turbvel)^2 * OCEAN.rho_a; 

OCEAN.ustar_i = sqrt((OCEAN.tau*OCEAN.rho_a/OCEAN.rho) * OCEAN.CD_i); 
OCEAN.ustar_o = sqrt((OCEAN.tau*OCEAN.rho_a/OCEAN.rho) * OCEAN.CD_o);
OCEAN.ustar =   sqrt((OCEAN.tau*OCEAN.rho_a/OCEAN.rho) * (FSTD.conc * OCEAN.CD_i + (1-FSTD.conc) * OCEAN.CD_o));

% Here we calculate what the temperature of the ocean surface will be
% according to the balance of surface heat fluxes
% The surface experiences three types of heat forcings
% 1) Radiative heating at the surface
% 2) Turbulent heating from the mixed layer
% 3) Exchange with the sea ice

% Formula for the ocean saturation specific humidity
OCEAN.qsat = @(T) .622 * OCEAN.pv(T) ./ (OCEAN.P_a - .378 * OCEAN.pv(T));

% The heating from the ocean is parameterized as linear turbulent exchange
% from a water at OCEAN.T to the surface. Positive means warming. 
F_surf_ml = @(T) OCEAN.rho * OCEAN.cw * OCEAN.ustar_o * (OCEAN.T - T);

% The flux from the mixed layer through to the sea ice is the same,
% assuming the sea ice is at its freezing temp. This occurs in the "lead
% region" of area Al. 
F_surf_ic = @(T) OCEAN.rho * OCEAN.cw * OCEAN.ustar_o * (OCEAN.Tfrz - T);

SH_out = @(T) OCEAN.do_SH * OCEAN.rho_a * OCEAN.c_a * OCEAN.CD_o * OCEAN.U_a * (T - OCEAN.T_a);

if OCEAN.presc_evap
    
    LH_out = @(T) OCEAN.do_LH * OCEAN.Evap_presc * OPTS.rho_water * OCEAN.L_v;

else
    
    LH_out = @(T) OCEAN.do_LH * OCEAN.rho_a * OCEAN.L_v * OCEAN.CD_o * OCEAN.U_a * (OCEAN.qsat(T) - OCEAN.q_a);
    
end

LW_out = @(T) OCEAN.epsilon * OPTS.sigma * (T + 273.14).^4;
LW_in = OCEAN.epsilon * OCEAN.LW;
SW_in = (1 - OCEAN.alpha) * ( 1 - OCEAN.Io) * OCEAN.SW;

%%
DUMMY.lead_width = OPTS.r_p; 

% Calculate the area of the lead region and open water region
% [OCEAN.Al,OCEAN.Ao] = calc_lead_area(FSTD,DUMMY);

% Calculate the side area shared b/w the sea ice and the surface layer
OCEAN.Aside = calc_side_area(FSTD,OPTS); 


% The surface heat flux budget over ocean. 
F_surf_at =  @(T) ... % Positive implies warming
    LW_in + ... % Longwave absorption at surface 
    SW_in - ... % Shortwave absorption at surface
    LH_out(T) - ... % Latent heat flux from surface
    SH_out(T) - ... % Sensible heat flux out of surface
    LW_out(T);    % Longwave heat flux out of surface
    
phi = (1-FSTD.conc);

hf_balance = @(T) ... % Net heat flux balance for surface layer. Negative means cooling
    phi*F_surf_at(T) ... % Surface radiative flux balance.
    + phi*F_surf_ml(T) ... % Heating from mixed layer turbulence.
    + (OCEAN.Aside) * F_surf_ic(T); % Heating going to ice

try

OCEAN.T_s = fzero(hf_balance,OCEAN.T);

catch err
   
    disp(err);
    
end

% Now we have the temperature of the surface ocean. If it is permissible,
% the open water heat flux that leads to pancake formation is 0. 
OCEAN.Q_o = 0; 

% On the other hand, if it is too cold, we need to freeze new floes. 
if OCEAN.T_s < OCEAN.Tfrz

    % The open water heat flux is the balance of fluxes at the freezing
    % point
    OCEAN.Q_o = phi * hf_balance(OCEAN.Tfrz); 

    % The surface temperature is now the freezing temperature    
    OCEAN.T_s = OCEAN.Tfrz; 
    
end

% All fluxes are calculated at the new surface temperature
OCEAN.Q_LH = LH_out(OCEAN.T_s); 
OCEAN.q = OCEAN.qsat(OCEAN.T_s); 
OCEAN.Q_SH = SH_out(OCEAN.T_s); 
OCEAN.Q_LW_out = LW_out(OCEAN.T_s); 
OCEAN.Q_SW_in = SW_in; 
OCEAN.Q_LW_in = LW_in;


%% 
% Now we convert from the fluxes we have evaluated, which are in a region
% of only open water, and compute what they mean in terms of the sea ice
% model. 

% This is the exchange of heat from the mixed layer below to the surface
% layer. Per unit of grid. 
OCEAN.Q_surf_ml = phi * F_surf_ml(OCEAN.T_s); 

% The lead heat flux is the exchange between the surface and the lead
% region. Per unit of grid. 
OCEAN.Q_lead = OCEAN.Aside * F_surf_ic(OCEAN.T_s); 

% This is the radiative/atmospheric exchange per unit of grid cell. The
% amount of heat coming from the atmosphere. 
OCEAN.Q_surf_at = phi * F_surf_at(OCEAN.T_s); 

% The exchange between the sea ice and the mixed layer. Per unit of grid
OCEAN.Q_mi = FSTD.conc * OCEAN.rho * OCEAN.cw * OCEAN.ch * OCEAN.ustar_i * (OCEAN.T - OCEAN.Tfrz); 