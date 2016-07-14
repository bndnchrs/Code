% This routine updates the ocean component of the model

dens = OCEAN.EOS(OCEAN.T,OCEAN.S);

prefac = OCEAN.cp_w * OCEAN.H * dens;


if ~THERMO.DO
    
    error('Ocean heating on but not thermodynamics package. We have to quit')
    
    % If not, we can get the heat flux
    % from somewhere else. Not included yet so we error out. 
    
end

%% Calculate pancake growth if cooling

% We want to calculate whether we will form pancakes or not. This replaces
% the calculation of THERMO.pancakes in Thermo_Timestep.m

% If we do, we first calculate how much heat would be required to cool the
% water to its freezing point (OCEAN.q_to_frz) in a time dt_sub. 

OCEAN.q_to_frz= (prefac /OPTS.dt_sub) * (OCEAN.Tfrz - OCEAN.T);

% It should always be negative, sicne Tfrz < T. If it is not, then we
% assume it is all basically frozen already. 

if OCEAN.q_to_frz > 0
    OCEAN.q_to_frz = 0; 
end

% If the heat flux that cools the mixed layer is greater in magnitude than
% this required value, the amount larger than this leads to the formation
% of ice pancakes. 

if THERMO.Q_o <= OCEAN.q_to_frz
    OCEAN.panQ = THERMO.Q_o - OCEAN.q_to_frz;
    OCEAN.Q_o = OCEAN.q_to_frz;
else
    OCEAN.panQ = 0; 
    OCEAN.Q_open = THERMO.Q_o;
end

% Ocean.Q is now the heating applied to the ocean

% There is a heat flux due to restoring to the initial temperature

OCEAN.Q_rest = -(prefac/OCEAN.lambda_rest) * (OCEAN.T - OCEAN.T_rest);

%% Calculate Temperature Tendency

% The time rate of change of ocean temperature is now the sum of three heat
% fluxes

% OCEAN.Q_open - The heating from 
% OCEAN.Q_rest - the restoring heat flux to the deep or lower latitude
% waters
% OCEAN.Qoi - The heat flux exchanged between the ocean and the sea ice

OCEAN.dTdt = (1/prefac) * (OCEAN.Q_open + OCEAN.Q_rest);

% If there is some heat flux to the pancakes, we make them out of it. 

if OCEAN.panQ < 0 
    % This is assuming they are all the same size and thickness. 
    OCEAN.pancakes = -OCEAN.panQ / (OPTS.L_f * OPTS.rho_ice * OPTS.h_p);
else
    OCEAN.pancakes = 0; 
end

% Pancake Growth array - same as in Thermo_Timestep.m
OCEAN.pancake_growth = 0*FSTD.meshR;
OCEAN.pancake_growth(THERMO.panloc_r,THERMO.panloc_h) ... 
    = OCEAN.pancakes*OPTS.dt_sub;

% In case we just have a single thickness category, these pancakes go there
OCEAN.dV_max_pancake = sum(OCEAN.pancake_growth(:,end)/OPTS.dt_sub)*OPTS.h_p;
OCEAN.V_max_in =  OCEAN.dV_max_pancake;

% The total dFSD/dt from the ocean. This is usually zero unless there are
% pancakes
OCEAN.diff = (1/OPTS.dt_sub) * (OCEAN.pancake_growth);
OCEAN.opening = -sum(OCEAN.diff(:));

%% Calculate Salinity Tendency
% This must happen as the last item in the code FSTD_Timestep. This is
% because the update to salinity comes from the change in the total ice
% volume and therefore needs FSTD.diff

% We calculate changes in salinity from changes in ice volume. if pancakes
% get formed, we need to add them in here. This is volume per unit time. 
FSTD.dV_ice = sum_FSTD(FSTD.diff + OCEAN.diff,FSTD.Hmid,0);

% Time rate of change of salinity
OCEAN.dSdt = OCEAN.S * (OPTS.rho_ice/dens) * (1 / OCEAN.H) * FSTD.dV_ice;