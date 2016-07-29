%% FD_timestep_thermo

% This routing calculates the tendency due to one thermodynamic timestep
% Determine what portion of the heat flux goes where over open water

[FSTD.latSA,THERMO.Q_open,THERMO.Q_o,OCEAN.Q_oi,THERMO.Q_lead,THERMO.Q_lat,THERMO.Q_bas] = partition_heat_flux(FSTD,OPTS,THERMO,OCEAN,EXFORC);

% This returns the following partitioning of a net heat flux EXFORC.Q_oc
% That is the total heating per square meter of water.

% Q_o - The heat flux impinging directly onto open water
% Q_lead - The heat flux that affects the development of sea ice
% Q_oi - The heat flux exchanged between the ocean and the ice at the base
% All of these are in units of W/m^2 of grid. So they are net heat fluxes.

% Q_lead is partitioned into:
% Q_lat - the heat flux that leads to the freezing/melting of ice laterally
% Q_bas - the heat flux that leads to the freezing/melting of ice
% Both are again W/square meter of grid.

% The total energy entering the system is equal to Q_lead + Q_o, which has
% units W/m^2 of grid.

%% Do Horizontal Growth Rate from Lead Heat FLux

% Horizontal floe growth rate
% We take the total heat flux Q_lat, and divide by the surface area over
% which it is applied (FSTD.latSA). This is a zero-layer model where all
% the heat it transformed into ice.

THERMO.drdt = -THERMO.Q_lat/(OPTS.rho_ice*OPTS.L_f*FSTD.latSA);

% Since this growth also increases area, we have to include this.

THERMO.edgegrowth = (2./FSTD.meshRmid).*FSTD.psi*THERMO.drdt;

%% Do Pancake Growth Rate from open water heat flux
if THERMO.Q_o < 0 
    
    % If we are freezing, and we aren't simulating the ocean, we just add
    % ice at the surface of the ocean in a simple way.
    
    THERMO.pancakes = -THERMO.Q_o / (OPTS.L_f * OPTS.rho_ice * FSTD.Hmid(THERMO.panloc_h) * FSTD.dA(THERMO.panloc_r,THERMO.panloc_h));
    
else
    
    % If we are heating, or the ocean is on, we don't do this quite yet,
    % and instead we do it within the OCEAN code.
    THERMO.pancakes = 0;
    
end

% Just make a matrix so we can add this to the other tendencies.

THERMO.pancake_growth = 0*FSTD.meshR;

% Put it in the locations labeled by panloc_r,panloc_h. If these aren't set
% in initialization, they default to the smallest size/thickness.

THERMO.pancake_growth(THERMO.panloc_r,THERMO.panloc_h) ...
    = THERMO.pancakes;

%% Now Do Vertical Thermodynamics and determine the time rate of change of ice thickness

% Using Semtner Thermo, simple 1-D model using net heat flux to the ice
% without longwave. This is calculated locally at each ice floe.
if THERMO.dosemtner
    
    % dhdt_surf - The time rate of change of ice thickness for each thickness
    % category due to surface melting
    % dhdt_base - The time rate of change of ice thickness for each thickness
    % category dut to basal melting
    % Q_cond - The conductive heat flux through ice at each thickness category
    % Q_vert - The net heat flux at the ice surface for each thickness category
    % T_ice - The internal ice temperature calculated for each thickness
    % category
    
    THERMO = semtner_1D_thermo_petty(FSTD,OPTS,THERMO,OCEAN,EXFORC);
    
    % This allows us to calculate the long-wave outgoing heat flux for each ice
    % thickness category
    THERMO.Q_lw = OPTS.sigma * (THERMO.T_ice + 273.14).^4;
    
    % Now we have the net heat flux at the ice surface. This should
    % be zero unless there is surface melting
    EXFORC.Q_ic = THERMO.surf_HF;
    
else % Just conserve energy
    
    % If we are fixing the heat flux to be constant
    if THERMO.fixQ
        %%
        % Then we know the energy entering the system is fixed. There won't be
        % a surface heat flux, so these terms should subtract to zero.
        % Here, EXFORC.Q_oc is the heat flux per square meter of grid.
        THERMO.Q_vert = EXFORC.Q_oc - THERMO.Q_lat - THERMO.Q_bas - THERMO.Q_o;
        
        % Time rate of change of surface melting should be zero in this case.
        THERMO.dhdt_surf = 0*[FSTD.H] - THERMO.Q_vert / (FSTD.conc * OPTS.L_f * OPTS.rho_ice);
        
    else
        
        THERMO.Q_vert = 0;
        THERMO.dhdt_surf = 0*FSTD.H;
        
    end
    
    %%
    % Just set these values to be zero since we aren't worrying about a
    % conductive heat flux/internal ice temperature
    THERMO.Q_cond =  0*FSTD.H;
    THERMO.T_ice =  0*FSTD.H;
    
    % In which case we assume no net heating at the ice surface
    THERMO.Q_lw = 0;
    EXFORC.Q_ic = 0;
    
    % The time rate of change of thickness is just the base heat flux (net,
    % so units W/m^2 of grid) divided by the area it is applied to - the
    % basal area = the concentration and L_f rho.
    THERMO.dhdt_base = -THERMO.Q_bas / (FSTD.conc * OPTS.L_f * OPTS.rho_ice);
    
end
% The net time rate of change of ice thickness is the sum of surface and
% basal components

THERMO.dhdt = (THERMO.dhdt_surf + THERMO.dhdt_base);

%% Now we actually compute the tendencies in each size-thickness category


% If there is no concentration, we can't change the thickness, so we set it
% to zero.

if FSTD.conc == 0
    THERMO.dhdt = zeros(size(THERMO.dhdt));
end

% drdt, dhdt set up advective tendencies in size and thickness space. We
% turn these into matrices of size of FSTD.psi.
% Advective Tendencies

v_r = repmat(THERMO.drdt,size(FSTD.psi));
v_h = repmat(THERMO.dhdt,[length(FSTD.Rint),1]);

% This is the simple upwind advection scheme for doing this. It conserves
% volume and area as in Hibler (1980). It does this by calculating terms
% like df(r) = v_r(r) * f(r) where r represent the edges of each box and
% f(r) is the spectrum (units 1/m^2) and not psi itself.

if THERMO.drdt == 0
 %   disp('zero drdt')
end

[THERMO.adv_tend,THERMO.meltoutR,THERMO.meltoutH] = advect2_upwind(FSTD.psi,FSTD.dA,FSTD.Rint,FSTD.H,FSTD.dR,FSTD.dH,v_r,v_h,THERMO.allow_adv_loss_R,THERMO.allow_adv_loss_H);

% Make sure we aren't advecting and dividing by infinity somehow

if(isnan(sum(THERMO.adv_tend(:))))
    OPTS.dt_sub
    FSTD.eflag = 1;
    disp('isnan advtend');
    
end

% The time rate of change is each term. Each one is unitless since we have
% incorporated time into them, so we divide by the temporary timestep

THERMO.diff = THERMO.adv_tend + THERMO.pancake_growth + THERMO.edgegrowth;

% The opening is the sum of the difference, naturally
THERMO.opening = -sum(THERMO.diff(:).*FSTD.dA(:));


%% Handle volume losses at the thickest floe size category
% There will be no advective tendency at the thickest ice class due to this
% ice becoming thicker since we don't have a thicker ice thickness
% category.

% We want to calculate the volume of ice that is added into the thickest
% ice category. There are four sources:

% 1: Advection (when freezing) of thinner ice up.
% 2: Growth/Melting of existing ice.
% 3: Edge growth of the thickest ice.
% 4: Pancakes that form of the thickest size. This is probably never gonna
% happen.

% We must keep track of this seperately because it is a flux from thinner
% ice, equal to something like dh/dt * c(j-1), and not dhdt * c(j).
% When the ice is melting, we don't have such a.

% This part of the volume change will only occur if there is other ice that
% becomes part of this category. Otherwise the change in volume is
% accounted for by the terms dV_max_basal and dV_max_edge which incorporate
% volume changes to existing thickest ice.

%if ~isempty(FSTD.H) && sum(THERMO.adv_tend(:,end))>0
    THERMO.dV_max_adv = sum(THERMO.adv_tend(:,end).*FSTD.dA(:,end)*FSTD.H_max) ;
%else
    % There is no advective flux between thickness classes as there are no
    % thickness classes
%    THERMO.dV_max_adv = 0;
% end


%%
% Added volume in freezing/melting directly to the thickess ice class is
% equal to dhdt * psi(r,h_max).
% if length(FSTD.H) > 1
%     % We compute the time rate of change of this volume
if THERMO.dhdt(:,end) > 0
    
    THERMO.dV_max_basal = sum(FSTD.psi(:,end).*FSTD.dA(:,end)*THERMO.dhdt(:,end));

else
    % This is accounted for in the advection scheme, which removes the ice
    % from the domain.
    THERMO.dV_max_basal = 0;
    %
end

% The ice added is accounted for by the edge growth term
THERMO.dV_max_edge = 2*FSTD.H_max*sum(FSTD.psi(:,end).*FSTD.dA(:,end)./FSTD.Rmid')*THERMO.drdt;

% If there is pancake growth into the largest floe size category, account
% for that. This is usually untrue.
THERMO.dV_max_pancake = sum(THERMO.pancake_growth(:,end).*FSTD.dA(:,end)/OPTS.dt_sub)*OPTS.h_p;

% There may also be loss due to advection out in size space from the
% smallest size category
THERMO.dV_max_meltout = -0*THERMO.meltoutR(end) * FSTD.dA(1,end) * FSTD.H_max;

% The increase in volume in the largest category is the sum
THERMO.V_max_in = THERMO.dV_max_basal + THERMO.dV_max_adv + THERMO.dV_max_edge + THERMO.dV_max_pancake + THERMO.dV_max_meltout;

% Independent of sign
THERMO.V_max_out = 0;
