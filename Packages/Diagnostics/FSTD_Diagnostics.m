%% FD_Diagnostics
% This file, executed every global timestep, computes diagnostics for
% saving. It is only called when DIAG.DO

diag_ind = FSTD.i+1; % When FSTD.i = 1, this is the first timestep. However, DIAG(1) is the initial value. So FSTD.i=1 corresponds to DIAG(2);

%% Evaluate Regular Diagnostics

if FSTD.DO
    
    DIAG.FSTD.psi(:,:,diag_ind) = FSTD.psi; % The full FSTD
    DIAG.FSTD.diff(:,:,diag_ind) = FSTD.diff; % The total change per timestep for all components
    DIAG.FSTD.time(diag_ind) = FSTD.time_now; % The model time
    DIAG.FSTD.numSC(diag_ind) = OPTS.numSC; % The total number of sub-intervals in each timestep
    DIAG.FSTD.conc(diag_ind) = integrate_FSTD(FSTD.psi,FSTD.one,FSTD.dA,0); % The ice concentration at each timestep
    DIAG.FSTD.Rmeanarea(diag_ind) = integrate_FSTD(FSTD.psi,FSTD.Rmid',FSTD.dA,1); % The mean floe size (area-weighted,normalized)
    DIAG.FSTD.Rmeannum(diag_ind) = integrate_FSTD(FSTD.NumberDist,FSTD.Rmid',FSTD.dA,1); % The mean floe size (number-weighted,normalized)
    DIAG.FSTD.Hmean(diag_ind) = integrate_FSTD(FSTD.psi,FSTD.Hmid,FSTD.dA,1); % The mean ice thickness (area-weighted, normalized)
    DIAG.FSTD.Vtot(diag_ind) = integrate_FSTD(FSTD.psi,FSTD.Hmid,FSTD.dA,0); % The total ice volume (area-weighted)
    DIAG.FSTD.Hmax(diag_ind) = FSTD.H_max; % The mean ice thickness (area-weighted, normalized)
    DIAG.FSTD.Amax(diag_ind) = integrate_FSTD(FSTD.psi(:,end),1,FSTD.dA(:,end),0); % Maximum Ice Thickness Category Area
    
end

if diag_ind > 1 % Only do this once we've started
    
    %% Evaluate Thermodynamic Diagnostics
    
    if THERMO.DO
        
        DIAG.THERMO.Q_lead(diag_ind) = THERMO.Q_lead;
        DIAG.THERMO.Q_lat(diag_ind) = THERMO.Q_lat;
        DIAG.THERMO.Q_o(diag_ind) = THERMO.Q_o;
        DIAG.THERMO.Q_open(diag_ind) = THERMO.Q_open;        
        DIAG.THERMO.Q_bas(diag_ind) = THERMO.Q_bas;
        DIAG.THERMO.dV(diag_ind) = sum_FSTD(THERMO.diff,FSTD.Hmid,0) + ...
            + THERMO.dV_max_basal;
        DIAG.THERMO.dVmax_basal(diag_ind) = THERMO.dV_max_basal;
        DIAG.THERMO.drdt(diag_ind) = THERMO.drdt;
        DIAG.THERMO.dc_adv(diag_ind) = sum(THERMO.adv_tend(:));
        DIAG.THERMO.dc_pan(diag_ind) = sum(THERMO.pancakes(:));
        DIAG.THERMO.dc_edge(diag_ind) = sum(THERMO.edgegrowth(:));
        DIAG.THERMO.dc_tot(diag_ind) = sum(THERMO.diff(:));
        
        DIAG.THERMO.dhdt(:,diag_ind) = THERMO.dhdt;  
        DIAG.THERMO.Tice(:,diag_ind) = THERMO.T_ice; 
        DIAG.THERMO.Q_cond(:,diag_ind) = THERMO.Q_cond; 
        
        DIAG.EXFORC.Q_oc(diag_ind) = EXFORC.Q_oc; 
        DIAG.EXFORC.Q_ic_noLW = EXFORC.Q_ic_noLW;
        DIAG.THERMO.diffnet(diag_ind) = sum(abs(THERMO.diff(:)));  
        
    end
    
    %% Evaluate Mechanical Diagnostics
    if MECH.DO && MECH.mag~=0 % Require something to happen (mag ~= 0) in order to run MECH_timestep, so don't make diagnostics otherwise
        
        DIAG.MECH.mag(diag_ind) = MECH.mag; 
        DIAG.MECH.epsI(diag_ind) = MECH.eps_I; 
        DIAG.MECH.epsII(diag_ind) = MECH.eps_II; 
        DIAG.MECH.diffnet(diag_ind) = sum(abs(MECH.diff(:)));  
        
        
    end
    
    if ADVECT.DO
        
        DIAG.ADVECT.diffnet(diag_ind) = sum(abs(ADVECT.diff(:))); 
        
    end
    
    if WAVES.DO && EXFORC.stormy(FSTD.i) == 1

        DIAG.WAVES.Omega(:,:,diag_ind) = WAVES.Omega; 
        DIAG.WAVES.tau(diag_ind) = WAVES.tau;
        DIAG.WAVES.In(:,:,diag_ind) = WAVES.In; % The full FSTD
        DIAG.WAVES.Out(:,:,diag_ind) = WAVES.Out; % The total change per timestep for all components
        DIAG.WAVES.diffnet(diag_ind) = sum(abs(WAVES.diff(:)));  
        
    end
    
    
    if OCEAN.DO && THERMO.DO
        
        DIAG.OCEAN.T(diag_ind) = OCEAN.T;
        DIAG.OCEAN.S(diag_ind) = OCEAN.S;
        DIAG.OCEAN.pancakes(diag_ind) = sum(OCEAN.pancakes(:));
        DIAG.OCEAN.Q_open(diag_ind) = OCEAN.Q_open;
        DIAG.OCEAN.Qrest(diag_ind) = OCEAN.Q_rest;
        DIAG.OCEAN.Q_to_ice(diag_ind) = OCEAN.Q_oi;
        
    end
    
end
