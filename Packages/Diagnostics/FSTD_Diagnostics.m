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
    DIAG.FSTD.conc(diag_ind) = sum_FSTD(FSTD.psi,FSTD.one,0); % The ice concentration at each timestep
    DIAG.FSTD.Rmeanarea(diag_ind) = sum_FSTD(FSTD.psi,FSTD.Rmid',1); % The mean floe size (area-weighted,normalized)
    DIAG.FSTD.Rmeannum(diag_ind) = sum_FSTD(FSTD.NumberDist,FSTD.Rmid',1); % The mean floe size (number-weighted,normalized)
    DIAG.FSTD.Hmean(diag_ind) = sum_FSTD(FSTD.psi,FSTD.Hmid,1); % The mean ice thickness (area-weighted, normalized)
    DIAG.FSTD.Vtot(diag_ind) = sum_FSTD(FSTD.psi,FSTD.Hmid,0); % The total ice volume (area-weighted)
    DIAG.FSTD.Hmax(diag_ind) = FSTD.H_max; % The mean ice thickness (area-weighted, normalized)
    DIAG.FSTD.Amax(diag_ind) = sum(FSTD.psi(:,end)); % Maximum Ice Thickness Category Area
    
end

if diag_ind > 1 % Only do this once we've started
    
    %% Evaluate Thermodynamic Diagnostics
    
    if THERMO.DO
        
        DIAG.THERMO.Q_lead(diag_ind) = THERMO.Q_lead;
        DIAG.THERMO.Q_lat(diag_ind) = THERMO.Q_lat;
        DIAG.THERMO.Q_o(diag_ind) = THERMO.Q_o;
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
        
        
    end
    
    %% Evaluate Mechanical Diagnostics
    if MECH.DO && MECH.mag~=0 % Require something to happen (mag ~= 0) in order to run MECH_timestep, so don't make diagnostics otherwise
        
        
    end
    
    
    
    if WAVES.DO && EXFORC.stormy(ind) == 1
        
    end
    
    
    if OCEAN.DO && THERMO.DO
        
        
    end
    
end
