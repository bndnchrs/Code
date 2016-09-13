%% FD_Diagnostics
% This file, executed every global timestep, computes diagnostics for
% saving. It is only called when DIAG.DO

diag_ind = FSTD.i+1; % When FSTD.i = 1, this is the first timestep. However, DIAG(1) is the initial value. So FSTD.i=1 corresponds to DIAG(2);

%% Evaluate Regular Diagnostics

if FSTD.DO
    
    DIAG.FSTD.psi(:,:,diag_ind) = FSTD.psi; % The full FSTD
    DIAG.FSTD.diff(:,:,diag_ind) = FSTD.diff; % The total change per timestep for all components
    DIAG.FSTD.dA(:,:,diag_ind) = FSTD.dA; % The full FSTD
    DIAG.FSTD.diff_FSD(:,diag_ind) = sum(FSTD.diff .* FSTD.dA,2);
    DIAG.FSTD.diff_ITD(:,diag_ind) = sum(FSTD.diff .* FSTD.dA,1);
    DIAG.FSTD.Hmid(:,diag_ind) = FSTD.Hmid;
    DIAG.FSTD.time(diag_ind) = FSTD.time_now; % The model time
    DIAG.FSTD.numSC(diag_ind) = OPTS.numSC; % The total number of sub-intervals in each timestep
    DIAG.FSTD.conc(diag_ind) = integrate_FSTD(FSTD.psi,FSTD.one,FSTD.dA,0); % The ice concentration at each timestep
    DIAG.FSTD.Rmeanarea(diag_ind) = integrate_FSTD(FSTD.psi,FSTD.Rmid',FSTD.dA,1); % The mean floe size (area-weighted,normalized)
    DIAG.FSTD.Rmeannum(diag_ind) = integrate_FSTD(FSTD.NumberDist,FSTD.Rmid',FSTD.dA,1); % The mean floe size (number-weighted,normalized)
    %    DIAG.FSTD.Rmeannum(diag_ind) = integrate_FSTD(FSTD.psi,FSTD.meshRmid./(pi * FSTD.meshRmid.^2),FSTD.dA,1); % The mean floe size (number-weighted,normalized)
    DIAG.FSTD.Hmean(diag_ind) = integrate_FSTD(FSTD.psi,FSTD.Hmid,FSTD.dA,1); % The mean ice thickness (area-weighted, normalized)
    DIAG.FSTD.Vtot(diag_ind) = integrate_FSTD(FSTD.psi,FSTD.Hmid,FSTD.dA,0); % The total ice volume (area-weighted)
    DIAG.FSTD.Ntot(diag_ind) = integrate_FSTD(FSTD.NumberDist,1,FSTD.dA,0); % The total ice volume (area-weighted)
    DIAG.FSTD.Hmax(diag_ind) = FSTD.H_max; % The mean ice thickness (area-weighted, normalized)
    DIAG.FSTD.Amax(diag_ind) = integrate_FSTD(FSTD.psi(:,end),1,FSTD.dA(:,end),0); % Maximum Ice Thickness Category Area
    DIAG.FSTD.resid_adjust(diag_ind) = integrate_FSTD(FSTD.resid_adjust,1,FSTD.dA,0);
    
    
end

if diag_ind > 1 % Only do this once we've started
    
    %% Evaluate Thermodynamic Diagnostics
    
    if THERMO.DO
        
        DIAG.THERMO.Q_lead(diag_ind) = THERMO.Q_lead;
        DIAG.THERMO.Q_lat(diag_ind) = THERMO.Q_lat;
        DIAG.THERMO.Q_o(diag_ind) = THERMO.Q_o;
        DIAG.THERMO.Q_open(diag_ind) = THERMO.Q_open;
        DIAG.THERMO.Q_bas(diag_ind) = THERMO.Q_bas;
        DIAG.THERMO.Q_vert(:,diag_ind) = THERMO.Q_vert;
        DIAG.THERMO.dV(diag_ind) = sum(THERMO.diff(:).*FSTD.dA(:).*FSTD.meshHmid(:)) + ...
            + THERMO.dV_max_basal;
        DIAG.THERMO.dVmax_basal(diag_ind) = THERMO.dV_max_basal;
        DIAG.THERMO.drdt(diag_ind) = THERMO.drdt;
        DIAG.THERMO.dc_adv(diag_ind) = sum(THERMO.adv_tend(:).*FSTD.dA(:));
        DIAG.THERMO.dc_pan(diag_ind) = sum(THERMO.pancakes(:).*FSTD.dA(:));
        DIAG.THERMO.dc_edge(diag_ind) = sum(THERMO.edgegrowth(:).*FSTD.dA(:));
        DIAG.THERMO.dc_tot(diag_ind) = sum(THERMO.diff(:).*FSTD.dA(:));
        
        DIAG.THERMO.dhdt(:,diag_ind) = THERMO.dhdt;
        DIAG.THERMO.Tice(:,diag_ind) = THERMO.T_ice;
        DIAG.THERMO.Q_cond(:,diag_ind) = THERMO.Q_cond;
        
        % DIAG.EXFORC.Q_oc(diag_ind) = EXFORC.Q_oc;
        % DIAG.EXFORC.Q_ic(:,diag_ind) = EXFORC.Q_ic;
        % DIAG.EXFORC.Q_ic_noLW(diag_ind) = EXFORC.Q_ic_noLW;
        DIAG.THERMO.diffnet(diag_ind) = sum(abs(THERMO.diff(:)));
        DIAG.THERMO.diff_FSD(:,diag_ind) = sum(THERMO.diff .* FSTD.dA,2);
        DIAG.THERMO.diff_ITD(:,diag_ind) = sum(THERMO.diff .* FSTD.dA,1);
        
    end
    
    if THERMO.DO && THERMO.dosemtner
        DIAG.THERMO.QSW(diag_ind) = OCEAN.SW*(1-OPTS.alpha_ic);
        DIAG.THERMO.QLW(:,diag_ind) = THERMO.Q_lw;
    end
    
    
    if THERMO.DO && THERMO.mergefloes && THERMO.drdt >= 0 && FSTD.conc <= 1
        
        DIAG.THERMO.alphamerge(diag_ind) = THERMO.alphamerge;
        DIAG.THERMO.mn1(diag_ind) = THERMO.Mn1;
        DIAG.THERMO.diffmerge(diag_ind) = sum(abs(THERMO.diff_merge(:).*FSTD.dA(:)));
        DIAG.THERMO.diff_FSD_merge(:,diag_ind) = sum(THERMO.diff_merge .* FSTD.dA,2);
        
    end
    
    %% Evaluate Mechanical Diagnostics
    if MECH.DO && MECH.mag~=0 % Require something to happen (mag ~= 0) in order to run MECH_timestep, so don't make diagnostics otherwise
        
        DIAG.MECH.mag(diag_ind) = MECH.mag;
        DIAG.MECH.epsI(diag_ind) = MECH.eps_I;
        DIAG.MECH.epsII(diag_ind) = MECH.eps_II;
        DIAG.MECH.diffnet(diag_ind) = sum(abs(MECH.diff(:).*FSTD.dA(:)));
        DIAG.MECH.diff_FSD(:,diag_ind) = sum(MECH.diff .* FSTD.dA,2);
        DIAG.MECH.diff_ITD(:,diag_ind) = sum(MECH.diff .* FSTD.dA,1);
        
        
    end
    
    if ADVECT.DO
        
        DIAG.ADVECT.diffnet(diag_ind) = sum(abs(ADVECT.diff(:).*FSTD.dA(:)));
        DIAG.ADVECT.dc_adv(diag_ind) = sum(ADVECT.diff(:).*FSTD.dA(:));
        DIAG.ADVECT.diff_FSD(:,diag_ind) = sum(ADVECT.diff .* FSTD.dA,2);
        DIAG.ADVECT.diff_ITD(:,diag_ind) = sum(ADVECT.diff .* FSTD.dA,1);
        
    end
    
    if WAVES.DO && EXFORC.stormy(FSTD.i) == 1
        
        DIAG.WAVES.Omega(:,:,diag_ind) = WAVES.Omega;
        DIAG.WAVES.tau(diag_ind) = WAVES.tau;
        DIAG.WAVES.In(:,:,diag_ind) = WAVES.In; % The full FSTD
        DIAG.WAVES.Out(:,:,diag_ind) = WAVES.Out; % The total change per timestep for all components
        DIAG.WAVES.diffnet(diag_ind) = sum(abs(WAVES.diff(:).*FSTD.dA(:)));
        DIAG.WAVES.diff_FSD(:,diag_ind) = sum(WAVES.diff .* FSTD.dA,2);
        DIAG.WAVES.diff_ITD(:,diag_ind) = sum(WAVES.diff .* FSTD.dA,1);
        
        
    end
    
    
    if OCEAN.DO && THERMO.DO
        
        DIAG.OCEAN.T(diag_ind) = OCEAN.T;
        DIAG.OCEAN.S(diag_ind) = OCEAN.S;
        DIAG.OCEAN.rho(diag_ind) = OCEAN.rho; 
        DIAG.OCEAN.W(diag_ind) = OCEAN.w;
        DIAG.OCEAN.q(diag_ind) = OCEAN.q; 
        
        if isfield('OCEAN','dopetty') && OCEAN.dopetty
            % we are doing this power budget thing
            DIAG.OCEAN.P_W(diag_ind) = OCEAN.Power_Wave;
            DIAG.OCEAN.P_B(diag_ind) = OCEAN.Power_Buoy;
            DIAG.OCEAN.P_E(diag_ind) = OCEAN.Power_Entrain;
            DIAG.OCEAN.db(diag_ind) = OCEAN.deltaB;
            
        end
        
        DIAG.OCEAN.QLH(diag_ind) = OCEAN.Q_LH;
        DIAG.OCEAN.QSH(diag_ind) = OCEAN.Q_SH;
        DIAG.OCEAN.QLW_out(diag_ind) = OCEAN.Q_LW_out;
        DIAG.OCEAN.QSW_in(diag_ind) = OCEAN.Q_SW_in;
        DIAG.OCEAN.QLW_in(diag_ind) = OCEAN.Q_LW_in;
        
        DIAG.OCEAN.Qlead(diag_ind) = OCEAN.Q_lead;
        DIAG.OCEAN.Qsurf_at(diag_ind) = OCEAN.Q_surf_at;
        DIAG.OCEAN.Qsurf_ml(diag_ind) = OCEAN.Q_surf_ml;
        
        DIAG.OCEAN.Q_mi(diag_ind) = OCEAN.Q_mi;
        DIAG.OCEAN.Q_ml_out(diag_ind) = OCEAN.Q_ml_out;
        DIAG.OCEAN.S_ml_out(diag_ind) = OCEAN.S_ml_out;
        DIAG.OCEAN.S_ml_ice(diag_ind) = OCEAN.S_ml_ice;
        DIAG.OCEAN.S_ml_precip(diag_ind) = OCEAN.S_ml_precip;
        DIAG.OCEAN.S_ml_evap(diag_ind) = OCEAN.S_ml_evap;
        
        
        DIAG.OCEAN.T_s(diag_ind) = OCEAN.T_s;
        
        DIAG.OCEAN.H_ml(diag_ind) = OCEAN.H_ml;
        
        DIAG.OCEAN.Evap(diag_ind) = OCEAN.Evap;
        
        DIAG.OCEAN.dSdt(diag_ind) = OCEAN.dSdt;
        DIAG.OCEAN.dTdt(diag_ind) = OCEAN.dTdt;
        
        DIAG.OCEAN.T_a(diag_ind) = OCEAN.T_a;
        
        DIAG.OCEAN.ustar(diag_ind) = OCEAN.ustar;
        DIAG.OCEAN.wturb(diag_ind) = OCEAN.w_turb;
        
        DIAG.OCEAN.Q_base_mix(diag_ind) = OCEAN.Q_base_mix;
        DIAG.OCEAN.S_base_mix(diag_ind) = OCEAN.S_base_mix;

        DIAG.OCEAN.Q_ml_SW(diag_ind) = OCEAN.Q_ml_SW; 
        DIAG.OCEAN.Aside(diag_ind) = OCEAN.Aside; 
        DIAG.OCEAN.dV_thermo(diag_ind) = OCEAN.dV_ice; 
        DIAG.OCEAN.kappa_turb(diag_ind) = OCEAN.kappa_calc; 
    end
    
end
