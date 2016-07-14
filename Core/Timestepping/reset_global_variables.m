%% Reset_global_variables
% This code handles the resetting of counting variables (like the number of
% internal timesteps) that must be reset at each large-scale timestep. It
% also resets matrices which are used each timestep
% Timestep Markers
OPTS.dt_sub =OPTS.dt;
OPTS.dt_temp = OPTS.dt_sub;
OPTS.numSC = 0;

%% Reset the calculation of gamma_ridge/raft

if MECH.DO
    
    MECH.gamma_raft = calc_gamma_raft_FD(FSTD.Hmid,FSTD.meshHmid,MECH.H_raft);
    MECH.gamma_ridge = 1 -  MECH.gamma_raft;
    
    if MECH.rafting == 1 && MECH.ridging == 0
        
        MECH.gamma_ridge = 0*MECH.gamma_ridge;
        
        MECH.gamma_raft = 0*MECH.gamma_raft + 1;
        
    else if MECH.rafting == 0 && MECH.ridging == 1
            
            MECH.gamma_ridge = 0*MECH.gamma_ridge + 1;
            
            MECH.gamma_raft = 0*MECH.gamma_raft;
        end
        
    end
    
end

%% Re-set concentration, etc

% Reset the open water and concentration
FSTD.conc = integrate_FSTD(FSTD.psi,FSTD.one,FSTD.dA,0);
FSTD.openwater = 1 - FSTD.conc;

FSTD.NumberDist = FSTD.psi./(pi*FSTD.meshRmid.^2);
