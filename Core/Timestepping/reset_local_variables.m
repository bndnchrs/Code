%% reset_local_variables
% This code handles the resetting of matrices which are updated during each
% sub-timestep.

if FSTD.DO
   
   FSTD.diff = 0*FSTD.diff; 
   FSTD.opening = 0*FSTD.opening; 
   FSTD.V_max = integrate_FSTD(FSTD.psi(:,end),FSTD.meshHmid(:,end),FSTD.dA(:,end),0);
   FSTD.NumberDist = FSTD.psi./(pi*FSTD.meshRmid.^2);
   FSTD.V_max_in = 0;
   FSTD.V_max_out= 0; 
   
end

%%

if THERMO.DO
    
    if THERMO.mergefloes
        THERMO.gain_merge = 0*FSTD.psi;
        THERMO.loss_merge = 0*FSTD.psi; 
    end
    
end

if MECH.DO
    
    % Large-Scale In and Out matrices
    MECH.In = 0*MECH.In;
    MECH.Out = 0*MECH.Out;
    
    % Variables related to the changes from single processes
    MECH.In_ridge = MECH.In;
    MECH.In_raft = MECH.In;
    MECH.Out_ridge = MECH.Out;
    MECH.Out_raft = MECH.Out;
    MECH.Out_max_ridge = zeros(length(FSTD.Rmid),1);
    MECH.Out_max_raft = MECH.Out_max_ridge;
    MECH.In_max_ridge = MECH.Out_max_ridge;
    MECH.In_max_raft = MECH.Out_max_ridge;
    MECH.In_max = MECH.In_max_raft;
    MECH.Out_max = MECH.Out_max_raft;
    
end

if WAVES.DO
    
    WAVES.In = 0*WAVES.In;
    WAVES.Out = 0*WAVES.Out; 
    
end




