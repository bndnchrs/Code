%% update_global_variables
% This routine updates the variables that change each large-scale timestep.
% This includes the making of diagnostic data
FSTD.A_max = sum(FSTD.psi(:,end));

if FSTD.A_max > eps
    FSTD.H_max = FSTD.V_max / FSTD.A_max;
else
    FSTD.H_max = FSTD.H_max_i;
end

if isnan(FSTD.A_max)
    FSTD.H_max = FSTD.H_max_i; 
end
