%% OCEAN_FORCING

if OCEAN.DO && (~isfield(THERMO,'fixed_Q') || ~THERMO.fixQ)
    
    % need to make sure that we aren't fixing the heat flux
    
    OCEAN.U_a = EXFORC.UATM(FSTD.i);
    
    OCEAN.q_a = EXFORC.QATM(FSTD.i);
    
    OCEAN.SW = EXFORC.QSW(FSTD.i);
    
    OCEAN.LW = EXFORC.QLW(FSTD.i);
    
    OCEAN.T_a = EXFORC.TATM(FSTD.i);
    
    OCEAN.P_a = EXFORC.PATM(FSTD.i); 
    
    OCEAN.Precip = EXFORC.PRECIP(FSTD.i); 
    
end
