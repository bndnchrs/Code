% Ocean Init
% this function initializes the ocean component of the FSTD model


if ~isfield(OCEAN,'Hml')
    
    OCEAN.Hml = 100; % Mixed Layer Depth of the ocean
    
end

if ~isfield(OCEAN,'S')
    
    OCEAN.S = 33; % Mixed Layer Salinity in psu
    
end

if ~isfield(OCEAN,'lambda_ll')
    % The restoring timescale to deep ocean values
    OCEAN.lambda_ll = 7*86400; 

end

if ~isfield(OCEAN,'T')
    
    OCEAN.T = -1; % Mixed Layer Temp
    
end


if ~isfield(OCEAN,'Tfrz')
    
    OCEAN.Tfrz = -1.8; % Freezing Temperature of seawater
    
end

if ~isfield(OCEAN,'ustar_oceice')
    
    OCEAN.ustar_oceice = 1e-4; % m/s. Ice-ocean roughness velocity
    
end

if ~isfield(OCEAN,'cp_w')
    
    OCEAN.cp_w = 4185; % Specific Heat Capacity of Water (J/kg K)
    
end

if ~isfield(OCEAN,'alpha_T')
    
    OCEAN.alpha_T = 2e-4; % 1/deg C Thermal expansion coeff of seawater
    
end

if ~isfield(OCEAN,'beta_S')
    
    OCEAN.beta_S = 7.4e-4; % 1/psu Haline contraction coeff of seawater
    
end

if ~isfield(OCEAN,'T_0')
    OCEAN.T_0 = 0; % Deg C Reference Temperature for linear EOS
end

if ~isfield(OCEAN,'S_0')
    OCEAN.S_0 = 32; % ppt Reference Salinity for linear EOS
end

if ~isfield(OCEAN,'EOS')
    % This function is the linear equation of state. This is the default. 
    OCEAN.EOS = @(T,S) OPTS.rho_water * (1 - OCEAN.alpha_T * (T - OCEAN.T_0) ...
        + OCEAN.beta_S * ( S - OCEAN.S_0));  
end

if ~isfield(OCEAN,'oi_hf')
    % We calculate a heat flux between the ocean and the ice using a
    % turbulent velocity exchange parameterization (with u^*). If we don't 
    % want to do this, this flag should be set to zero. Default is to 
    % actually have this flux
    OCEAN.oi_hf = 1; 
end

if ~isfield(OCEAN,'taui')
    
    OCEAN.taui = .5*86400; % Relaxation timescale of the ice temperature
    
end

if ~isfield(OCEAN,'StrainInvar') && ~ADVECT.DO
    % The ocean strain rate tensor, 2 by nt long, for doing advection by
    % the ocean currents. If we don't have advection we can just let this
    % be equal to zero so there is no convergence.  
    OCEAN.StrainInvar = zeros(2,OPTS.nt); 

end

if ~isfield(OCEAN,'Uvel') && ADVECT.DO
    % The zonal velocity matrix. The situation may be diagrammed like
    % pack ice ----- MIZ ---- open water
    % Therefore the zonal velocity directs floes from the pack to the MIZ
    % to the open water
    OCEAN.Uvel = zeros(1,OPTS.nt);  
end

if ~isfield(OCEAN,'Vvel') && ADVECT.DO
    % The meridional velocity matrix. The situation may be diagrammed like
    % pack ice ----- MIZ ---- open water
    % Therefore the meridional velocity doesn't direct floes in or out of
    % the MIZ
    OCEAN.Vvel = zeros(1,OPTS.nt);  
end
